#[! /usr/bin/env python3
# coding: utf-8

'''
A module that analyzes a Hayes-and-Wilson-formatted LearningData.txt file (given a Features.txt file) to identify clusters of segments that have the distributions of complex segments.

creates a new Features.txt file and a new LearningData.txt file

checks the data against features and vice versa

To see more, for running from the command line:

$ python3 compseg.py --help

Or, import the module into python and invoke help there:

>>> import compseg
>>> help(compseg)

'''

#standard modules
import os 
import sys 
from re import sub as resub
import shutil


#custom modules

try:
    import pynatclasses as pnc
    import natclasscounter as nc
    import datachecker as dc
    import plot_insep as pins
    import messages as msg
except ModuleNotFoundError:
    import compseg.code.natclasscounter as nc
    import compseg.code.pynatclasses as pnc
    import compseg.code.datachecker as dc
    import compseg.code.plot_insep as pins
    import compseg.code.messages as msg


##########################################
#get original segments-features structure
##########################################

def old_feats(featfilepath, **kwargs):
    '''
    given a path to features.txt, returns
    (featurenames, segmentlines)
    '''
    return pnc.read_feat_file(featfilepath, **kwargs)


def get_bidir_transprobs(featpath, datapath, **kwargs):
    '''
    calculates the likelihood of any two segments preceding and following each other
    either as bidirectional transitional probabilities or using the inseparability measure
    returns a dictionary of seqs and numbers referred to in other functions as "clusters"
    '''
    return nc.bidir_prob_wrapper(featpath, datapath, **kwargs)

def find_diff_feats(clusters, featpath, threshold=1, **kwargs):
    '''
    finds features that two segs in a bigram share, and features that two segs in a bigram differ on. 
    takes as input a dictionary of bidirectional transitional probabilities.
    needs a path to the features file, Features.txt
    the threshold value determines which ngrams are considered. by default, only clusters whose bidirectional trans prob is equal or greater than 1 will be analyzed.
    '''
    segdic = pnc.segs_to_feats(pnc.read_feat_file(featpath, **kwargs), **kwargs)
    ngramlist = [x for x in clusters if clusters[x]>=threshold]
    differ = {}.fromkeys(ngramlist, [])
    same = {}.fromkeys(ngramlist, [])
    for ngram in ngramlist:
        segs = ngram.split(" ")
        differ[ngram] = list(set(segdic[segs[0]]) ^ set(segdic[segs[1]]))
        same[ngram] = list(set(segdic[segs[0]]) & set(segdic[segs[1]]))
    return (same, differ, segdic)

def analyze_cons_feats(same, differ, segdic):
    '''
    applies phonological logic to figure out what features the complex segments should be.
    affricates become -son, -cont, +strid
    prenasalized stops are +nas, -son and have other features of the 2nd half of sequence
    labiovelars have both place features (lab and velar or dorsal)

    returns a dictionary of new segments, along with their feature value vectors. that is, given:
    same == {'x y' : [+feat, -otherfeat], 'z q': ...}
    differ == {'x y': [+thirdfeat, -thirdfeat, ] }
    segdic: {'x': [+feat, +otherfeat],}
    it will return:
    {'xy': [+feat, -otherfeat, +thirdfeat ] }
    the choice of "thirdfeat" value is phonological.
    '''
    compsegs = {}
    for ngram in same:
        seg2_feats = segdic[ngram.split(" ")[1]]
        seg1_feats = segdic[ngram.split(" ")[0]]
        cseg = ngram.replace(" ", "")
        compsegs[cseg] = same[ngram]
        #affricates default to -cont, +strid
        if any([thing.startswith('-son') for thing in same[ngram]]):
            if any([thing.startswith('-cont') for thing in differ[ngram]]):
                compsegs[cseg].append([x for x in differ[ngram] if x.startswith('-cont')][0])
            if any([thing.startswith('+strid') for thing in differ[ngram]]):
                compsegs[cseg].append([x for x in differ[ngram] if x.startswith('+strid')][0])
            if any([thing.startswith('+del') for thing in differ[ngram]]):
                compsegs[cseg].append([x for x in differ[ngram] if x.startswith('+del')][0])
        #prenas segs default to +nas, -son -- though this allows for prenasalized sonorants in principle
        if any([thing.startswith('+nas') for thing in differ[ngram]]):
            for x in ['+nas', '-son']:
                compsegs[cseg].extend([y for y in differ[ngram] if y.startswith(x)])
        #labiovelars are both labial and velar
        if not any([thing.startswith('-cons') for thing in differ[ngram]]):
            if any([thing.startswith('+lab') for thing in differ[ngram]]):
                compsegs[cseg].append([x for x in differ[ngram] if x.startswith("+lab")][0])
            if any([thing.startswith('+vel') for thing in differ[ngram]]):
                compsegs[cseg].append([x for x in differ[ngram] if x.startswith("+vel")][0])
            elif any([thing.startswith('+dor') for thing in differ[ngram]]):
                compsegs[cseg].append([x for x in differ[ngram] if x.startswith("+dor")][0])
            for ngram in differ:
                feats_to_fix = set([x.lstrip('+-') for x in differ[ngram]])
                curr_feats = set([x.lstrip('+-') for x in compsegs[cseg]])
                for feat in feats_to_fix:
                    if not feat in curr_feats:
                        if ('+'+feat in compsegs[cseg]) or ('-'+feat in compsegs[cseg]):
                            continue
                        else:
                            if '+'+feat in seg2_feats:
                                compsegs[cseg].append('+'+feat)
                            elif '-'+feat in seg2_feats:
                                compsegs[cseg].append('-'+feat)

        #glides become secondary articulations; first seg's place and manner kept (i.e. all other feats)
        elif any([thing.startswith('-cons') for thing in differ[ngram]]):
            if '+round' in seg2_feats:
                compsegs[cseg].append('+round') 
            elif '-back' in seg2_feats:
                compsegs[cseg].append('-back')
            elif '-round' in seg2_feats:
                compseg[cseg].append('-round')
        #now take all other segs from seg1
            for ngram in differ:
                feats_to_fix = set([x.lstrip('+-') for x in differ[ngram]])
                curr_feats = set([x.lstrip('+-') for x in compsegs[cseg]])
                for feat in feats_to_fix:
                    if not feat in curr_feats:
                        if ('+'+feat in compsegs[cseg]) or ('-'+feat in compsegs[cseg]):
                            continue
                        else:
                            if '+'+feat in seg1_feats:
                                compsegs[cseg].append('+'+feat)
                            elif '-'+feat in seg1_feats:
                                compsegs[cseg].append('-'+feat)
        for thing in compsegs:
            compsegs[thing] = sorted(list(set(compsegs[thing])))
    return compsegs


def analyze_vocoid_feats(same, differ, segdic):
    '''
    this function takes in two dictionaries listing features in common and different features; it detects contradictory features and adds features to specify the new complex segs as diphthongs
    returns a dictionary of new segments, along with their feature value vectors. that is, given:
    same == {'x y' : [+feat, -otherfeat], 'z q': ...}
    differ == {'x y': [+thirdfeat, -thirdfeat, ] }
    segdic: {'x': [+feat, +otherfeat],}
    it will return:
    {'xy': [+feat, -otherfeat, +thirdfeat ] }
    the choice of "thirdfeat" value is based on vowel height; the features of the lower (more sonorous) vowel are preserved
    the vowel features this function can handle are very conservative, SPE-style.
    '''
    compsegs = {}
    #for glides, override their features in favor of vowels
    for ngram in same:
        seg2_feats = segdic[ngram.split(" ")[1]]
        seg1_feats = segdic[ngram.split(" ")[0]]
        cseg = ngram.replace(" ", "")
        compsegs[cseg] = same[ngram]
        #glides first: the vowel half wins
        if any([thing.startswith('-syll') for thing in differ[ngram]]):
            if any([thing.startswith('-syll') for thing in seg1_feats]):
                feats_to_keep = seg2_feats
            else:
                feats_to_keep = seg1_feats
        #vowel-vowel diphthongs: the lower ones win
        elif '+low' in differ[ngram]:
            if '+low' in seg1_feats:
                feats_to_keep = seg1_feats
            elif '+low' in seg2_feats:
                feats_to_keep = seg2_feats
        elif '+high' in differ[ngram]:
            if '+high' in seg1_feats:
                feats_to_keep = seg2_feats
            elif '+high' in seg2_feats:
                feats_to_keep = seg2_feats
        #and otherwise, keep the first one's feats (assume left-dominant dipthongs. this is arbitrary)
        else:
            feats_to_keep = seg1_feats 
        #resolve differences in favor of one value
        for ngram in differ:
                feats_to_fix = set([x.lstrip('+-') for x in differ[ngram]])
                curr_feats = set([x.lstrip('+-') for x in compsegs[cseg]])
                for feat in feats_to_fix:
                    if not feat in curr_feats:
                        if ('+'+feat in compsegs[cseg]) or ('-'+feat in compsegs[cseg]):
                            continue
                        else:
                            if '+'+feat in feats_to_keep:
                                compsegs[cseg].append('+'+feat)
                            elif '-'+feat in feats_to_keep:
                                compsegs[cseg].append('-'+feat)
    for thing in compsegs:
        compsegs[thing] = sorted(list(set(compsegs[thing])))
    return compsegs




def get_new_segs(featpath, compdic, **kwargs):
    '''
    takes in a feature file and a dictionary of clusters to be converted to complex segs.
    wrapper function.
    {'xy' : [-feat, +feat, ...]}
    '''
    fd = find_diff_feats(compdic, featpath, **kwargs)
    if 'vowels' in kwargs and kwargs['vowels']==True:
        return analyze_vocoid_feats(fd[0], fd[1], fd[2])
    else:
        return analyze_cons_feats(fd[0], fd[1], fd[2])

def make_new_feats(featpath, compdic, **kwargs): 
    '''
    gets a path to a feature file, and a dictionary of new complex segments and features that define them.
    returns a tuple with a list of feature names and a list of segments and their feature values, in order, ready to write to a feature file. 
    '''
    featlist = old_feats(featpath, **kwargs)
    segdic = pnc.segs_to_feats(featlist, **kwargs)
    if 'vowels' in kwargs:
        vowels=kwargs['vowels']
    else:
        vowels=False
    if vowels and not 'diph' in featlist[0]:
        featlist.append('diph')
    for seg in compdic:
        segdic[seg] = compdic[seg]
        #first get the missing (zero) feature values
        for feat in featlist[0]:
            if not ('+'+feat in segdic[seg]) or ('-'+feat in segdic[seg]):
                segdic[seg].append('0'+feat)
    #then, create an entry in feature dictionary for each segment (old and new)
    for seg in segdic:
        if vowels:
            segdic[seg].append('0diph')
        #first, fill in zero feat values for all the old segments
        for feat in featlist[0]:
            if not any([x for x in segdic[seg] if x.endswith(feat)]):
                segdic[seg].append('0'+feat)
    #append a list of all the new segments, and their feature values, to the tuple of features
    for seg in compdic:
        its_line = [seg]
        for pos in range(len(featlist[0])):
            its_line.append([x[0] for x in segdic[seg] if x[1:]==featlist[0][pos]][0])
        featlist[1].append(its_line)
        if vowels:
            featlist[1].append('+diph')
    return featlist

def write_new_feats(outpath, featlist):
    with open(outpath, 'w', encoding='utf-8') as f:
        f.write("\t" + "\t".join(featlist[0])+'\n')
        for line in featlist[1]:
            f.write('\t'.join(line) + '\n')


def write_new_ld_file(clusters, oldpath, newpath, threshold=1, **kwargs):
    '''
    clusters is created by get_new_segs above
    oldpath and newpath are the locations of LearningData.txt files, old and new 
    threshold defaults to 1 (inseparability measure)
    bigram clusters are sorted by their inseparability value;
    thus, if something that eventually becomes a trigram or tetragram has two-way parts in the current inseparability table, the bigram that is higher on the list will be replaced first.
    '''
    clustlist = sorted([x for x in clusters if clusters[x]>=threshold], key=clusters.get, reverse=True)
    with open(oldpath, 'r', encoding='utf-8') as f:
        with open(newpath, 'w', encoding='utf-8') as out:
            for line in f:
                word = line.strip()
                if '\t' in line:
                    word = line.split('\t')[0]
                    rest = line.split('\t')[1:]
                for clust in clustlist: 
                    x = r'(^|\s)'+(clust)+r'(\s|$)'
                    y = r'\1'+''.join(clust.split(" "))+r'\2'
                    word = resub(x, y, word)
                if '\t' in line:
                    out.write(f"word\trest")
                else:
                    out.write(word +'\n')
    msg.env_render(message=f"\n\nWrote modified learning data to {newpath.split('simulation')[1]}", **kwargs)
                            

def check_new_segs(newpath, oldfeats, newfeats, **kwargs): 
    '''
    given an interim data file, checks to see which segs from old feats are no longer in the new data.
    also checks which segs from the new feature file are not in the new data, because that can happen if the complex segment is really a trigram:
    old data:
    d z
    n d
    m b
    ... (but really, both 'nd' and 'ndzh' should be complex segs--see fijian)
    '''
    ldsegs = dc.collectLDSegs(newpath)
    oldfeatsegs = sorted(pnc.segs_to_feats(oldfeats, **kwargs).keys())
    newfeatsegs = sorted(pnc.segs_to_feats(newfeats, **kwargs).keys())
    missing_segs = sorted([x for x in ldsegs if not x in newfeatsegs])
    extra_segs = sorted(list(set(newfeatsegs) - set(ldsegs)))
    if missing_segs:
        msg.env_render(message=f"\n\nThe feature file is missing the following segments: \n{' '.join(missing_segs)}", **kwargs)
        return (missing_segs, 'missing')
    if extra_segs:
        kwargs['message']=f"\n\nThe following segments are in the feature file but are not in the data file, and will be removed from feature file:\n{' '.join(extra_segs)}"
        msg.env_render(**kwargs)
        return (extra_segs, 'extra')
    else:
        kwargs['message']="\n\nYour segments are all defined in the feature file."
        msg.env_render(**kwargs)
        return (False, False) 

def write_feats_w_check(newfeatpath, datapath, oldfeats, newfeats, skipcheck=False, **kwargs): 
    '''
    checks the data file against a new feature list (newfeats). this is a collection of featlines.
    if any segs are missing, raises an error and exits
    if any segs are extra, writes all but the extra segs to newfeatpath. 
    '''
    segs=None
    if not skipcheck:
        segs = check_new_segs(datapath, oldfeats, newfeats, **kwargs)
    with open(newfeatpath, 'w', encoding='utf-8') as f:
        f.write('\t'+'\t'.join(newfeats[0])+'\n')
        for line in newfeats[1]:
            if segs and segs[1]=='extra':
                if line[0] in segs[0]:
                    pass
                else:
                    f.write('\t'.join(line)+'\n')
            else:
                f.write('\t'.join(line)+'\n')
    return(newfeats) #just a path

def complexify(**kwargs):
    '''
    the algorithm.
    feats is a full path to Features.txt or another appropriately formatted feature file
    ld is a full path to a LearningData.txt or other learning data file, with space-separated words
    outdir is a writable directory where the learner will save results.
    the function creates a subfolder inside this directory, called 'simulation', and creates new versions of learning data and features. NOTE: any existing simulation directories will be deleted without warning.
    threshold is the cutoff for the inseparability measure. clusters above the threshold get converted into complex segments. defaults to 1.
    '''
    if not 'vowels' in kwargs:
        vowels=False
    else:
        vowels=kwargs['vowels']
    feats=kwargs.get('feats')
    ld=kwargs.get('ld')
    outdir=kwargs.get('outdir')
    threshold = kwargs.get('threshold')
    alpha = kwargs.get('alpha')
    ofpth = os.path.join(outdir, 'simulation_report.txt')
    kwargs['outfilepath']=ofpth
    msg.env_render(message="\nSearching for complex segments.", **kwargs)
    msg.env_render(message=f"\nInseparability threshold: {threshold}\nAlpha level for Fisher's Exact Test: {alpha}", **kwargs)
    oldfeats = pnc.read_feat_file(feats, outfilepath=ofpth) #before starting the recursion, first version of features is arg passed to the command
    msg.env_render(message="\nChecking feature file...\n", **kwargs)
    if (vowels==False) and (not any([feat.startswith('syll') for feat in oldfeats[0]])):
        message="\nYou need to have a 'syll(abic)' feature in your feature file. The learner needs a list of consonants, [-syll], to get started.\n"
        x = f"Your features are: [{','.join(oldfeats[0])}]"
        msg.env_render(message=message, **kwargs)
        if __name__=='__main__':
            raise SystemExit
        else:
            return message
    elif vowels and (not any([feat.startswith('cons') for feat in oldfeats[0]])):
        message = "\nYou need to have a 'cons(onantal)' feature in your feature file. The learner needs a list of vocoids, [-cons], to get started.\n"
        msg.env_render(message=message, **kwargs)
        if __name__=='__main__':
            raise SystemExit
    elif pnc.check_feats(oldfeats, **kwargs):
        if os.path.isdir(os.path.join(outdir, 'simulation')):
            shutil.rmtree(os.path.join(outdir, 'simulation'))
        os.mkdir(os.path.join(outdir, 'simulation'))
        step = 1
        while step:
            wdir = os.path.join(outdir, 'simulation', 'iteration'+str(step))
            os.mkdir(wdir)
            #get numbers from feats and ld files
            temp = get_bidir_transprobs(feats, ld, **kwargs)
            if temp[0]=={}:
                shutil.rmtree(wdir)
            else:
                nc.write_insep(temp, os.path.join(wdir, 'inseparability.txt'))
            clusters = {k:v for (k,v) in temp[0].items() if temp[0][k]>=threshold}
            #fisher's test check to avoid unifying clusters that are too infrequent
            for c in clusters.copy():
                counts = temp[1][c].split('\t')
                if float(counts[5]) > alpha:
                    del clusters[c]
            if not clusters:
                msg.env_render(message=f'\nNo complex segments identified in {os.path.split(ld)[1]}. That is the final version of your learning data.', **kwargs)
                msg.env_render(message="\n\nSimulation Finished", **kwargs)
                step = 0
            else:
                nc.write_insep(temp, os.path.join(wdir, 'inseparability.txt'))
                msg.tab_render(d=clusters, message=f'\nFound complex segments on iteration {step}:\n', **kwargs)
                newfeats = make_new_feats(feats, get_new_segs(feats, clusters, **kwargs)) 
                msg.env_render(message="\nChecking learner-generated feature file...\n", **kwargs)
                if not pnc.check_feats(newfeats, **kwargs):
                    msg.env_render(message = msg.messages['badfeatswarning'], **kwargs)
                write_new_ld_file(clusters, ld, os.path.join(wdir, 'LearningData.txt'), **kwargs)
                write_feats_w_check(os.path.join(wdir, 'Features.txt'), os.path.join(wdir, 'LearningData.txt'), oldfeats, newfeats, **kwargs)
                msg.env_render(message=f"\nExamining data from iteration {step}.\n", **kwargs)
                ld= os.path.join(wdir, 'LearningData.txt')
                feats = os.path.join(wdir, 'Features.txt') 
                step+=1
        shutil.move(ofpth, os.path.join(outdir, 'simulation', 'simulation_report.txt'))
        return None
    else: #failed feature check on first pass, cannot proceed
        if __name__=="__main__":
            raise SystemExit
        else:
            if os.path.isfile(os.path.join(outdir, 'simulation_report.txt')):
                with open(os.path.join(outdir, 'simulation_report.txt'), 'r', encoding='utf-8') as f:
                    errors = f.read().replace('\n', "<br>")
            return "Your feature file does not allow segments to be distinguished from each other. Perhaps try again with <a href='media/generic/Features.txt'>this generic feature file</a>?.<br>Here is how far the learner got:<br>"+errors
    

def plot_insep(simpath, threshold=1, takefirst=15, show=False, ftype='pdf'):
    '''
    searches contents of insepath for inseparability.txt files, and produces plots of them.
    saves each in the same location as the inseparability.txt file it represents.
    '''
    try:
        pins.plot_all_vert(simpath, threshold, takefirst, show, ftype)
        msg.env_render(message=f'\n {ftype.upper()} Plot generated for {simpath}')
    except FileNotFoundError:
        msg.env_render(message=f'\n Check that {simpath} exists and can be written to')

#############################################
##########command line options ##############
#############################################

if __name__=='__main__':
    import argparse
    parser= argparse.ArgumentParser(description=msg.messages['help'])
    parser.add_argument('--ld', help="full path to learning data file")
    parser.add_argument('--feats', help="full path to feature file")
    parser.add_argument('--outdir', help='full path to the location of output files. Warning: any folder called "simulation" in that location will be overwritten without a prompt.')
    parser.add_argument('--vowels', help='makes the learner count vocoids rather than consonants', type=bool, default=False)
    parser.add_argument('--language', help='if only this argument is specified, the learner will look for an appropriately named folder within "data" (located at the same level as "code") and will run the simulation on the learning data and features files inside that folder. For example, python compseg.py --lang=english/celex/broad runs the learner on the learning data and features inside ../data/english/celex/broad')
    parser.add_argument('--threshold', help='threshold value for inseparability', nargs='?', const=1.0, type=float, default=1.0)
    parser.add_argument('--alpha', help="alpha value for Fisher's Exact Test", nargs='?', const=0.05, type=float, default=0.05)
    args=parser.parse_args()
    kwargs = vars(args)
    if args.language:
        print(kwargs['language'])
        lgpath = os.path.join(os.path.dirname(os.getcwd()), 'data', args.language)
        kwargs['ld'] = os.path.join(lgpath, 'LearningData.txt')
        kwargs['feats'] = os.path.join(lgpath, 'Features.txt')
        kwargs['outdir']=lgpath
        simpath = os.path.join(lgpath, 'simulation')
        try: 
            complexify(**kwargs)
            plot_insep(simpath, ftype='png')
            plot_insep(simpath, ftype='pdf')
        except FileNotFoundError:
            msg.env_render(message=f'\nCould not locate the Learning Data or Features or output path at data/{language}')
            raise
    else:
        try:
            complexify(**kwargs)
            plot_insep(os.path.join(args.outdir, 'simulation'), ftype='png')
        except:
            msg.env_render(message=f"attempting to plot: {os.path.join(sys.argv[3], 'simulation')} but something went wrong. Are the simulation files at that location?")
            raise
