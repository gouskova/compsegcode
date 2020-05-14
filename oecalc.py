#!/usr/bin/env python3


'''
Supplies a function for calculating Observed/Expected values over pairwise combinations of segments.

Input format for the data file:

p a t a
p i k u b e
s a mb u k i
k a tʃ o
...

Input format for the segment list:

"seg1 seg2 seg3"

(e.g., "a e i o u" or "ph t tʃ")

The OE() function will print a table of pairwise OE calculations in Terminal or save it to a file in the location you specify; see doc entry for "OE" for details.

Example use:


from the command line (e.g., bash):

$ oecalc yourdatafile.txt 'p t k b d g' local


you can import the module into an interactive python3 shell:


start python3 however you always do it, import oecalc, and call help(oecalc)

'''

# LICENSE: Released under the FreeBSD license. 
# http://www.freebsd.org/copyright/freebsd-license.html




from itertools import product
import re
from nltk import ngrams

# make a dictionary of pairs of segs and collect their O/E values into it


def OEcalc(**pars):
    '''
    arguments:
    *filepath* is a path to the file you want to calculate O/E over.
    formatting: same as the input to Hayes and Wilson's UCLA Phonotactic Learner; segments separated by spaces.

    p a t a
    p i k u b e
    s a mb u k i
    p t a k u

    *segs* is the list of symbols you want to evaluate for Observed/Expected. separate it by spaces and surround by quotes, "e i o"
    O/E is calculated as follows:
    Expected: N(S1) * N(S2)/ N of all pairs
    Observed: N(S1S2)
    the function will return unrounded OE, as well as a value rounded to the parameter given by the "rounded" argument.
    Defaults to 2, so an O/E value of 1.3432 will be printed as 1.34.
    *local* is boolean and determines whether segment pairs are adjacent (e.g., "p a" in "p a t i") or nonlocal (as in "p t" in "p a t i") 
    '''
    filepath = pars['filepath']
    segs = pars['segs']
    local = pars['local']
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            words =f.readlines()
    except FileNotFoundError:
        print("please make sure there is a file to read at the location.")
        pass
    seglist = segs.strip().split(' ')
    wordlist = [x.strip().split(' ') for x in words]
    pairs = {}.fromkeys(' '.join(list(x)) for x in product(seglist, repeat=2)) #creates a dictionary with S1,S2 pairs from seglist, every possible combination
    segs = {} #reusing variable because that's a good thing
    for seg in seglist:
        segs[seg+'1']=0 #initializing positional frequency counts at zero; segs and pairs are counted separately in diff dictionaries
        segs[seg+'2']=0
    for x in pairs:
        pairs[x] = {'observed':0, 'expected':0}
    paircount = 0
    for word in wordlist:
        if not local:
            wrd = [x for x in word if x in seglist]
        else:
            wrd = word
        if len(wrd)<2:
                continue
        else:
                segpairsinword = [wrd[x]+' ' + wrd[x+1] for x in range(0, len(wrd)-1)] #a list of 2-seg pairs in the word
                paircount += len(segpairsinword)
                for pair in segpairsinword:
                    if pair in pairs:
                        pairs[pair]['observed']+=1
                        segs[pair.split(' ')[0]+'1'] += 1 #count actual observed freq of seg 1 in pair
                        segs[pair.split(' ')[1]+'2'] += 1 #count actual observed freq of seg 2 in pair
    for pair in pairs:
            pr = pair.split(" ")
            seg1,seg2 = pr[0],pr[1]
            pairs[pair]['expected'] = segs[seg1+'1']*segs[seg2+'2']/paircount
            try:
                pairs[pair]['OE'] = pairs[pair]['observed']/pairs[pair]['expected']
                pairs[pair]['OErnd']=round(pairs[pair]['OE'],3)
            except ZeroDivisionError:
                pairs[pair]['OE'] = 'NA'
                pairs[pair]['OErnd']='NA'
    return(pairs)


# make the dictionary printable   

def makeOETable(**pars):
    '''
    arranges the O/E values and sorts them into a table for display.
    the input should be the dictionary that's returned by pairOEcalc()
    output is a list of tab-separated lines that, if printed, looks like a table
    '''
    filepath = pars['filepath']
    segs = pars['segs']
    local = pars['local']
    if 'rounded' in pars:
        rounded = pars['rounded']
    else:
        rounded=True
    pairsdic = OEcalc(**pars)
    seglist = segs.strip().split(' ')
    header = '\t'+'\t'.join(seglist)
    rows = [header]
    for seg in seglist:
        row = [seg]
        for otherseg in seglist:
            pair = seg+" "+otherseg
            if rounded:
                row.append(str(pairsdic[pair]['OErnd']))
            else:
                row.append(str(pairsdic[pair]['OE']))
        outrow = '\t'.join(row)
        rows.append(outrow)
    return(rows)

# printable observed and expected values (as opposed to ratio) for each pair

def makeCountTable(**pars):
        '''
        print how often each segment in a pair was observed and how often it was expected

        '''
        filepath = pars['filepath']
        segs = pars['segs']
        local = pars['local']
        pairsdic=OEcalc(**pars)
        seglist = segs.strip().split(' ')
        header = '\t'+'\t'.join(seglist)
        rows = [header]
        for seg in seglist:
                row1 = [seg+' observed']
                row2 = [seg+' expected']
                for otherseg in seglist:
                        pair = seg+' '+otherseg
                        row1.append(str(pairsdic[pair]['observed']))
                        row2.append(str(round(pairsdic[pair]['expected'], 2)))
                outrow1 = '\t'.join(row1)
                outrow2 = '\t'.join(row2)
                rows.append(outrow1)
                rows.append(outrow2)
        return(rows)


#saves OE table to file
def writeOE(filepath, segs, outfilepath, local, ounded=True):
        '''
        filepath is where your wordlist is located.
        segs is the list of segments, space-separated, over which you want O/E calculated
        outfilepath is for recording the results
        if rounded is True, the output will be rounded to 3 decimals (otherwise its Python3s default for integers, which is unreadably long)
        '''
        pairsdic = OEcalc(filepath, segs, local)
        out = makeOETable(filepath, segs,local, rounded)
        f = open(outfilepath, 'w', encoding='utf-8')
        for row in out:
            f.write(row+'\n')
        f.close()

def trigramOEcalc(filepath, trigram, projection, verbose=False):
    '''
    arguments:
    *filepath* is a path to the file you want to calculate O/E over.
    formatting: same as the input to Hayes and Wilson's UCLA Phonotactic Learner; segments separated by spaces.

    p a t a
    p i k u b e
    s a mb u k i
    p t a k u

    The trigram argument is a specific non-local sequence that you want to evaluate. It must be parameterized against a specific projection. That is,
    say you want to see how often the sequence "a e i" occurs in a word. You presumably want to look at the nonlocal trigram of these vowels, but not include in your counts examples where other vowels intervene--that is, you want to count "p a t e k i" but not "p a t e u k i". In order to make this happen, add "a e i o u" as the projection that you are counting on.

    O/E of the trigram is calculated as follows:
    Expected: N(S1) * N(S2) * N(s3) / N of all trigrams
    Observed: N(S1S2S3)
    the function will return unrounded OE, as well as a value rounded to the parameter given by the "rounded" argument.
    '''
    #some prep for counting
    projection = projection.strip().split(' ')
    trigram = trigram.strip().split(' ')
    observed_counts=0
    expec_dic={'seg1':0, 'seg2':0, 'seg3':0}
    n_of_all_trigrams = 0
    with open(filepath, 'r', encoding='utf-8') as f:
        for word in f:
            projsegs = [x for x in word.strip().split(' ') if x in projection]
            if len(projsegs)<3:
                continue
            if len(projsegs)==3 and projsegs!=trigram:
                n_of_all_trigrams+=1
                if projsegs[0] == trigram[0]:
                    expec_dic['seg1']+=1
                if projsegs[1] == trigram[1]:
                    expec_dic['seg2']+=1
                if projsegs[2] == trigram[2]:
                    expec_dic['seg3']+=1
            else:
                for subtrig in ngrams(projsegs, 3):
                    n_of_all_trigrams+=1
                    if list(subtrig) == trigram:
                        observed_counts+=1
                    if subtrig[0]==trigram[0]:
                        expec_dic['seg1']+=1
                    if subtrig[1]==trigram[1]:
                        expec_dic['seg2']+=1
                    if subtrig[2]==trigram[2]:
                        expec_dic['seg3']+=1
    if verbose:
        print('observed counts for trigram %s : %s' % (' '.join(trigram), observed_counts))
        print('number of all trigrams in wordlist: ' + str(n_of_all_trigrams))
        print('positional frequencies for each segment:\n %s : %s \n %s : %s \n %s : %s' % (trigram[0], expec_dic['seg1'], trigram[1], expec_dic['seg2'], trigram[2], expec_dic['seg3']))
    expected = expec_dic['seg1']*expec_dic['seg2']*expec_dic['seg3']/n_of_all_trigrams
    try:
        return round(observed_counts/expected,6)
    except ZeroDivisionError:
        print("The O/E is not undefined (division by zero)")
        
           

#if you want to run it from command line--just to print table of OE values to screen
if __name__ == "__main__":
    import sys
    try:
        filepath = [x for x in sys.argv if x.endswith('.txt')][0]
        segs = [x for x in sys.argv if " " in x][0]
        if 'local' in sys.argv:
            local = True
        elif 'nonlocal' in sys.argv:
            local = False
        else:
            print('Please specify whether you want to look at local or nonlocal co-occurrence of your segments. Enter either "local" or "nonlocal."') 
        if 'raw' in sys.argv:
            table = makeCountTable(filepath, segs, local)
        else: 
            table = makeOETable(filepath, segs, local)
            if len(table)==1:
                print('something went wrong')
        for row in table:
            print(row)
    except:
        print("You entered: ")
        print(sys.argv)
        print("you need to supply three arguments to oecalc: the file where your words are (ending in .txt), the segments you want to calculate O/E statistics for, encased in quotation marks (e.g., 'a e i o u'), and a local/nonlocal choice.\nfor example, enter the following at the bash prompt:\n$ oecalc /home/yourname/directory/yourfile.txt 'a e i o u' local")
