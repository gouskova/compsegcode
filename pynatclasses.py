#!/usr/bin/env python3
# coding: utf-8 

'''
a module that identifies all unique natural classes defined by the feature combinations in a table of phonological features.

dependencies are standard modules (itertools, sys, os)

The main functions:

FEAT LIST TO SEGS LOOKUP
	feats_to_segs_wrapper(featlist, featfilepath)	

NAT CLASS DICTIONARY
	wrap_classes(featfilepath)

GETTING ALL THE CONSONANTS
	get_consonants(featfilepath)

Look through the rest for more

'''

import os
import itertools
import sys


try:
    import messages as msg
except ModuleNotFoundError:
    import compseg.code.messages as msg



def read_feat_file(featfilepath, **kwargs):
	'''
	the input argument is a path to features.txt, tab-formatted according to Hayes and Wilson rules
	returns a vector of feature names, and a vector of segments plus their feature values in order
	featnames: [syll, cons, son, dor, ...]
	seglines: [k, -, -, -, +, ...] 
	'''
	try:
		with open(featfilepath, 'r', encoding='utf-8') as f:
			feats = f.readlines()
			featnames = feats[0].strip().split('\t')
			seglines = [x.strip().split('\t') for x in feats[1:]]
			return (featnames,seglines)
	except FileNotFoundError:
		kwargs['message']=f'could not open {featfilepath}'
		msg.env_render(**kwargs)
		raise SystemExit

def segs_to_feats(featlines, **kwargs):
	'''
	takes as input a tuple: (featnames, seglines) see read_feat_file
	returns a dictionary with segment name keys and +feat -feat lists as values
	{k: [-syll, -cons, -son, +dor, ...], p: [], etc}
	'''
	featnames = featlines[0]
	seglines = featlines[1]
	segdict = {}
	for line in seglines:
		segdict[line[0]]=[]
		for feat in featnames:
			featvalue = line[featlines[0].index(feat)+1]
			if not featvalue in ['+', '-', '0']:
				kwargs['message']='Your feature file is malformed. Feature values have to be "+", "-", or "0".'
				msg.env_render(**kwargs)
			if not featvalue=='0':
				segdict[line[0]].append(featvalue+feat)
	return segdict


def make_feat_vectors(featlines):
	'''
	takes as input a tuple: (featnames, seglines) see read_feat_file
	returns a dictionary of feature values along with segments that have those feature values.
	e.g.,
{-syll: [k, t, p, w, j, n,...]}
	'''
	featdict = {}
	for feat in featlines[0]:
		featdict["+"+feat]=[]
		featdict["-"+feat]=[]
		featindex = featlines[0].index(feat)
		for line in featlines[1]:
			seg = line[0]
			featvalue = line[featindex+1]
			if featvalue=="+":
				featdict['+'+feat].append(seg)
			elif featvalue=="-":
				featdict["-"+feat].append(seg)
	#drop feature values not associated with seg lists (e.g., -dorsal if dorsal is privative)
	finaldic = {}
	for x in featdict:
		if featdict[x]:
			finaldic[x] = featdict[x] 
	return finaldic 



def check_feats(featlines, **kwargs):
    '''
    returns seg and feat value if the feature specifications of one seg are a proper subset of the other. when this holds, the first seg cannot be uniquely identified using its features, so the user should be told.
    '''
    segdict = segs_to_feats(featlines)
    problemsegs = []
    for seg, otherseg in itertools.combinations(segdict.keys(), 2):
            if set(segdict[seg]).issubset(set(segdict[otherseg])):
                problemsegs.append((seg, otherseg))
    if not problemsegs:
         kwargs['message'] = "\nThe new feature file is well-formed. all the segments can be uniquely identified.\n"
         msg.env_render(**kwargs)
         return True
    else:
         kwargs['message'] = "\nThe new feature file does not allow certain segments to be distinguished from each other:\n\n"
         msg.env_render(**kwargs)
         for x in problemsegs:
             kwargs['message']=f'\n{x[0]}  has a subset of the features of {x[1]}'
             msg.env_render(**kwargs)
         return False 
   
    

def get_nat_classes(featlines, **kwargs):
	'''
	takes as input a tuple: (featnames, seglines) see read_feat_file
	returns a dictionary of natural classes as keys, and segments as values
	{-son,+cons: [p, t, n, l, ...]}
	'''
	segdict = segs_to_feats(featlines)
	featdict = make_feat_vectors(featlines)
	#find all pairs of segs that share feature values:
	natclasslist = []
	for seg in segdict:
		for otherseg in segdict:
			overlap = set(segdict[seg]) & set(segdict[otherseg])
			if len(overlap)>0 and not overlap in natclasslist:
				natclasslist.append(list(overlap))
	natclassdic = {}
	#compile lists of segments that natural classes expand to:
	for cl in natclasslist:
		clname = ','.join(sorted(cl))
		natclassdic[clname]=[]
		#for every seg, check if its feature values are in that natural class description
		for seg in segdict:
			if set(cl).issubset(set(segdict[seg])): 
				natclassdic[clname].append(seg)
	kwargs['message']=f'\nNumber of natural classes: {len(natclassdic)}'
	msg.env_render(**kwargs)
	return natclassdic

def feat_to_seg_lookup(feats, featdict):
    '''
    given a list of features (feats) and a dictionary of them (featdict), returns all the segs that have that feature combo
    '''
    segs = []
    for feat in feats:
        segs.append(set(featdict[feat]))
    return sorted(list(set.intersection(*segs)))



def feats_to_segs_wrapper(feats, featfilepath):
	'''
	takes as args a list of features (in any order), and a path to the feature file
        returns a list of segments
	example:
	pynatclasses.feats_to_segs_wrapper(['-syll', '-son', '+cont'], /home/path/to/feats.txt')
	this is just a wrapper for feat_to_seg_lookup, and it adds the extra step of opening the feature file (so this only happens once per feat line)
	'''
	featdict = make_feat_vectors(read_feat_file(featfilepath))
	return feat_to_seg_lookup(feats, featdict)

def powerset(thing):
	'''
	uses itertools recipe for powerset. see python docs
	'''
	x = list(thing)
	return itertools.chain.from_iterable(itertools.combinations(x, r) for r in range(len(x)+1))

def which_bigger(featdict, feat1, feat2):
	'''
	given a dictionary mapping features to segment lists, and a couple of feature names, finds out which feature covers a bigger natural class. if they are tied, it will return a boolean False
	'''
	if len(featdict[feat1])>len(featdict[feat2]):
		return feat1
	elif len(featdict[feat2])>len(featdict[feat1]):
		return feat2
	else:
		return False

def avg_cl_size(featdict, featlist):
	'''
	calculates average number of segs that a class refers to
	'''
	lenths = [len(featdict[feat]) for feat in featlist]
	try:
		return sum(lenths)/len(featlist)
	except ZeroDivisionError:	
		return 1 


def find_shortest_descriptions(natclassdic, featlines):
	'''
	takes features from a verbose nat class dictionary, looks to see if there is a shorter subset of features that could describe same class, returns a less verbose dictionary of natural classes 
	'''
	finaldic = {}
	featdict = make_feat_vectors(featlines)
	for cl in sorted(natclassdic):
		segs = sorted(natclassdic[cl])
		feats = sorted(cl.split(','))
		if len(feats)==1:
			finaldic[cl]=segs
		else:
			#takes every sub-combination of features in the verbose dic
			#if any of them picks out an equiv. set of segs, takes the shortest 
			equiv_classes = []
			for featuple in powerset(feats):
				if (not len(featuple)==0) and (feat_to_seg_lookup(list(featuple), featdict) == segs):
					equiv_classes.append(sorted(list(featuple)))
			if len(equiv_classes)==1:
				finaldic[cl] = segs
			else:
				lenths = [len(x) for x in equiv_classes]
				shortest=[x for x in equiv_classes if len(x)==min(lenths)]
				if len(shortest)==1:
					finaldic[','.join(shortest[0])]=segs
				elif len(shortest)>1:
					cl_sizes = [avg_cl_size(featdict, cl) for cl in equiv_classes]
					#reward classes for using "bigger" features
					generalest = [x for x in equiv_classes if avg_cl_size(featdict, x) == min(cl_sizes)]
					#and then give up, just give the first if there are ties
					finaldic[','.join(generalest[0])]=segs
	return finaldic

def wrap_classes(featfilepath):
	'''
	takes in a full path to a feature file, Features.txt
        returns a dictionary of natural classes, and the segs they contain:
        natclassdict = {'-son,-cont': ['p', 't','k']...}
	'''
	featlines = read_feat_file(featfilepath)
	return find_shortest_descriptions(get_nat_classes(featlines), featlines)


def get_consonants(featfilepath, **kwargs):
	'''
	given a full path to a features.txt file, returns a list of all the symbols that are specified as -syll or -syllabic. Those feature names are special.
        kwargs are passed on to messages module for error handling
	'''
	featlines = read_feat_file(featfilepath)
	feats = make_feat_vectors(featlines)
	if '-syll' in feats:
		return feats['-syll']
	elif '-syllabic' in feats:
		return feats['-syllabic']
	else:
		kwargs['message']=f'\nThe feature file {featfilepath.split("simulation")[1]} does not have a column for -syll or -syllabic. The learner needs this feature to separate consonants from vowels. Fix this and try again.'
		msg.env_render(**kwargs)

def get_vowels(featfilepath, **kwargs):
    '''
    first argument is a path to a features.txt file. returns a list of vowel and glide symbols. kwargs are passed to essages module for `error handling
    '''
    featlines = read_feat_file(featfilepath)
    feats = make_feat_vectors(featlines)
    if '+syll' in feats:
        return feats['+syll']
    elif '+syllabic' in feats:
        return feats['+syllabic']
    else:
        kwargs['message']=f'\nThe feature file {featfilepath.split("simulation")[1]} does not have a column for +syll or +syllabic. The learner needs this feature to separate vowels from true consonants. Fix this and try again.'
        msg.env_render(**kwargs)



def get_vocoids(featfilepath, **kwargs):
    '''
    first argument is a path to a features.txt file. returns a list of vowel and glide symbols. kwargs are passed to essages module for `error handling
    '''
    featlines = read_feat_file(featfilepath)
    feats = make_feat_vectors(featlines)
    if '-cons' in feats:
        return feats['-cons']
    elif '-consonantal' in feats:
        return feats['-consonantal']
    else:
        kwargs['message']=f'\nThe feature file {featfilepath.split("simulation")[1]} does not have a column for -cons or -consonantal. The learner needs this feature to separate vocoids from true consonants. Fix this and try again.'
        msg.env_render(**kwargs)


def outwrite_classes(featfilepath, outpath):
	'''
	gets natural classes from the file at featfilepath (needs to be a full path, /home/you/etc/Features.txt), and saves natural classes to a file at outpath.
	'''
	import os
	classdic = wrap_classes(featfilepath)
	if not os.path.isfile(outpath):
		with open(outpath, 'w', encoding='utf-8') as f:
			for cl in sorted(classdic):
				f.write(cl + '\t' + ','.join(sorted(classdic[cl]))+'\n')
	else:
		overwrite = input("File " + outpath + " already exists. Overwrite? [y/n] ")
		if overwrite=='y':
			with open(outpath, 'w', encoding='utf-8') as f:
				for cl in sorted(classdic):
					f.write(cl + '\t' + ','.join(sorted(classdic[cl]))+'\n')
		else:
			raise SystemExit

def seglist_to_feats(seglist, segdict):
	'''
	returns features *not* shared by a list of features.
	'''
	contrastfeats = set()
	for seg in seglist:
		contrastfeats = contrastfeats.union({x for x in segdict[seg] if not x.startswith('0')})
	return sorted(list({x.lstrip('+-') for x in contrastfeats}))


def make_custom_proj(feats, featfilepath, outpath, **kwargs):
	'''
	feats: some feature(s) defining a natural class. e.g., '+son' or '-son,-cont'. If more than 1, must be a comma-separated string. no spaces
	featfilepath: path to feature file from which to read natural classes.
	outpath: where to put projections.txt.
	this function writes a projections file in the format used by Wilson's MaxEnt learner:
	projname    feats_defining_class    feats_visible_on_proj   ngrams
	default proj always included.
	this is a stand-alone function, it opens the feature file rather than be fed pre-read lines
	'''
	featl = [x.strip() for x in feats.split(',')]
	kwargs['message']=f'\n{featl}'
	msg.env_render(**kwargs)
	thesegs = feats_to_segs_wrapper(feats, featfilepath)
	feats_to_project = seglist_to_feats(thesegs, read_feat_file(featfilepath)).append('wb')
	with open(outpath, 'w', encoding='utf-8') as f:
		f.write('\t'.join(['default', 'any', 'all', '3']))
		f.write('\t'.join([feats, ''.join(featl), ','.join(feats_to_project), '2', '3']))
     


if __name__ == '__main__':
	import sys
	HelpString = '\n\nThis utility finds natural classes in a feature file formatted according to Hayes and Wilson (2009, Linguistic Inquiry) conventions. Basic usage: \n\n$ python3 pynatclasses.py /home/full/path/to/file/Features.txt /home/full/path/to/output.txt\n\n You can also get all the consonants from a feature file from a command line call: \n $ python3 pynatclasses.py /home/full/path/to/Features.txt cons\n\n This last option requires there being a -syll or -syllabic feature in the file.\n\n\n To see other options, import it into python and try help(pynatclasses)'
	CLError = '\n\nPlease provide the name of a feature file and a place to save the natural classes to. \n\nFor example: "python3 pynatclasses.py /home/you/Desktop/features.txt /home/you/Desktop/natclasses.txt"\n\n'
	basepath = os.path.dirname(os.path.dirname(os.getcwd()))
	if "help" in sys.argv:
		msg.env_render(message=HelpString)
	elif 'check' in sys.argv:
		feats = os.path.join(basepath, 'data', sys.argv[1], 'Features.txt')
		featlines = read_feat_file(feats)
		check_feats(featlines)
	elif 'cus' in sys.argv:
		feats = os.path.join(basepath, 'data', sys.argv[1], 'Features.txt')
		make_custom_proj(sys.argv[2], feats, '/home/maria/Desktop/projections.txt')
	elif not "cons" in sys.argv:
		try: 
			feats = sys.argv[1]
			outfile = sys.argv[2]
			outwrite_classes(feats, outfile)
		except IndexError:
			msg.env_render(message=CLError)
	else:
		try:
			get_consonants(sys.argv[1])
		except IndexError:
			msg.env_render(message="Please provide a full path to the feature file, like this: \n $ python3 pynatclasses.py /home/full/path/to/features.txt cons")

