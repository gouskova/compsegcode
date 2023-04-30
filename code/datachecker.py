#!/usr/bin/env python3
#-*- coding: utf-8 -*-

'''
this is a utility function for checking data for the UCLAPL learner (command line version)
it goes through a file called LearningData.txt, and checks it against your Features.txt file 
if any segments turn up that aren't in both files, the function alerts you.

usage:
    datachecker.findOrphans('/home/you/WordListFile.txt', '/home/you/Featfile.txt')
'''
import pynatclasses
import messages as msg

def collectLDSegs(somepath):
	segs = []
	try:
		with open(somepath, 'r', encoding='utf-8') as ldatafile:
			for line in ldatafile:
				wordsegs = line.strip('\n').split()
				for seg in wordsegs:
					if not seg in segs:
						segs.append(seg.strip())
	except FileNotFoundError:
		msg.env_render(message=f"\n\n\nNo file at {somepath} \n\n\n\n")
	return segs

def collectFSegs(somepath):
	return pynatclasses.segs_to_feats(pynatclasses.read_feat_file(somepath)).keys()

def findOrphans(learningdata, featurefile, verbose=False):
	learningdatasegs = collectLDSegs(learningdata)
	featurefilesegs = collectFSegs(featurefile)
	for x in ['+', '-', '*', '|']:
		if x in featurefilesegs:
			msg.env_render(message=f"Please do not use {x} in your segment list. Choose a symbol that is not used in regular expressions. Your feature file is incompatible with the gain-based version of the MaxEnt Phonotactic Learner.")
	if verbose:
		msg.env_render(message=f"\nThe segments in your data file are: \n {','.join(sorted(learningdatasegs))}")
		msg.env_render(message=f"\nThe segments in your Features.txt file are: \n {','.join(sorted(featurefilesegs))}")
	orphans = list(set(learningdatasegs) - set(featurefilesegs))
	if len(orphans)==0:
		msg.env_render(message=f'\nAll the segments are defined in the feature file.')
	else:
		if verbose:
			msg.env_render(message=f"\nyour orphan segments are: \n {','.join(orphans)}")
			return orphans
		else:
			#pass
			return orphans



if __name__ == '__main__':
	import sys
	import os
	helpmessage = 'please provide the full locations of the learning data file and the feature file as follows: \n$ python3 datachecker.py /home/me/Desktop/LearningData.txt /home/me/Desktop/Features.txt. Alternatively:\n\n $ python3 datachecker.py russian/wds_t_s_t_sh dirc\n\n this will check the LearningData.txt against the Features.txt file.'
	if 'dirc' in sys.argv:
		basepath = os.path.dirname(os.getcwd())
		lgname = sys.argv[1]
		learningdata = os.path.join(basepath, 'data', lgname, 'LearningData.txt')
		featurefile = os.path.join(basepath, 'data', lgname, 'Features.txt')
		findOrphans(learningdata, featurefile, verbose=True)
	elif len(sys.argv)>1:
		learningdata=sys.argv[1]
		featurefile = sys.argv[2]
		msg.env_render(message=f'learning data: {learningdata}')
		msg.env_render(message=f'"feature file: {featurefile}')
		findOrphans(learningdata, featurefile, verbose=True)
	else:
		try:
			findOrphans(learningdata, featurefile, verbose=True)
		except FileNotFoundError:
			msg.env_render(message=helpmessage)
