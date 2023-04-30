#!/usr/bin/env python3
# coding: utf-8

#standard
import os
import sys
import itertools
from datetime import datetime

#special/external
import scipy
from nltk import ngrams

import pynatclasses as pnc
import messages as msg

'''
During the review process, the question arose why bigram probabilities were calculated over segments rather than over natural class sequences. This module explores the natural class bigram option.
'''

#todo:
    # make a dictionary of all possible bigrams of nat class subsets of -syll
    # from inseparability.txt, calculate unigram probs for each nat class in NC1 NC2 positions.
    # then calculate bigrams same way

def get_unigram_counts_per_class(**kwargs):
    '''
    requires a wordlist, a list of consonants, and a natural class dictionary (e.g., {"-syll,-cont": ['p', 't', 'k']...})
    returns a dictionary with segment counts for each natural class
    '''
    wdlist = kwargs['data']
    cons = kwargs['conslist'] #pnc.get_consonants(featfilepath)
    natclassdic = kwargs['natclassdic'] # pnc.wrap_classes(featfilepath)
    outdic={}
    for cl in natclassdic:
        if set(natclassdic[cl]).issubset(set(cons)):
            outdic[cl]=0        
    for wd in wdlist:
        wd = wd.strip('\n').split(' ')
        for cl in outdic:
            outdic[cl] += len([x for x in wd if x in natclassdic[cl]])
    return outdic


def make_natclass_bigrams(**kwargs):
    '''
    this is for consonants only!
    '''
    consdic = kwargs['consdic'] # get_unigram_counts_per_class()
    outlist = [' '.join(x) for x in itertools.product(consdic.keys(), repeat=2)]
    print('\n\n\nNumber of consonantal natural class bigrams\t' + str(len(outlist)) +'\n\n\n')
    return outlist

def count_natclass_bigrams(**kwargs):
    '''
    gets a nat class bigram list and a bigram count (from inseparability.txt)
    returns counts for each nat class bigram by looking up its seg sequences in the list
    '''
    natclassdic = kwargs['natclassdic']
    bigrams = {}.fromkeys(kwargs['bigrams'], 0) #make_natclass_bigrams() output
    with open(kwargs['inseppath'], 'r', encoding='utf-8') as f:
        lines = f.readlines()[1:]
    for line in lines:
        line = line.strip('\n').split('\t')
        seg1 = line[0].split(' ')[0]
        seg2 = line[0].split(' ')[1]
        for bigram in bigrams:
            cl1 = bigram.split(" ")[0]
            cl2 = bigram.split(" ")[1]
            if seg1 in natclassdic[cl1] and seg2 in natclassdic[cl2]:
                bigrams[bigram]+=int(line[2])
    return {k: v for k, v in bigrams.items() if v is not 0}
             
    
def forward_tp(clustercount, unigramcount):
    tpd = {}
    for ngram in [x for x in clustercount if not clustercount[x]==0]:
        cons = ngram.split(" ")[0]
        if sum(clustercount.values())==0 or sum(unigramcount.values())==0:
            tpd[ngram] = 0
        else:
            prob_ngram = clustercount[ngram]/float(sum(clustercount.values()))
            prob_cons = unigramcount[cons]/float(sum(unigramcount.values()))
            tpd[ngram] = prob_ngram/prob_cons
    return tpd

def backward_tp(clustercount, unigramcount):
    tpd = {}
    for ngram in [x for x in clustercount if not clustercount[x]==0]:
        cons = ngram.split(" ")[1]
        if sum(clustercount.values())==0 or sum(unigramcount.values())==0:
            tpd[ngram]=0
        else:
            prob_ngram = clustercount[ngram]/float(sum(clustercount.values()))
            prob_cons = unigramcount[cons]/float(sum(unigramcount.values()))
            tpd[ngram] = prob_ngram/prob_cons
    return tpd

# stats calculation

def findPValue(ngram,clustercount,allclustercounts):
    '''
    this does a fisher's exact tests with the goal of determining whether or not the observed number of a particular cluster is significantly different than 0.  it takes two arguments: the number of a particular kind of cluster (produced by count_clusters()) and the unigramcount (produced by uni_counts)).
    '''
    pValue = scipy.stats.fisher_exact([[clustercount,allclustercounts-clustercount],[0,allclustercounts]])
    return pValue[1]


def insep(clustercount, unigramcount):
    '''
    returns a dictionary of bidirectional inseparability measures: given a dictionary of ngrams xy, calculates prob of x being followed by y of all things (forward prob), and of y being preceded by x of all things (backward prob). needs a dictionary of clusters/ngram counts in a wordlist, and a dictionary of unigram counts. (produced by count_clusters() and uni_counts() respectively)
    '''
    forwards = forward_tp(clustercount, unigramcount)
    print('calculated forward tp')
    backwards = backward_tp(clustercount, unigramcount)
    print('calculated backward tp')
    tpd = {}
    countd = {}
    pvals = {} 
    for ngram in [x for x in clustercount if not clustercount[x]==0]:
        if ngram in forwards and ngram in backwards:
            tpd[ngram] = forwards[ngram]*backwards[ngram]
            pvals[ngram] = findPValue(ngram,clustercount[ngram],sum(clustercount.values()))
            outstring = [ngram, round(tpd[ngram], 4), clustercount[ngram], unigramcount[ngram.split(' ')[0]], unigramcount[ngram.split(' ')[1]], round(pvals[ngram], 3)] 
            countd[ngram] = '\t'.join(str(x) for x in outstring) 
        else:
            tpd[ngram]=0
    return (tpd, countd)


def bidir_prob_wrapper(**kwargs):
    '''
    a wrapper function for insep. featpath leads to Features.txt, and datapath leads to LearningData.txt
    '''
    clustercount = kwargs['bigramcounts']
    unigramcount = kwargs['unigramcounts']
    bidic = insep(clustercount, unigramcount)
    return bidic


def write_insep(**kwargs):
    d= kwargs['bidir']
    outpath = kwargs['outpath']
    natclassdic = kwargs['natclassdic']
    with open(outpath, 'w', encoding='utf-8') as f:
        f.write("\t".join(['segs', "ngram", "insep", "N(C1C2)", "N(C1)", "N(C2)", "p(C1C2)"])+'\n')
        for k in sorted(d[0], key=d[0].get, reverse=True):
            x = k.split(' ')
            s1 = '|'.join(natclassdic[x[0]])
            s2 = '|'.join(natclassdic[x[1]])
            f.write(f'[{s1}] [{s2}]\t{d[1][k]}\n')


if __name__ == '__main__':
    import sys, os
    stime = datetime.now()
    basepath = os.path.dirname(os.getcwd())
    #datapath = 'data/turkish/tell'
    #datapath = 'data/quechua/roots'
    #datapath = 'data/hebrew'
    #datapath = 'data/russian/zaliznjak'
    #datapath = 'data/english/celex/broad'
    #datapath = 'data/latin/long'
    #datapath = 'data/greek'
    #datapath = 'data/fijian'
    #datapath = 'data/navajo/stems'
    datapath = 'data/ngbaka'
    #datapath = 'data/english/celex/narrow'
    #datapath = 'data/french/childes_cds'
    #datapath = 'data/polish/polex/narrow'
    print(datapath)
    outpath = os.path.join(basepath, datapath, 'natclass_bigram_insep.txt')
    featfilepath = os.path.join(basepath, datapath, 'Features.txt')
    wordlist = os.path.join(basepath, datapath, 'LearningData.txt')
    inseppath = os.path.join(basepath, datapath, 'simulation/iteration1/inseparability.txt')
    conslist = pnc.get_consonants(featfilepath)
    natclassdic = pnc.wrap_classes(featfilepath)
    with open(wordlist, 'r', encoding='utf-8') as f:
        unigramcounts = get_unigram_counts_per_class(data=f, conslist=conslist, natclassdic=natclassdic) 
    bigrams = make_natclass_bigrams(consdic=unigramcounts)
    print('doing bigram counts now')
    bigramcounts = count_natclass_bigrams(natclassdic = natclassdic, bigrams = bigrams, inseppath = inseppath)
    x = bidir_prob_wrapper(bigramcounts=bigramcounts, unigramcounts=unigramcounts)
    write_insep(bidir=x, outpath=outpath, natclassdic=natclassdic)
    print(datetime.now()-stime)

