#!/usr/bin/env python3
# coding: utf-8

#standard
import os
import sys
import itertools

#special/external
import scipy
from nltk import ngrams

try:
    import pynatclasses as pnc
    import messages as msg
except ModuleNotFoundError:
    import compseg.code.pynatclasses as pnc
    import compseg.code.messages as msg
'''
the module extracts sequences of consonants from a corpus, and calculates the inseparability measures of certain clusters from probabilities of individual Cs and CCs in a phonological corpus. 
'''

def list_clusters(conslist, gramsize): 
    '''
    does not open anything
    compiles a list of all possible clusters of length gramsize,
    given the consonants in conslist. 
    '''
    exhaust_clust= [' '.join(x) for x in itertools.product(conslist, repeat=gramsize)]
    return exhaust_clust

def count_clusters(clustlist, datapath, gramsize):
    '''
    walks through each word in learning data and counts up its ngrams
    datapath is path to LearningData.txt, H&Wilson format
    clustlist is dictionary returned by list_clusters above
    gramsize is the length of ngrams (usu. betw 2 and 3)
    '''
    clustercount = {}.fromkeys(clustlist, 0)
    with open(datapath, 'r', encoding='utf-8') as f:
        for line in f:
            if '\t' in line:
                line = line.split('\t')[0]
            grams = line.strip().split(" ")
            word = [' '.join(x) for x in ngrams(grams, gramsize)]
            for cluster in clustercount:
                clustercount[cluster]+=len([x for x in word if x==cluster])
    return clustercount


def uni_counts(conslist, datapath):
    '''
    counts up frequencies of each consonant in data
    '''
    unigramcount={}.fromkeys(conslist, 0)
    with open(datapath, 'r', encoding='utf-8') as f:
        for word in f:
            word = word.strip().split(" ")
            for c in unigramcount:
                unigramcount[c] += len([x for x in word if x == c])
    return unigramcount

def first_counts(conslist, clustercount):
    '''
    counts up how often each consonant occurs as C1 in a bigram
    '''
    count={}.fromkeys(conslist, 0)
    for c in count:
        the_counts = [clustercount[ngram] for ngram in clustercount if ngram.split(" ")[0] == c]
        count[c] += sum(the_counts)
    return count

def last_counts(conslist, clustercount):
    '''
    counts up how often each consonant occurs as C2 in a bigram
    '''
    count={}.fromkeys(conslist, 0)
    for c in count:
        the_counts = [clustercount[ngram] for ngram in clustercount if ngram.split(" ")[1] == c]
        count[c] += sum(the_counts)
    return count


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
    backwards = backward_tp(clustercount, unigramcount)
    tpd = {}
    countd = {}
    pvals = {} 
    for ngram in [x for x in clustercount if not clustercount[x]==0]:
        if ngram in forwards and ngram in backwards:
            tpd[ngram] = forwards[ngram]*backwards[ngram]
            pvals[ngram] = findPValue(ngram,clustercount[ngram],sum(clustercount.values()))
            outstring = [ngram, round(tpd[ngram], 2), clustercount[ngram], unigramcount[ngram.split(' ')[0]], unigramcount[ngram.split(' ')[1]], round(pvals[ngram], 3)] 
            countd[ngram] = '\t'.join(str(x) for x in outstring) 
        else:
            tpd[ngram]=0
    return (tpd, countd)


def bidir_prob_wrapper(featpath, datapath, **kwargs):
    '''
    a wrapper function for insep. featpath leads to Features.txt, and datapath leads to LearningData.txt
    '''
    if 'vowels' in kwargs:
        msg.env_render(message='\nGetting vocoids...', **kwargs)
        conslist = pnc.get_vowels(featpath, **kwargs)
    else:
        msg.env_render(message="\nGetting consonants...", **kwargs)
        conslist = pnc.get_consonants(featpath, **kwargs)
    msg.env_render(message="\nGetting clusters...", **kwargs)
    clustlist = list_clusters(conslist, 2)
    msg.env_render(message="\nCounting clusters...", **kwargs)
    clustercount = count_clusters(clustlist, datapath, 2)
    msg.env_render(message="\nCalculating probabilities...", **kwargs)
    d= uni_counts(conslist, datapath)
    bidic = insep(clustercount, d)
    return bidic


def write_insep(d, outpath):
    with open(outpath, 'w', encoding='utf-8') as f:
        f.write("\t".join(["ngram", "insep", "N(C1C2)", "N(C1)", "N(C2)", "p(C1C2)"])+'\n')
        for k in sorted(d[0], key=d[0].get, reverse=True):
            f.write(d[1][k]+'\n')

