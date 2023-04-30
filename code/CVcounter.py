#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import regex as re

'''
a module for counting various CV patterns in a .txt file. The file should have one word per line, with segments separated by spaces. Needs a list of vowels, "vs", to function.
'''

def countCVC(**kwargs):
    '''
    path is a path to a LearningData.txt file (or similar)
    '''
    ld = kwargs['ld']
    searchseqs = kwargs['searchseqs']
    vs = kwargs['vs']
    if ',' in searchseqs:
        seqs = [x.strip() for x in searchseqs.split(",")]
    elif "\n" in searchseqs:
        seqs = [x.strip() for x in searchseqs.split("\n")]
    if type(vs) == str:
        vowels = vs.split(' ')
    elif type(vs)==list:
        vowels = vs
#    print('the vowels are: ' + ' '.join(vowels))
    CVdic = {}.fromkeys(seqs, 0)
    with open(ld, 'r', encoding='utf-8') as f:
        for line in f:
            word = line.rstrip("\n").split(" ")
            CVword = []
            for i in range(len(word)):
                if word[i] in vowels:
                    CVword.append("V")
                else:
                    CVword.append("C")
            CVword = ' '.join(CVword)
            for seq in CVdic:
                #to count overlapped instances of CVC in every word
                #special module (not overwrought regex)
                matches = re.findall(seq, CVword, overlapped=True)
                n = len(matches)
                CVdic[seq]+=n
    return CVdic
