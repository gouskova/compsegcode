#!/usr/bin/env python3


'''
returns a simple count of the number of times a given segment occurs in a file.

file has to contain one word per line, with spaces between segments. like this:

p a t a
n a ts u
e k w a

etc.

Usage:

$python3 rawcounts.py /home/you/filename.txt 'a'


will print "4" for the input given above.

'''


def rawcount(path, seg):
    count = 0
    with open(path, 'r', encoding='utf-8') as f:
        for word in f:
            segs = word.rstrip('\n').split(' ')
            for x in segs:
                if x == seg:
                    count+=1
    print(count)
    return count


if __name__ == "__main__":
    import sys
    rawcount(sys.argv[1],sys.argv[2])
