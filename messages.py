'''
messages to be printed while the compseg learner is running.
'''

messages = {
        'badfeatswarning': '\n\nWarning: some of the segments can no longer be defined using combinations of existing features (e.g.: "k" [+dorsal], "kp" [+dorsal, +labial]. No constraint can refer to just "k" without also referring to [kp]). You will have to fix this by hand later by editing your final feature file yourself.',
        'help': "\n\nThe learner identifies consonant sequences that have the distribution of complex segments. If you specify the 'language' argument, as in \n\n $ python3 compseg.py --language=english/celex/broad\n\n The learner will look for 'english/celex/broad' inside the 'data' subfolder, at the same level as 'code'. The 'simulation' output folder will then be at the same level as learning data and features files. n\n The second way to run the learner is to supply full paths to the learning data, features, and an output directory:\n\n $ python3 compseg.py --ld=/path/to/learningdata.txt --feats=/path/to/features.txt --outdir=/path/to/output\n\n \n\n Threshold and alpha arguments, both optional, modify the quantitative parameters of the learner.",
        }



def env_render(**kwargs):
    message = kwargs.get('message')
    if __name__=='__main__':
        print(message)
    else:
        print(message)
        outfilepath = kwargs.get('outfilepath')
        if outfilepath:
            with open(outfilepath, 'a', encoding='utf-8') as f:
                f.write(message)


def tab_render(**kwargs):
    dic = kwargs.get('d')
    sort = kwargs.get('sort', True)
    env_render(message=kwargs.get('message'), outfilepath=kwargs.get('outfilepath'))
    for key in sorted(dic, key=dic.get, reverse=True):
        env_render(message=f'{key}\t{dic[key]}\n', outfilepath=kwargs.get('outfilepath'))

