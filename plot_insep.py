'''
produces plots for inseparability values, from inseparability.txt files
requires matplotlib and dependencies (tkinter)

this module is used in compseg via a wrapper function, plot_insep.

i removed some of the earlier legacy options (plotting individual inseparability.txt files) to simplify maintenance
'''


import os


import matplotlib.pyplot as plt
import numpy as np

try:
    import messages as msg
except ModuleNotFoundError:
    import compseg.code.messages as msg

def get_paths(simpath):
    '''
    returns a dictionary with iterations as keys, and values are full paths to all the inseparability.txt files inside simpath
    input arg is a full path to the simulation directory
    '''
    outdic = {}
    for tup in os.walk(simpath): 
        if 'inseparability.txt' in tup[2] and tup[1]==[]:
            iteration = int(os.path.split(tup[0])[1].replace('iteration', ''))
            insepath = os.path.join(simpath, 'iteration'+str(iteration), 'inseparability.txt')
            outdic[str(iteration)] = insepath 
    return outdic


def get_plot_title(simpath):
    '''
    checks if the path is a full path, and finds the Language + version substring of it
    '''
    if os.path.isabs(simpath):
        #return ' '.join([x.capitalize() for x in os.path.split(simpath)[0].split('data')[1].split(os.path.sep)])
        return 'Simulation plot'
    else:
        return ' '.join([x.capitalize() for x in simpath.split(os.path.sep)])

def plot_all_vert(simpath, threshold, takefirst, show, ftype):
    '''
    given a path to a data/language folder, plots all the simulations in a single figure, sharing the x-axes
    the clusters in this version of the plotting function appear on the y-axis.
    '''
    iterations = get_paths(simpath)
    q = len(iterations)
    if q<=3:
        dims=(1, q)
        if q==1:
            size = (4, 8)
        else:
            size = (7, 4)
        stack = False
    else:
        if q % 2 == 0:
            dims = (2, int(q/2))
        else:
            dims = (2, int(q/2)+1)
        size = (4*q, 3*q)
        stack = True
    fig, axs = plt.subplots(dims[0], dims[1], figsize = size, sharex= True)
    if len(iterations)>3:
        axs=axs.flatten()
    libfont = {'fontname': 'Linux Libertine', 'size': 'x-large'}
    if 'mbay' in simpath:
        libfont = {'fontname': 'Linux Libertine', 'size': 'large'}
    maxval = 0
    for iteration in iterations:
        with open(iterations[iteration], 'r', encoding='utf-8') as f:
            clusters = []
            ins_vals = [] #inseparability values
            p_vals = [] #p > 0.05
            lines = f.readlines()
            for line in lines[1:]:
                vals = line.strip().split('\t')
                clusters.append(vals[0])
                ins_vals.append(float(vals[1]))
                if float(vals[5])>=0.05:
                    p_vals.append('kX')
                else:
                    p_vals.append('bo')
            if takefirst:
                clusters = clusters[:takefirst]
                ins_vals = ins_vals[:takefirst]
                p_vals = p_vals[:takefirst]
            clusters = clusters[::-1] #reverse order
            ins_vals=ins_vals[::-1]
            p_vals = p_vals[::-1]
            theplot = axs
            if len(iterations)>1:
                try:
                    theplot = axs[int(iteration)-1]
                except IndexError:
                    theplot = axs[int(iteration)-2]
            theplot.set_xlabel("Inseparability", **libfont)
            if iteration==1:
                theplot.set_ylabel("Clusters", **libfont)
            if max(ins_vals)>maxval:
                maxval = max(ins_vals)
            for i in range(0,len(ins_vals)):
                theplot.plot(ins_vals[i], clusters[i], p_vals[i])
            theplot.set_title('Iteration %s' % iteration, **libfont)
            theplot.axvline(threshold, color='r')
            plt.setp(theplot.yaxis.get_majorticklabels(), **libfont)
    ptit= get_plot_title(simpath) 
    if stack:
        plt.tight_layout()
        plt.suptitle("%s\n\n" % (ptit), **libfont, va='baseline')
    else:
        plt.suptitle("%s\n\n" % (ptit), **libfont)
    if maxval<=1:
        plt.xlim(0,1.1)
    else:
        plt.xlim(0)
    plt.subplots_adjust(left=0.20, bottom=0.20, wspace=0.4)
    if show:
        plt.show()
    fig.savefig(os.path.join(simpath, '.'.join(['insep_plots', ftype])))
    #because matplotlib does not take out the trash
    fig.clf()
    plt.close(fig)



if __name__=='__main__':
    import sys
    if 'help' in sys.argv:
        msg.env_render(message="Looks inside the language/simulation directory and plots inseparability values for the top 15 clusters in each iteration. The individual plots will be placed at the same level as inseparability.txt files that inspired them.\n\nUsage:\n\n$ python3 plot_insep.py languagename\n\n. The 'languagename' argument is a full path to the location of the 'simulation' folder that contains the inseparability.txt file. You can also plot simulations for languages in the 'data' folder.") 
    else:
        try:
            plot_all_vert(sys.argv[1], threshold=1, takefirst=15, show=False, ftype='pdf')
            plot_all_vert(sys.argv[1], threshold=1, takefirst=15, show=False, ftype='png')
        except FileNotFoundError:
            simpath = os.path.join(os.path.dirname(os.getcwd()), 'data', sys.argv[1], 'simulation')
            plot_all_vert(simpath, threshold=1, takefirst=15, show=False, ftype = 'pdf')
            plot_all_vert(simpath, threshold=1, takefirst=15, show=False, ftype = 'png')
        msg.env_render(message=f'saved figures for {sys.argv[1]}"')
    

