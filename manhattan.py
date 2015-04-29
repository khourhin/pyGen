#! /usr/bin/env python

# this is not properly working yet >>> make a check with other program
# Based on manhattan.py by Brent Pedersen
# https://github.com/brentp/bio-playground/blob/master/plots/manhattan-plot.py

#currently the size of the complete genome will not exact on the plot
# xs are incremented with the last position of the previous chromosomes, but
# this will be the last SNP position and not the last position on the complete
#chromosome

import basics_fasta as bf

#-------------------------------------------------------------------------------
def get_chro_len(g_fas):

    seq_d = bf.fasta_2_dict(g_fas)
    len_d = dict()

    for chro in seq_d:
        len_d[chro] = len(seq_d[chro])

    return sum(len_d.values())

#-------------------------------------------------------------------------------
def gen_data(out_weir_fst):

    with open(out_weir_fst, "r") as f:
        f.next()
        for line in f:
            chro, pos, fst = line.split()
            if fst != "-nan":
                yield chro, int(pos), float(fst)

#-------------------------------------------------------------------------------
def prep_for_plot(data, out_graph):
    # for fixing no display error on server
    import matplotlib as mpl
    mpl.use('Agg')

    from matplotlib import pyplot as plt
    from itertools import groupby, cycle

    
    colors = "br"
    colors = cycle(colors)
    tick_sp = cycle(["", "    ", "        "])
    last_pos = 0
    xs = []
    ys = []
    cs = []
    xs_by_chr = []
    
    for chro, dlist in groupby( data, lambda x: x[0] ):
        
        color = colors.next()
        dlist = list(dlist)
        xs_tmp = [ x[1] + last_pos for x in dlist ]
        xs.extend(xs_tmp)
        ys.extend([ x[2] for x in dlist ])
        cs.extend( [color] * len(dlist) )
        print "Chromosome: %s" % chro
        print "min: %e" % min(xs_tmp)
        print "max: %e" % max(xs_tmp)
        last_pos = xs_tmp[-1]

        # for ticking chr names
        sp = tick_sp.next()
        xs_by_chr.append( (sp + chro, (min(xs_tmp) + max(xs_tmp)) / 2 ))

    # Plotting
    plt.close()
    f = plt.figure()
    ax = f.add_axes((0.1, 0.09, 0.88, 0.85))
    ax.scatter(xs,ys,c=cs, s=2, lw=0)
    plt.xlim(0, xs[-1])
    plt.ylim(-0.2,1.1)
    plt.xticks([x[1] for x in xs_by_chr], [x[0] for x in xs_by_chr],
               rotation=-90, size=7)
    plt.savefig( out_graph + ".png", format="png" )
   
#-------------------------------------------------------------------------------
# MAIN
#-------------------------------------------------------------------------------

#USAGE : 
# $1: weir fst file
# $2: name for graph_out

if __name__ == "__main__":
    import sys
    #get_chro_len(sys.argv[2])
    data = gen_data(sys.argv[1])
    prep_for_plot(data, sys.argv[2])
