#-------------------------------------------------------------------------------
def plot_hist(len_d, outfile):

       # for fixing no display error on server                           
       import matplotlib as mpl
       mpl.use('Agg')

       import matplotlib.pyplot as plt

       plt.hist(len_d.values(),bins=50)

       plt.suptitle("Contigs lengths distribution") # the title          
       plt.xlabel("length (in bp)")
       plt.ylabel("Number of contigs")
       plt.savefig( outfile, format="png" )

