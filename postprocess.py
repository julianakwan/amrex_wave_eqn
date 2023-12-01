import yt
from yt.frontends.boxlib.data_structures import AMReXDataset
from matplotlib import rc_context
from matplotlib.animation import FuncAnimation
import sys
import glob
import numpy as np


def plot_single_slice(file_dir, stepno):
    #Directory to AMReX outputs

    suffix = ""
    
    try:
        input_filename = glob.glob(file_dir+"plt{:0>5d}.old.*".format(stepno))
        f = open(input_filename+"/Header", 'r')
        suffix = "plt{:0>5d}.old.*"
        print("Here at try\n");
    except:
        f = open(file_dir+"plt{:0>5d}/Header".format(stepno), 'r')
        suffix = "plt{:0>5d}*"
        print("Here at except\n");
    finally:
        input_filename = glob.glob(file_dir+suffix.format(stepno))
        print("Here at finally\n");
        print(input_filename)
        f.close()


    input_filename = input_filename[0]


    #Output directory for plots
    output_filename = file_dir+"plt{:0>5d}".format(stepno)
    
    ds = AMReXDataset(input_filename)

    #can query values using 
    #field_name = 'u'
    #coords = [0,0,0]
    #ds.find_field_values_at_point(field_name, coords) 

    axes = ["x", "y", "z"]

    
    for i in range(len(axes)):
        slc = yt.SlicePlot(ds, axes[i], "u")
        slc.save(output_filename+"_Slice_{:s}_u.png".format(axes[i]))

    return



def time_series_movie(file_dir):

    ts = yt.load(file_dir+"plt00???.old.*")

    plot = yt.SlicePlot(ts[0], "z", "u")
    plot.set_zlim("u", -5e-3, 5e-3)
    
    fig = plot.plots["u"].figure



    def animate(i):
        ds = ts[i]
        plot._switch_ds(ds)


    animation = FuncAnimation(fig, animate, frames=len(ts))


    
    with rc_context({"mathtext.fontset":"stix"}):
        animation.save(file_dir+"test_scaling_animation.gif", writer="pillow")

    return

if __name__ == '__main__':
    file_dir = "/home/dc-kwan1/rds/rds-dirac-dp002/dc-kwan1/AMReX/wave/sphere_average_down/"



    noutputs = int(sys.argv[1])
    
    for i in range(0, noutputs):
        plot_single_slice(file_dir, i)

    
