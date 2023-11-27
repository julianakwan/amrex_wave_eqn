import yt
from yt.frontends.boxlib.data_structures import AMReXDataset
from matplotlib import rc_context
from matplotlib.animation import FuncAnimation
import sys
import glob


def plot_single_slice(file_dir, stepno):
    #Directory to AMReX outputs
    input_filename = glob.glob(file_dir+"plt{:0>5d}.old.*".format(stepno))


    #Output directory for plots
    output_filename = file_dir+"plt{:0>5d}".format(stepno)
    
    ds = AMReXDataset(input_filename[0])

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

    ts = yt.load(file_dir+"plt00???")

    plot = yt.SlicePlot(ts[0], "z", "u")

    fig = plot.plots["u"].figure



    def animate(i):
        ds = ts[i]
        plot._switch_ds(ds)


    animation = FuncAnimation(fig, animate, frames=len(ts))


    
    with rc_context({"mathtext.fontset":"stix"}):
        animation.save("animation.gif", writer="pillow")

    return

if __name__ == '__main__':
    file_dir = "/home/dc-kwan1/rds/rds-dirac-dp002/dc-kwan1/AMReX/wave/orig_4th_order/"



    noutputs =1
    
    for i in range(0, noutputs):
        plot_single_slice(file_dir, i)

    
