from isolatesubhalo import *
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_agg import FigureCanvasAgg
from gaepsi.cosmology import agn, Cosmology

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("snapid", type=int)
    args = parser.parse_args()
    halosfrMain(args) 


def halosfrMain(args):
    snapid = args.snapid
    outputdir = ROOT + '/subhalos/%03d' % snapid
    subdump = numpy.fromfile(outputdir + '/subhalotab.raw', subdtype)
    
    bhmdot = numpy.fromfile(outputdir + '/5/bhmdot.raw', dtype)
    
