from isolatesubhalo import *
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_agg import FigureCanvasAgg
from gaepsi.cosmology import agn, Cosmology

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("snap", type=SnapDir)
    args = parser.parse_args()
    halosfrMain(args) 

def halosfrMain(args):
    snap = args.snap
    g = snap.readsubhalo()

    for ptype, field in [
            (0, 'sfr'), 
            (5, 'bhmdot'), 
            (5, 'bhmass')]:
        dtype = extradtype[field]
        try:
            wrong_file_or_die(snap.filename('subhalo', field),
                dtype.itemsize * len(g))
        except: 
            continue
        input = snap.load(ptype, field, g)
        target = numpy.memmap(snap.filename('subhalo', field),
            shape=len(g), dtype=dtype, mode='w+')
        for i in range(len(g)):
            target[i] = input[i].sum()
           
        target.flush()

if __name__ == '__main__':
    main()

