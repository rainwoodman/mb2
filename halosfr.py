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
    bhmdot = snap.load(5, 'bhmdot', g)
    bhmass = snap.load(5, 'bhmass', g)
    sfr = snap.load(0, 'sfr', g)

    if os.path.isfile(snap.subhalodir + '/subhalosfr.raw'):
        size = os.path.getsize(snap.subhalodir + '/subhalosfr.raw')
        N = len(g)
        if size == 4 * N:
            raise Exception("File seems to be right. quit!")
        else:
            print 'ovewrite', size, 4 , N
    halosfr = numpy.memmap(snap.subhalodir + '/subhalosfr.raw', shape=len(g),
            dtype='f4', mode='w+')
    halomdot = numpy.memmap(snap.subhalodir + '/subhalobhmdot.raw',
            shape=len(g), dtype='f4', mode='w+')
    halobhmass = numpy.memmap(snap.subhalodir + '/subhalobhmass.raw',
            shape=len(g), dtype='f4', mode='w+')
    for i in range(len(sfr)):
        halosfr[i] = sfr[i].sum()
        halomdot[i] = bhmdot[i].sum()
        halobhmass[i] = bhmass[i].sum()
        print i, len(sfr)
    halosfr.flush()
    halomdot.flush()


if __name__ == '__main__':
    main()

