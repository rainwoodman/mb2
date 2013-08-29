from isolatesubhalo import *
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_agg import FigureCanvasAgg
from gaepsi.cosmology import agn, Cosmology
from sys import stdout
def main():
    """ this code finds the most massive bh and most luminous bh
        in a subhalo and save the list to a binary file
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("snap", type=SnapDir)
    args = parser.parse_args()
    mostbh(args) 

def mostbh(args):
    snap = args.snap
    g = snap.readsubhalo()
    bhpos = snap.load(5, 'pos', g)
    bhvel = snap.load(5, 'vel', g)
    bhid = snap.load(5, 'id', g)
    bhmdot = snap.load(5, 'bhmdot', g)
    bhmass = snap.load(5, 'bhmass', g)

    wrong_file_or_die(snap.subhalodir + '/subhalobhmassive.raw', len(g) *
            bhdtype.itemsize)

    massive = numpy.memmap(snap.subhalodir + '/subhalobhmassive.raw',
            mode='w+', shape=len(g), dtype=bhdtype)
    luminous = numpy.memmap(snap.subhalodir + '/subhalobhluminous.raw',
            mode='w+', shape=len(g), dtype=bhdtype)

    gfields = [
            'mass', 'len', 'pos', 'vel', 'vdisp', 'vcirc', 'rcirc', 'grnr',
            'massbytype', 'lenbytype']
    filtered = (g['lenbytype'][:, 5] > 0).nonzero()[0]
    for i in filtered:
        im = bhmass[i].argmax()
        il = bhmdot[i].argmax()
        massive[i]['pos'] = bhpos[i][im]
        massive[i]['vel'] = bhvel[i][im]
        massive[i]['id'] = bhid[i][im]
        massive[i]['bhmdot'] = bhmdot[i][im]
        massive[i]['bhmass'] = bhmass[i][im]
        luminous[i]['pos'] = bhpos[i][il]
        luminous[i]['vel'] = bhvel[i][il]
        luminous[i]['id'] = bhid[i][il]
        luminous[i]['bhmdot'] = bhmdot[i][il]
        luminous[i]['bhmass'] = bhmass[i][il]
    filtered = (g['lenbytype'][:, 5] == 0).nonzero()[0]
    print snap.snapid, 'no bh', len(filtered)
    massive['bhmass'][filtered] = numpy.nan
    massive['bhmdot'][filtered] = numpy.nan
    luminous['bhmass'][filtered] = numpy.nan
    luminous['bhmdot'][filtered] = numpy.nan
    massive.flush()
    luminous.flush()

if __name__ == '__main__':
    main()

