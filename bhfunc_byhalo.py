from isolatesubhalo import *
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_agg import FigureCanvasAgg
from gaepsi.cosmology import agn, Cosmology

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('snap', type=SnapDir)
    args = parser.parse_args()
    bhfuncMain(args) 

def bhfuncMain(args):
    snap = args.snap
    g = snap.readsubhalo()
    U = snap.U
    g = g[g['lenbytype'][:, 5] > 0]
    print len(g)
    boxsize = snap.C['boxsize']
    bhmass = numpy.log10(snap.load(5, 'bhmass') * 1e10)
    bhmdot = snap.load(5, 'bhmass')
    bhlum = numpy.log10(agn.bhLbol(U, bhmdot))
    
    groupmass = numpy.log10(g['mass'] * 1e10)
    print groupmass.max(), groupmass.min()
    if False:
        bhmass = numpy.array([i.max() if len(i) > 0 else 0 for i in bhmass])
        bhlum= numpy.array([i.max() if len(i) > 0 else 0 for i in bhlum])
    else:
        bhmass = bhmass
        bhlum = bhlum
        groupmass = numpy.repeat(groupmass, g['lenbytype'][:, 5])
    figure = Figure((3, 3), dpi=200)
    ax = figure.gca()
    ax.plot(groupmass, bhmass, '. ', label='z = %0.2g' % snap.C['redshift'], color='k')
    canvas = FigureCanvasAgg(figure)
    ax.set_xlabel(r'$\log\, h M_\mathsf{Halo} / M_\odot$', fontsize='x-large')
    ax.set_ylabel(r'$\log\, h M_\mathsf{BH} / M_\odot$', fontsize='x-large')
    ax.set_ylim(6, 10 + 0.2)
    ax.set_yticks([7, 8, 9, 10])
    ax.set_xlim(10 - 0.2, numpy.nanmax(groupmass) + 0.2)
    ax.set_xticks([9, 10, 11, 12, 13, 14])
    ax.set_title('z = %0.2g' % snap.C['redshift'], position=(0.2, 0.8))
    figure.tight_layout()
    figure.savefig('bhmass-halomass-scatter-%03d.pdf' % snap.snapid)
    figure.savefig('bhmass-halomass-scatter-%03d.png' % snap.snapid)

    figure = Figure((3, 3), dpi=200)
    ax = figure.gca()
    ax.plot(groupmass, bhlum, '. ', label='z = %0.2g' % snap.C['redshift'], color='k')
    canvas = FigureCanvasAgg(figure)
    ax.set_xlabel(r'$\log\, h M_\mathsf{Halo} / M_\odot$', fontsize='x-large')
    ax.set_ylabel(r'$\log\, L_\mathsf{BH} / L_\odot$', fontsize='x-large')
    ax.set_xticks([9, 10, 11, 12, 13, 14])
    ax.set_yticks([9, 10, 11, 12, 13, 14])
    ax.set_ylim(9.5, 12.5)
    ax.set_xlim(10 - 0.2, numpy.nanmax(groupmass) + 0.2)
    ax.set_title('z = %0.2g' % snap.C['redshift'], position=(0.2, 0.8))
    figure.tight_layout()
    figure.savefig('bhlum-halomass-scatter-%03d.pdf' % snap.snapid)
    figure.savefig('bhlum-halomass-scatter-%03d.png' % snap.snapid)
    
    gmbins = numpy.linspace(7, 14, 6, endpoint=True)
    dig = numpy.digitize(groupmass, gmbins)
    print bhmass.max(), bhmass.min()
    print bhlum.max(), bhlum.min()
    bmbins = numpy.linspace(5, 10, 11, endpoint=True)
    blbins = numpy.linspace(7, 14, 30, endpoint=True)
    figure = Figure(dpi=200, figsize=(12, 4))

    print bhmass.shape
    print groupmass.shape
    print dig.shape
    ax = figure.add_subplot(133)
    h, b = numpy.histogram(groupmass, gmbins)
    h = h / numpy.diff(gmbins)
    h /= (boxsize / 1000) ** 3
    bc = 0.5 * (b[1:] + b[:-1])
    ax.plot(bc, h, 'x-', label='Halo mass func')
    ax.set_xticks(b)
    ax.set_xlabel(r'$M \,\,h^{-1}M_\odot$')
    ax.set_ylabel(r'$dN / V d log M $')
    ax.set_yscale('log')
    
    for bin in range(1, len(gmbins) + 1):
        sel = dig == bin
        ax = figure.add_subplot(131)
        h, junk = numpy.histogram(bhmass[sel], 
                bins=bmbins)
        h = h / numpy.diff(bmbins)
        h /= (boxsize / 1000)** 3
        bc = 0.5 * (bmbins[1:] + bmbins[:-1])
        ax.plot(bc, h, 'x-', label='%g' %gmbins[bin - 1])

        ax = figure.add_subplot(132)
        h, jlulnk = numpy.histogram(bhlum[sel], 
                bins=blbins)
        h = h / numpy.diff(blbins)
        h /= (boxsize / 1000)** 3
        bc = 0.5 * (blbins[1:] + blbins[:-1])
        ax.plot(bc, h, 'x-', label='%g' % gmbins[bin - 1])

    ax = figure.add_subplot(131)
    ax.set_ylabel(r'$dN / V d \log M $')
    ax.set_xlabel(r'log $h^{-1}M_\odot$')
    ax.set_yscale('log')
    ax.legend()
    ax = figure.add_subplot(132)
    ax.set_ylabel(r'$dN / V d \log L $')
    ax.set_yscale('log')
    ax.set_xlabel(r'log $L_\odot$')
    ax.legend()
    figure.tight_layout()
    canvas = FigureCanvasAgg(figure)
    figure.savefig("bhfunc-%03d.pdf" % snap.snapid)



if __name__ == '__main__':
    main()
