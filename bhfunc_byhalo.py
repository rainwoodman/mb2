from isolatesubhalo import *
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_agg import FigureCanvasAgg
from gaepsi.cosmology import agn, Cosmology

def main():
    parser = makeparser()
    args = parser.parse_args()
    bhfuncMain(args) 

def bhfuncMain(args):
    snapid = args.snapid
    snapname = args.snapname
    outputdir = args.outputdir
    groupdump = numpy.fromfile(outputdir + '/grouphalotab.raw', groupdtype)
    groupdump = numpy.fromfile(outputdir + '/subhalotab.raw', subdtype)
    U = args.U
    boxsize = args.C['boxsize']
    bhmass = numpy.log10(args.open(5, 'bhmass') * 1e10)
    bhmdot = args.open(5, 'bhmass')
    bhlum = numpy.log10(agn.bhLbol(U, bhmdot))
    groupid = numpy.repeat(numpy.arange(len(groupdump)),
            groupdump['lenbytype'][:, 5])
    groupmass = numpy.log10(groupdump['mass'][groupid] * 1e10)

    figure = Figure((4, 3), dpi=200)
    ax = figure.gca()
    ax.plot(groupmass, bhmass, '. ', args.C['redshift'])
    canvas = FigureCanvasAgg(figure)
    ax.set_xlabel(r'$\log M_\mathsf{Halo} / h^{-1}M_\odot$')
    ax.set_ylabel(r'$\log M_\mathsf{BH} / h^{-1}M_\odot$')
    ax.set_ylim(6, 9)
    ax.set_xlim(10, 12)
    ax.legend()
    figure.tight_layout()
    figure.savefig('bhmass-halomass-scatter-%03d.pdf' % snapid)

    figure = Figure((4, 3), dpi=200)
    ax = figure.gca()
    ax.plot(groupmass, bhlum, '. ', label='z = %0.2g' % args.C['redshift'])
    canvas = FigureCanvasAgg(figure)
    ax.set_xlabel(r'$\log M_\mathsf{Halo} / h^{-1}M_\odot$')
    ax.set_ylabel(r'$\log L_\mathsf{BH} / L_\odot$')
    ax.set_ylim(9, 12)
    ax.set_xlim(10, 12)
    ax.legend()
    figure.tight_layout()
    figure.savefig('bhlum-halomass-scatter-%03d.pdf' % snapid)
    
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
    figure.savefig("bhfunc-%03d.pdf" % snapid)



if __name__ == '__main__':
    main()
