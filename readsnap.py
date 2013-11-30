import numpy

class Lazy(object):
    def __init__(self, calculate_function):
        self._calculate = calculate_function

    def __get__(self, obj, _=None):
        if obj is None:
            return self
        value = self._calculate(obj)
        setattr(obj, self._calculate.func_name, value)
        return value

class SnapFile(numpy.memmap):
    def _getblock(self, N, dtype, ptr):
        dtype = numpy.dtype(dtype)
        if isinstance(N, tuple) and len(N) == 0:
            N = 1
            t = True
        else: 
            t = False
        bs1 = self[ptr:ptr+4].view(dtype='i4')[0]
        assert bs1 == dtype.itemsize * N
        ptr += 4
        block = self[ptr:ptr+bs1].view(dtype=dtype.base)
        block = block.reshape([-1] + list(dtype.shape))
        if t: block = block[0]
        ptr += bs1
        bs2 = self[ptr:ptr+4].view(dtype='i4')[0]
        ptr += 4
        assert bs2 == bs1
        return block, ptr

    def __new__(cls, filename):
        self = numpy.memmap.__new__(cls, filename, mode='r', dtype='u1')
        ptr = 0
        self.header, ptr = self._getblock((), 
            dtype=numpy.dtype([
              ('N', ('u4', 6)),
              ('mass', ('f8', 6)),
              ('time', 'f8'),
              ('redshift', 'f8'),
              ('flag_sfr', 'i4'),
              ('flag_feedback', 'i4'),
              ('Ntot_low', ('u4', 6)),
              ('flag_cool', 'i4'),
              ('Nfiles', 'i4'),
              ('boxsize', 'f8'),
              ('OmegaM', 'f8'),
              ('OmegaL', 'f8'),
              ('h', 'f8'),
              ('flag_sft', 'i4'),
              ('flag_met', 'i4'),
              ('Ntot_high', ('u4', 6)),
              ('flag_entropy', 'i4'),
              ('flag_double', 'i4'),
              ('flag_ic_info', 'i4'),
              ('flag_lpt_scalingfactor', 'i4'),
              ('flag_pressure_entropy', 'i1'),
              ('Ndims', 'i1'),
              ('densitykerneltype', 'i1'),
              ('unused', ('u1', 45))]), ptr=ptr)
        N = self.header['N']
        ftype = [numpy.dtype(numpy.float32), numpy.dtype(numpy.float64)][self.header['flag_double']]
        self.pos, ptr = self._getblock(N.sum(), (ftype, 3), ptr)
        self.vel, ptr = self._getblock(N.sum(), (ftype, 3), ptr)
        bs = self[ptr:ptr+4].view(dtype='i4')[0]
        if bs / N.sum() == 8:
            idtype = numpy.dtype('u8')
        else:
            idtype = numpy.dtype('u4')
        self.id, ptr = self._getblock(N.sum(), idtype, ptr)
        Nmass = N * (self.header['mass'] == 0)
        mass, ptr = self._getblock(Nmass.sum(), ftype, ptr)
        self._rawmass = mass
        self.N = N
        return self
    @Lazy
    def ptype(self):
        return numpy.repeat(range(6), self.N)
    @Lazy
    def mass(self):
        """ mass of particles. correctly fill the entries if mass in header is
        nonzero"""
        N = self.N
        Nmass = N * (self.header['mass'] == 0)
        mass = self._rawmass
        start = numpy.concatenate(([0], Nmass.cumsum()))
        end = Nmass.cumsum()
        fullmass = numpy.concatenate(
                [
                    mass[start[i]:end[i]] if Nmass[i] > 0 else
                    numpy.repeat(self.header['mass'][i], N[i]) for i in
                    range(6)])
        return fullmass 
