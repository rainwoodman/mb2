import numpy
import os.path

subdtype = numpy.dtype([
    ('mass', 'f4'), 
    ('len', 'i4'), 
    ('pos', ('f4', 3)), 
    ('vel', ('f4', 3)),   
    ('vdisp', 'f4'),  
    ('vcirc', 'f4'),  
    ('rcirc', 'f4'),  
    ('parent', 'i4'),  
    ('massbytype', ('f4', 6)),  
    ('lenbytype', ('u4',6)), 
    ('unused', 'f4'),    
    ('groupid', 'u4'), 
   ])
groupdtype = numpy.dtype([
    ('mass', 'f4'), 
    ('len', 'i4'), 
    ('pos', ('f4', 3)), 
    ('vel', ('f4', 3)), 
    ('nhalo', 'i4'),   
    ('massbytype', ('f4', 6)),  
    ('lenbytype', ('u4',6)), 
   ])

pdtype = numpy.dtype([('pos', ('f4', 3)),
    ('vel', ('f4', 3)),
    ('mass', 'f4'),
    ('id', 'u8'),
    ('type', 'u1')
    ])

bhdtype = numpy.dtype([
    ('pos', ('f4', 3)),
    ('vel', ('f4', 3)),
    ('id', 'u8'),
    ('bhmass', 'f8'),
    ('bhmdot', 'f8')])

# extra properties of a subhalo / group
extradtype = numpy.dtype([
    ('bhmassive', bhdtype),
    ('bhluminous', bhdtype),
    ('sfr', 'f4'),
    ('bhmass', 'f4'),
    ('bhmdot', 'f4')])

class SnapDir(object):
    def __init__(self, snapid, ROOT):
        """ the default schema is a fake one, generated from header.txt
            subclass of the writeable version of SnapDir use a real schema.

            this is to avoid pulling in gaepsi in this public code
        """
        snapid = int(snapid)
        self.snapid = snapid
        self.subhalodir = ROOT + '/subhalos/%03d' % snapid
        self.subhalofile = self.subhalodir + '/subhalotab.raw'
        self.groupfile = self.subhalodir + '/grouphalotab.raw'
        self.headerfile = self.subhalodir + '/header.txt'
        try:
            for line in file(self.headerfile, 'r'):
                if line.startswith('flag_double(header) = 0'):
                   flag_double = False
                if line.startswith('flag_double(header) = 1'):
                   flag_double = True
            class fakeschema:
                def __getitem__(self, index, flag_double=flag_double):
                    if flag_double: return numpy.dtype('f8')
                    return numpy.dtype('f4')
            self.schema = fakeschema()
        except IOError:
            self.schema = None

    def readsubhalo(self, g=None):
        rt = numpy.memmap(self.subhalofile, mode='r', dtype=subdtype)
        if g is not None:
            rt = packarray(rt, g['nhalo'] + 1)
        return rt

    def readgroup(self):
        return numpy.memmap(self.groupfile, mode='r', dtype=groupdtype)

    def load(self, type, comp, g=None):
        """ this will read in property comp of partile type type.
            if g is given (either from readsubhalo or readgroup)
            the returned array A will be chunked so that
            A[0] is the property of particles in the first group
            A[1] is the property of particles in the second group
            so on.
            
            snapdir = Snapdir(18, './')
            g = snapdir.readsubhalo()
            v = snapdir.load(4, 'vel', g)
            m = snapdir.load(4, 'mass', g)
           
            mv2 = m * v ** 2
            mv2 = array([i.sum() for i in mv2])
           
            to readin extra properties of subhalos, use
            snapdir.load('subhalo', 'bhmassive')
            (for the most massive bh in the subhalo)
            or
            snapdir.load('subhalo', 'bhluminous')
            (for the most luminous bh in the subhalo)
            snapdir.load('subhalo', 'sfr')
            (for the sfr of the subhalo)

            this in princple also works for a group, though
            the group properties are not assembled yet.

            see extradtype for a list
        """
        if isinstance(type, basestring):
            dtype = extradtype[comp]
        else:
            if comp in pdtype.fields:
                dtype = pdtype[comp]
            else:
                dtype=self.schema[comp].dtype
        size = os.path.getsize(self.filename(type, comp))
        if size == 0:
            rt = numpy.fromfile(self.filename(type, comp),
                    dtype=dtype)
        else:
            rt = numpy.memmap(self.filename(type, comp),
                            mode='r',
                            dtype=dtype)
        if g is None:
            return rt
        rt = packarray(rt, g['lenbytype'][:, type])
        return rt
    def open(self, type, comp, mode='r'):
        return file(self.filename(type, comp), mode=mode)
    def filename(self, type, comp):
        """ the file name of a type/comp """
        if isinstance(type, basestring):
            return self.subhalodir + '/%s%s.raw' % (type, comp)
        else:
            return self.subhalodir + '/%d/%s.raw' % (type, comp)

class packarray(numpy.ndarray):
  """ A packarray packs/copies several arrays into the same memory chunk.

      It feels like a list of arrays, but because the memory chunk is continuous,
      
      arithmatic operations are easier to use(via packarray)
  """
  def __new__(cls, array, start=None, end=None):
    """ if end is none, start contains the sizes. 
        if start is also none, array is a list of arrays to concatenate
    """
    self = array.view(type=cls)
    if end is None and start is None:
      start = numpy.array([len(arr) for arr in array], dtype='intp')
      array = numpy.concatenate(array)
    if end is None:
      sizes = start
      self.start = numpy.zeros(shape=len(sizes), dtype='intp')
      self.end = numpy.zeros(shape=len(sizes), dtype='intp')
      self.end[:] = sizes.cumsum()
      self.start[1:] = self.end[:-1]
    else:
      self.start = start
      self.end = end
    self.A = array
    return self
  @classmethod
  def adapt(cls, source, template):
    """ adapt source to a packarray according to the layout of template """
    if not isinstance(template, packarray):
      raise TypeError('template must be a packarray')
    return cls(source, template.start, template.end)

  def __repr__(self):
    return 'packarray: %s, start=%s, end=%s' % \
          (repr(self.A), 
           repr(self.start), repr(self.end))
  def __str__(self):
    return repr(self)

  def copy(self):
    return packarray(self.A.copy(), self.start, self.end)

  def compress(self, mask):
    count = self.end - self.start
    realmask = numpy.repeat(mask, count)
    return packarray(self.A[realmask], self.start[mask], self.end[mask])

  def __getitem__(self, index):
    if isinstance(index, basestring):
      return packarray(self.A[index], self.end - self.start)

    if isinstance(index, slice) :
      start, end, step = index.indices(len(self))
      if step == 1:
        return packarray(self.A[self.start[start]:self.end[end]],
            self.start[start:end] - self.start[start],
            self.end[start:end] - self.start[start])

    if isinstance(index, (list, numpy.ndarray)):
      return packarray(self.A, self.start[index], self.end[index])

    if numpy.isscalar(index):
      start, end = self.start[index], self.end[index]
      if end > start: return self.A[start:end]
      else: return numpy.empty(0, dtype=self.A.dtype)
    raise IndexError('unsupported index type %s' % type(index))

  def __len__(self):
    return len(self.start)

  def __iter__(self):
    for i in range(len(self.start)):
      yield self[i]

  def __reduce__(self):
    return packarray, (self.A, self.end - self.start)

  def __array_wrap__(self, outarr, context=None):
    return packarray.adapt(outarr.view(numpy.ndarray), self)
