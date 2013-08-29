import numpy

subdtype = numpy.dtype([
    ('mass', 'f4'), 
    ('len', 'i4'), 
    ('pos', ('f4', 3)), 
    ('vel', ('f4', 3)),   
    ('vdisp', 'f4'),  
    ('vcirc', 'f4'),  
    ('rcirc', 'f4'),  
    ('grnr', 'i4'),  
    ('massbytype', ('f4', 6)),  
    ('lenbytype', ('u4',6)), 
    ('unused', 'f4'),    
    ('unused2', 'u4'), 
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


class SnapDir(object):
    def __init__(self, snapid, ROOT):
        snapid = int(snapid)
        self.snapid = snapid
        self.subhalodir = ROOT + '/subhalos/%03d' % snapid
        self.subhalofile = self.subhalodir + '/subhalotab.raw'
        self.groupfile = self.subhalodir + '/grouphalotab.raw'

    def readsubhalo(self):
        return numpy.fromfile(self.subhalofile, dtype=subdtype)

    def readgroup(self):
        return numpy.fromfile(self.groupfile, dtype=groupdtype)

    def load(self, type, comp, g=None):
        if comp in pdtype.fields:
            dtype = pdtype[comp]
        else:
            dtype=self.schema[comp].dtype
        rt = numpy.fromfile(self.filename(type, comp),
                            dtype=dtype)
        if g is None:
            return rt
        rt = packarray(rt, g['lenbytype'][:, type])
        return rt
    def open(self, type, comp, mode='r'):
        return file(self.filename(type, comp), mode=mode)
    def filename(self, type, comp):
        """ the file name of a type/comp """
        return self.subhalodir + '/%d/%s.raw' % (type, comp)

class packarray(numpy.ndarray):
  """ A packarray packs/copies several arrays into the same memory chunk.

      It feels like a list of arrays, but because the memory chunk is continuous,
      
      arithmatic operations are easier to use(via packarray.A)
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
#    return numpy.ndarray.__array_wrap__(self.view(numpy.ndarray), outarr, context)
