"""
A package for interfacing to uvfits.  Given that there is no true standard uvfits, this is a bit
ill-defined.  This is designed to work with the MWA style uvfits that come out of cotter. 
"""

__version__ = '0.1.1'

import numpy as n
from astropy.io import fits

def echo(uv, p, d): return p, d

str2pol = {
    'I' :  1,   # Stokes Paremeters
    'Q' :  2,
    'U' :  3,
    'V' :  4,
    'rr': -1,   # Circular Polarizations
    'll': -2,
    'rl': -3,
    'lr': -4,
    'xx': -5,   # Linear Polarizations
    'yy': -6,
    'xy': -7,
    'yx': -8,
}

pol2str = {}
for k in str2pol: pol2str[str2pol[k]] = k

#primary table
itemtable1 = [
    'SIMPLE',
]

#antenna table
itemtable2 = [
    'XTENSION',
]

#               __ _ _       
#              / _(_) |      
#  _   ___   _| |_ _| |_ ___ 
# | | | \ \ / /  _| | __/ __|
# | |_| |\ V /| | | | |_\__ \
#  \__,_| \_/ |_| |_|\__|___/ 

class UV():
    """Top-level interface to a uvfits data set."""
    def __init__(self, filename, status='old', **kwargs):
        """Open a uvfits file.  status can be ('old','new','append').""" 
        assert(status in ['old', 'new', 'append'])
        #uvfits.UV.__init__(self, filename, status) #XXX needed?
        self.hdulist = fits.open(filename,mmap=True)
        self.status = status
        self.it = n.nditer(self.hdulist[0].data)
        if self.status == 'old':
            #XXX see how relevant keeping NANTS, NTIMES, NBLS is
            self.nants = self.hdulist[1].header['NAXIS2']
            self.nbls = (self.nants * (self.nants-1))/2 #generically true, or can bls be flagged?
            self.ntimes = self.hdulist[0].header['GCOUNT']/self.nbls
            self.npols = self.hdulist[0].header['NAXIS3'] #XXX how generic?
            self.pol = pol2str[self.hdulist[0].header['CRVAL3'] + self.hdulist[0].header['CDELT3']]
            self.iterindex = 0
            self.polindex = 0
    def items(self):
        """Return a dictionary of available tables and their header items."""
        tables = {}
        for hdu in self.hdulist:
            tables[hdu.name] = hdu.header.keys()
        return tables
    def _rdhd(self, name):
        """Provide read access to header items."""
        for hdu in self.hdulist:
            if name != 'pol':
                try: rv = hdu.header[name]
                except(KeyError): continue
            else:
                rv = self.pol
        return rv 
    def _wrhd(self, table, name, val):
        """Provide write access to header items."""
        if name != 'pol':
            self.hdulist[table].header[name] = val
        else:
            self.pol = val
    def __getitem__(self, name):
        """Allow access to variables and header items via uv[name]."""
        try: 
            rv = self._rdhd(name) 
            return rv
        except(UnboundLocalError): 
            print 'Variable does not exist'
            return None
    def __setitem__(self, table, name, val):
        """Allow setting variables and header items via uv[name] = val."""
        self._wrhd(table,name,val)
    def select(self, name, n1, n2, include=1): #XXX need change
        """Choose which data are returned by read().  
            name    This can be: 'decimate','time','antennae','visibility',
                    'uvrange','pointing','amplitude','window','or','dra',
                    'ddec','uvnrange','increment','ra','dec','and', 'clear',
                    'on','polarization','shadow','auto','dazim','delev'
            n1,n2   Generally this is the range of values to select. For
                    'antennae', this is the two antennae pair to select
                    (indexed from 0); a -1 indicates 'all antennae'.
                    For 'decimate', n1 is every Nth integration to use, and
                    n2 is which integration within a block of N to use.
                    For 'shadow', a zero indicates use 'antdiam' variable.
                    For 'on','window','polarization','increment','shadow' only
                    p1 is used.
                    For 'and','or','clear','auto' p1 and p2 are ignored.
            include If true, the data is selected. If false, the data is
                    discarded. Ignored for 'and','or','clear'."""
        if name == 'antennae':
            jarr = (self.hdulist[0].data['BASELINE'] % 256 - 1).astype(int)
            iarr = ((self.hdulist[0].data['BASELINE']  - (jarr+1))/256 - 1).astype(int)
            if n1 != -1: iinds = n.where(iarr == n1)
            else: iinds = n.arange(len(iarr))
            if n2 != -1: jinds = n.where(jarr == n2)
            else: jinds = n.arange(len(jarr))
            inds = n.intersect1d(iinds,jinds)
            self.hdulist[0].data = self.hdulist[0].data[inds]
        elif name == 'time':
            self.hdulist[0].data = self.hdulist[0].data[n1*self.nbls:n2*self.nbls]
        elif name == 'polarization':
            print 'Polarization select not yet functional' #XXX may be hard with uvfits
        else:
            print '%s is not a valid option for the uvfits module' % name           
        self.it = n.nditer(self.hdulist[0].data)
    def read(self, raw=False):
        """Return the next data record. 'raw' causes data and flags 
        to be returned seperately."""
        if self.polindex < self.npols - 1: rec = self.it.value
        else: rec = self.it.next()
        u,v,w,bl,t = rec['UU'], rec['VV'], rec['WW'], rec['BASELINE'], rec['DATE'] + self.hdulist[0].header['PZERO4']
        i,j = bl2ij(bl) #XXX test with non-PAPER uvfits
        preamble = (n.array([u,v,w]).astype(float),t,(i,j))
        datpol = n.rollaxis(rec['DATA'].squeeze(),0,3)[self.polindex]
        self.pol = pol2str[self.hdulist[0].header['CRVAL3'] + self.hdulist[0].header['CDELT3']*self.polindex]
        data = datpol[0]+datpol[1]*1j
        flags = n.logical_not(-1*datpol[2])
        if self.polindex < self.npols - 1: 
            self.polindex += 1
        else: 
            self.polindex = 0
        if raw: 
            return preamble, data, flags
        else:
            data = n.ma.array(data, mask=flags)
            return preamble, data
    def all(self, raw=False):
        if self.hdulist[0].header['PTYPE4'] != 'DATE': print 'WARNING:, PZERO4 is not a date field; times are likely to be incorrect.' #XXX better way to do this?
        while True:
            try: yield self.read(raw=raw)
            except(StopIteration): return
    def write(self, preamble, data, flags=None):
        """Write the next data record.  data must be a complex, masked
        array.  preamble must be (uvw, t, (i,j)), where uvw is an array of 
        u,v,w, t is the Julian date, and (i,j) is an antenna pair."""
        if data is None: return
        if not flags is None: flags = n.logical_not(flags)
        elif len(data.mask.shape) == 0:
            flags = n.ones(data.shape)
            data = data.unmask()
        else:
            flags = n.logical_not(data.mask)
            #data = data.filled(0)
            data = data.data
        self.raw_write(preamble,data.astype(n.complex64),flags.astype(n.int32))
    def init_from_uv(self, uv, override={}, exclude=[]):
        """Initialize header items and variables from another UV.  Those in 
        override will be overwritten by override[k], and tracking will be 
        turned off (meaning they will not be updated in pipe()).  Those in
        exclude are omitted completely."""
        for k in uv.items():
            if k in exclude: continue
            elif k in override: self._wrhd(k, override[k])
            else: self._wrhd(k, uv[k])
        self.vartable = {}
        for k in uv.vars():
            if k in exclude: continue
            # I don't understand why reading 'corr' segfaults miriad,
            # but it does.  This is a cludgy work-around.
            elif k == 'corr': continue
            elif k in override:
                self.vartable[k] = uv.vartable[k]
                self._wrvr(k, uv.vartable[k], override[k])
            else:
                self.vartable[k] = uv.vartable[k]
                self._wrvr(k, uv.vartable[k], uv[k])
                uv.trackvr(k, 'c') # Set to copy when copyvr() called 
    def pipe(self, uv, mfunc=echo, append2hist='', raw=False):
        """Pipe in data from another UV through the function
        mfunc(uv,preamble,data), which should return (preamble,data).  If 
        mfunc is not provided, the dataset will just be cloned, and if the 
        returned data is None, it will be omitted.  The string 'append2hist' 
        will be appended to history."""
        self._wrhd('history', self['history'] + append2hist)
        # Pipe all data through mfunc
        if raw:
            for p,d,f in uv.all(raw=raw):
                np, nd, nf = mfunc(uv, p, d, f)
                self.copyvr(uv)
                self.write(np, nd, nf)
        else:
            for p, d in uv.all():
                np, nd = mfunc(uv, p, d)
                self.copyvr(uv)
                self.write(np, nd)
    def add_var(self, name, type):
        """Add a variable of the specified type to a UV file."""
        self.vartable[name] = type
    def rewind(): #XXX needs to be defined
        print "You haven't written this code yet."

def bl2ij(bl):
    j = int(bl % 256 - 1)
    i = int((bl - (j+1))/256 - 1)
    return i,j

def ij2bl(i, j):
    if i > j: i,j = j,i
    return int(256*(i+1) + (j+1))

