#! /usr/bin/env python
"""
Rotate zenith UV data to a particular source.  

Author: Aaron Parsons
"""

import aipy as a, numpy as n, sys, os, optparse

o = optparse.OptionParser()
o.set_usage('phs2src.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, loc=True, src=True)
o.add_option('--setphs', dest='setphs', action='store_true',
    help='Instead of rotating phase, assign a phase corresponding to the specified source.')
opts,args = o.parse_args(sys.argv[1:])

# Parse command-line options
uv = a.miriad.UV(args[0])
aa = a.loc.get_aa(opts.loc, uv['sdf'], uv['sfreq'], uv['nchan'])
srclist,cutoff = a.scripting.parse_srcs(opts.src)
src = a.loc.get_catalog(opts.loc, srclist, cutoff).values()[0]
del(uv)

# A pipe to use for phasing to a source
curtime,cache,top = None, {}, None
def phs(uv, p, d, f):
    global curtime, cache, top
    uvw, t, (i,j) = p
    if curtime != t:
        curtime = t
        aa.set_jultime(t)
        src.compute(aa)
    if i == j: return p, d, f
    try:
        if opts.setphs: d = aa.unphs2src(n.abs(d), src, i, j)
        else: d = aa.phs2src(d, src, i, j)
    except(a.ant.PointingError): d *= 0
    return p, d, f

# Process data
for filename in args:
    print filename
    uvofile = filename + '.' + opts.src
    if os.path.exists(uvofile):
        print 'File exists: skipping'
        continue
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(uvofile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=phs, raw=True)

