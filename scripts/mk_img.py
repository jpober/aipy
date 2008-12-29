#! /usr/bin/env python
"""
This is a general-purpose script for making images from MIRIAD UV files.  Data
(optionally selected for baseline, channel) are read from the file, phased
to a provided position, normalized for passband/primary beam effects, gridded
to a UV matrix, imaged, and optionally deconvolved by a corresponding PSF to
produce a clean image.

Author: Aaron Parsons
"""

import aipy as a, numpy as n, sys, optparse, ephem, pyfits

o = optparse.OptionParser()
o.set_usage('plot_img.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True, loc2=True,
    src=True, dec=True)
o.add_option('-r', '--residual', dest='residual', action='store_true',
    help='Display only the residual in the 4th panel (otherwise display sum of clean image and residual).')
o.add_option('-d', '--deconv', dest='deconv', default='mem',
    help='Attempt to deconvolve the dirty image by the dirty beam using the specified deconvolver (none,mem,lsq,cln,ann).')
o.add_option('-o', '--outfile', dest='outfile', default='img.fits',
    help='FITS file to save deconvolved data to.')
o.add_option('--var', dest='var', type='float', default=.6,
    help='Starting guess for variance in maximum entropy fit (defaults to variance of dirty image.')
o.add_option('--tol', dest='tol', type='float', default=1e-6,
    help='Tolerance for successful deconvolution.  For annealing, interpreted as cooling speed.')
o.add_option('--maxiter', dest='maxiter', type='int', default=200,
    help='Number of allowable iterations per deconvolve attempt.')
o.add_option('-u', '--uniform', dest='uniform', type='float', default=0,
    help="Use uniform (rather than natural) weighting for uv bins that have a weight above the specified fraction of the maximum weighting.")
o.add_option('--skip_amp', dest='skip_amp', action='store_true',
    help='Do not use amplitude information to normalize visibilities.')
o.add_option('--skip_bm', dest='skip_bm', action='store_true',
    help='Do not weight visibilities by the strength of the primary beam.')
o.add_option('--size', dest='size', type='int', default=200,
    help='Size of maximum UV baseline.')
o.add_option('--res', dest='res', type='float', default=0.5,
    help='Resolution of UV matrix.')
o.add_option('--no_w', dest='no_w', action='store_true',
    help="Don't use W projection.")
o.add_option('--dyn_rng', dest='dyn_rng', type='float', default=3.,
    help="Dynamic range in color of image (log).")
o.add_option('--phs', dest='phs', action='store_true',
    help="If plotting the UV matrix, show phases.")
o.add_option('--buf_thresh', dest='buf_thresh', default=1.8e6, type='float',
    help='Maximum amount of data to buffer before gridding.  Excessive gridding takes performance hit, but if buffer exceeds memory available... ouch.')
opts, args = o.parse_args(sys.argv[1:])

# Parse command-line options
locs = a.scripting.files_to_locs(opts.loc, args, sys.argv)
aas = {}
for L in locs:
    uv = a.miriad.UV(locs[L][0])
    chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    aa = a.loc.get_aa(L, uv['sdf'], uv['sfreq'], uv['nchan'])
    aa.select_chans(chans)
    aas[L] = aa
    afreqs = aa.ants[0].beam.afreqs
    cfreq = n.average(afreqs)
src = a.scripting.parse_srcs(opts.src, force_cat=True)
s = src.values()[0]
if opts.no_w: im = a.img.Img(opts.size, opts.res, mf_order=0)
else: im = a.img.ImgW(opts.size, opts.res, mf_order=0)
DIM = int(opts.size/opts.res)
del(uv)

# Gather data
uvw, dat, wgt = [], [], []
cnt, curtime = 0, None
for L in locs:
  aa = aas[L]
  # Read each file
  for filename in locs[L]:
    sys.stdout.write('.'); sys.stdout.flush()
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    # Read all data from each file
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        if curtime != t:
            curtime = t
            cnt = (cnt + 1) % opts.decimate
            if cnt == 0:
                aa.set_jultime(t)
                src.compute(aa)
                s_eq = src.get_crds('eq', ncrd=3)
                aa.sim_cache(s_eq)
        if cnt != 0: continue
        d,f = d.take(chans), f.take(chans)
        if not opts.skip_amp: d /= aa.passband(i,j)
        try:
            # Throws PointingError if not up:
            d = aa.phs2src(d, s, i, j)
            xyz = aa.gen_uvw(i,j,src=s)
            if not opts.skip_bm:
                # Calculate beam strength for weighting purposes
                w = aa.bm_response(i,j,pol=opts.pol).squeeze()
                # For optimal SNR, down-weight data that is already attenuated 
                # by beam  by another factor of the beam response (modifying 
                # weight accordingly).
                d *= w; w *= w
            else: w = n.ones(d.shape, dtype=n.float)
        except(a.ant.PointingError): continue
        valid = n.logical_not(f)
        d = d.compress(valid)
        if len(d) == 0: continue
        dat.append(d)
        uvw.append(xyz.compress(valid, axis=0))
        wgt.append(w.compress(valid))
        # If data buffer is full, grid data
        if len(dat) * len(chans) > opts.buf_thresh:
            sys.stdout.write('|'); sys.stdout.flush()
            dat = n.concatenate(dat)
            uvw = n.concatenate(uvw); uvw.shape = (uvw.size / 3, 3)
            wgt = n.concatenate(wgt).flatten()
            # Grid data into UV matrix
            uvw,dat,wgt = im.append_hermitian(uvw,dat,wgt)
            im.put(uvw, dat, wgt)
            uvw, dat, wgt = [], [], []

# Grid remaining data into UV matrix
if len(uvw) == 0: raise ValueError('No data to plot')
sys.stdout.write('|\n'); sys.stdout.flush()
dat = n.concatenate(dat)
uvw = n.concatenate(uvw); uvw.shape = (uvw.size / 3, 3)
wgt = n.concatenate(wgt).flatten()
uvw,dat,wgt = im.append_hermitian(uvw,dat,wgt)
im.put(uvw, dat, wgt)
if opts.uniform > 0: im.uniform_wgt(thresh=opts.uniform)

# Form dirty images/beams
im_img = im.image((DIM/2, DIM/2))
bm_img = im.bm_image(term=0)
bm_gain = a.img.beam_gain(bm_img)

# Deconvolve
if opts.deconv == 'mem':
    cl_img,info = a.deconv.maxent_findvar(im_img, bm_img, f_var0=opts.var,
        maxiter=opts.maxiter, verbose=True, tol=opts.tol, maxiterok=True)
elif opts.deconv == 'lsq':
    cl_img,info = a.deconv.lsq(im_img, bm_img, 
        maxiter=opts.maxiter, verbose=True, tol=opts.tol)
elif opts.deconv == 'cln':
    cl_img,info = a.deconv.clean(im_img, bm_img, 
        maxiter=opts.maxiter, verbose=True, tol=opts.tol)
elif opts.deconv == 'ann':
    cl_img,info = a.deconv.anneal(im_img, bm_img, maxiter=opts.maxiter, 
        cooling=lambda i,x: opts.tol*(1-n.cos(i/50.))*(x**2), verbose=True)

print 'Gain of dirty beam:', bm_gain
bm_img = im.bm_image((DIM/2,DIM/2), term=0)

if opts.deconv != 'none':
    rs_img = info['res'] / bm_gain
    if not opts.residual: rs_img += cl_img
else: rs_img = im_img

print 'Saving data to', opts.outfile
rs_img.shape = rs_img.shape + (1,1)
cen = ephem.Equatorial(s.ra, s.dec, epoch=aa.epoch)
cen = ephem.Equatorial(cen, epoch=ephem.J2000)
L,M = im.get_LM()
a.img.to_fits(opts.outfile, rs_img, clobber=True,
    object=opts.src, obs_date=str(aa.date), 
    ra=cen.ra*a.img.rad2deg, dec=cen.dec*a.img.rad2deg, epoch=2000.,
    d_ra=L[-1,-1]*a.img.rad2deg, d_dec=M[-1,-1]*a.img.rad2deg,
    freq=n.average(aa.ants[0].beam.afreqs), 
)
