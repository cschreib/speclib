ob = 'A'
arm = 'NIR'
epoch = 'XSHOO.2018-02-19T00:33:41.063-NIR'
target = '7329'
project_dir = '/media/cschreib/leiden/backup/reductions/xshooter/reduced/7329/'

root_name = project_dir+'reduced/'+ob+'/stacked_'+target+'_'+epoch
profile_file = root_name+'_profile.fits'
f = mrdfits(profile_file, 1, hdr, /silent)
e = mrdfits(profile_file, 2, /silent)

npix = n_elements(f)
nmodel = long(sxpar(hdr, 'NSRC'))
m = fltarr(nmodel, npix)
for i=0, nmodel-1 do begin
    m[i,*] = mrdfits(profile_file, 3+i, /silent)
endfor

tmp = mrdfits(root_name+'_spec2d.fits', 1, hdr, /silent)
aspix = float(sxpar(hdr, 'CDELT2'))

x = (indgen(npix) - npix/2)*aspix
plot, x, f, xtit='slit offset [arcsec]'
errplot, x, f-e, f+e

for i=0, nmodel-1 do begin
    oplot, x, m[i,*], col='ff'x
endfor

end
