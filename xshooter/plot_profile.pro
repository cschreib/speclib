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

plot, x, f, xtit='slit offset [arcsec]', ytit='flux [normalized]', $
    yr=[min(f-e), max(f+e)], /nodata, charsize=2

errplot, x, f-e, f+e, col='cccccc'x, thick=2
oplot, x, f, col='888888'x, thick=3

colors = ['ff'x, 'ff0000'x, 'ff00'x, 'ff00ff'x, 'ffff'x]
ncol = n_elements(colors)
for i=0, nmodel-1 do begin
    oplot, x, m[i,*], col=colors[i mod ncol], thick=3
endfor

end
