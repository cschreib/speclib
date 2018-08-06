ob = 'A'
arm = 'NIR'
project_dir = '/home/cschreib/reductions/xshooter/7329/'
speclib_dir = '/home/cschreib/code/speclib/xshooter/'

r = mrdfits(project_dir+'reduced/'+ob+'/telluric_'+arm+'_intermediate.fits',1,/silent)
plot, r.lam, r.flx, charsize=2, xr=[0.75,0.8], /xs
oplot, r.blam, r.bflx, psym=symcat(16), col='ff'x
oplot, u.lam, u.flux*!y.crange[1], col='ffff'x

readcol, speclib_dir+'telluric_regions.dat', l1, l2

for i=0, n_elements(l1)-1 do begin
    oplot, [l1[i], l2[i]], !y.crange[1]*[1,1], thick=5, $
        col='ff'x*(i mod 2 eq 0) + 'ff00'x*(i mod 2 eq 1), /noclip
endfor

plot, r.lam, r.baseline/r.flx

end
