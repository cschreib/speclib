spawn, './get.sh'

r = mrdfits('best_fit_slit_star.fits',1,/silent)
plot, r.lam, r.flx, charsize=2, xr=[0.75,0.8], /xs
oplot, r.blam, r.bflx, psym=symcat(16), col='ff'x
oplot, u.lam, u.flux*!y.crange[1], col='ffff'x

; minm = 8
; sm = 8
; mask = smooth(float(u.flux lt 0.95), sm) lt 1.0/float(2*sm)

; iin = 0
; in = 0
; for i=0, n_elements(mask)-1 do begin
;     if ~in and mask[i] then begin
;         iin = i
;         in = 1
;     endif
;     if in and ~mask[i] then begin
;         in = 0
;         if i-iin lt minm then begin
;             mask[iin:i] = 0
;         endif; else begin
;         ;     if (u.lam)[i] gt 0.965 and (u.lam)[i] lt 0.985 then begin
;         ;         print, (u.lam)[iin], (u.lam)[i-1]
;         ;     endif
;         ; endelse
;     endif
; endfor

; curve_fill, u.lam, mask*0.5*!y.crange[1], col='ffff'x

readcol, 'telluric_regions.dat', l1, l2

for i=0, n_elements(l1)-1 do begin
    oplot, [l1[i], l2[i]], !y.crange[1]*[1,1], thick=5, col='ff'x*(i mod 2 eq 0) + 'ff00'x*(i mod 2 eq 1), /noclip
endfor

plot, r.lam, r.baseline/r.flx

end
