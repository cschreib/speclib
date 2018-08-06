masks = ['cosmos', 'uds_k1', 'uds_k2', 'aegis']
getk  = [1,        0,        1,        1]
geth  = [1,        1,        1,        0]

; suffix = ''

for m=0, n_elements(masks)-1 do begin
    suffix = '_s3'
    ; if masks[m] eq 'uds_k1' then suffix = '_s9'

    master_file = masks[m]+'/reduced/telluric/template_correction'+suffix+'.fits'
    fmaster = mrdfits(master_file, 1, /silent)
    lmaster = get_lambda(master_file)

    print, masks[m]
    for b=0, 1 do begin
        if b eq 0 and geth[m] then begin
            xr = [1.42,1.85]
            xxr = [1.5,1.6]
            band = 'H'
        endif else if b eq 1 and getk[m] then begin
            xr = [1.9,2.45]
            xxr = [2.1,2.3]
            band = 'K'
        endif else continue

        print, masks[m]+'-'+band
        plot, lmaster, fmaster, xr=xr, charsize=2, yr=[-1,5], /xs, /ys

        oplotline, 0, 0, col='ff'x
        stop

        oplot, lmaster, fmaster, col='ff00'x
        files = file_search(masks[m]+'/reduced/telluric/*_gauss_correction'+suffix+'_fit.fits')
        for f=0, n_elements(files)-1 do begin
            r = mrdfits(files[f], 1, /silent)
            idl = where(lmaster ge xxr[0] and lmaster le xxr[1])
            oplot, lmaster, r/median(abs(r[idl]))
        endfor
        stop
    endfor
endfor

end
