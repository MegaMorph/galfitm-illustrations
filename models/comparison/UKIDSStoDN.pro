PRO UKIDSStoDN

band = ['Y', 'J', 'H', 'K']

;; gain = 4.5
;; exptime = [45.732, 48.157, 45.360, 45.424]  ; seconds
;; sky = [466.873, 580.646, 6822.956, 5991.750]  ; counts/pixel

for n=0, 3 do begin
b = band[n]
framename= b+'.fits'

img = mrdfits(framename,1,hdr, /dscale)
sky = sxpar(hdr, 'SKYLEVEL')
gain = sxpar(hdr, 'GAIN')
nelec= (img - sky)*gain

;; write out
outname = b+'-e.fits'
mwrfits, nelec, outname

endfor
end
