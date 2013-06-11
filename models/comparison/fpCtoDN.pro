PRO fpCtoDN

gain_bands = [1.71, 3.85, 4.73, 4.94, 4.62]
band = ['u','g','r','i','z']

for n=0, 4 do begin
b = band[n]
framename= b+'.fits'

;; 1. read in the FITS image from HDU0; the resulting image will be
;;    sky-subtracted as well as calibrated in nanomaggies/pixel
img= mrdfits(framename,0,hdr, /dscale)
nrowc= (size(img,/dim))[1]

;; 2. read in sky, and interpolate to full image size; this returns a
;;    sky image the same size as the frame image, in units of counts
sky= mrdfits(framename,2)
simg= interpolate(sky.allsky, sky.xinterp, sky.yinterp, /grid)

;; 3. read in calibration, and expand to full image size; this returns
;;    a calibration image the same size as the frame image, in units of
;;    nanomaggies per count
calib= mrdfits(framename,1)
cimg= calib#replicate(1.,nrowc)

;; In counts:
dn= img/cimg ;+simg

;; Get gain
for m=0, 4 do begin
   if (b eq band[m]) then begin
      gain = gain_bands[m]
      print, m, band[m], gain
   endif
endfor

;; In photo-electrons:
nelec= dn*gain

;; write out
outname = b+'-e.fits'
mwrfits, nelec, outname

endfor
end
