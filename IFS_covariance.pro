pro IFS_covariance


;following the methods in Greco & Brandt (2016)
;psi is the correlation matrix, which is a NxN matrix where N = the number of spectral channels
;for SPHERE-IFS this is a 39x39 matrix

;Kevin Wagner - Steward Observatory

psi=fltarr(39,39)

image=readfits('/Users/kevin/Desktop/Sco1b-fit/HD131399_ifs_klip_med_cube.fits',hdr)
;image=readfits('/Users/kevin/Desktop/Sco1b-fit/data_cube_psf.fits')
image=readfits('/Users/kevin/Desktop/Sco1b-fit/HD131399_ifs_sdiklip_xyl.fits',hdr)
wavelengths=readfits('/Users/kevin/Data/VLT/SPHERE/Scorpion/HD131399-working/HD131399-merged/IFS_YJH/products/data_wavelength.fits', wavelengthhd)


;image[where(finite(image) ne 1)]=0.

;image=sqrt(image*image)
;image=abs(image)

writefits,'~/Desktop/image.fits',image
;Measuring PSI_ij

for i=0, 38 do begin

for j=0, 38 do begin

print, 'Working on i=',i,' j=',j

rho=0.83;	radius in arcsec

wedge=0

lambdai=wavelengths[i]*1.0E-6
lambdaj=wavelengths[j]*1.0E-6

pxscale=0.00746 ;arcsec per pixel

apers=0

	inner_rad=(rho/0.00746) - 2.
	outer_rad=(rho/0.00746) + 2.
	
FWHMi=lambdai / 8.2 * 206265. / pxscale
aperturei=0.75*FWHMi ;aperture radius, e.g. 0.75xFWHM = 1.5FWHM diameter
napersi=2.*!PI*(rho/pxscale) / (aperturei*2.)	;how many apertures can fit into a circle of this rho?
dthetai=2.*!PI / napersi

FWHMj=lambdaj / 8.2 * 206265. / pxscale
aperturej=0.75*FWHMj ;aperture radius, e.g. 0.75xFWHM = 1.5FWHM diameter
napersj=2.*!PI*(rho/pxscale) / (aperturej*2.)	;how many apertures can fit into a circle of this rho?
dthetaj=2.*!PI / napersj


inoise=[] ;makes an array to start counting noise
jnoise=[]
ijnoise=[]
xt=163. & yt=24.

for xx=0.,269. do begin
	for yy=0.,269. do begin
		if sqrt ( (xx-xt)^2. + (yy-yt)^2. ) lt 5. then image[xx,yy,*]=!values.f_nan
	endfor
endfor

xc=135. & yc=135.

print, 'working with',napersi,napersj




theta=atan((yt-yc)/(xt-xc))

ntot=fix(napersi)-1

if apers then begin
for nn=1, ntot do begin;fix(napersi)-1 do begin	;loops through the noise apertures
	xn=xc+rho*Cos(theta+ (nn*dthetai) )
	yn=yc+rho*Sin(theta+ (nn*dthetai) )
	
	
	
	;if wedge then xn=xc+rho*(nn/fix(napersi)-1)*Cos(theta);+ (nn*dthetai) )
	;if wedge then yn=yc+rho*(nn/fix(napersi)-1)*Sin(theta);+ (nn*dthetai) )
	
	aper, image[*,*,i], xn, yn, noisei, noiseerr, sky, skyerr, 1.75*32., aperturei,[aperturei,aperturei+2],[-99E99,99E99], /exact, /flux, /silent, setskyval=0;,SETSKYVAL=0.
	inoise=[inoise,noisei]
endfor
endif else begin
	newimg=image
	rad=4.
	inoise=[]
	for xx=0,269 do begin
		for yy=0, 269 do begin
			if sqrt ( (xx-135.)^2. + (yy-135.)^2. ) gt inner_rad or $
				sqrt ( (xx-135.)^2. + (yy-135.)^2. ) lt outer_rad then inoise=[inoise,newimg[xx,yy,i]]
		endfor
	endfor
	
	;newimg[where(newimg eq 0.)]=!values.f_nan
		;writefits,'~/Desktop/img.fits',newimg

	;inoise=newimg[*,*,i]
	;print, inoise
endelse

;inoise=inoise[where(inoise gt 0)]

inoisesq_ex=mean(inoise*inoise,/nan)
inoise_ex=mean(inoise,/nan)

if apers then begin
for nn=1, ntot do begin;fix(napersj)-1 do begin	;loops through the noise apertures
	xn=xc+rho*Cos(theta+ (nn*dthetaj) )
	yn=yc+rho*Sin(theta+ (nn*dthetaj) )
	
	
	;if wedge then xn=xc+rho*(nn/fix(napersi)-1)*Cos(theta);+ (nn*dthetai) )
	;if wedge then yn=yc+rho*(nn/fix(napersi)-1)*Sin(theta);+ (nn*dthetai) )

	aper, image[*,*,j], xn, yn, noisej, noiseerr, sky, skyerr, 1.75*32., aperturej,[aperturej,aperturej+2],[-99E99,99E99], /exact, /flux, /silent, setskyval=0;,SETSKYVAL=0.
	jnoise=[jnoise,noisej]
	
	;ijnoise=[ijnoise,noisei*noisej]
endfor
endif else begin
	newimg=image
	rad=4.
	jnoise=[]
	for xx=0,269 do begin
		for yy=0, 269 do begin
			if sqrt ( (xx-135.)^2. + (yy-135.)^2. ) gt inner_rad or $
				sqrt ( (xx-135.)^2. + (yy-135.)^2. ) lt outer_rad then jnoise=[jnoise,newimg[xx,yy,j]]
		endfor
	endfor
endelse

;jnoise=jnoise[where(jnoise gt 0)]
;expectation values
jnoisesq_ex=mean(jnoise*jnoise,/nan)
jnoise_ex=mean(jnoise,/nan)

ijnoise_ex=mean(inoise*jnoise,/nan)

print, ijnoise_ex, inoisesq_ex, jnoisesq_ex

psi[i,j]= ijnoise_ex / sqrt(inoisesq_ex*jnoisesq_ex)

print, 'Psi_ij=',psi[i,j]

print, ijnoise_ex
print, sqrt(inoisesq_ex*jnoisesq_ex)

endfor

plot, psi[i,*]

endfor

print, psi

print, size(psi)

save,filename='~/Desktop/covariances.sav',psi



end