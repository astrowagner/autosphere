pro contrastcurve,name

;Kevin Wagner - Steward Observatory

;A contrast curve generator for the Scorpion survey

;USER INPUTS
separations=[0.1,0.125,0.15,0.2,0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0,3.0,5.0];, 3.0, 5.0]
;separations=[2.0, 3.0, 5.0]
npa=10.
pas=findgen(npa)*360./npa
threshold=5. 
pxscale=0.01225
contrast_guess=1e-4;0.0001 ;guess at the starting point for the contrast curve


;-------------------- END USER INPUT -------------------------------

;save start time
contstarttime=systime(/JULIAN)

;set up arrays
results=separations & results[*]=0.
contrasts=separations & contrasts[*]=contrast_guess
curves=[]

;start a loop over position angle
for m=0,n_elements(pas)-1 do begin

PA=separations & PA[*]=pas[m]*!DTOR ;set all same PA for now

;start a loop over separation
for n=0,n_elements(separations)-1 do begin
	print, 'Working on test separation',n,'/',n_elements(separations)-1
	result=99E99 ;start significance off at infinity
	runs=0
	while result gt 1.05*threshold or result lt 0.95*threshold do begin
		rho=separations[n] & theta=PA[n] & contrast=contrasts[n]
		if n ge 1 and finite(contrast) eq 0 then contrast=contrasts[n-1]
		print, 'Testing rho, theta, contrast = ',rho,theta/!DTOR,contrast
		runs=runs+1

		;set boundaries to remove the line (NEAR SPECIFIC)
		;if rho lt 0.6 then begin 
		;	aa=5. & bb=5. & ab=5. & ba=5. 
		;	if pas[m] lt 200. then ab=10.
		;endif else begin 
		;	aa=10. & bb=10.  & ab=10. & ba=10. 
		;endelse
		;aa=0. & bb=0.  & ab=0. & ba=0. ;fix to zero


		;perform reduction
		;name='S37'
		;reduce_near,rho=rho,theta=theta,contrast=contrast,/block_burn,aa=aa,bb=bb,ab=ab,ba=ba]

		contrasti=contrast ;store original for later

		if rho lt 0.15 then contrast=contrast/0.9 ;correct for coro transmission
		
		autosphere,name,/ird,/inj,/skip,/klip,rho=rho,theta=theta,contrast=contrast		


		contrast=contrasti ;restore
		
		;read image
		image=readfits(strcompress('/Users/kevinwagner/Data/SPHERE/'+name+'/IRDIS/processed/'+name+'_left_klip_inj.fits',/rem))

		;what is rho, theta on the image?
		xs=(size(image))(1) & ys=(size(image))(2)
		xhs=xs/2. & yhs=ys/2. ;half-sizes
		print, xhs, yhs
		print, rho
		xp=xhs + ( (rho/pxscale) * Cos(theta)  ) & yp=yhs + ( (rho/pxscale) * Sin(theta)  )
		print, 'Xp,Yp = ',xp,yp

		;measure signal
		signal=max(image[xp-1:xp+1,yp-1:yp+1])

		nimage=image ;make a copy

		;block known companions


		if name eq 'SA29' then begin
			xc=518. & yc=504.	
			boxehs=12. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif

		if name eq 'SA10' then begin
			xc=548. & yc=470.	
			boxehs=20. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif
	
		if name eq 'SA7' then begin
			xc=365. & yc=499.	
			boxehs=32. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif

		if name eq 'SA20' then begin
			xc=530. & yc=468.	
			boxehs=20. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif
		if name eq 'SA24' then begin
			xc=516. & yc=516.	
			boxehs=25. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif
		if name eq 'SA17' then begin
			xc=483. & yc=657.	
			boxehs=111. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif

		if name eq 'S5' or name eq 'S5_2' then begin
			xc=748. & yc=592.	
			boxehs=180. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan

			xc=506. & yc=514.	
			boxehs=8. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif

		if name eq 'S10' then begin
			xc=587. & yc=614.	
			boxehs=10. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif


		if name eq 'S25_2' then begin
			xc=507. & yc=503.	
			boxehs=12. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif

		if name eq 'S11' then begin
			xc=666. & yc=655.	
			boxehs=150. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif

		if name eq 'SA13' then begin
			xc=478. & yc=454.	
			boxehs=20. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif


		if name eq 'S6' or name eq 'S6_2' or name eq 'S6_3' then begin
			xc=449. & yc=456.	
			boxehs=60. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan

			xc=454. & yc=575.	
			boxehs=15. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif
		if name eq 'S8' or name eq 'S8_2' then begin
			xc=508. & yc=485.	
			boxehs=20. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan

			xc=459. & yc=660.	
			boxehs=23. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan

			xc=751. & yc=573.	
			boxehs=23. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan


			xc=200. & yc=770.	
			boxehs=23. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan

			xc=912. & yc=596.	
			boxehs=23. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan


			xc=660. & yc=586.	
			boxehs=10. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif
		if name eq 'S12_2' then begin
			xc=509. & yc=526.	
			boxehs=10. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif
		if name eq 'S14' or name eq 'S14_2' then begin
			xc=556. & yc=530.	
			boxehs=50. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif

		if name eq 'S20' or name eq 'S20_2' then begin
			xc=635. & yc=588.	
			boxehs=65. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif

		if name eq 'S22' or name eq 'S22_2' or name eq 'S22_3' or name eq 'S22_4' then begin
			xc=686. & yc=319.	
			boxehs=120. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif

		if name eq 'S23_2' then begin
			xc=767. & yc=558.	
			boxehs=178. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif


		if name eq 'S23' or name eq 'S23_3' then begin
			xc=258. & yc=469.	
			boxehs=178. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif

		if name eq 'S25' or name eq 'S25_2' then begin
			xc=274.5 & yc=863.	
			boxehs=200. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:1023]=!values.f_nan
		endif

		if name eq 'S31' then begin
			xc=110. & yc=275.	
			boxehs=150. ;source exclusion half-size
			nimage[0:xc+boxehs-1,0:yc+boxehs-1]=!values.f_nan
		endif

		if name eq 'S30' then begin
			xc=525. & yc=492.	
			boxehs=10. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif
		if name eq 'S31' then begin
			xc=534. & yc=500.	
			boxehs=20. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif

		if name eq 'S32' or name eq 'S32_2' then begin
			xc=406. & yc=459.	
			boxehs=25. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif

		if name eq 'S33' then begin
			xc=518. & yc=528.	
			boxehs=15. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan

			xc=491. & yc=554.	
			boxehs=15. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif

		if name eq 'S37_3' then begin
			xc=462. & yc=535.	
			boxehs=25. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan

			xc=375. & yc=465.	
			boxehs=25. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif
		
		if name eq 'S38' then begin
			xc=677. & yc=512.	
			boxehs=118. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif

		if name eq 'S40' then begin
			xc=505. & yc=507.	
			boxehs=10. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif
		if name eq 'S41' then begin
			xc=510. & yc=508.	
			boxehs=10. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
			xc=347. & yc=521.	
			boxehs=15. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
			xc=667. & yc=707.	
			boxehs=15. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif
		if name eq 'S43' or name eq 'S43_2'   then begin
			xc=724. & yc=293.	
			boxehs=52. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
			xc=471. & yc=921.	
			boxehs=61. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan

			xc=635. & yc=119.	
			boxehs=61. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan

			xc=145. & yc=438.	
			boxehs=100. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
			xc=500. & yc=919.	
			boxehs=200. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:min([1023,yc+boxehs-1])]=!values.f_nan
		endif
		
		if name eq 'S45' or name eq 'S45_2' then begin
			xc=571. & yc=615.	
			boxehs=100. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
			xc=304. & yc=464.	
			boxehs=100. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif


		if name eq 'S44' or name eq 'S44_2' then begin
			xc=462. & yc=640.	
			boxehs=100. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif
		if name eq 'S46' then begin
			xc=613. & yc=295.	
			boxehs=132. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif
		if name eq 'S48' or name eq 'S48_2'  or name eq 'S48_3' then begin
			xc=611. & yc=549.	
			boxehs=26. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif
		if name eq 'S56' then begin
			xc=495. & yc=498.	
			boxehs=20. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif
		if name eq 'S50' then begin
			xc=527. & yc=497.	
			boxehs=10. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif
		if name eq 'S52' then begin
			xc=513. & yc=361.	
			boxehs=131. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan

			xc=699. & yc=270.	
			boxehs=29. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif
		if name eq 'S53' then begin
			xc=682. & yc=557.	
			boxehs=90. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif
		if name eq 'S54' then begin
			xc=674. & yc=554.	
			boxehs=15. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif
		if name eq 'S57' then begin
			xc=490. & yc=520.	
			boxehs=10. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif

		if name eq 'S61' then begin
			xc=493. & yc=522.	
			boxehs=25. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
			xc=887. & yc=456.	
			boxehs=102. ;source exclusion half-size
			nimage[xc-boxehs:1023,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif
		if name eq 'S63' or name eq 'S63_2' then begin
			xc=405. & yc=726.	
			boxehs=90. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan

			xc=577. & yc=398.	
			boxehs=20. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan

		xc=133. & yc=366.	
			boxehs=20. ;source exclusion half-size
			nimage[xc-boxehs:xc+boxehs-1,yc-boxehs:yc+boxehs-1]=!values.f_nan
		endif

		;measure noise

		;writefits,'~/Desktop/noise_test.fits',nimage

		boxehs=7. ;source exclusion half-size
		if name eq 'S37' then boxehs=12
		nimage[xp-boxehs:xp+boxehs-1,yp-boxehs:yp+boxehs-1]=!values.f_nan
		noise=stddev(nimage,/nan)
		result=signal/noise


	fwhm=(2.1E-6)/8.2 *206265 / 0.01225
	print, fwhm

		;measure signal in an aperture
		aper_rad=fwhm/2. ;& skyradi=4. & skyrado=4.
		zimage=image
		zimage[where(finite(zimage) ne 1)]=0.
		aper, zimage, xp,yp, signal,err,sky,skyerr,1,aper_rad,[0,0],[-99E99,99E99],/flux,/exact,/silent,SETSKYVAL = 0;bkg
		

	;copied from find_sources
	;start at source, progress by each pixel (to 2PiR) and calculate flux in each aperture
	ntheta=float( round( ( 2.0*!PI*float(rho/pxscale)/(aper_rad*2.0) )-1 ) ) ;about one angular bin per resolution element

		nzimage=nimage
		nzimage[where(finite(nzimage) ne 1)]=0.

	if ntheta lt 6 then begin stitheta=1 & stendtheta=ntheta-1 & endif else begin & stitheta=2 & stendtheta=ntheta-2 & endelse
	for itheta=stitheta,stendtheta do begin
		ttheta=theta+(float(itheta)/float(ntheta))*2.*!PI;+(!PI/31.)
		xtest=xhs+(float(rho/pxscale)*cos(ttheta)) & ytest=yhs+(float(rho/pxscale)*sin(ttheta))
		;print, size(frame)
		aper, nzimage, xtest, ytest, flux,fluxerr,sky,skyerr,1.75,aper_rad,[0,0],[-99E99,99E99],/silent,/flux,SETSKYVAL=0,/exact

		if itheta eq stitheta then fluxarr=flux else fluxarr=[fluxarr,flux]
		if itheta eq stitheta then xarr=xtest else xarr=[xarr,xtest]
		if itheta eq stitheta then yarr=ytest else yarr=[yarr,ytest]
	endfor ;;itheta for
		noise=stddev(fluxarr,/nan)

		signal=signal-mean(fluxarr,/nan)
		noise=noise*sqrt(1.+(1./n_elements(fluxarr))) ;small sample statistics correction: Mawet+2014
			result=signal/noise

		if result lt 0 then result=0.5 ;if no signal, double next the contrast guess next time
		print, 'Source recovered with significance = ',result


			
		print, 'Source recovered with significance = ',result

;		if runs gt 5 then begin & result=threshold & contrasts[n] = !values.f_nan & endif
		if runs gt 7 and abs(1.-(result/threshold) ) gt 0.1 then begin & result=threshold & contrasts[n] = !values.f_nan & endif
		if runs gt 7 and abs(1.-(result/threshold) ) le 0.1 then result=threshold 
		
		;if result is not within 5% of the threshold, revise the injected contrast
		;if rho gt 0.4 then mincr=0.0005
		if rho lt 0.4 then mincr=0.001
		if rho lt 0.3 then mincr=0.01
		if rho lt 0.2 then mincr=0.1
		if rho lt 0.15 then mincr=0.25
		if rho lt 0.125 then mincr=1.
		if abs(1.-(result/threshold) ) gt 0.05 then contrasts[n:n_elements(separations)-1]=min([mincr,contrast/((result/threshold))])	
		;if result lt 0.95*threshold then contrasts[n:n_elements(separations)-1]=contrast/((result/threshold)) 

		
		results[n]=result ;store results for diagnostics
		print, contrasts
		print, results
		print, separations
		;hak

		;start plotting things
		cgps_open,strcompress('/Users/kevinwagner/Data/SPHERE/'+name+'/IRDIS/processed/'+name+'_contrast_curve.ps',/rem),xsize=7.25, ysize=5, /inches, portrait=1

			;open a new plot if it is the first position angle
			if m eq 0 then cgplot,separations,contrasts,xtitle='Separation (arcsec)',ytitle='Contrast',/ylog,yrange=[1E-7,1E-1],xrange=[0.08,6],title='result = '+string(sigfig(result,3))+' '+'Floor = '+string(sigfig(min(contrasts,/nan),3)),charsize=1,/xlog
		
			;open a new plot for the first PA and then oplot for the remaining
			if m gt 0 then begin
				cgplot,separations,curves[*,0],xtitle='Separation (arcsec)',ytitle='Contrast',/ylog,yrange=[1E-7,1E-1],xrange=[0.08,6],title='result = '+string(sigfig(result,3))+' '+'Floor = '+string(sigfig(min(contrasts,/nan),3)),charsize=1,color='gray',/xlog
				if m gt 1 then for l=1,m-1 do cgoplot,separations,curves[*,l],color='gray'

				cgoplot, separations, contrasts,color='gray' ;now plot the current contrasts
			endif

			;cgoplot,[0,2],[threshold*0.67E-6,threshold*0.67E-6],linestyle=2
	

		;plot mean contrast curve from the data reductions

		if m gt 1 then begin
			medcurve=mean(curves,dim=2,/nan)
			cgoplot,separations,medcurve,thick=4
		endif




		cgps_close,/pdf
		writefits,strcompress('/Users/kevinwagner/Data/SPHERE/'+name+'/IRDIS/processed/'+name+'_contrastcurve_testimage.fits',/rem),image
		writefits,strcompress('/Users/kevinwagner/Data/SPHERE/'+name+'/IRDIS/processed/'+name+'_contrastcurve_noiseimage.fits',/rem),nimage

		;output images once the threshold has been found
		if result ge 0.95*threshold and result le 1.05*threshold then $
			if n eq 0 and m eq 0 then imgs=image else imgs=[ [[imgs]], [[image]] ]
		if n gt 0 or m gt 0 then writefits,strcompress('/Users/kevinwagner/Data/SPHERE/'+name+'/IRDIS/processed/'+name+'_contrastcurve_images.fits',/rem),imgs
	endwhile
endfor ;n

curves=[ [curves], [contrasts] ]

save,filename=strcompress('/Users/kevinwagner/Data/SPHERE/'+name+'/IRDIS/processed/'+name+'_contrast_output.sav',/rem),curves,pas,separations,threshold

endfor ;m



print, 'Completed contrast curve generation in ',(systime(/JULIAN)-contstarttime)*86400./60.,' minutes.'

end
