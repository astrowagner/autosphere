pro autosphere, object_names, skip=skip, fresh=fresh, single=single, no_interupt=no_interupt,suffix=suffix,klip=klip,inject=inject,ifs_only=ifs_only,irdis_only=irdis_only, second_cen=second_cen,no_flat_ird=no_flat_ird,no_dark_ird=no_dark_ird,rho=rho,theta=theta,contrast=contrast;,waffle=waffle

;--------[   Keywords :   ]-----------
; SKIP: binary keyword to skip to the post-processing steps (irdis and IFS)
; SINGLE: binary keyword to enable a HAK after each object that is reducted. Useful if you want to
;         stop and check results after each run before going on to the next.
; NO_INTERRUPT: bianry keyword that if enabled will not interupt if calibration files have different dates than science
; SUFFIX: will add a suffix to the files (expects a string input)
; KLIP: binary keyword that will run the klip algorith on IRDIS data
; INEJCT: binary keyword that will inject fake planets
; IFS_ONLY: binary keyword to only process IFS
; IRDIS_ONLY: binary keyword to only process IRDIS
; SECOND_CEN: will use the second center frame in IFS. For IRDIS, this functionality is set in reduce_irdis.
; NO_FLAT_IRD: will skip the flat fielding for irdis
; NO_DARK_IRD: will skip dark correction for irdis
;
;
;--------[  Written by K. Wagner (UA)   ]-----------
;
;--------[   Description:   ]-----------
; This script reads in a folder of fits files downloaded from the ESO archive (including science and
; calibrations from IRDIS and or IFS), or simply just a download script, and then does all the rest of 
; the steps for IRDIS (DBI and CI modes) and IFS (YJ and YJH modes). 
;
;
;--------[   Instructions:   ]-----------
;
; The user simply needs to set the reduction_path (defined just below) and then use as
; IDL> autosphere, 'object-name',<keywords>


reduction_path='/Users/kevinwagner/Data/SPHERE/'

;-----------------------------------------------------------------------------------------


for iii=0, n_elements(object_names)-1 do begin


cgcleanup	;clean plot windows from last time, if they exist

object_name=object_names[iii]
print, '-------------[ Beginning SPHERE reduction on ',object_name,' ]--------------'

speaking=1
if not keyword_set(skip) then skip=0
if not keyword_set(irdis_only) then irdis_only=0
if not keyword_set(klip) then klip=0
if not keyword_set(ifs_only) then ifs_only=0
if keyword_set(irdis_only) then process_ifs=0 else process_ifs=1

;following switches are for IRDIS data


;fresh=1	;will delete the processed folder when beginning...
		;made the above into a keyword 2017-11-14
if keyword_set(fresh) then fresh=1 else fresh=0

if keyword_set(single) then single=1 else single=0	
;runs HAK between runs (must hit any key to continue, or ESC to end reduction train). 

;object_name='S47'

;will go straight to the final post-processing
	if keyword_set(skip) then skip_start=1 else skip_start=0

	if skip_start then fresh=0 ;if skipping start then fresh makes no sense





input_directory=strcompress(reduction_path+object_name+'/',/rem)

;setting the path to the raw data (extracted .Z files)
download_subdir=strcompress(input_directory+'data_with_raw_calibs/',/rem)

dark_comb_imtype='median' ;'median' or 'mean'

if keyword_set(no_dark_ird) then sub_dark=0	else sub_dark=1 ;switch to turn dark subtraction on/off
	
	combine_all_darks=0	;will scale all darks by DIT and median combine into master dark
				;poor performance compared to the auto dark selection

if keyword_set(no_flat_ird) then div_flat=0 else div_flat=1	;switch to turn flat division on/off


remove_bads=1 ;very, very slow but worth it
	badmap_loc=reduction_path+'calib/SPHERE_badpix-map.fits'
	find_bads=0	;look for additional bad pixels?
	bad_thresh=5	;DN above which is to be considered a bad pixel in dark frame
	neg_bad_thresh=-5


reduce_irdis_sizes=1 	;re-zip raw data and delete intermediate products to save space?
	reduce_ifs_sizes=reduce_irdis_sizes	
	if ifs_only then reduce_irdis_sizes=0
	if keyword_set(irdis_only) then reduce_ifs_sizes=0
post_process_irdis=1	;run the reduce_irdis script?
	auto_cen=1	;automatically center IRDIS data
	if ifs_only then post_process_irdis=0


unzip=1	;will first uncompress the raw data, may or may not speed up overall script runtime. Fastest is to unzip in UNIX and delete the *.Z files. 




;-----------------------------------------------------------------------------------------

if skip_start then starttime=systime(/JULIAN)
if not skip_start then begin

	;has the script been run before? If not, is there downloaded data? If so, download it.
	files=file_search(input_directory,'*',count=cnt_total)
	dl_file=file_search(input_directory,'*.sh',count=cnt_dl)


	if cnt_dl eq 1 and cnt_total eq 1 then begin

		dlstarttime=systime(/JULIAN)

		cd,'/'
;

		spawn, 'chmod +x '+dl_file
		
		result = STRPOS( dl_file, '/dow' , /REVERSE_SEARCH )
		length=strlen(dl_file)
		dl_file=strmid(dl_file,length-1-result-1,/rev)
		cd,input_directory
			spawn,  './'+string(dl_file) 

		print, 'Completed download from ESO archive on ',object_name,' in ',(systime(/JULIAN)-dlstarttime)*86400./60.,' minutes.'

	endif

	;start clock on program run time after the download step is finished.
	starttime=systime(/JULIAN)



print, 'Step 0: Cleaning downloaded files...'
	files=file_search(download_subdir,'M.SPHERE*',count=cnt)
	if cnt gt 0 then file_delete,files
	files=file_search(download_subdir,'*.txt',count=cnt)
	if cnt gt 0 then file_delete,files
	files=file_search(download_subdir,'*.xml',count=cnt)
	if cnt gt 0 then file_delete,files
	
	if unzip then begin
		files=file_search(download_subdir,'*.fits.Z',count=cnt)
		for ii=0, cnt-1 do begin 
			print, 'Unzipping file ',ii+1,'/',cnt
			a=readfits(files[ii],hd) 
			writefits, strcompress(files[ii]+'.fits',/rem),a,hd 
		endfor;file_gunzip, files[ii]
		if cnt gt 0 then file_delete, files
	endif

print, 'Step 1: Making subdirectories... '
	if fresh then file_delete,input_directory+'IRDIS/calib',$
				input_directory+'IFS/processed',$
				input_directory+'IFS/products',$
				input_directory+'IFS/calib',$
				input_directory+'IFS/interim',$
				input_directory+'IRDIS/processed',$
				input_directory+'IRDIS/flux',$
				input_directory+'IRDIS/center',$
				/quiet,/recursive

	;clearing old reduction interim files, if present
	testfile=file_search(input_directory+'IRDIS/products/',count=testcnt)
	if testcnt ge 1 then file_delete,input_directory+'IRDIS/products/',/recursive

		file_mkdir,	input_directory+'IRDIS/raw',$
				input_directory+'IRDIS/calib',$
				input_directory+'IRDIS/flux',$
				input_directory+'IRDIS/center',$
				input_directory+'IRDIS/products',$
				input_directory+'IRDIS/processed',$
				input_directory+'IFS/raw',$			
				input_directory+'IFS/extra',$
				input_directory+'IFS/products',$
				input_directory+'IFS/processed'


	if not ifs_only then begin
print, 'Looking for previously zipped data...'
	files=file_search(input_directory+'IRDIS/','raw.zip',count=cnt)
	if cnt gt 0 then file_unzip,files
	endif

print, '	Done.'
	
;-----------------------------------------------------------------------------------------
	
	
print, 'Step 2: Searching for .fits files... Done.'
	files=file_search(download_subdir,'*.fits.Z',count=cnt)	;first search for compressed files
	if cnt le 0 then files=file_search(download_subdir,'*.fits',count=cnt)

	
print, '	Found ', cnt, ' fits files.'
if cnt gt 0 then begin

;-----------------------------------------------------------------------------------------


print, 'Step 3: Looping through fits headers and sorting files...'
for ii=0, cnt-1 do begin
	frame=readfits(files[ii],hdr,/silent)
	
	
	detector=''
	imtype=''
	date=''
	exptime=''
	ndit=''
	dbfilt=''
	ndfilt=''
	
	;extracting info...
	;want to know what detector, imtype, date, DIT, ND filter, DB filter
	detector=string(esopar(hdr,'HIERARCH ESO DET NAME '))
	imtype=string(esopar(hdr,'HIERARCH ESO DPR TYPE '))
	date=fxpar(hdr,'DATE')
	object=string(fxpar(hdr,'OBJECT'))
	date_full=date
	date=strmid(date,0,10)	;only keep YYYY-MM-DD
	exptime=fxpar(hdr,'EXPTIME')
	ndit=esopar(hdr,'HIERARCH ESO DET NDIT ')
	dbfilt=esopar(hdr,'ESO INS1 OPTI2 NAME ')
	ndfilt = esopar(hdr,'ESO INS4 COMB IND')
	if isa(ndfilt, 'string') eq 0 then ndfilt = esopar(hdr,'ESO INS4 FILT2 NAME')  

	if ii eq 0 then begin
		detectors=detector
		dates=date
		full_dates=date_full
		imtypes=imtype
		exptimes=exptime
		ndfilts=ndfilt
		dbfilts=dbfilt
		ndits=ndit
		objects=object
	endif else begin
		detectors=[detectors,detector]
		dates=[dates,date]
		imtypes=[imtypes,imtype]
		exptimes=[exptimes,exptime]
		ndfilts=[ndfilts,ndfilt]
		dbfilts=[dbfilts,dbfilt]
		ndits=[ndits,ndit]
		full_dates=[full_dates,date_full]
		objects=[objects,object]
	endelse

	if speaking then begin
		print, strcompress('	Frame '+string(ii)+':')
		print, '	Date       = ',date
		print, '	Object     = ',object
		print, '	Detector   = ',detector
		print, '	Frame imtype = ',imtype
		print, '	Exptime    = ',exptime
		print, '	NDITs    = ',ndit
		print, '	ND filter  = ',ndfilt
		print, '	DB Filter  = ',dbfilt
	endif
endfor
print, '	Done.'

;separating dates into other forms
	yyyys=strmid(dates,0,4)
	mms=strmid(dates,5,2)
	dds=strmid(dates,9,2)
	
;-----------------------------------------------------------------------------------------


	
;separate IRDIS and IFS... 
print, 'Step 4: Separating IFS frames. Done with IFS step. Proceed to AVs scripts.'
	if total(where(strcompress(detectors,/rem) eq 'IFS')) gt 0 then file_move,$
		files(where(strcompress(detectors,/rem) eq 'IFS')),input_directory+'IFS/raw'

;-----------------------------------------------------------------------------------------


print, 'Step 5: Separating IRDIS frames.'
	if total(where(strcompress(detectors,/rem) eq 'IRDIS')) gt 0 then file_move,$
	files(where(strcompress(detectors,/rem) eq 'IRDIS')),input_directory+'IRDIS/raw'

endif else print, 'No files found in download directory... assuming they were already sorted, proceeding to step 6.'



;-----------------------------------------------------------------------------------------

if not ifs_only then begin

print, 'Step 6: Looping through IRDIS fits headers and sorting files...'
	files=file_search(input_directory+'IRDIS/raw','*.fits.Z',count=cnt)	;first find compressed files

	if cnt le 0 then files=file_search(input_directory+'IRDIS/raw','*.fits',count=cnt)

for ii=0, cnt-1 do begin
	frame=readfits(files[ii],hdr,/silent)
	
	;extracting info...
	;want to know what detector, imtype, date, DIT, ND filter, DB filter
	detector=esopar(hdr,'HIERARCH ESO DET NAME')
	imtype=String(esopar(hdr,'HIERARCH ESO DPR TYPE'))
	date=fxpar(hdr,'DATE')
	date_full=date
	date=strmid(date,0,10)	;only keep YYYY-MM-DD
	exptime=float(fxpar(hdr,'EXPTIME'))
	ndit=esopar(hdr,'HIERARCH ESO DET NDIT')
	dbfilt=esopar(hdr,'ESO INS1 OPTI2 NAME')
	object=fxpar(hdr,'OBJECT')
	ndfilt = esopar(hdr,'ESO INS4 COMB IND')
	if isa(ndfilt, 'string') eq 0 then ndfilt = esopar(hdr,'ESO INS4 FILT2 NAME')  

	if ii eq 0 then begin
		detectors=detector
		dates=date
		imtypes=imtype
		exptimes=exptime
		ndfilts=ndfilt
		dbfilts=dbfilt
		ndits=ndit
		full_dates=date_full
		objects=object
	endif else begin
		detectors=[detectors,detector]
		dates=[dates,date]
		imtypes=[imtypes,imtype]
		exptimes=[exptimes,exptime]
		ndfilts=[ndfilts,ndfilt]
		dbfilts=[dbfilts,dbfilt]
		ndits=[ndits,ndit]
		full_dates=[full_dates,date_full]
		objects=[objects,object]
	endelse

	if speaking then begin
		print, strcompress('	Frame '+string(ii)+':')
		print, '	Date       = ',date
		print, '	Object     = ',object
		print, '	Detector   = ',detector
		print, '	Frame imtype = ',imtype
		print, '	Exptime    = ',exptime
		print, '	NDITs    = ',ndit
		print, '	ND filter  = ',ndfilt
		print, '	DB Filter  = ',dbfilt
	endif
endfor
print, '	Done.'	

dates_full=full_dates

;-----------------------------------------------------------------------------------------


print, 'Step 7: Collecting object files...'


	obj_inds=where(strcompress(imtypes,/rem) eq 'OBJECT' and strcompress(objects,/rem) ne 'OBJECT')
	;if keyword_set(waffle) then obj_inds=where(strcompress(imtypes,/rem) eq 'CENTER'); and strcompress(objects,/rem) ne 'OBJECT')

	;if no object files, assume it is all center files
	if obj_inds[0] eq -1 then obj_inds=where(strcompress(imtypes,/rem) eq 'OBJECT,CENTER'); and strcompress(objects,/rem) ne 'OBJECT')
	
	if obj_inds[0] eq -1 then begin
		print, 'Warning! No object files found!'
		return
	endif
	print, '	Found ',n_elements(obj_inds),' object frames.'
	;check if the dates are the same
	obj_dates=dates[obj_inds]
		obj_dates_full=full_dates[obj_inds]

		date_check=total( where( obj_dates ne obj_dates[0] ) )

	print, '	Date of first exposure: ',obj_dates[0]
	print, '	Date of last exposure: ',obj_dates[n_elements(obj_inds)-1]

	if date_check gt 0 then begin
		print, '	Object dates differ!!!'
		;print, obj_dates
		print, '	If dates are only different because the observation went through midnight then proceed, otherwise manually separate the files and re-run this script.'
		;hak
	endif
	
		obj_files=files[obj_inds]

	
	obj_exptimes=exptimes[obj_inds]
		exp_check=total( where( obj_exptimes ne obj_exptimes[0] ) )
		
	if exp_check gt 0 then begin
		print, 'Object exposure times differ. Check input files'
		;print, obj_dates
		
		print, obj_exptimes
		
		print, 'Continuing with most prevalent exposure:'
	
		distfreq = Histogram(obj_exptimes, MIN=Min(obj_exptimes))
  		 maxfreq= Max(distfreq)
   		mode = (Where(distfreq EQ maxfreq) + Min(obj_exptimes))
   	
   		
   		obj_exptimes=obj_exptimes[where(obj_exptimes eq mode[0])]
   		obj_dates=obj_dates[where(obj_exptimes eq mode[0])]
   		obj_files=obj_files[where(obj_exptimes eq mode[0])]
		obj_dates_full=obj_dates_full[where(obj_exptimes eq mode[0])]
		;hak
		
		obj_exptime=mode[0]
		
	endif else obj_exptime=obj_exptimes[0]
	
	print, '	Object exposure times = ',obj_exptime
	
	
;-----------------------------------------------------------------------------------------


print, 'Step 8: Collecting star center files...'

	cen_inds=where(strcompress(imtypes,/rem) eq 'OBJECT,CENTER' )
	if cen_inds[0] eq -1 then begin
		print, 'Warning! No star center files found!'
		;return

		;in this case use the object frames themselves
		;typically their absence is due to a non-coronagraphic sequence
		cen_files=obj_files[0]
		cen_inds=obj_inds[0]

		no_cen_files_check=1
	endif else no_cen_files_check=0
	print, '	Found ',n_elements(cen_inds),' center frames.'
	;check if the dates are the same
	cen_dates=dates[cen_inds]
		cen_dates_full=full_dates[cen_inds]

		date_check=total( where( cen_dates ne cen_dates[0] ) )

	print, '	Date of first exposure: ',cen_dates[0]
	print, '	Date of last exposure: ',cen_dates[n_elements(cen_inds)-1]

	if date_check gt 0 then begin
		print, '	Center dates differ!!!'
		;print, cen_dates
		print, '	If dates are only different because the observation went through midnight then proceed, otherwise manually separate the files and re-run this script.'
		if not keyword_set(no_interupt) then hak
	endif
	
	cen_exptimes=exptimes[cen_inds]
		exp_check=total( where( cen_exptimes ne cen_exptimes[0] ) )
		
	cen_files=files[cen_inds]

		
	if exp_check gt 0 then begin
		print, 'Star center exposure times differ. Check input files'
		;print, cen_dates
		
		
		print, cen_exptimes
		
		print, 'Continuing with most prevalent exposure:'
	
		distfreq = Histogram(cen_exptimes, MIN=Min(cen_exptimes))
  		 maxfreq= Max(distfreq)
   		mode =Where(distfreq EQ maxfreq) + Min(cen_exptimes)
   		Print, mode[0]
   		
   		cen_exptimes=cen_exptimes[where(cen_exptimes eq mode[0])]
   		cen_dates=cen_dates[where(cen_exptimes eq mode[0])]
   		cen_files=cen_files[where(cen_exptimes eq mode[0])]
		cen_dates_full=cen_dates_full[where(cen_exptimes eq mode[0])]
		;hak
		
		cen_exptime=mode[0]
		
		
		;hak
	endif else cen_exptime=cen_exptimes[0]
	
	print, '	Star center exposure times = ',cen_exptime
	
		
;-----------------------------------------------------------------------------------------

		
print, 'Step 9: Collecting flux calibration files...'


	flux_inds=where(strcompress(imtypes,/rem) eq 'OBJECT,FLUX' )
	if flux_inds[0] eq -1 then begin
		print, 'Warning! No flux files found!'
		;return


		;in this case use also the object files
		;typically their absence is due to a non-coronagraphic sequence

		flux_inds=obj_inds[0]
		flux_files=obj_files[0]
	endif

print, '	Found ',n_elements(flux_inds),' flux calibration frames.'
	;check if the dates are the same
	flux_dates=dates[flux_inds]
	flux_dates_full=full_dates[flux_inds]
		date_check=total( where( flux_dates ne flux_dates[0] ) )

	print, '	Date of first exposure: ',flux_dates[0]
	print, '	Date of last exposure: ',flux_dates[n_elements(flux_inds)-1]

	if date_check gt 0 then begin
		print, '	Flux calibration dates differ!!!'
		;print, flux_dates
		print, '	If dates are only different because the observation went through midnight then proceed, otherwise manually separate the files and re-run this script.'
		;hak
	endif
	
	flux_exptimes=exptimes[flux_inds]
		exp_check=total( where( flux_exptimes ne flux_exptimes[0] ) )
		
		
	flux_files=files[flux_inds]
	
		
	if exp_check gt 0 then begin
		print, 'Flux calibration exposure times differ. Check input files'
		;print, flux_dates
		print, flux_exptimes
		
		print, 'Continuing with most prevalent exposure:'
	
		distfreq = Histogram(flux_exptimes, MIN=Min(flux_exptimes))
  		 maxfreq= Max(distfreq)
   		mode =Where(distfreq EQ maxfreq) + Min(flux_exptimes)
   		Print, mode[0]
   		
   		flux_exptimes=flux_exptimes[where(flux_exptimes eq mode[0])]
   		flux_dates=flux_dates[where(flux_exptimes eq mode[0])]
   		flux_files=flux_files[where(flux_exptimes eq mode[0])]
		flux_dates_full=flux_dates_full[where(flux_exptimes eq mode[0])]
		;hak
		
		flux_exptime=mode[0]
		
	endif else flux_exptime=flux_exptimes[0]
	
	print, '	Flux calibration exposure times = ',flux_exptime

	
;-----------------------------------------------------------------------------------------


print, 'Step 10: Collecting dark frames...'
	dark_inds=where(strcompress(imtypes,/rem) eq 'DARK' )


	if dark_inds[0] eq -1 then begin
		print, 'Warning! No dark files found!'
		return
	endif

	print, '	Found ',n_elements(dark_inds),' dark frames.'
	
	

	dark_dates=dates[dark_inds]
	dark_exptimes=exptimes[dark_inds]
	print, '	Darks are from dates: ',dark_dates
	print, '	With exposure times :', dark_exptimes
	
	dark_files=files[dark_inds]

obj_darks=dark_files[where(float(dark_exptimes) eq float(obj_exptime))]
	
	
	print, 'Found ',n_elements(obj_darks),' object darks.'
	print, obj_darks
	cen_darks=dark_files[where(dark_exptimes eq cen_exptime)]
	print, 'Found ',n_elements(cen_darks),' center darks.'
	print, cen_darks
	flux_darks=dark_files[where(dark_exptimes eq flux_exptime)]
	print, 'Found ',n_elements(flux_darks),' flux darks.'
	print, flux_darks
	
	
	obj_dark_date=dark_dates[where(dark_exptimes eq obj_exptime)]
	cen_dark_date=dark_dates[where(dark_exptimes eq cen_exptime)]
	flux_dark_date=dark_dates[where(dark_exptimes eq flux_exptime)]
	
	;limiting to dates that were actually observed
	
	if n_elements(where(dark_dates eq obj_dates and dark_exptimes eq obj_exptime)) gt 1 then begin 
		obj_darks=dark_files[where(dark_dates eq obj_dates and dark_exptimes eq obj_exptime)]
		print, 'Using object dark from date of observation only.'
		obj_dark_date=dark_dates[where(dark_dates eq obj_dates and dark_exptimes eq obj_exptime)]
		endif
		
		if n_elements(where(dark_dates eq flux_dates and dark_exptimes eq flux_exptime)) gt 1 then begin 
		flux_darks=dark_files[where(dark_dates eq flux_dates and dark_exptimes eq flux_exptime)]
		print, 'Using flux dark from date of observation only.'
		flux_dark_date=dark_dates[where(dark_dates eq flux_dates and dark_exptimes eq flux_exptime)]
		endif
		
if n_elements(where(dark_dates eq cen_dates and dark_exptimes eq cen_exptime)) gt 1 then begin 
		cen_darks=dark_files[where(dark_dates eq cen_dates and dark_exptimes eq cen_exptime)]
		
		cen_dark_date=dark_dates[where(dark_dates eq cen_dates and dark_exptimes eq cen_exptime)]
		print, 'Using center dark from date of observation only.'
		endif
		

	
	if n_elements(obj_darks) le 0 or  n_elements(cen_darks) le 0 or n_elements(flux_darks) le 0 then begin
		print, 'Dark files not found for object, flux, and center sequences! Halting sequence.'
		return
	endif

print, size(obj_dark_date)

if n_elements( size(obj_dark_date) ) ge 4 then obj_dark_date=obj_dark_date[0]
if n_elements( size(cen_dark_date) ) ge 4 then cen_dark_date=cen_dark_date[0]
if n_elements( size(flux_dark_date) ) ge 4 then flux_dark_date=flux_dark_date[0]	

print, obj_dates[0], obj_dark_date[0]

	if obj_dark_date ne obj_dates[0] or cen_dark_date ne obj_dates[0] or flux_dark_date ne obj_dates[0] then 	print, 'Warning - Science sequence date does not correspond to date of calibration files. Proceeding anyways, but consider using the calibration files from the correct date, if they exist.'
	
	
;-----------------------------------------------------------------------------------------


;begin copy space

	for ii=0, n_elements(obj_darks)-1 do begin
	print, obj_darks
		frame=readfits(obj_darks[ii],hdr,/silent)
		if ii eq 0 then obj_dark_cube=frame else obj_dark_cube=[ [[obj_dark_cube]], [[frame]] ]		
	endfor
	
if n_elements(size(obj_dark_cube)) gt 5 then begin
	print, 'More than one object dark frame found... combining.'

	if dark_comb_imtype eq 'median' then begin
		if (size(obj_dark_cube))(3) mod 2 eq 0 then obj_dark=median(obj_dark_cube,dim=3,/even) else obj_dark=median(obj_dark_cube,dim=3)
	
	endif
	
	if dark_comb_imtype eq 'mean' then obj_dark=mean(obj_dark_cube,dim=3)
endif else obj_dark=obj_dark_cube

print, 'Writing calibration file...'

writefits,input_directory+'IRDIS/calib/obj_dark.fits',obj_dark,hdr
writefits,input_directory+'IRDIS/calib/obj_dark-prebad.fits',obj_dark,hdr

;end copy space

for ii=0, n_elements(cen_darks)-1 do begin
		frame=readfits(cen_darks[ii],hdr,/silent)
		if ii eq 0 then cen_dark_cube=frame else cen_dark_cube=[ [[cen_dark_cube]], [[frame]] ]		
	endfor
	
	
if n_elements(size(cen_dark_cube)) gt 5 then begin
	print, 'More than one star center dark frame found... combining.'
	if dark_comb_imtype eq 'median' then begin
		if (size(cen_dark_cube))(3) mod 2 eq 0 then cen_dark=median(cen_dark_cube,dim=3,/even) else cen_dark=median(cen_dark_cube,dim=3)
	
	endif
	
	if dark_comb_imtype eq 'mean' then cen_dark=mean(cen_dark_cube,dim=3)
endif else cen_dark=cen_dark_cube

print, 'Writing calibration file...'

writefits,input_directory+'IRDIS/calib/cen_dark.fits',cen_dark,hdr
writefits,input_directory+'IRDIS/calib/cen_dark-prebad.fits',cen_dark,hdr

for ii=0, n_elements(flux_darks)-1 do begin
		frame=readfits(flux_darks[ii],hdr,/silent)
		if ii eq 0 then flux_dark_cube=frame else flux_dark_cube=[ [[flux_dark_cube]], [[frame]] ]		
	endfor
	
if n_elements(size(flux_dark_cube)) gt 5 then begin
	print, 'More than one flux calibration dark frame found... combining.'

	if dark_comb_imtype eq 'median' then begin
		if (size(flux_dark_cube))(3) mod 2 eq 0 then flux_dark=median(flux_dark_cube,dim=3,/even) else flux_dark=median(flux_dark_cube,dim=3)
	
	endif
	
	if dark_comb_imtype eq 'mean' then flux_dark=mean(flux_dark_cube,dim=3)
endif else flux_dark=flux_dark_cube

print, 'Writing calibration file...'

writefits,input_directory+'IRDIS/calib/flux_dark.fits',flux_dark,hdr
writefits,input_directory+'IRDIS/calib/flux_dark-prebad.fits',flux_dark,hdr


;-----------------------------------------------------------------------------------------



;use all available darks and scale by exptime
if combine_all_darks then begin
	for ii=0,n_elements(dark_files)-1 do begin
		frame=readfits(dark_files[ii])
		frame=frame/dark_exptimes[ii]
		if ii eq 0 then dark_cube=frame else dark_cube=[ [[dark_cube]], [[frame]] ]
	endfor
	
	medarr,dark_cube,master_dark

	obj_dark=master_dark*obj_exptime
	flux_dark=master_dark*flux_exptime
	cen_dark=master_dark*cen_exptime

	writefits,input_directory+'IRDIS/calib/obj_dark.fits',obj_dark,hdr
	writefits,input_directory+'IRDIS/calib/obj_dark-prebad.fits',obj_dark,hdr
	
	writefits,input_directory+'IRDIS/calib/cen_dark.fits',cen_dark,hdr
	writefits,input_directory+'IRDIS/calib/cen_dark-prebad.fits',cen_dark,hdr

	writefits,input_directory+'IRDIS/calib/flux_dark.fits',flux_dark,hdr
	writefits,input_directory+'IRDIS/calib/flux_dark-prebad.fits',flux_dark,hdr
	


endif

	

	


if remove_bads then begin

;print, 'Step 10a: Correcting bad pixels in dark frame...'



badmap=readfits(badmap_loc)



if find_bads then begin


;badmap[*]=1
badsig=3.
cs=3
smallcube=obj_dark;readfits(cen_files[0])


 print, 'Step 10a... Finding bad pixels... may take a while....'
  for ix=cs, (size(badmap))(1)-1-cs do begin
     for iy=cs, (size(badmap))(2)-1-cs do begin
           box = smallcube[ix-cs:ix+cs,iy-cs:iy+cs,*]
           this =  box[cs,cs,*]  ; this pixel

         ; if mean(this) gt 50 or mean(this) lt -10 then badmap[ix,iy]=0
                     if mean(this) gt bad_thresh or mean(this) lt neg_bad_thresh then badmap[ix,iy]=0

           box[cs,cs,*] = !values.f_nan
           box[cs-1,cs+1,*] = !values.f_nan
           box[cs-1,cs,*] = !values.f_nan
           box[cs-1,cs-1,*] = !values.f_nan
           
           box[cs+1,cs-1,*] = !values.f_nan
           box[cs+1,cs,*] = !values.f_nan
           box[cs+1,cs+1,*] = !values.f_nan
           
           box[cs,cs+1,*] = !values.f_nan
           box[cs,cs-1,*] = !values.f_nan



           if mean(this) gt (mean(box,/nan)+badsig*stddev(box,/nan)) or mean(this) lt (mean(box,/nan)-badsig*stddev(box,/nan)) then badmap[ix,iy] = 0 

        endfor
   ;  print, ' Left Column : ', ix
     endfor
     
     ;repeat for negative (dead) pixels on a frame with some actual light
     if 1 eq 1 then begin
     smallcube=readfits(obj_files[0])
    print, 'Step 10aa... Finding dead pixels... may take a minute....'
    for ix=cs, (size(badmap))(1)-1-cs do begin
     for iy=cs, (size(badmap))(2)-1-cs do begin
           box = smallcube[ix-cs:ix+cs,iy-cs:iy+cs,0]
           this =  box[cs,cs]  ; this pixel

         ; if mean(this) gt 50 or mean(this) lt -10 then badmap[ix,iy]=0
                     if  this lt 0 then badmap[ix,iy]=0

			
           box[cs,cs,*] = !values.f_nan
           box[cs-1,cs+1,*] = !values.f_nan
           box[cs-1,cs,*] = !values.f_nan
           box[cs-1,cs-1,*] = !values.f_nan
           
           box[cs+1,cs-1,*] = !values.f_nan
           box[cs+1,cs,*] = !values.f_nan
           box[cs+1,cs+1,*] = !values.f_nan
           
           box[cs,cs+1,*] = !values.f_nan
           box[cs,cs-1,*] = !values.f_nan

           if mean(this) lt (mean(box,/nan)-badsig*stddev(box,/nan)) then badmap[ix,iy] = 0 

        endfor
   ;  print, ' Left Column : ', ix
     endfor  
     endif
 
endif ;find_bads  

writefits,strcompress(input_directory+'IRDIS/calib/bad_pixel_mask.fits',/rem),badmap


	
;-----------------------------------------------------------------------------------------



print, 'Step 10b: Removing bad pixels from flux dark.'
remove_bad_pixels, badmap,input_directory+'IRDIS/calib/flux_dark.fits';flux_dark
print, 'Step 10b: Removing bad pixels from object dark.'
remove_bad_pixels, badmap,input_directory+'IRDIS/calib/obj_dark.fits';obj_dark
print, 'Step 10b: Removing bad pixels from center dark.'
remove_bad_pixels, badmap,input_directory+'IRDIS/calib/cen_dark.fits';cen_dark

flux_dark=readfits(input_directory+'IRDIS/calib/flux_dark.fits')
obj_dark=readfits(input_directory+'IRDIS/calib/obj_dark.fits')
cen_dark=readfits(input_directory+'IRDIS/calib/cen_dark.fits')


endif ;remove bads if



;-----------------------------------------------------------------------------------------


if sub_dark then print, 'Step 11: Subtracting dark current from object files...'


	
	for ii=0, n_elements(obj_files)-1 do begin
		frame=readfits(obj_files[ii],hdr,/silent)
		if n_elements(size(frame)) gt 5 then begin $
			if sub_dark then for jj=0, (size(frame))(3)-1 do frame[*,*,jj]=frame[*,*,jj]-obj_dark & $
			endif else frame=frame-obj_dark
		writefits,strcompress(input_directory+'IRDIS/products/'+'IRDIS_OBJECT_'+string(obj_dates_full[ii])+'_.fits',/rem),frame,hdr
	endfor
	
;-----------------------------------------------------------------------------------------

	
if sub_dark then print, 'Step 12: Subtracting dark current from star center files...'
	
	
	for ii=0, n_elements(cen_files)-1 do begin
		frame=readfits(cen_files[ii],hdr,/silent)
		if n_elements(size(frame)) gt 5 then begin $
			if sub_dark then for jj=0, (size(frame))(3)-1 do frame[*,*,jj]=frame[*,*,jj]-cen_dark ; $
			endif else frame=frame-cen_dark
		writefits,strcompress(input_directory+'IRDIS/center/'+'IRDIS_CENTER_'+string(cen_dates_full[ii])+'_.fits',/rem),frame,hdr
	endfor
	
;-----------------------------------------------------------------------------------------


if sub_dark then print, 'Step 13: Subtracting dark current from flux calibration files...'
	

	
	for ii=0, n_elements(flux_files)-1 do begin
		frame=readfits(flux_files[ii],hdr,/silent)
		if n_elements(size(frame)) gt 5 then begin $
			if sub_dark then for jj=0, (size(frame))(3)-1 do frame[*,*,jj]=frame[*,*,jj]-flux_dark; & $
			endif else frame=frame-flux_dark
		writefits,strcompress(input_directory+'IRDIS/flux/'+'IRDIS_FLUX_'+string(flux_dates_full[ii])+'_.fits',/rem),frame,hdr
	endfor
	
;-----------------------------------------------------------------------------------------

	
print, 'Step 14: Collecting flat frames...'
	flat_inds=where(strcompress(imtypes,/rem) eq 'FLAT,LAMP' )
	flat_exptimes=exptimes[flat_inds]
	
	;Are there two one second flats? If so the first is probably bad, so let's ignore it.
	if flat_exptimes[0] eq flat_exptimes[1] then flat_inds=flat_inds[1:n_elements(flat_inds)-1]
	
	print, '	Found ',n_elements(flat_inds),' flat frames.'
	if n_elements(flat_inds) gt 5 then flat_inds=flat_inds[0:4]  ;limit to five flats since ESO often includes extras
	flat_dates=dates[flat_inds]
	flat_exptimes=exptimes[flat_inds]
	print, '	Flats are from dates: ',flat_dates
	print, '	With exposure times :', flat_exptimes

	flat_files=files[flat_inds]

	print, 'Flat Files:'
	print, flat_files
	;hak

	obj_flat=fltarr(2048,1024)
	flux_flat=fltarr(2048,1024)
	cen_flat=fltarr(2048,1024)
	

	;first generate a cube of flats
	for ii=0,n_elements(flat_inds)-1 do begin	;loop through flat images
		flat=readfits(flat_files[ii],flathd)	
		if n_elements(size(flat)) gt 5 then medarr, flat, flat	;if cube then median combine

		
		;subtract dark from flat frame
		scale=flat_exptimes[ii]/flux_exptime
		flat=flat-scale*flux_dark

		if ii eq 0 then print, 'Prepping flat frame by subtracting dark current and correcting bad pixels.'
  		fixpix,flat,badmap,outimg
		flat=outimg

		if ii eq 0 then flat_cube=flat else flat_cube=[ [[flat_cube]], [[flat]] ]			
	endfor

	if not keyword_set(no_flat_ird) then no_flat_ird=0

	print, 'Determining flat with second order polynomial... may take a few minutes.'
	for xx=0,2047 do begin	;loop through x pixels
		for yy=0,1023 do begin	;loop through y pixels
			for ii=0,n_elements(flat_inds)-1 do begin	;loop through flat images
				pxval=flat_cube[xx,yy,ii]
				if ii eq 0 then pxvals=pxval else pxvals=[pxvals,pxval]
			endfor
			
			;now fit a polynomial to the values as a function of exptime
			
		
			if not no_flat_ird then begin
			result = poly_fit(flat_exptimes, pxvals , 2,MEASURE_ERRORS=sqrt(pxvals),SIGMA=sigma) 
			; Print the coefficients:
			;PRINT, 'Coefficients: ', result
			;print, 'Pixel values: ', pxvals
			;print, 'Fit values:', poly(flat_exptimes,result)
			;print, 'Residual / pixel value = ', (poly(flat_exptimes,result)-pxvals)/pxvals
			;print, 'Exposure times: ',flat_exptimes

			pxval_obj=poly(obj_exptime,result)
			pxval_flux=poly(flux_exptime,result)
			pxval_cen=poly(cen_exptime,result)
			endif else begin ;if not using a flat, just set to 1
			pxval_obj=1.
			pxval_flux=1.
			pxval_cen=1.
			endelse
			obj_flat[xx,yy]=pxval_obj
			flux_flat[xx,yy]=pxval_flux
			cen_flat[xx,yy]=pxval_cen
			
			;hak
		endfor
	endfor


		
		;ensure flats have median of unity in illuminated regions to not bias fluxes
		print, 'Dividing flats by median of illuminated region.'
		frame=obj_flat
			;frame[0:50,*]=!values.f_nan
			;frame[935:1080,*]=!values.f_nan
			;frame[1960:2047,*]=!values.f_nan
			;frame[*,0:35]=!values.f_nan
			;frame[*,1010:1023]=!values.f_nan

		;mmm, frame, obj_flat_scale
		;obj_flat_scale=median(frame)
		mmm,frame[100:900,100:900],flat_scalel
		mmm,frame[1124:1940,100:900],flat_scaler
		frame[0:1024,*]=frame[0:1024,*]/flat_scalel
		frame[1024:2047,*]=frame[1024:2047,*]/flat_scaler
		obj_flat=frame
		;print, 'Object flat median:',obj_flat_scale
		;obj_flat=obj_flat/obj_flat_scale

		frame=cen_flat
			mmm,frame[100:900,100:900],flat_scalel
		mmm,frame[1124:1940,100:900],flat_scaler
		frame[0:1024,*]=frame[0:1024,*]/flat_scalel
		frame[1024:2047,*]=frame[1024:2047,*]/flat_scaler
		cen_flat=frame

		frame=flux_flat
			mmm,frame[100:900,100:900],flat_scalel
		mmm,frame[1124:1940,100:900],flat_scaler
		frame[0:1024,*]=frame[0:1024,*]/flat_scalel
		frame[1024:2047,*]=frame[1024:2047,*]/flat_scaler
		flux_flat=frame

		obj_flat[where(finite(obj_flat)) ne 1 ]=1.	;if no flat data for a particular pixel, no correction is made
		cen_flat[where(finite(cen_flat)) ne 1 ]=1.
		flux_flat[where(finite(flux_flat)) ne 1 ]=1.


		writefits,strcompress(input_directory+'IRDIS/calib/obj_flat.fits',/rem),obj_flat
		writefits,strcompress(input_directory+'IRDIS/calib/cen_flat.fits',/rem),cen_flat
		writefits,strcompress(input_directory+'IRDIS/calib/flux_flat.fits',/rem),flux_flat



	print,'Done. Masking non-illuminated detector regions...'


	print, 'Done. Calibrating data with flat field division...'

	obj_files=file_search(input_directory+'IRDIS/products/','*.fits')

	
	for ii=0, n_elements(obj_files)-1 do begin
		frame=readfits(obj_files[ii],hdr,/silent)

		if n_elements(size(frame)) gt 5 then begin 
			frame[0:50,*,*]=0.
			frame[935:1080,*,*]=0.
			frame[1960:2047,*,*]=0.
			frame[*,0:35,*]=0.
			frame[*,1010:1023,*]=0.
		endif else begin
			frame[0:50,*]=0.
			frame[935:1080,*]=0.
			frame[1960:2047,*]=0.
			frame[*,0:35]=0.
			frame[*,1010:1023]=0.

		endelse

		if n_elements(size(frame)) gt 5 then begin $
			if div_flat then for jj=0, (size(frame))(3)-1 do frame[*,*,jj]=frame[*,*,jj]/obj_flat & $
			endif else frame=frame/obj_flat
		frame[where(finite(frame) eq 0)]=0.
		writefits,strcompress(input_directory+'IRDIS/products/'+'IRDIS_OBJECT_'+string(obj_dates_full[ii])+'_.fits',/rem),frame,hdr
	endfor
	


	

	
	
	for ii=0, n_elements(cen_files)-1 do begin
		frame=readfits(cen_files[ii],hdr,/silent)

		if n_elements(size(frame)) gt 5 then begin 
			frame[0:50,*,*]=0.
			frame[935:1080,*,*]=0.
			frame[1960:2047,*,*]=0.
			frame[*,0:35,*]=0.
			frame[*,1010:1023,*]=0.
		endif else begin
			frame[0:50,*]=0.
			frame[935:1080,*]=0.
			frame[1960:2047,*]=0.
			frame[*,0:35]=0.
			frame[*,1010:1023]=0.

		endelse

		if n_elements(size(frame)) gt 5 then begin $
			if div_flat then for jj=0, (size(frame))(3)-1 do frame[*,*,jj]=frame[*,*,jj]/cen_flat ; $
			endif else frame=frame/cen_flat
		frame[where(finite(frame) eq 0)]=0.
		writefits,strcompress(input_directory+'IRDIS/center/'+'IRDIS_CENTER_'+string(cen_dates_full[ii])+'_.fits',/rem),frame,hdr
	endfor
	




	

	
	for ii=0, n_elements(flux_files)-1 do begin
		frame=readfits(flux_files[ii],hdr,/silent)

		if n_elements(size(frame)) gt 5 then begin 
			frame[0:50,*,*]=0.
			frame[935:1080,*,*]=0.
			frame[1960:2047,*,*]=0.
			frame[*,0:35,*]=0.
			frame[*,1010:1023,*]=0.
		endif else begin
			frame[0:50,*]=0.
			frame[935:1080,*]=0.
			frame[1960:2047,*]=0.
			frame[*,0:35]=0.
			frame[*,1010:1023]=0.

		endelse

		if n_elements(size(frame)) gt 5 then begin $
			if div_flat then for jj=0, (size(frame))(3)-1 do frame[*,*,jj]=frame[*,*,jj]/flux_flat; & $
			endif else frame=frame/flux_flat
		frame[where(finite(frame) eq 0)]=0.
		writefits,strcompress(input_directory+'IRDIS/flux/'+'IRDIS_FLUX_'+string(flux_dates_full[ii])+'_.fits',/rem),frame,hdr
	endfor


	if div_flat eq 0 then print, 'Skipping flat field division...'



	
	;skipping flat fielding for now...
	
	;print, 'Skipping flat fielding... Account for flat field variations in final uncertainty.'
	

;-----------------------------------------------------------------------------------------

if remove_bads  then begin


print, 'Step 16: Removing bad pixels from flux calibration data.'
remove_bad_pixels, badmap,strcompress(input_directory+'IRDIS/flux/',/rem)

;-----------------------------------------------------------------------------------------


print, 'Step 17: Removing bad pixels from star center data.'
remove_bad_pixels, badmap,strcompress(input_directory+'IRDIS/center/',/rem)

;-----------------------------------------------------------------------------------------


print, 'Step 18: Removing bad pixels from object data.'
remove_bad_pixels, badmap,strcompress(input_directory+'IRDIS/products/',/rem)

endif else print, 'Skipping bad pixel correction.'	

;-----------------------------------------------------------------------------------------


print, 'Step 19: Splitting left and right.'

;leftmask=readfits('/Users/kevinwagner/Data/VLT/SPHERE/calib/irdis_left_mask.fits')
;rightmask=readfits('/Users/kevinwagner/Data/VLT/SPHERE/calib/irdis_right_mask.fits')

	files=file_search(input_directory+'IRDIS/flux/','*_.fits',count=cnt)
	for ii=0,cnt-1 do begin
		frame=readfits(files[ii],hdr,/silent)
		if n_elements(size(frame)) gt 5 then left = frame[0:1023,*,*] else left=frame[0:1023,*]
		if n_elements(size(frame)) gt 5 then right = frame[1024:2047,*,*] else right=frame[1024:2047,*]
		
		writefits,strcompress(input_directory+'IRDIS/flux/'+'IRDIS_FLUX_'+string(flux_dates_full[ii])+'_left.fits',/rem),left,hdr
							writefits,strcompress(input_directory+'IRDIS/flux/'+'IRDIS_FLUX_'+string(flux_dates_full[ii])+'_right.fits',/rem),right,hdr
							
							

	endfor
	
	file_delete,files
	
	files=file_search(input_directory+'IRDIS/center/','*_.fits',count=cnt)
	for ii=0,cnt-1 do begin
		frame=readfits(files[ii],hdr,/silent)
		if n_elements(size(frame)) gt 5 then left = frame[0:1023,*,*] else left=frame[0:1023,*]
		if n_elements(size(frame)) gt 5 then right = frame[1024:2047,*,*] else right=frame[1024:2047,*]
		writefits,strcompress(input_directory+'IRDIS/center/'+'IRDIS_CENTER_'+string(cen_dates_full[ii])+'_left.fits',/rem),left,hdr
							writefits,strcompress(input_directory+'IRDIS/center/'+'IRDIS_CENTER_'+string(cen_dates_full[ii])+'_right.fits',/rem),right,hdr

	endfor
	
		file_delete,files

	
	files=file_search(input_directory+'IRDIS/products/','*_.fits',count=cnt)
	for ii=0,cnt-1 do begin
		frame=readfits(files[ii],hdr,/silent)
		if n_elements(size(frame)) gt 5 then left = frame[0:1023,*,*] else left=frame[0:1023,*]
		if n_elements(size(frame)) gt 5 then right = frame[1024:2047,*,*] else right=frame[1024:2047,*]
		writefits,strcompress(input_directory+'IRDIS/products/'+'IRDIS_OBJECT_'+string(obj_dates_full[ii])+'_left.fits',/rem),left,hdr
							writefits,strcompress(input_directory+'IRDIS/products/'+'IRDIS_OBJECT_'+string(obj_dates_full[ii])+'_right.fits',/rem),right,hdr

	endfor
	
		file_delete,files

endif ;ifs_only if

endif ;skip initial steps if

;delete the download folder if it is still there
if skip_start then begin

	test_file=file_search(input_directory+'/data_with_raw_calibs/','*',count=testcnt)
	if testcnt ge 1 then print, 'Found extra downloaded files that are not needed... deleting to save space.'
	if testcnt ge 1 then file_delete, input_directory+'/data_with_raw_calibs/',/recursive
endif

;-----------------------------------------------------------------------------------------

cgcleanup

if not skip and irdis_only then if no_cen_files_check then auto_cen=0

if  keyword_set(inject) then begin proc_orig = 0 & proc_inj=1 & endif else begin $
				   proc_orig = 1 & proc_inj=0 & endelse

;if neither keyword is set then process the original data only
if not keyword_set(proc_inj) and not keyword_set(proc_orig) then proc_orig=1

if post_process_irdis then begin


if keyword_set(second_cen) then begin

if keyword_set(klip) then begin

if keyword_set(suffix) then begin
	if skip_start then begin
	print, 'Running reduce_irdis.pro with injections'

	if proc_inj then begin	
	if auto_cen then reduce_irdis,object_name,reduction_path=reduction_path, /inj,/autocen,/skip_start,suffix=suffix,/klip,/second_cen,rho=rho,theta=theta,contrast=contrast else $
		 reduce_irdis,object_name,reduction_path=reduction_path, /inj,/skip_start,suffix=suffix,/klip ,/second_cen,rho=rho,theta=theta,contrast=contrast
	endif

	if proc_orig then begin
	print, 'Running reduce_irdis.pro '	
	if auto_cen then  reduce_irdis,object_name,reduction_path=reduction_path, /autocen,/skip_start,suffix=suffix,/klip ,/second_cen else $
		reduce_irdis,object_name,reduction_path=reduction_path, /skip_start,suffix=suffix,/klip ,/second_cen
	endif	
	endif else begin


	if proc_inj then begin
	print, 'Running reduce_irdis.pro with injections'	
	if auto_cen then reduce_irdis,object_name,reduction_path=reduction_path, /inj,/autocen,suffix=suffix,/klip  ,/second_cen,rho=rho,theta=theta,contrast=contrast else $
		 reduce_irdis,object_name,reduction_path=reduction_path, /inj,suffix=suffix,/klip,/second_cen,rho=rho,theta=theta,contrast=contrast
	endif
	if proc_orig then begin
	print, 'Running reduce_irdis.pro '	
	if auto_cen then  reduce_irdis,object_name,reduction_path=reduction_path, /autocen,suffix=suffix,/klip  ,/second_cen else $
		reduce_irdis,object_name,reduction_path=reduction_path,suffix=suffix,/klip ,/second_cen
	endif
	endelse	

endif else begin

	if skip_start then begin
	if proc_inj then begin
	print, 'Running reduce_irdis.pro with injections'	
	if auto_cen then reduce_irdis,object_name,reduction_path=reduction_path, /inj,/autocen,/skip_start ,/klip,/second_cen,rho=rho,theta=theta,contrast=contrast else $
		 reduce_irdis,object_name,reduction_path=reduction_path, /inj,/skip_start,/klip ,/second_cen,rho=rho,theta=theta,contrast=contrast
	endif
	if proc_orig then begin
	print, 'Running reduce_irdis.pro '	
	if auto_cen then  reduce_irdis,object_name,reduction_path=reduction_path, /autocen,/skip_start,/klip,/second_cen  else $
		reduce_irdis,object_name,reduction_path=reduction_path, /skip_start,/klip ,/second_cen
	endif
	endif else begin

	if proc_inj then begin
	print, 'Running reduce_irdis.pro with injections'	
	if auto_cen then reduce_irdis,object_name,reduction_path=reduction_path, /inj,/autocen,/klip,/second_cen ,rho=rho,theta=theta,contrast=contrast else $
		 reduce_irdis,object_name,reduction_path=reduction_path, /inj,/klip,/second_cen,rho=rho,theta=theta,contrast=contrast 
	endif
	if proc_orig then begin
	print, 'Running reduce_irdis.pro '	
	if auto_cen then  reduce_irdis,object_name,reduction_path=reduction_path, /autocen, /klip,/second_cen  else $
		reduce_irdis,object_name,reduction_path=reduction_path,/klip ,/second_cen
	endif
	endelse	

endelse	

endif else begin

if keyword_set(suffix) then begin
	if skip_start then begin
	print, 'Running reduce_irdis.pro with injections'

	if proc_inj then begin	
	if auto_cen then reduce_irdis,object_name,reduction_path=reduction_path, /inj,/autocen,/skip_start,suffix=suffix,/second_cen,rho=rho,theta=theta,contrast=contrast else $
		 reduce_irdis,object_name,reduction_path=reduction_path, /inj,/skip_start,suffix=suffix,/second_cen,rho=rho,theta=theta,contrast=contrast
	endif

	if proc_orig then begin
	print, 'Running reduce_irdis.pro '	
	if auto_cen then  reduce_irdis,object_name,reduction_path=reduction_path, /autocen,/skip_start,suffix=suffix,/second_cen else $
		reduce_irdis,object_name,reduction_path=reduction_path, /skip_start,suffix=suffix,/second_cen
	endif	
	endif else begin


	if proc_inj then begin
	print, 'Running reduce_irdis.pro with injections'	
	if auto_cen then reduce_irdis,object_name,reduction_path=reduction_path, /inj,/autocen,suffix=suffix,/second_cen,rho=rho,theta=theta,contrast=contrast else $
		 reduce_irdis,object_name,reduction_path=reduction_path, /inj,suffix=suffix,/second_cen,rho=rho,theta=theta,contrast=contrast
	endif
	if proc_orig then begin
	print, 'Running reduce_irdis.pro '	
	if auto_cen then  reduce_irdis,object_name,reduction_path=reduction_path, /autocen, suffix=suffix,/second_cen else $
		reduce_irdis,object_name,reduction_path=reduction_path,suffix=suffix,/second_cen
	endif
	endelse	

endif else begin

	if skip_start then begin
	if proc_inj then begin
	print, 'Running reduce_irdis.pro with injections'	
	if auto_cen then reduce_irdis,object_name,reduction_path=reduction_path, /inj,/autocen,/skip_start,/second_cen,rho=rho,theta=theta,contrast=contrast else $
		 reduce_irdis,object_name,reduction_path=reduction_path, /inj,/skip_start,/second_cen,rho=rho,theta=theta,contrast=contrast
	endif
	if proc_orig then begin
	print, 'Running reduce_irdis.pro '	
	if auto_cen then  reduce_irdis,object_name,reduction_path=reduction_path, /autocen,/skip_start,/second_cen else $
		reduce_irdis,object_name,reduction_path=reduction_path, /skip_start,/second_cen
	endif
	endif else begin

	if proc_inj then begin
	print, 'Running reduce_irdis.pro with injections'	
	if auto_cen then reduce_irdis,object_name,reduction_path=reduction_path, /inj,/autocen,/second_cen,rho=rho,theta=theta,contrast=contrast else $
		 reduce_irdis,object_name,reduction_path=reduction_path, /inj,/second_cen,rho=rho,theta=theta,contrast=contrast
	endif
	if proc_orig then begin
	print, 'Running reduce_irdis.pro '	
	if auto_cen then  reduce_irdis,object_name,reduction_path=reduction_path, /autocen,/second_cen else $
		reduce_irdis,object_name,reduction_path=reduction_path,/second_cen
	endif
	endelse	

endelse	

endelse

endif

if not keyword_set(second_cen) then begin

if keyword_set(klip) then begin

if keyword_set(suffix) then begin
	if skip_start then begin
	print, 'Running reduce_irdis.pro with injections'

	if proc_inj then begin	
	if auto_cen then reduce_irdis,object_name,reduction_path=reduction_path, /inj,/autocen,/skip_start,suffix=suffix,/klip,rho=rho,theta=theta,contrast=contrast else $
		 reduce_irdis,object_name,reduction_path=reduction_path, /inj,/skip_start,suffix=suffix,/klip ,rho=rho,theta=theta,contrast=contrast
	endif

	if proc_orig then begin
	print, 'Running reduce_irdis.pro '	
	if auto_cen then  reduce_irdis,object_name,reduction_path=reduction_path, /autocen,/skip_start,suffix=suffix,/klip  else $
		reduce_irdis,object_name,reduction_path=reduction_path, /skip_start,suffix=suffix,/klip 
	endif	
	endif else begin


	if proc_inj then begin
	print, 'Running reduce_irdis.pro with injections'	
	if auto_cen then reduce_irdis,object_name,reduction_path=reduction_path, /inj,/autocen,suffix=suffix,/klip,rho=rho,theta=theta,contrast=contrast  else $
		 reduce_irdis,object_name,reduction_path=reduction_path, /inj,suffix=suffix,/klip ,rho=rho,theta=theta,contrast=contrast
	endif
	if proc_orig then begin
	print, 'Running reduce_irdis.pro '	
	if auto_cen then  reduce_irdis,object_name,reduction_path=reduction_path, /autocen,suffix=suffix,/klip  else $
		reduce_irdis,object_name,reduction_path=reduction_path,suffix=suffix,/klip 
	endif
	endelse	

endif else begin

	if skip_start then begin
	if proc_inj then begin
	print, 'Running reduce_irdis.pro with injections'	
	if auto_cen then reduce_irdis,object_name,reduction_path=reduction_path, /inj,/autocen,/skip_start ,/klip,rho=rho,theta=theta,contrast=contrast else $
		 reduce_irdis,object_name,reduction_path=reduction_path, /inj,/skip_start,/klip ,rho=rho,theta=theta,contrast=contrast
	endif
	if proc_orig then begin
	print, 'Running reduce_irdis.pro '	
	if auto_cen then  reduce_irdis,object_name,reduction_path=reduction_path, /autocen,/skip_start,/klip  else $
		reduce_irdis,object_name,reduction_path=reduction_path, /skip_start,/klip 
	endif
	endif else begin

	if proc_inj then begin
	print, 'Running reduce_irdis.pro with injections'	
	if auto_cen then reduce_irdis,object_name,reduction_path=reduction_path, /inj,/autocen,/klip,rho=rho,theta=theta,contrast=contrast  else $
		 reduce_irdis,object_name,reduction_path=reduction_path, /inj,/klip ,rho=rho,theta=theta,contrast=contrast
	endif
	if proc_orig then begin
	print, 'Running reduce_irdis.pro '	
	if auto_cen then  reduce_irdis,object_name,reduction_path=reduction_path, /autocen, /klip  else $
		reduce_irdis,object_name,reduction_path=reduction_path,/klip 
	endif
	endelse	

endelse	

endif else begin

if keyword_set(suffix) then begin
	if skip_start then begin
	print, 'Running reduce_irdis.pro with injections'

	if proc_inj then begin	
	if auto_cen then reduce_irdis,object_name,reduction_path=reduction_path, /inj,/autocen,/skip_start,suffix=suffix,rho=rho,theta=theta,contrast=contrast else $
		 reduce_irdis,object_name,reduction_path=reduction_path, /inj,/skip_start,suffix=suffix,rho=rho,theta=theta,contrast=contrast
	endif

	if proc_orig then begin
	print, 'Running reduce_irdis.pro '	
	if auto_cen then  reduce_irdis,object_name,reduction_path=reduction_path, /autocen,/skip_start,suffix=suffix else $
		reduce_irdis,object_name,reduction_path=reduction_path, /skip_start,suffix=suffix
	endif	
	endif else begin


	if proc_inj then begin
	print, 'Running reduce_irdis.pro with injections'	
	if auto_cen then reduce_irdis,object_name,reduction_path=reduction_path, /inj,/autocen,suffix=suffix,rho=rho,theta=theta,contrast=contrast else $
		 reduce_irdis,object_name,reduction_path=reduction_path, /inj,suffix=suffix,rho=rho,theta=theta,contrast=contrast
	endif
	if proc_orig then begin
	print, 'Running reduce_irdis.pro '	
	if auto_cen then  reduce_irdis,object_name,reduction_path=reduction_path, /autocen, suffix=suffix else $
		reduce_irdis,object_name,reduction_path=reduction_path,suffix=suffix
	endif
	endelse	

endif else begin

	if skip_start then begin
	if proc_inj then begin
	print, 'Running reduce_irdis.pro with injections'	
	if auto_cen then reduce_irdis,object_name,reduction_path=reduction_path, /inj,/autocen,/skip_start,rho=rho,theta=theta,contrast=contrast else $
		 reduce_irdis,object_name,reduction_path=reduction_path, /inj,/skip_start,rho=rho,theta=theta,contrast=contrast
	endif
	if proc_orig then begin
	print, 'Running reduce_irdis.pro '	
	if auto_cen then  reduce_irdis,object_name,reduction_path=reduction_path, /autocen,/skip_start else $
		reduce_irdis,object_name,reduction_path=reduction_path, /skip_start
	endif
	endif else begin

	if proc_inj then begin
	print, 'Running reduce_irdis.pro with injections'	
	if auto_cen then reduce_irdis,object_name,reduction_path=reduction_path, /inj,/autocen,rho=rho,theta=theta,contrast=contrast else $
		 reduce_irdis,object_name,reduction_path=reduction_path, /inj,rho=rho,theta=theta,contrast=contrast
	endif
	if proc_orig then begin
	print, 'Running reduce_irdis.pro '	
	if auto_cen then  reduce_irdis,object_name,reduction_path=reduction_path, /autocen else $
		reduce_irdis,object_name,reduction_path=reduction_path
	endif
	endelse	

endelse	

endelse

endif

	;if proc_inj then print, 'Deleting data cube with injected planets (this may be easily regenerated).'
	;if proc_inj then file_delete,input_directory+'IRDIS/processed/'+object_name+'_left_cube_inject.fits
	;if proc_inj then file_delete,input_directory+'IRDIS/processed/'+object_name+'_right_cube_inject.fits'
	;if klip then file_delete,input_directory+'IRDIS/processed/'+object_name+'_left_cube_klip.fits'
	;if klip then file_delete,input_directory+'IRDIS/processed/'+object_name+'_right_cube_klip.fits'

	
endif



;-----------------------------------------------------------------------------------------

if reduce_irdis_sizes and not skip_start then begin
	print, 'Zipping raw data and deleting old files - may take a minute...'
	file_zip, input_directory+'IRDIS/raw/', input_directory+'IRDIS/raw.zip'
	file_delete,input_directory+'IRDIS/raw/',/recursive
	print, 'Deleting intermediate data products...'
	file_delete,input_directory+'IRDIS/products/',/recursive
	print, 'Done.'	
endif else print, 'Skipping file size reduction...'
;-----------------------------------------------------------------------------------------

cgcleanup

if process_ifs then begin

if not skip_start then begin

print, 'Starting IFS reduction...'
;-----------------------------------------------------------------------------------------
;  IFS data is currently just split on the basis of whether or not it has the IFS tag. There
;  are probably too many files for the data_reduction_ifs.sh script to know what to do with.
;  Thus, first step of this script is to sort the files further, keeping only the calibrations
;  with the correct dates, etc. similar to the IRDIS method above. 	
;-----------------------------------------------------------------------------------------



ROOT=strcompress(input_directory+'IFS/',/rem)
;MODE='YJH'



;first need to read in the raw IFS files and separate them into raw and raw-extra


print, 'Searching for .fits files... Done.'
	files=file_search(input_directory+'IFS/raw/','*.fits',count=cnt)	;first search for compressed files
	if cnt lt 1 then begin
	print, 'Searching for zip files... Done.'
	files=file_search(input_directory+'IFS/raw.zip',count=cnt)	;first search for compressed files
	if cnt ge 1 then print, 'Unzipping...'	
	if cnt ge 1 then file_unzip,files
	files=file_search(input_directory+'IFS/raw/','*.fits',count=cnt)	;first search for compressed files
	endif

	
print, '	Found ', cnt, ' fits files.'


;-----------------------------------------------------------------------------------------


print, 'Looping through fits headers and sorting files...'
for ii=0, cnt-1 do begin
	frame=readfits(files[ii],hdr,/silent)
	
	
	detector=''
	imtype=''
	date=''
	exptime=''
	ndit=''
	dbfilt=''
	ndfilt=''
	
	;extracting info...
	;want to know what detector, imtype, date, DIT, ND filter, DB filter
	detector=string(esopar(hdr,'HIERARCH ESO DET NAME '))
	imtype=esopar(hdr,'HIERARCH ESO DPR TYPE ')
	imcomb=string(esopar(hdr,'HIERARCH ESO INS2 COMB IFS '))
	date=fxpar(hdr,'DATE')
	object=string(fxpar(hdr,'OBJECT'))
	date_full=date
	date=strmid(date,0,10)	;only keep YYYY-MM-DD
	exptime=fxpar(hdr,'EXPTIME')
	ndit=esopar(hdr,'HIERARCH ESO DET NDIT ')
	dbfilt=esopar(hdr,'ESO INS1 OPTI2 NAME ')
	ndfilt = esopar(hdr,'ESO INS4 COMB IND')
	if isa(ndfilt, 'string') eq 0 then ndfilt = esopar(hdr,'ESO INS4 FILT2 NAME')  

	if ii eq 0 then begin
		detectors=detector
		dates=date
		full_dates=date_full
		imtypes=imtype
		imcombs=imcomb
		exptimes=exptime
		ndfilts=ndfilt
		dbfilts=dbfilt
		ndits=ndit
		objects=object
	endif else begin
		detectors=[detectors,detector]
		dates=[dates,date]
		imtypes=[imtypes,imtype]
		imcombs=[imcombs,imcomb]
		exptimes=[exptimes,exptime]
		ndfilts=[ndfilts,ndfilt]
		dbfilts=[dbfilts,dbfilt]
		ndits=[ndits,ndit]
		full_dates=[full_dates,date_full]
		objects=[objects,object]
	endelse

	if speaking then begin
		print, strcompress('	Frame '+string(ii)+':')
		print, '	Date       = ',date
		print, '	Object     = ',object
		print, '	Detector   = ',detector
		print, '	Frame imtype = ',imtype
		print, '	Frame imcomb = ',imcomb
		print, '	Exptime    = ',exptime
		print, '	NDITs    = ',ndit
		print, '	ND filter  = ',ndfilt
		print, '	DB Filter  = ',dbfilt
	endif
endfor
print, '	Done.'


if n_elements(where( strcompress(imcombs,/rem) eq 'OBS_YJ' )) gt 1 then yj=1 else yj=0


;see if darks are missing and copy them if needed
ifsdarks=where(imcombs eq strcompress('CAL_DARK',/rem) and exptimes eq 1.6507260)

if n_elements(ifsdarks) lt 2 then begin

	fs=file_search(input_directory+'IFS/raw/SPHER.2016-06-11T13:21:03.270.fits',COUNT=COUNT)

	if COUNT lt 1 then file_copy,'/Users/kevinwagner/Data/SPHERE/SPHERE IFS darks/1.65s/SPHER.2016-06-11T13:21:03.270.fits' ,input_directory+'IFS/raw/SPHER.2016-06-11T13:21:03.270.fits'

	fs=file_search(input_directory+'IFS/raw/SPHER.2016-06-11T13:22:28.145.fits',COUNT=COUNT)

	if COUNT lt 1 then file_copy,'/Users/kevinwagner/Data/SPHERE/SPHERE IFS darks/1.65s/SPHER.2016-06-11T13:22:28.145.fits' ,input_directory+'IFS/raw/SPHER.2016-06-11T13:22:28.145.fits'

endif


if yj then mode='YJ' else mode='YJH'

;separating dates into other forms
	yyyys=strmid(dates,0,4)
	mms=strmid(dates,5,2)
	dds=strmid(dates,9,2)

;script will fail if it encounters more than two flats, so here we choose to keep only the first two and move the rest
;can make this more intelligent by selecting based on dates instead

white_flats=where(imcombs eq strcompress('CAL_BB_2_'+mode,/rem))
print, files[white_flats]
print, n_elements(white_flats)




if n_elements(white_flats) gt 2 then file_move, files[white_flats[2:n_elements(white_flats)-1]], input_directory+'IFS/extra/'


nb1_flats=where(imcombs eq strcompress('CAL_NB1_1_'+mode,/rem))
print, files[nb1_flats]
print, n_elements(nb1_flats)

if n_elements(nb1_flats) gt 2 then file_move, files[nb1_flats[2:n_elements(nb1_flats)-1]], input_directory+'IFS/extra/'

nb2_flats=where(imcombs eq strcompress('CAL_NB2_1_'+mode,/rem))
print, files[nb2_flats]
print, n_elements(nb2_flats)

if n_elements(nb2_flats) gt 2 then file_move, files[nb2_flats[2:n_elements(nb2_flats)-1]], input_directory+'IFS/extra/'

nb3_flats=where(imcombs eq strcompress('CAL_NB3_1_'+mode,/rem))
print, files[nb3_flats]
print, n_elements(nb3_flats)

if n_elements(nb3_flats) gt 2 then file_move, files[nb3_flats[2:n_elements(nb3_flats)-1]], input_directory+'IFS/extra/'


nb4_flats=where(imcombs eq strcompress('CAL_NB4_2_'+mode,/rem))
print, files[nb4_flats]
print, n_elements(nb4_flats)

if n_elements(nb4_flats) gt 2 then file_move, files[nb4_flats[2:n_elements(nb4_flats)-1]], input_directory+'IFS/extra/'

specpos=where(imtypes eq strcompress('SPECPOS,LAMP',/rem))
print, files[specpos]
print, n_elements(specpos)

if n_elements(specpos) ge 2 then file_move, files[specpos[1:n_elements(specpos)-1]], input_directory+'IFS/extra/'


wavcals=where(imtypes eq strcompress('WAVE,LAMP',/rem))
print, files[wavcals]
print, n_elements(wavcals)

if n_elements(wavcals) ge 2 then file_move, files[wavcals[1:n_elements(wavcals)-1]], input_directory+'IFS/extra/'

save,filename=input_directory+'IFS/products/mode.sav',mode, yj

;print, 'Running command: ./idl/sphere-tools/SPHERE-legacy-master/ifs_reduction/data_reduction_ifs_new.sh '+ root+' '+mode
;spawn,'./idl/sphere-tools/SPHERE-legacy-master/ifs_reduction/data_reduction_ifs_new.sh '+root+' '+ mode



;if not keyword_set(skip_AV_pre) then begin
findpro, 'data_reduction_ifs', dirlist=dirlist
path_to_AV='Users/kevinwagner/idl/sphere-tools/SPHERE-legacy-master/ifs_reduction/';dirlist
print, 'Running bash command: ',strcompress('.'+path_to_AV+'data_reduction_ifs.sh',/rem)+' '+root+' '+ mode
cd, '/'
spawn, strcompress( './'+path_to_AV+'data_reduction_ifs.sh',/rem)+' '+root+' '+ mode
;endif




;hak
;finishing pre-processing
;if not keyword_set(skip_AV_post) then begin
if keyword_set(second_cen) then data_reduction_ifs, root,/second_cen else data_reduction_ifs, root
;endif

;clean display
cgcleanup

endif ;skip_start if

restore,filename=input_directory+'IFS/products/mode.sav'


if keyword_set(suffix) then begin

if klip then begin
;PSF subtraction and final processing with fake planet injections
	if proc_inj and yj then reduce_ifs,object_name,/inj,/yj,reduction_path=reduction_path,/klip,suffix=suffix
	if proc_inj and not yj then reduce_ifs,object_name,/inj,reduction_path=reduction_path,/klip,suffix=suffix

	;PSF subtraction and final processing
	if proc_orig and yj then reduce_ifs,object_name,/yj,reduction_path=reduction_path,/klip,suffix=suffix
	if proc_orig and not yj then reduce_ifs,object_name,reduction_path=reduction_path,/klip,suffix=suffix


endif else begin

;PSF subtraction and final processing with fake planet injections
	if proc_inj and yj then reduce_ifs,object_name,/inj,/yj,reduction_path=reduction_path,suffix=suffix
	if proc_inj and not yj then reduce_ifs,object_name,/inj,reduction_path=reduction_path,suffix=suffix

	;PSF subtraction and final processing
	if proc_orig and yj then reduce_ifs,object_name,/yj,reduction_path=reduction_path,suffix=suffix
	if proc_orig and not yj then reduce_ifs,object_name,reduction_path=reduction_path,suffix=suffix

endelse

endif 

if not keyword_set(suffix) then begin
	
if klip then begin
;PSF subtraction and final processing with fake planet injections
	if proc_inj and yj then reduce_ifs,object_name,/inj,/yj,reduction_path=reduction_path,/klip
	if proc_inj and not yj then reduce_ifs,object_name,/inj,reduction_path=reduction_path,/klip

	;PSF subtraction and final processing
	if proc_orig and yj then reduce_ifs,object_name,/yj,reduction_path=reduction_path,/klip
	if proc_orig and not yj then reduce_ifs,object_name,reduction_path=reduction_path,/klip


endif else begin

;PSF subtraction and final processing with fake planet injections
	if proc_inj and yj then reduce_ifs,object_name,/inj,/yj,reduction_path=reduction_path
	if proc_inj and not yj then reduce_ifs,object_name,/inj,reduction_path=reduction_path

	;PSF subtraction and final processing
	if proc_orig and yj then reduce_ifs,object_name,/yj,reduction_path=reduction_path
	if proc_orig and not yj then reduce_ifs,object_name,reduction_path=reduction_path

endelse

endif


test_file=file_search(input_directory+'IFS/interm/','*',count=testcnt)
	if testcnt lt 1 then reduce_ifs_sizes=0
if reduce_ifs_sizes then begin
	print, 'Zipping raw data and deleting unnecessary files - may take several minutes.'
	file_zip, input_directory+'IFS/raw/', input_directory+'IFS/raw.zip'
	file_delete,input_directory+'IFS/raw/',/recursive
	file_zip, input_directory+'IFS/extra/', input_directory+'IFS/extra.zip'
	file_delete,input_directory+'IFS/extra/',/recursive
	print, 'Deleting intermediate data products...'
	file_delete,input_directory+'IFS/interm/',$
		input_directory+'IFS/calib/',/recursive

	preproc_files=file_search(input_directory+'IFS/products/','*.fits',count=cnt)
	for ii=0,cnt-1 do if strpos(preproc_files[ii],'preproc') ne -1 then file_delete,preproc_files[ii]
	file_delete,input_directory+'IFS/spectra_positions_distortion.fits'
	file_delete,input_directory+'IFS/ifs_instrument_flat.fits'
	;file_delete,input_directory+'IFS/ifs_science_dr.fits'	
	file_delete,input_directory+'IFS/input_image2.fits'
	file_delete,input_directory+'IFS/input_image_row.fits'
	file_delete,input_directory+'IFS/esorex.log'

print, 'Done.'	
endif else print, 'Skipping file size reduction...'



endif ;process_ifs if

print, 'Completed SPHERE reduction on ',object_name,' in ',(systime(/JULIAN)-starttime)*86400./60.,' minutes.'


if single then hak

cgcleanup

endfor


end


