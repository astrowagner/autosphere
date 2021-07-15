
pro reduce_ifs,obj,inj=inj,yj=yj, suffix=suffix,reduction_path=reduction_path, klip=klip

;Version: 2017/11/29

;This routine performs PSF subtraction on SPHERE/IFS data products that are 
;pre-processed using Arthur Vigan's pipeline available at http://astro.vigan.fr/

;Author: Kevin Wagner - University of Arizona - kwagner@as.arizona.edu
;
;INPUTS:
; OBJ: object name, that the routine expects to find a folder for in the reduction path
; REDUCTION_PATH: self explanatory; where your object folders live
;OPTIONAL INPUTS:
; INJ: binary keyword for whether to inject fake planets (defined below)
; SUFFIX: extra suffix to be added to files
; YJ: binary, enable YJ mode? 
; KLIP: binary, enable KLIP ADI+SDI reduction of the data?
; 
;Outputs: None - all work is done on files in the observation directory. 
;
;Modify the lines at the top of this routine to fit your reduction needs. 
;

notify_adi=0
notify_adiklip=0
notify_sdiklip=0
notify_sdiklip_ltfilt=0	;will apply a T-dwarf filter to boost SNR

;Establish paths:
if keyword_set(reduction_path) then root=reduction_path+obj+'/IFS/' else root='/Users/kevinwagner/Data/VLT/SPHERE/'+obj+'/IFS'
klipfolder = 'processed/'	;will store processed files in root+klipfolder

if not keyword_set(yj) then yj=0 	;set to 1 for YJ mode on IFS, or to 0 for YJH mode

if not keyword_set(suffix) then suffix=''


;Reject bad frames: (NOTE: IDL notation is 0 = first frame) 
bads=[]			;manual bad frames
cross_thresh=0.8	;for automatic bad frame recognition set to a value less than 1. Set to 0.0 to keep all frames.

field_stab_clean=0	;useful for tight binaries (<1"), otherwise the cleaning does not work

comb_type ='nw-mean'	;'mean', 'median', 'nw-mean' (Bottom et al. 2017)

;What type of KLIP reduction?
if keyword_set(klip) then adi=1	else adi=0	;run adi?
if keyword_set(klip) then sdi=1	else sdi=0	;run sdi? (will find old data cube if adi has already run)

;Reference Differential Imaging (RDI) with KLIP?
rdi=0	;Will prompt for a data_cube_coro.psf to be dragged into terminal (you will need to put in '' manually). 
	;Overrides ADI- and SDI-klip.
if rdi then begin adi=0 & sdi=0 & endif

rdi_adi=0
rdi_sdi=0	;run sdi on an rdi cube

fresh_sdi=0 		;this uses the raw data_cube_coro in the SDI. Bypasses ADI.


auto_bin=1 ;automatically bin the cube to ~50 temporal elements before ADI

classical_adi = 1 ;will output classical ADI results before running KLIP (useful for quick looks and some disks). 

hyp=0	;hyperklip a particular section?
	pxscale=0.00746 ;arcsec/pixel
	rho=[0.5]/pxscale;	spot location input in a radius in arcsec, then converted to pixels
	phi=[90.] & phi=(phi)*!DTOR;
	spot_radius=7.;	spot-width

annmode=0 & if hyp  then annmode=0
	annmode_inout=[10,40]	;will process only an annulus with these settings (in n_ang segments)
	
	
debug=0			;1 for on; 0 for off. Outputs extra files for development purposes.
				;Leaving this on can eat up a LOT of hard drive space.

destripe_ifs=1	;really useful if doing a fresh_sdi reduction in which the field rotation is also minimal. 
	destripe_iter=3	;will process up to destripe_angle_<destripe_iter>
	destripe_angle=100.
	destripe_angle_2=101.
	destripe_angle_3=100.5
	destripe_angle_4=10.
	destripe_angle_5=99.5
	destripe_angle_6=9.5
	range_sz=140.		;range from center to destripe
	destripe_level=0.0	



;Include forward modelling of fake planets?

if keyword_set(inj) then addplanets=1 else addplanets=0
if addplanets then suffix=strcompress('_inj'+suffix,/rem)


;--------------- uniform grid of planets, can easily define other architectures

theta=148.; -  0.1* runs 
rhop=0.2; + 0.001*runs	;right 0.825, left: 0.830	

planet_r=[rhop,rhop,rhop,rhop,rhop,rhop,rhop,rhop]
planet_theta=[theta,theta+180.,theta+90.,theta-90.,theta+45.,theta+180.+45.,theta+90.+45.,theta-90.+45.]

rplanet=[planet_r,planet_r+0.1,planet_r+0.2,planet_r+0.3,planet_r+0.4,planet_r+0.5,planet_r+0.6]
offang=65.
tplanet=[planet_theta,planet_theta-offang,planet_theta-2.*offang,planet_theta-3.*offang,planet_theta-4.*offang,planet_theta-5.*offang,planet_theta-6.*offang]

rplanet=rplanet/0.00746
tplanet=tplanet*!DTOR

nplanets=n_elements(rplanet)
contrast=fltarr(nplanets)
contrast[*]=5.0E-4

;----------------

identify=0 ;idenfity sources?

;General reduction parameters:
szz=140.			;half size of processed frame - do not change!!
k_adiklip=7;17.		;number of KL basis vectors to retain (7)
k_sdiklip=7;15.		;number of KL basis vectors to retain (7)
k_rdiklip=7;15.		;number of KL basis vectors to retain (7)

wr = 14. 		;Width of annuli in pixels (14) (12 recently)
wr_sdi = 14. 		;Width of annuli in pixels (14) (12 recently)
wr_rdi = 14. 		;Width of annuli in pixels (14) (12 recently)
nrings = fix(szz/wr) 	;Number of annuli (15)
nrings_rdi = fix(szz/wr_rdi) 	;Number of annuli (15)
nrings_sdi = fix(szz/wr_sdi) 	;Number of annuli (15)
n_ang = 6.		;Number of segments (6)
ANGSEP=0.5	;Exclude frames from ADI KLIP whose angular separation are within ANGSEP x FWHM of the target frame
anglemax=360.	;or greater than anglemax from the target frame. 
SDISEP=1.5		;Exclude frames from SDI KLIP that are within SDISEP x FWHM of the target frame
filter=11. 		;high pass filter width (set to 0 or 1 to disable) 
filter_low=0	;smoothing width (set to 0 or 1 to disable)

radmask=1		;masks the region exterior to maskrad and interior to innerrad, useful for full frame SDI
maskrad=1.2/0.00746	;pixels away from center
innerrad=0.08/0.00746

savesplit=0 	;Turn to 1 to also save first and second half of data separately.
				;Useful for determining artifacts / speckles.

;--------------------[         Begin object specific input:         ]-----------------------------;




if obj eq 'S25_2'  or obj eq 'S25' then begin

	annmode=1; & if hyp  then annmode=0
	annmode_inout=[7,50]
	wr=30

	ANGSEP=0.5	;Exclude frames from ADI KLIP whose angular separation are within ANGSEP x FWHM of the target frame
	anglemax=360.	;or greater than anglemax from the target frame. 
	filter=9. 		;high pass filter width (set to 0 or 1 to disable) 
	n_ang=4.

	k_adiklip=2;17.		;number of KL basis vectors to retain (7)
	k_sdiklip=4

	rhop=0.125; + 0.001*runs	;right 0.825, left: 0.830	
	if obj eq 'S25' then rhop=0.15

	planet_r=[rhop,rhop,rhop]
	planet_theta=[-90.,90.,180.]+90.-35.;-130.
	if obj eq 'S25' then planet_theta=[-90.,90.,180.]+90.-35.-20.
	rplanet=planet_r
	tplanet=planet_theta;
	rplanet=rplanet/0.00746
	tplanet=tplanet*!DTOR

	nplanets=n_elements(rplanet)
	contrast=fltarr(nplanets)
	contrast[*]=0.8E-3
	
	cross_thresh=0.95

	classical_adi=1
	radmask=0	

	adi=1
	sdi=1	
endif




if obj eq 'HD131488' then begin

	filter=0

endif

	if obj eq 'MWC758' then begin 
		annmode=1 & annmode_inout=[5,110] 
		filter=11.
		n_ang=1

		
		theta=300.; -  0.1* runs 
		rhop=0.6; + 0.001*runs	;right 0.825, left: 0.830	

		planet_r=[rhop,rhop+0.1]
		planet_theta=[theta,theta+180.]

		rplanet=[planet_r]
		tplanet=[planet_theta]

		rplanet=rplanet/0.00746
		tplanet=tplanet*!DTOR

		nplanets=n_elements(rplanet)
		contrast=fltarr(nplanets)
		contrast[*]=1.0E-5

	endif

if obj eq 'S7' or obj eq 'S7_2' then begin

	annmode=1 & annmode_inout=[0,30]

endif

if obj eq 'S26' then bads=[1:5]-1	;burn in effect
if obj eq 'S34' then bads=[1:10]-1	;burn in effect
if obj eq 'S35' then bads=[1:10]-1	;burn in effect
if obj eq 'S38' then bads=[48:64]-1	
if obj eq 'S44' then bads=[46,[53:59]]-1	
if obj eq 'S8_2' then bads=[1:10]-1	;burn in effect
if obj eq 'S7_2' then bads=[1:21]-1	;burn in effect
if obj eq 'S22_11' then bads=[1:11]-1	;burn in effect

if obj eq 'S13_2' then cross_thresh=0.6
if obj eq 'S2' then cross_thresh=0.6
if obj eq 'S14' then cross_thresh=0.
if obj eq 'S20' then cross_thresh=0.95
if obj eq 'S33' then cross_thresh=0.8
if obj eq 'S34' then cross_thresh=0.0	;includes all frames
if obj eq 'S38' then cross_thresh=0.7
if obj eq 'HIP79860_C' then cross_thresh=0.	;include all

if obj eq 'S8' then field_stab_clean=1	;binary
if obj eq 'S8_2' then field_stab_clean=1	;binary
if obj eq 'S14' then field_stab_clean=1	;binary
if obj eq 'S14_2' then field_stab_clean=1	;binary
if obj eq 'S14_C2' then field_stab_clean=1	;binary
if obj eq 'S30' or obj eq 'HIP67036_C' or obj eq 'HIP67036_X' then field_stab_clean=1	;binary
if obj eq 'S33' then field_stab_clean=1	;binary
if obj eq 'S34' then field_stab_clean=1	;binary
if obj eq 'S38' then field_stab_clean=1	;binary
if obj eq 'S50' then field_stab_clean=1	;binary
if obj eq 'S61' then field_stab_clean=1	;binary
if obj eq 'S12' then field_stab_clean=1	;binary
if obj eq 'S12_2' then field_stab_clean=1	;binary
;if obj eq 'S14' then comb_type='median' ;for bright binaries a median combination can work better
if obj eq 'S20' then comb_type='median' ;AO errors leave lots of speckles. Median works better here.
if obj eq 'S21' then comb_type='median' ;AO errors leave lots of speckles. Median works better here.
	;leaving off for consistent reduction

if obj eq 'S46' then fresh_sdi=1





if obj eq 'S45' then begin

	annmode=0; & if hyp  then annmode=0
	annmode_inout=[7,40]

	ANGSEP=0.25	;Exclude frames from ADI KLIP whose angular separation are within ANGSEP x FWHM of the target frame
	anglemax=20.	;or greater than anglemax from the target frame. 
	filter=11. 		;high pass filter width (set to 0 or 1 to disable) 
	n_ang=6.
	theta=148.; -  0.1* runs 
	rhop=0.2; + 0.001*runs	;right 0.825, left: 0.830	

	planet_r=[rhop,rhop,rhop]
	planet_theta=[90.,270.,0.]-130.

	rplanet=planet_r;[planet_r,planet_r+0.1,planet_r+0.2,planet_r+0.3,planet_r+0.4,planet_r+0.5,planet_r+0.6]
	offang=65.
	tplanet=planet_theta;[planet_theta,planet_theta-offang,planet_theta-2.*offang,planet_theta-3.*offang,planet_theta-4.*offang,planet_theta-5.*offang,planet_theta-6.*offang]

	rplanet=rplanet/0.00746
	tplanet=tplanet*!DTOR

	nplanets=n_elements(rplanet)
	contrast=fltarr(nplanets)
	contrast[*]=0.75E-4


endif

if obj eq '51Eri' then begin
	angsep=0.5
	sdisep=1.
	cross_thresh=0.7
	filter=11.
	destripe_ifs=1
	annmode=1 & annmode_inout=[0.35,0.55]/0.00746
	;anglemax=20.	;or greater than anglemax from the target frame. 

	k_adiklip=7;17.		;number of KL basis vectors to retain (7)
	k_sdiklip=7;15.		;number of KL basis vectors to retain (7)

	notify_adi=0
	notify_adiklip=0
	notify_sdiklip=1 & notify_sdiklkp_ltfilt=0

	auto_bin=0

	default_bin=12




	theta=148.; -  0.1* runs 
	rhop=0.45; + 0.001*runs	;right 0.825, left: 0.830	

	planet_r=[rhop,rhop,rhop]
	planet_theta=[90.,270.,0.]-100.

	rplanet=planet_r;[planet_r,planet_r+0.1,planet_r+0.2,planet_r+0.3,planet_r+0.4,planet_r+0.5,planet_r+0.6]
	;offang=65.
	tplanet=planet_theta;[planet_theta,planet_theta-offang,planet_theta-2.*offang,planet_theta-3.*offang,planet_theta-4.*offang,planet_theta-5.*offang,planet_theta-6.*offang]

	rplanet=rplanet/0.00746
	tplanet=tplanet*!DTOR

	nplanets=n_elements(rplanet)
	contrast=fltarr(nplanets)
	contrast[*]=0.12E-5

endif

if obj eq '51Eri_3' then begin
	angsep=0.5
	sdisep=1.
	cross_thresh=0.7
	filter=11.
	destripe_ifs=1
	annmode=1 & annmode_inout=[0.35,0.55]/0.00746
	;anglemax=20.	;or greater than anglemax from the target frame. 

	k_adiklip=7;17.		;number of KL basis vectors to retain (7)
	k_sdiklip=7;15.		;number of KL basis vectors to retain (7)

	notify_adi=0
	notify_adiklip=0
	notify_sdiklip=1

	auto_bin=1
		default_bin=1


	theta=148.; -  0.1* runs 
	rhop=0.45; + 0.001*runs	;right 0.825, left: 0.830	

	planet_r=[rhop,rhop,rhop]
	planet_theta=[90.,270.,0.]-100.

	rplanet=planet_r;[planet_r,planet_r+0.1,planet_r+0.2,planet_r+0.3,planet_r+0.4,planet_r+0.5,planet_r+0.6]
	;offang=65.
	tplanet=planet_theta;[planet_theta,planet_theta-offang,planet_theta-2.*offang,planet_theta-3.*offang,planet_theta-4.*offang,planet_theta-5.*offang,planet_theta-6.*offang]

	rplanet=rplanet/0.00746
	tplanet=tplanet*!DTOR

	nplanets=n_elements(rplanet)
	contrast=fltarr(nplanets)
	contrast[*]=0.12E-5

endif


if not annmode then annmode_inout=[0.,139.]	;set to min and max radial range to be used to cut off edges in find_sources 

;--------------------[ END User input, begin reduction script:      ]-----------------------------;


path=root+'/products/'	;where are the files from AV's scripts?

scicube=readfits(path+'data_cube_coro.fits', scihd)	;load the science cube

if rdi_adi  then scicube=readfits(root+klipfolder+obj+'_ifs_rdiklip_xyln.fits', scihd)  

if rdi_adi  then scicube[where(finite(scicube) eq 0)]=0

if fresh_sdi then begin adi=0 & classical_adi=0 & endif

info=readfits(path+'data_info.fits', infohd)



if rdi  then begin
	refcube='empty'
	read, refcube, PROMPT='Enter reference cube path (data_cube_coro.fits): '
	refcube=readfits(strcompress(string(refcube),/remove_all),refhead)
endif

print,'Size of data cube:',size(scicube)
psfcube=readfits(path+'data_cube_psf.fits', psfhd)
medpsf=median(psfcube,dim=3)
if n_elements(size(medpsf)) gt 5 then medarr,medpsf,medpsf
writefits,path+'data_cube_psf_med.fits',medpsf
xpsf=(size(psfcube))(1)
psfcubezero=psfcube

;prepping PSF for fake planet injections (if enabled)
sky=0.
for x=0,xpsf-1 do begin
	for y=0,xpsf-1 do begin
		if (x-xpsf/2.)^2. + (y-xpsf/2.)^2. gt 10 then psfcubezero[x,y,*]=0
		;this step makes it so that noise in the PSF isn't added to the science region

		if radmask  then begin
		if Sqrt(  (x-(xpsf/2.) )^2. + (y-(xpsf/2.) )^2.  ) gt maskrad $
						then scicube[x,y,*,*]=sky
		if Sqrt(  (x-(xpsf/2.) )^2. + (y-(xpsf/2.) )^2.  ) lt innerrad $
						then scicube[x,y,*,*]=sky
		
		endif
	endfor
endfor



centers=readfits(path+'data_centers.fits', centershd)
scaling=readfits(path+'data_scaling.fits', scalinghd)

print, 'Scaling factors:', scaling

wavelength=readfits(path+'data_wavelength.fits', wavelengthhd)
wave=wavelength
print,wave
;extracting center of frames in each wavelength slice
xc=centers[*,1]
yc=centers[*,0]

xc[*]=139.5
yc[*]=139.5	;these are for Arthur Vigan's output, which is already centered. 


;making a mask identical in size to a single frame and reading sizes
mask=fltarr(2*szz,2*szz)	;making a mask the same size as the science frame
sizes=size(mask)
width=sizes[1]
hwidth=float(width)/2.-1.
specn=sizes[3]
mask=fltarr(width,width) ;make the mask 2-dimensional
mask[*]=0.

;Setting up things for KLIP:
anglemask = jhrangmask(mask)
;Get the dimensions of a single frame
x = float((size(anglemask))(1))
y = float((size(anglemask))(2))
x = min([x,y])
distmask = mask  ; use the first image in cube to create mask
N = x  ; Define the size of N x N output array
dist_circle, distmask, N, [x/2., y/2.]  ; create distance mask
				

;reading angles
info = mrdfits(path+'data_info.fits',1)   ;; read data
;plot,info[*].time,info[*].pa         ;; plot pa values
pa = info[*].pa  
pupoff = info[*].pupoff  

print, pupoff 
;hak

;pupoff is 33.64 after AV routines, should be 33.76 to be consistent with Maire et al. 2016
;IRDIS offset: -135.99+/-0.1
;IRDIS -> IFS offset: 100.45+/-0.1
;TN=-1.75+/-0.08
;thus we add 0.12 to the offset
pupoff=pupoff+0.12

;note this isn't exact for astrometry before July 2016 (see Maire 2016)

print, 'Angles:',pa
angles=pa


;binning cube

;bins to about a size of 20
if auto_bin then begin
	bin=1
	if n_elements(pa) gt 40 then bin=2
	if n_elements(pa) gt 80 then bin=4
	if n_elements(pa) gt 120 then bin=6
	if n_elements(pa) gt 160 then bin=8
	if n_elements(pa) gt 200 then bin=10
endif else bin = default_bin

ncubesize=fix(float(n_elements(angles))/float(bin))
ncube=fltarr(280,280,39,ncubesize)
print, 'New cube size will be: ',size(ncube)

if bin gt 1 then begin
	for kk=0,38 do begin
	binnum=0
	;st=1
	for ii=0.,(size(scicube))(4)-1. do begin
		print, fix(ii+1.) mod fix(bin)
		if fix(ii+1.) mod fix(bin) eq 1 then binned=scicube[*,*,kk,ii] else binned=[ [[binned]], [[ scicube[*,*,kk,ii] ]] ]
		if fix(ii+1.) mod fix(bin) eq 0 then begin
		

			 binned=mean(binned,dim=3)
			

			print, kk, binnum
			ncube[*,*,kk,binnum]=binned
			binnum=binnum+1
		endif
	endfor	
	;print, size(binned_cube)
	;kscicube=binned_cube
	endfor


	st=1
	for ii=0,(size(scicube))(4)-1. do begin
		if fix(ii+1.) mod fix(bin) eq 1 then binned_angle=angles[ii] else binned_angle=binned_angle+angles[ii]
		if fix(ii+1.) mod fix(bin) eq 0 then begin
			if st eq 1 then begin binned_angles = binned_angle/bin & st=0 & endif else binned_angles=[binned_angles,binned_angle/bin]
		endif
	endfor	
	angles=binned_angles
scicube=ncube
endif



;end binning code



;performing auto bad frame recognition here
;first calculate median of science cube in wavelength
for ii=0,n_elements(angles)-1 do begin
	yjhframe=median(reform(scicube[*,*,*,ii]),dim=3)
	if field_stab_clean then yjhframe=rot(yjhframe,-pa[ii]-pupoff[ii],/interp)
	if ii eq 0 then medframes=yjhframe else medframes=[ [[medframes]],[[yjhframe]] ]

endfor
;then compute median of those frames

medarr, medframes, pup_median
;now check frames against the median
for ii=0,n_elements(angles)-1 do begin
	yjhframe=median(reform(scicube[*,*,*,ii]),dim=3)
	if field_stab_clean then yjhframe=rot(yjhframe,-pa[ii]-pupoff[ii],/interp)
	maxcor=max(crosscorr(yjhframe,pup_median))
			print, 'Frame ',ii,' has max correlation value of ', maxcor, ' with pupil median.'
		if ii eq 0 then corrs=maxcor else corrs=[corrs,maxcor]
endfor

;keep those above cross_thresh percentile (defined out of unity), e.g. 0.8=80%

	;lcube=lcube[*,*,where(corrs ge cross_thresh)]
	;rcube=rcube[*,*,where(corrs ge cross_thresh)]
;scicube=scicube[*,*,*,where(corrs ge cross_thresh)]
	

	print, 'Cleaned ', n_elements(where(corrs lt cross_thresh)) ,' frames.'
	print, 'Bad frames:', where(corrs lt cross_thresh)

if  n_elements(where(corrs lt cross_thresh)) ge 0.5*n_elements(corrs) then begin
		print, '!!! Warning! Autoclean has rejected over 50% of frames !!! 
		print, 'You will porbably want to go back to the reduce_ifs.pro routine and change the cross_thresh parameter to a lower value (0 will accept all frames).'
		print, 'If this is a binary star, this may have happened because of comparing single frames (with a bright star) to a median-combined pupil containing a smeared image of the star. Try enabling field_stab_clean in this case.'
		print, 'It is recommended to apply the /skip keyword next run to reduce tiem spet on initial steps.'
		print, 'Hit any key to continue anyways, or ESC to stop the reduction.'
		hak
	endif

goods=fltarr((size(angles))(1))
goods[*]=1
goods[bads]=0

goods[where(corrs lt cross_thresh)]=0

;print,'offset:',pupoff

angles=angles[where(goods )]
pa=angles
pupoff=pupoff[where(goods )]
scicube=scicube[*,*,*,where(goods )]

cubesize=size(scicube)
cnt=cubesize[4]
print, size(angles)
print, size(scicube)
print, 'Proceeding with ',cnt,' frames.'
fullcube=fltarr(2*szz,2*szz,39,cnt)
fullcubeps=fullcube	;pupil-stabilized cube to be used for SDI
fullcubeps[*]=0.


bigarr=fltarr(4.*szz,4.*szz,39,cnt)
bigpsfarr=fltarr(4.*szz,4.*szz,39)

print, size(bigpsfarr[szz:3.*szz-1,szz:3.*szz-1,*])
print, size(psfcubezero)

;use the first PSF, if there is more than one PSF cube in the file (i.e. a 4D cube)
if n_elements(psfcubezero) gt 6 then psfcubezero=reform(psfcubezero[*,*,*,0])
bigpsfarr[szz:3.*szz-1,szz:3.*szz-1,*]=reform(psfcubezero)

for k=0,38 do begin
	for ii=0, cnt-1 do begin	
		
		
	if addplanets  then begin



			;print, size(bigarr[szz:3.*szz-1,szz:3.*szz-1,k,ii])
			;print, size(scicube[*,*,k,ii])

			bigarr[szz:3.*szz-1,szz:3.*szz-1,k,ii]=scicube[*,*,k,ii]



 			
 			bigarr[*,*,k,ii]=rot(bigarr[*,*,k,ii],-pa[ii]-pupoff[ii],/INTERP)
 				
 				
 					xcc=140.
 					ycc=140.
 					if k eq 0 and ii eq 0 then print, $
 						'-------- Adding synthetic planets -------------'
 			
 				for jj=0,nplanets-1 do begin 

				print, 'Injecting planet # ', jj,k,ii, ' out of ', nplanets-1, 38, cnt-1				
 			 				
 				bigarr[*,*,k,ii]= bigarr[*,*,k,ii] + $
 					fshift(bigpsfarr[*,*,k], -rplanet[jj]*cos(tplanet[jj]),-rplanet[jj]*sin(tplanet[jj]))*contrast[jj]

				 print,'Injection coordinates: ',135.-rplanet[jj]*cos(tplanet[jj]),135.-rplanet[jj]*sin(tplanet[jj])
 				
 				endfor
				if ii eq 0 then print, 'Please wait...'
 			 	bigarr[*,*,k,ii]= rot(bigarr[*,*,k,ii],pa[ii]+pupoff[ii],/INTERP)
 			 	scicube[*,*,k,ii]=bigarr[szz:3.*szz-1,szz:3.*szz-1,k,ii]
				
 	endif ;addplanets if
		
		;spatial filtering:
		if filter_low gt 1 then scicube[*,*,k,ii]=smooth(scicube[*,*,k,ii],filter_low)
		if filter gt 1 then scicube[*,*,k,ii]=scicube[*,*,k,ii]-smooth(scicube[*,*,k,ii],filter)
		
	if destripe_ifs  then begin
	
		if k eq 0 and ii eq 0 then print, 'Destriping (may take some time...)'
		if k eq 0 and ii eq 0 then print, 'Cube sizes:', size(scicube)
		
		destripe_range=[xpsf/2.-range_sz:xpsf/2.+range_sz-1.]
		
 			if destripe_iter le 1 then $
 			scicube[destripe_range,destripe_range,k,ii]=destripe(scicube[destripe_range,destripe_range,k,ii],$
 			 destripe_angle, clip_level=destripe_level, /nodisp)
 			 
 			if destripe_iter le 2 then $
 			scicube[destripe_range,destripe_range,k,ii]=destripe(scicube[destripe_range,destripe_range,k,ii], $
 			destripe_angle_2, clip_level=destripe_level, /nodisp)
 			
 			if destripe_iter le 3 then $
 			scicube[destripe_range,destripe_range,k,ii]=destripe(scicube[destripe_range,destripe_range,k,ii], $
 			destripe_angle_3, clip_level=destripe_level, /nodisp)
 			
 			if destripe_iter le 4 then $
 			scicube[destripe_range,destripe_range,k,ii]=destripe(scicube[destripe_range,destripe_range,k,ii], $
 			destripe_angle_4, clip_level=destripe_level, /nodisp)
 			
 			if destripe_iter le 5 then $
 			scicube[destripe_range,destripe_range,k,ii]=destripe(scicube[destripe_range,destripe_range,k,ii], $
 			destripe_angle_5, clip_level=destripe_level, /nodisp)
 			
 			if destripe_iter le 6 then $
 			scicube[destripe_range,destripe_range,k,ii]=destripe(scicube[destripe_range,destripe_range,k,ii], $
 			destripe_angle_6, clip_level=destripe_level, /nodisp)
 	endif
 					
 					
 	;repeating for reference cube
 if rdi  then begin
 	;spatial filtering:
		if filter_low gt 1 then begin
			if ii le (size(refcube))(4)-1 then refcube[*,*,k,ii]=smooth(refcube[*,*,k,ii],filter_low)
		endif
		
		if filter gt 1 then begin
			if ii le (size(refcube))(4)-1 then refcube[*,*,k,ii]=refcube[*,*,k,ii]-smooth(refcube[*,*,k,ii],filter)
		endif
		
		
if destripe_ifs  then begin
	if k eq 0 and ii eq 0 then print, 'Destriping reference cube (may take some time...)'
	if k eq 0 and ii eq 0 then print, 'Cube sizes:', size(refcube)

destripe_range=[xpsf/2.-range_sz:xpsf/2.+range_sz-1.]

	if ii lt ((size(refcube))(4)) then begin
 			if destripe_iter le 1 then	refcube[destripe_range,destripe_range,k,ii]= $
 			destripe(refcube[destripe_range,destripe_range,k,ii], destripe_angle, clip_level=destripe_level, /nodisp)
 			if destripe_iter le 2 then $
 			refcube[destripe_range,destripe_range,k,ii]=destripe(refcube[destripe_range,destripe_range,k,ii], $
 			destripe_angle_2, clip_level=destripe_level, /nodisp)
 			
 			if destripe_iter le 3 then $
 			refcube[destripe_range,destripe_range,k,ii]=destripe(refcube[destripe_range,destripe_range,k,ii], $
 			destripe_angle_3, clip_level=destripe_level, /nodisp)
 			
 			if destripe_iter le 4 then refcube[destripe_range,destripe_range,k,ii]= $
 			destripe(refcube[destripe_range,destripe_range,k,ii], destripe_angle_4, clip_level=destripe_level, /nodisp)
 			
 			if destripe_iter le 5 then refcube[destripe_range,destripe_range,k,ii]= $
 			destripe(refcube[destripe_range,destripe_range,k,ii], destripe_angle_5, clip_level=destripe_level, /nodisp)
 			
 			if destripe_iter le 6 then refcube[destripe_range,destripe_range,k,ii]= $
 			destripe(refcube[destripe_range,destripe_range,k,ii], destripe_angle_6, clip_level=destripe_level, /nodisp)
 	endif
 endif
 	
endif ;rdi if

	endfor
		

endfor
		;derotates the frames and median combines for each band
		derotrawcube=scicube
		adicuberaw=scicube
		for k=0,38 do begin
			;make a median in wavelength k
			medarr,reform(adicuberaw[*,*,k,*]),pupilk
			for ii=0, cnt-1 do adicuberaw[*,*,k,ii]=adicuberaw[*,*,k,ii]-pupilk
			for ii=0, cnt-1 do adicuberaw[*,*,k,ii]=rot(adicuberaw[*,*,k,ii],-pa[ii]-pupoff[ii],/interp)
			for ii=0, cnt-1 do derotrawcube[*,*,k,ii]=rot(derotrawcube[*,*,k,ii],-pa[ii]-pupoff[ii],/interp)
			derotrawcubek=reform(derotrawcube[*,*,k,*])
			adicubek=reform(adicuberaw[*,*,k,*])

			if comb_type eq 'median' then medarr, derotrawcubek, kframemedraw
			if comb_type eq 'mean' then kframemedraw=mean(derotrawcubek,dim=3)
			if comb_type eq 'nw-mean' then kframemedraw=nw_ang_comb(derotrawcubek,pa+pupoff)
			if comb_type eq 'median' then medarr, adicubek, kadiframe
			if comb_type eq 'mean' then kadiframe=mean(adicubek,dim=3)
			if comb_type eq 'nw-mean' then kadiframe=nw_ang_comb(adicubek, pa+pupoff)			
			if k eq 0 then adisum=kadiframe else adisum=adisum+kadiframe
			if k eq 0 then adicube=kadiframe else adicube=[[[adicube]],[[kadiframe]]]
			

			if k eq 0 then sumraw=kframemedraw else sumraw=sumraw+kframemedraw
			if k eq 0 then sumrawcube=kframemedraw else sumrawcube=[[[sumrawcube]],[[kframemedraw]]]
		endfor
		sumraw=sumraw/39.
		
		writefits, STRCOMPRESS(root+''+klipfolder+''+obj+'_ifs_raw_xyl'+suffix+'.fits'), sumrawcube, scihd
		if yj then writefits, STRCOMPRESS(root+''+klipfolder+''+obj+'_ifs_raw_YJ'+suffix+'.fits'), sumraw, scihd else writefits, STRCOMPRESS(root+''+klipfolder+''+obj+'_ifs_raw_YJH'+suffix+'.fits'), sumraw, scihd 
		
		if classical_adi  then begin
		
		writefits, STRCOMPRESS(root+''+klipfolder+''+obj+'_ifs_adi_xyl'+suffix+'.fits'), adicube, scihd
		if yj then writefits, STRCOMPRESS(root+''+klipfolder+''+obj+'_ifs_adi_YJ'+suffix+'.fits'), adisum, scihd else writefits, STRCOMPRESS(root+''+klipfolder+''+obj+'_ifs_adi_YJH'+suffix+'.fits'), adisum, scihd

		;add further YJH if statements here once figuring out the ranges in YJ mode
		
		
		if yj then yadicube = adicube[*,*,2:15] else yadicube=adicube[*,*,0:8]
		
		;medarr, yadicube, yadi
		yadi=mean(yadicube,dimension=3)
		writefits, root+''+klipfolder+''+obj+'_ifs_adi_Y'+suffix+'.fits', yadi, head

		if yj then jadicube = adicube[*,*,17:38] else jadicube=adicube[*,*,10:21]
		;medarr, jadicube, jadi
		jadi=mean(jadicube,dimension=3)

		writefits, root+''+klipfolder+''+obj+'_ifs_adi_J'+suffix+'.fits', jadi, head
	
		if not yj then hadicube=adicube[*,*,28:36]
		;medarr, hadicube, hadi
		if not yj then hadi=mean(hadicube,dimension=3)

		if not yj then writefits, root+''+klipfolder+''+obj+'_ifs_adi_H'+suffix+'.fits', hadi, head
		
		;make a convolved cube
		conv_xylcube=adicube
		for ccc=0,38 do begin
	
			sz=280.-1.
			width=(wavelength[ccc]*1E-6) / (8.2) * 206265. / 0.00746
			print, 'PSF Width: ',width
			PSF = psf_Gaussian(NPIX=sz, FWHM=[width,width])
			PSFN = PSF/MAX(PSF)
			conv_xylcube[where(finite(conv_xylcube) ne 1)]=0.
			conv_xylcube[*,*,ccc] = convolve(conv_xylcube[*,*,ccc], PSFN)

		endfor

		writefits, root+''+klipfolder+''+obj+'_ifs_adi_xyl_conv'+suffix+'.fits', conv_xylcube, head

		if yj then yadicube=conv_xylcube[*,*,2:15] else yadicube=conv_xylcube[*,*,0:8]
		;medarr, yadicube, yadi
		yadi=mean(yadicube,dimension=3)
		writefits, root+''+klipfolder+''+obj+'_ifs_adi_Y_conv'+suffix+'.fits', yadi, head

		if yj then jadicube=conv_xylcube[*,*,17:38] else jadicube=conv_xylcube[*,*,10:21]
		;medarr, jadicube, jadi
		jadi=mean(jadicube,dimension=3)

		writefits, root+''+klipfolder+''+obj+'_ifs_adi_J_conv'+suffix+'.fits', jadi, head
	
		if not yj then hadicube=conv_xylcube[*,*,28:36]
		;medarr, hadicube, hadi
		if not yj then hadi=mean(hadicube,dimension=3)

		if not yj then writefits, root+''+klipfolder+''+obj+'_ifs_adi_H_conv'+suffix+'.fits', hadi, head


		;medarr,conv_xylcube, yjhadi
		yjhadi=mean(conv_xylcube,dimension=3)

		if not yj then writefits, root+''+klipfolder+''+obj+'_ifs_adi_YJH_conv'+suffix+'.fits', yjhadi, head else writefits, root+''+klipfolder+''+obj+'_ifs_adi_YJ_conv'+suffix+'.fits', yjhadi, head


cgcleanup
		print,strcompress(root+klipfolder+''+obj+'_ifs_adi_YJH'+suffix+'.fits',/rem)
		if identify and not yj then find_sources, strcompress(root+klipfolder+''+obj+'_ifs_adi_YJH'+suffix+'.fits',/rem),sigma=5.0, reference=root+'products/'+'data_cube_psf_med.fits',FWHM=((1.3*1.0E-6)/(8.4))*206265./0.00746,platescale=0.00746,clip_level=0.03,stretch='linear',correction_factor=1.0,inrad=annmode_inout[0]+4, outrad=annmode_inout[1]-4
		if identify and yj then find_sources, strcompress(root+klipfolder+''+obj+'_ifs_adi_YJ'+suffix+'.fits',/rem),sigma=5.0, reference=root+'products/'+'data_cube_psf_med.fits',FWHM=((1.3*1.0E-6)/(8.4))*206265./0.00746,platescale=0.00746,clip_level=0.03,stretch='linear',correction_factor=1.0,inrad=annmode_inout[0]+4, outrad=annmode_inout[1]-4




if notify_adi and not YJ then begin


body='Classical IFS ADI reduction with high pass filter width = '+string(filter)+' --- Message automatically generated on ' + String(systime())
notify_string='echo "'+body+'" | mail -s "ZONASPHERE '+ obj + ' IFS Data Reduction Results" kwagner@as.arizona.edu -A '+strcompress(root+klipfolder+obj+'_ifs_adi_YJH'+suffix+'.fits_source_map.png',/rem) + $
' -A '+strcompress(root+klipfolder+obj+'_ifs_adi_YJH'+suffix+'.fits_sources.txt',/rem) + $
' -A '+strcompress(root+klipfolder+obj+'_ifs_adi_YJH'+suffix+'.fits_radial_sensitivity.png',/rem) + $
' -A '+strcompress(root+klipfolder+obj+'_ifs_adi_YJH'+suffix+'.fits_radial_sensitivity_2.0_arcsec.png',/rem) 

spawn, notify_string


sources_rho=0. & sources_x=0. & sources_y=0. & candidates=!null ;predefining to keep things running
restore,root+klipfolder+obj+'_ifs_adi_YJH'+suffix+'.fits_sources.sav'

;select candidates based on contrast and separation

candidates = where(sources_rho le 2.0 and sources_rho ge 0.)

;candidates=indgen(7) 

if n_elements(candidates) eq 1 then if candidates eq -1 then candidates=!null
if sources_x[0] eq 0. or sources_y[0] eq 0. then candidates=!null


if candidates ne !NULL then begin



candidates_contrast= sources_contrast[candidates]
candidates_x=sources_x[candidates]
candidates_y=sources_y[candidates]
candidates_rho=sources_rho[candidates]
candidates_theta=sources_theta[candidates]
candidates_sigma=sources_sigma[candidates]


smalls=readfits(root+klipfolder+obj+'_ifs_adi_YJH'+suffix+'.fits_sources.fits')
smalls=smalls[*,*,candidates]

rho_string='' & sigma_string='' & contrast_string='' & attach_string = '' & theta_string='' & candidates_string=''

cgcleanup

for mmm=0,n_elements(candidates)-1 do begin

plotcolor='Snow'
	cand=smalls[*,*,mmm]
	cand[where(finite(cand) ne 1)]=0.
	rho_string=rho_string+string(candidates_rho[mmm])
	theta_string=theta_string+string(candidates_theta[mmm])
	sigma_string=sigma_string+string(candidates_sigma[mmm])
	contrast_string=contrast_string+string(candidates_contrast[mmm])
	candidates_string=candidates_string+string(candidates[mmm]+1)
	
	cgps_open,root+klipfolder+obj+'_ifs_adi_YJH'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png'


   cgLoadCT, 1;, NColors=100,BOTTOM=100;,BOTTOM=100
	cgDisplay,40,40,/device

	;cgImage,cand;,stretch=stretch,xrange=[0,(size(smalls))(2)],yrange=[0,(size(smalls))(2)],/keep_aspect_ratio,minvalue=min(cand),maxvalue=max(cand),mean=mean(cand);,position=[0,0,1,1]

cgImage,cand,xrange=[0,39],yrange=[0,39],/keep_aspect_ratio;=[0,0,1,1];,minvalue=min(cand),maxvalue=max(cand),mean=mean(cand);,position=[0,0,1,1]

	cgps_close,/png

	cgcleanup

	attach_string = attach_string+ ' -A '+strcompress(root+klipfolder+obj+'_ifs_adi_YJH'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png',/rem)
endfor


;for cc=0,n_elements(candidates)-1 do begin

	body='ALERT!'+string(n_elements(candidates))+' Candidate exoplanets discovered around '+obj+' in cADI IFS reduction with significance = '+sigma_string+' contrast = '+contrast_string+', rho='+rho_string+', theta='+theta_string+' source numbers ='+candidates_string+' --- Message automatically generated by ZONASPHERE on ' + String(systime())

	notify_string='echo "'+body+'" | mail -s "ALERT!!! '+string(n_elements(candidates)) +' Exoplanet Candidates Detected Around '+ obj + '!!!" kwagner@as.arizona.edu -A '+strcompress(root+klipfolder+obj+'_ifs_adi_YJH'+suffix+'.fits_source_map.png',/rem) +' -A '+strcompress(root+klipfolder+obj+'_ifs_adi_YJH'+suffix+'.fits_sources.txt',/rem) + ' -A '+strcompress(root+klipfolder+obj+'_ifs_adi_YJH'+suffix+'.fits_SNR_map.png',/rem)  + attach_string

	spawn, notify_string 

for mmm=0,n_elements(candidates)-1 do begin
	file = file_Search(strcompress(root+klipfolder+obj+'_ifs_adi_YJH'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png',/rem),count=cnt)
	if cnt gt 0 then file_delete,file
endfor


;endfor
endif
endif

if notify_adi and YJ then begin


body='Classical IFS ADI reduction with high pass filter width = '+string(filter)+' --- Message automatically generated on ' + String(systime())
notify_string='echo "'+body+'" | mail -s "ZONASPHERE '+ obj + ' IFS Data Reduction Results" kwagner@as.arizona.edu -A '+strcompress(root+klipfolder+obj+'_ifs_adi_YJ'+suffix+'.fits_source_map.png',/rem) + $
' -A '+strcompress(root+klipfolder+obj+'_ifs_adi_YJ'+suffix+'.fits_sources.txt',/rem) + $
' -A '+strcompress(root+klipfolder+obj+'_ifs_adi_YJ'+suffix+'.fits_radial_sensitivity.png',/rem) + $
' -A '+strcompress(root+klipfolder+obj+'_ifs_adi_YJ'+suffix+'.fits_radial_sensitivity_2.0_arcsec.png',/rem) 

spawn, notify_string


sources_rho=0. & sources_x=0. & sources_y=0. & candidates=!null ;predefining to keep things running
restore,root+klipfolder+obj+'_ifs_adi_YJ'+suffix+'.fits_sources.sav'

;select candidates based on contrast and separation

candidates = where(sources_rho le 2.0 and sources_rho ge 0.)

;candidates=indgen(7) 

if n_elements(candidates) eq 1 then if candidates eq -1 then candidates=!null
if sources_x[0] eq 0. or sources_y[0] eq 0. then candidates=!null


if candidates ne !NULL then begin



candidates_contrast= sources_contrast[candidates]
candidates_x=sources_x[candidates]
candidates_y=sources_y[candidates]
candidates_rho=sources_rho[candidates]
candidates_theta=sources_theta[candidates]
candidates_sigma=sources_sigma[candidates]


smalls=readfits(root+klipfolder+obj+'_ifs_adi_YJ'+suffix+'.fits_sources.fits')
smalls=smalls[*,*,candidates]

rho_string='' & sigma_string='' & contrast_string='' & attach_string = '' & theta_string='' & candidates_string=''

cgcleanup

for mmm=0,n_elements(candidates)-1 do begin

plotcolor='Snow'
	cand=smalls[*,*,mmm]
	cand[where(finite(cand) ne 1)]=0.
	rho_string=rho_string+string(candidates_rho[mmm])
	theta_string=theta_string+string(candidates_theta[mmm])
	sigma_string=sigma_string+string(candidates_sigma[mmm])
	contrast_string=contrast_string+string(candidates_contrast[mmm])
	candidates_string=candidates_string+string(candidates[mmm]+1)
	
	cgps_open,root+klipfolder+obj+'_ifs_adi_YJ'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png'


   cgLoadCT, 1;, NColors=100,BOTTOM=100;,BOTTOM=100
	cgDisplay,40,40,/device

	;cgImage,cand;,stretch=stretch,xrange=[0,(size(smalls))(2)],yrange=[0,(size(smalls))(2)],/keep_aspect_ratio,minvalue=min(cand),maxvalue=max(cand),mean=mean(cand);,position=[0,0,1,1]

cgImage,cand,xrange=[0,39],yrange=[0,39],/keep_aspect_ratio;=[0,0,1,1];,minvalue=min(cand),maxvalue=max(cand),mean=mean(cand);,position=[0,0,1,1]

	cgps_close,/png

	cgcleanup

	attach_string = attach_string+ ' -A '+strcompress(root+klipfolder+obj+'_ifs_adi_YJ'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png',/rem)
endfor


;for cc=0,n_elements(candidates)-1 do begin

	body='ALERT!'+string(n_elements(candidates))+' Candidate exoplanets discovered around '+obj+' in cADI IFS reduction with significance = '+sigma_string+' contrast = '+contrast_string+', rho='+rho_string+', theta='+theta_string+' source numbers ='+candidates_string+' --- Message automatically generated by ZONASPHERE on ' + String(systime())

	notify_string='echo "'+body+'" | mail -s "ALERT!!! '+string(n_elements(candidates)) +' Exoplanet Candidates Detected Around '+ obj + '!!!" kwagner@as.arizona.edu -A '+strcompress(root+klipfolder+obj+'_ifs_adi_YJ'+suffix+'.fits_source_map.png',/rem) +' -A '+strcompress(root+klipfolder+obj+'_ifs_adi_YJ'+suffix+'.fits_sources.txt',/rem) + ' -A '+strcompress(root+klipfolder+obj+'_ifs_adi_YJ'+suffix+'.fits_SNR_map.png',/rem)  + attach_string

	spawn, notify_string 

for mmm=0,n_elements(candidates)-1 do file_delete,strcompress(root+klipfolder+obj+'_ifs_adi_YJ'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png',/rem)



;endfor
endif
endif


		endif ;classical_adi if


		if adi  then begin
for k=0,38 do begin ;loops through the spectral positions
	split=cnt/2-1	;determines how to split data after klip (if desired)


	;making a cube of dits in the wavelength channel (stored as kcube)
	for nc=0,cnt-1 do begin

	;first we shift the frame so that the center appears exactly in between the middle
	;two pixels, and then we will work with an even number of pixels

	indxck=fix(xc[k])
	indyck=fix(yc[k])
	
	remxck=xc[k]-float(indxck)
	remyck=yc[k]-float(indyck)
	
	print, 'Centering on pixel grid by ',-remxck+0.5,-remyck+0.5
	scicube(*,*,k,nc)=fshift(scicube(*,*,k,nc),-remxck+0.5,-remyck+0.5)
	;now the center lies exactly at indxck+0.5,indyck+0.5 
	
	
	frame=scicube[indxck-szz+1:indxck+szz,indyck-szz+1:indyck+szz,k,nc]

		;UPDATING HEADER 
		FXADDPAR, scihd, 'Filter Width', filter,'High-pass filter pixel width.'
		FXADDPAR, scihd, 'ADI K_adiklip', k_adiklip,'ADI K_adiklip'
		FXADDPAR, scihd, 'SDI K_adiklip', k_sdiklip,'SDI K_adiklip'
		FXADDPAR, scihd, 'ADI ANNULI WIDTH', wr,'ADI ANNULAR WIDTH'
		FXADDPAR, scihd, 'SDI ANNULI WIDTH', wr_sdi,'SDI ANNULI WIDTH'
		FXADDPAR, scihd, 'ADI N RINGS', NRINGS,'ADI N RINGS'
		FXADDPAR, scihd, 'SDI N RINGS', NRINGS_SDI,'SDI NRINGS'
		FXADDPAR, scihd, 'N ANG SEG', N_ANG,'N ANG SEG'
		FXADDPAR, scihd, 'ADI ANGSEP', angsep,'ADI ANGSEP'
		FXADDPAR, scihd, 'SDI WAVESEP', SDISEP,'SDI WAVESEP'
		FXADDPAR, scihd, 'Planets Injected?', addplanets,'Were artifical planets injected?'
		FXADDPAR, scihd, 'Injected Contrast:', contrast[0],'Injected artificial planet contrast with star (first planet).'
		FXADDPAR, scihd, 'NPlanets', nplanets,'Number of artificial planet injections.'
		head=scihd
		
		for ij =0, nplanets-1 do FXADDPAR, scihd, strcompress('InjP'+string(ij+1)+'-Con'), $
				contrast[ij], strcompress('Planet '+string(ij+1)+' contrast')
		
		for ij =0, nplanets-1 do FXADDPAR, scihd, strcompress('InjP'+string(ij+1)+'-R'), $
				rplanet[ij], strcompress('Planet '+string(ij+1)+' rho (px)')
		
		for ij =0, nplanets-1 do FXADDPAR, scihd, strcompress('InjP'+string(ij+1)+'-T'), $
				tplanet[ij]/!DTOR, strcompress('Planet '+string(ij+1)+' theta (deg E of N)')

		derotframe=rot(frame,-pa[nc]-pupoff[nc],1.0,/INTERP)	;returns to pupil-stabilized 
	
	print, 'Frame size:',size(frame)
	if nc eq 0 then kcube=frame else kcube=[[[kcube]],[[frame]]]
	if nc eq 0 then derotkcube=derotframe else derotkcube=[[[derotkcube]],[[derotframe]]]


endfor ;done making cube

;splitting cube
if savesplit  then begin
	derotkcube1=derotkcube[*,*,0:split]
	derotkcube2=derotkcube[*,*,split:cnt-1]

	medarr,derotkcube1,korigmed1
	writefits, STRCOMPRESS(root+''+klipfolder+''+obj+'_ifs_raw_median1_k'+String(k+1)+'prelim.fits',/rem), korigmed1, scihd

	medarr,derotkcube2,korigmed2
	writefits, STRCOMPRESS(root+''+klipfolder+''+obj+'_ifs_raw_median2_k'+String(k+1)+'prelim.fits',/rem), korigmed2, scihd

	if k eq 0 then medcube1=korigmed1 else medcube1 = [  [[medcube1]],[[korigmed1]]  ]
	if k eq 0 then medcube2=korigmed2 else medcube2 = [  [[medcube2]],[[korigmed2]]  ]
endif ;savesplit if

print, 'Size of k cube:',size(kcube)

medarr,derotkcube,korigmed
writefits, strcompress(root+''+klipfolder+''+obj+'_ifs_adiklip_input_k'+String(k+1)+'prelim.fits',/rem), kcube ;will be deleted later
if k eq 0 then medcube=korigmed else medcube = [  [[medcube]],[[korigmed]]  ]

	;this loops through the frames in each particular band and performs the KLIP routine
	for ii=0, cnt-1 do begin
		;print, 'input size:', size(kcube)
		print, '--[ KLIPing Image ', ii, ' ]--'
		print, 'Object: ',obj
		print, 'Wavelength:', wavelength[k],'/',wavelength[38]
			
  		if hyp  then  klip = adiklip(kcube, k_adiklip, target=ii, $
  		 anglemask=anglemask, distmask=distmask, posang=-angles,  $
 		 	wl=wave[k],diam=8.,pixelscale=0.00746,angsep=angsep,anglemax=anglemax,obj=obj,nrings=nrings,wr=wr,$
  		  n_ang =n_ang,spot_radius=spot_radius,rho=rho,phi=phi,/hyper) 
  		 
  		if annmode  then klip = adiklip(kcube, k_adiklip, target=ii, anglemask=anglemask, distmask=distmask,$
  		  posang=-angles,  wl=wave[k],diam=8.,pixelscale=0.00746,angsep=angsep,anglemax=anglemax,$
  		  obj=obj,nrings=nrings,wr =wr, n_ang =n_ang,annmode_inout=annmode_inout)
  		  
  		  if annmode eq 0 and hyp eq 0 then klip = adiklip(kcube, k_adiklip, target=ii, anglemask=anglemask, 	$
  		  distmask=distmask, posang=-angles,  wl=wave[k],diam=8.,pixelscale=0.00746,angsep=angsep,anglemax=anglemax,$
  		  obj=obj,nrings=nrings,wr =wr, n_ang =n_ang)
  		 ;	print, 'output size:', size(klip)

		if ii eq 0 then klipcube=klip else klipcube=[[[klipcube]],[[klip]]]
	endfor
			
		;derotates the frames and median combines for each band
		
		print, size(fullcubeps)
		print, size(klipcube)
		
		for ii=0, cnt-1 do begin
			fullcubeps[*,*,k,ii]=klipcube[*,*,ii]	
 			klipcube[*,*,ii]=rot(klipcube[*,*,ii],-pa[ii]-pupoff[ii],/INTERP)
 			fullcube[*,*,k,ii]=klipcube[*,*,ii]
	
		endfor
		

		if comb_type eq 'median' then medarr, klipcube, kframemed
		if comb_type eq 'mean' then kframemed=mean(klipcube,dim=3)
		if comb_type eq 'nw-mean' then kframemed=nw_ang_comb(klipcube,pa+pupoff)

		
		
		if debug  then medarr, reform(fullcubeps[*,*,k,*]), kframemedps

		if debug  then writefits, root+klipfolder+'/DEBUG_adiklipmedps'+String(k)+'_ouput.fits', kframemedps


		if debug  and k eq 38 then writefits, root+klipfolder+'/DEBUG_adiklip'+String(k)+'_ouput_full.fits', fullcubeps
		if debug  and k eq 38 then writefits, root+klipfolder+'/DEBUG_adiklip'+String(k)+'_ouput_full-derot.fits', fullcube

	for ii=0, cnt-1 do if ii eq 0 then sumk=klipcube[*,*,ii]$
			else sumk=sumk+klipcube[*,*,ii]

	
		for ii=0, split do if ii eq 0 then sumk1=klipcube[*,*,ii]$
			else sumk1=sumk1+klipcube[*,*,ii]
	
		for ii=split,cnt-1. do if ii eq split then sumk2=klipcube[*,*,ii]$
			else sumk2=sumk2+klipcube[*,*,ii]
	
	
	klipcube1=klipcube[*,*,0:split]
	medarr,klipcube1,kframemed1
	klipcube2=klipcube[*,*,split:cnt-1.]
	medarr,klipcube2,kframemed2

	kframemed[where(finite(kframemed) eq 0)]=0.0
		sumk[where(finite(sumk) eq 0)]=0.0
		
			kframemed1[where(finite(kframemed1) eq 0)]=0.0
		sumk1[where(finite(sumk1) eq 0)]=0.0
		
			kframemed2[where(finite(kframemed2) eq 0)]=0.0
		sumk2[where(finite(sumk2) eq 0)]=0.0
				
	sz=280.-1.
	width=(wavelength[k]*1E-6) / (8.2) * 206265. / 0.00746
	print, 'PSF Width: ',width
	PSF = psf_Gaussian(NPIX=sz, FWHM=[width,width])
	PSFN = PSF/MAX(PSF)
	kframemed_conv=kframemed
	kframemed_conv[where(finite(kframemed_conv) ne 1)]=0.
	kframemed_conv = convolve(kframemed_conv, PSFN)

	
	writefits, STRCOMPRESS(root+''+klipfolder+''+obj+'_ifs_adiklip_median_k'+String(k+1)+'prelim.fits',/rem), kframemed, scihd
	
	
	writefits, STRCOMPRESS(root+''+klipfolder+''+obj+'_ifs_adiklip_median_conv_k'+String(k+1)+'prelim.fits',/rem), kframemed_conv, scihd
	
if savesplit  then 	writefits, STRCOMPRESS(root+''+klipfolder+''+obj+'_ifs_adiklip_median1_k'+String(k+1)+'prelim.fits',/rem), kframemed1, scihd
	
if savesplit  then 	writefits, STRCOMPRESS(root+''+klipfolder+''+obj+'_ifs_adiklip_median2_k'+String(k+1)+'prelim.fits',/rem), kframemed2, scihd
	

		if k eq 0 then newcube=fltarr(2*szz,2*szz,39)
	if k eq 0 then newcube[*]=0.
	if k eq 0 then sumcube=newcube
	
	if k eq 0 then newcube1=fltarr(2*szz,2*szz,39)
	if k eq 0 then newcube1[*]=0.
	if k eq 0 then sumcube1=newcube1
	
	if k eq 0 then newcube2=fltarr(2*szz,2*szz,39)
	if k eq 0 then newcube2[*]=0.
	if k eq 0 then sumcube2=newcube2
		
		
		aasz=2*szz
		kframemed=kframemed[0:aasz-1,0:aasz-1]
		sumk=sumk[0:aasz-1,0:aasz-1]
		kframemed1=kframemed1[0:aasz-1,0:aasz-1]
		sumk1=sumk1[0:aasz-1,0:aasz-1]
		kframemed2=kframemed2[0:aasz-1,0:aasz-1]
		sumk2=sumk2[0:aasz-1,0:aasz-1]
		

		print, size(newcube)
				print, size(kframemed)

	newcube[*,*,k]=kframemed
	sumcube[*,*,k]=sumk
	
	newcube1[*,*,k]=kframemed1
	sumcube1[*,*,k]=sumk1
	
	newcube2[*,*,k]=kframemed2
	sumcube2[*,*,k]=sumk2
	
	
	;generating Y image
	if k eq 8 then begin
	

	if yj then yklipcube=newcube[*,*,2:15] else yklipcube=newcube[*,*,0:8]
	sumy=mean(yklipcube,dim=3)
	writefits, root+''+klipfolder+''+obj+'_ifs_adiklip_Y'+suffix+'.fits', sumy

	endif
	


	;generating J image
	if k eq 21 then begin


	if yj then jklipcube=newcube[*,*,17:38] else jklipcube=newcube[*,*,10:21]
	sumj=mean(jklipcube,dim=3)
	writefits, root+''+klipfolder+''+obj+'_ifs_adiklip_J'+suffix+'.fits', sumj
	
	endif
endfor
;----------------ADIKLIP complete

;fullcube contains the raw ADI-processed cube, newcube contains the x,y,39 median-combined cube

;this is the cube of the 39 images, the individual frames may be deleted once this file is generated
writefits, root+''+klipfolder+''+obj+'_ifs_adiklip_xyln.fits', fullcube, scihd
writefits, root+''+klipfolder+''+obj+'_ifs_adiklip_xyln_pupst.fits', fullcubeps, scihd


writefits, root+''+klipfolder+''+obj+'_ifs_adiklip_xyl'+suffix+'.fits', newcube, scihd
;writefits, root+''+klipfolder+''+obj+'_ifs_adiklip_coadd_cube.fits', sumcube, scihd
;writefits, root+''+klipfolder+''+obj+'_ifs_original_cube.fits', medcube, scihd

if savesplit  then begin 
writefits, root+''+klipfolder+''+obj+'_ifs_adiklip_med_cube-1'+suffix+'.fits', newcube1, scihd
writefits, root+''+klipfolder+''+obj+'_ifs_adiklip_coadd_cube-1'+suffix+'.fits', sumcube1, scihd
writefits, root+''+klipfolder+''+obj+'_ifs_original_cube-1'+suffix+'.fits', medcube1, scihd
writefits, root+''+klipfolder+''+obj+'_ifs_adiklip_med_cube-2'+suffix+'.fits', newcube2, scihd
writefits, root+''+klipfolder+''+obj+'_ifs_adiklip_coadd_cube-2'+suffix+'.fits', sumcube2, scihd
writefits, root+''+klipfolder+''+obj+'_ifs_original_cube-2'+suffix+'.fits', medcube2, scihd
endif

;make a convolved cube
conv_xylcube=newcube
for ccc=0,38 do begin
	
	sz=280.-1.
	width=(wavelength[ccc]*1E-6) / (8.2) * 206265. / 0.00746
	print, 'PSF Width: ',width
	PSF = psf_Gaussian(NPIX=sz, FWHM=[width,width])
	PSFN = PSF/MAX(PSF)
	conv_xylcube[where(finite(conv_xylcube) ne 1)]=0.
	conv_xylcube[*,*,ccc] = convolve(conv_xylcube[*,*,ccc], PSFN)

endfor

writefits, root+''+klipfolder+''+obj+'_ifs_adiklip_xyl_conv'+suffix+'.fits', conv_xylcube, head

if yj then yadicube=conv_xylcube[*,*,2:15] else yadicube=conv_xylcube[*,*,0:8]
;medarr, yadicube, yadi
yadi=mean(yadicube,dimension=3)
writefits, root+''+klipfolder+''+obj+'_ifs_adiklip_Y_conv'+suffix+'.fits', yadi, head

if yj then jadicube=conv_xylcube[*,*,17:38] else jadicube=conv_xylcube[*,*,10:21]
;medarr, jadicube, jadi
jadi=mean(jadicube,dimension=3)

writefits, root+''+klipfolder+''+obj+'_ifs_adiklip_J_conv'+suffix+'.fits', jadi, head
	
if not yj then hadicube=conv_xylcube[*,*,28:36]
;medarr, hadicube, hadi
if not yj then hadi=mean(hadicube,dimension=3)

if not yj then writefits, root+''+klipfolder+''+obj+'_ifs_adiklip_H_conv'+suffix+'.fits', hadi, head

;medarr,conv_xylcube, yjhadi
yjhadi=mean(conv_xylcube,dimension=3)

if not yj then writefits, root+''+klipfolder+''+obj+'_ifs_adiklip_YJH_conv'+suffix+'.fits', yjhadi, head else writefits, root+''+klipfolder+''+obj+'_ifs_adiklip_YJ_conv'+suffix+'.fits', yjhadi, head

;stacking images H, YJH

start=28
for num=start,36 do begin
 if num eq start then sumh=newcube[*,*,num] else sumh=sumh+newcube[*,*,num]
endfor
sumh=sumh/9.
if not yj then writefits, root+''+klipfolder+''+obj+'_ifs_adiklip_H'+suffix+'.fits', sumh, head

;h first half
start=28
for num=start,36 do begin
 if num eq start then sumh1=newcube1[*,*,num] else sumh1=sumh1+newcube1[*,*,num]
endfor
sumh1=sumh1/9.
if savesplit and not yj then writefits, root+''+klipfolder+''+obj+'_ifs_adiklip_H-1'+suffix+'.fits', sumh1, head

;second half
start=28
for num=start,36 do begin
 if num eq start then sumh2=newcube2[*,*,num] else sumh2=sumh2+newcube2[*,*,num]
endfor
sumh2=sumh2/9.
if savesplit and not yj then writefits, root+''+klipfolder+''+obj+'_ifs_adiklip_H-2'+suffix+'.fits', sumh2, head
		
start=0
for num=start,38 do begin
 if num eq start then sumyjh=newcube[*,*,num] else sumyjh=sumyjh+newcube[*,*,num]
endfor
sumyjh=sumyjh/39.
if not yj then writefits, root+''+klipfolder+''+obj+'_ifs_adiklip_YJH'+suffix+'.fits', sumyjh, head  else writefits, root+''+klipfolder+''+obj+'_ifs_adiklip_YJ'+suffix+'.fits', sumyjh, head

;writefits, root+''+klipfolder+''+obj+'_ifs_adiklip_xyln.fits', fullcube, head

finalcube=fullcube
yjhcube=[]
for kz=0, (size(finalcube))(4)-1 do begin &	kzcube=reform(finalcube[*,*,*,kz]) & medarr,kzcube,kzmed & yjhcube=[[[yjhcube]],[[kzmed]]] & endfor	
	if not yj then writefits, root+''+klipfolder+''+obj+'_ifs_adiklip_YJH_cube'+suffix+'.fits', yjhcube, head else writefits, root+''+klipfolder+''+obj+'_ifs_adiklip_YJ_cube'+suffix+'.fits', yjhcube, head 

start=0
for num=start,38 do begin
 if num eq start then sumyjh1=newcube1[*,*,num] else sumyjh1=sumyjh1+newcube1[*,*,num]
endfor
sumyjh1=sumyjh1/39.
if savesplit  then writefits, root+''+klipfolder+''+obj+'_ifs_adiklip_BB-1'+suffix+'.fits', sumyjh1, head

start=0
for num=start,38 do begin
 if num eq start then sumyjh2=newcube2[*,*,num] else sumyjh2=sumyjh2+newcube2[*,*,num]
endfor
sumyjh2=sumyjh2/39.
if savesplit  then writefits, root+''+klipfolder+''+obj+'_ifs_adiklip_BB-2'+suffix+'.fits', sumyjh2, head

path=root+''+klipfolder+''
;oldfiles=FILE_SEARCH(path,'*noklip.fits',COUNT=cc)
oldfiles=FILE_SEARCH(path,'*prelim.fits',COUNT=cc)
FILE_DELETE, oldfiles


if identify and not yj then find_sources, strcompress(root+klipfolder+''+obj+'_ifs_adiklip_YJH'+suffix+'.fits',/rem),sigma=5.0, reference=root+'products/'+'data_cube_psf_med.fits',FWHM=((1.3*1.0E-6)/(8.4))*206265./0.00746,platescale=0.00746,clip_level=0.03,stretch='linear',correction_factor=1.0,inrad=annmode_inout[0]+4, outrad=annmode_inout[1]-4
if identify and yj then find_sources, strcompress(root+klipfolder+''+obj+'_ifs_adiklip_YJ'+suffix+'.fits',/rem),sigma=5.0, reference=root+'products/'+'data_cube_psf_med.fits',FWHM=((1.3*1.0E-6)/(8.4))*206265./0.00746,platescale=0.00746,clip_level=0.03,stretch='linear',correction_factor=1.0,inrad=annmode_inout[0]+4, outrad=annmode_inout[1]-4




if notify_adiklip and not YJ then begin


body='IFS adiklip reduction with high pass filter width = '+string(filter)+' --- Message automatically generated on ' + String(systime())
notify_string='echo "'+body+'" | mail -s "ZONASPHERE '+ obj + ' IFS Data Reduction Results" kwagner@as.arizona.edu -A '+strcompress(root+klipfolder+obj+'_ifs_adiklip_YJH'+suffix+'.fits_source_map.png',/rem) + $
' -A '+strcompress(root+klipfolder+obj+'_ifs_adiklip_YJH'+suffix+'.fits_sources.txt',/rem) + $
' -A '+strcompress(root+klipfolder+obj+'_ifs_adiklip_YJH'+suffix+'.fits_radial_sensitivity.png',/rem) + $
' -A '+strcompress(root+klipfolder+obj+'_ifs_adiklip_YJH'+suffix+'.fits_radial_sensitivity_2.0_arcsec.png',/rem) 

spawn, notify_string


sources_rho=0. & sources_x=0. & sources_y=0. & candidates=!null ;predefining to keep things running
restore,root+klipfolder+obj+'_ifs_adiklip_YJH'+suffix+'.fits_sources.sav'

;select candidates based on contrast and separation

candidates = where(sources_rho le 2.0 and sources_rho ge 0.)

;candidates=indgen(7) 

if n_elements(candidates) eq 1 then if candidates eq -1 then candidates=!null
if sources_x[0] eq 0. or sources_y[0] eq 0. then candidates=!null


if candidates ne !NULL then begin



candidates_contrast= sources_contrast[candidates]
candidates_x=sources_x[candidates]
candidates_y=sources_y[candidates]
candidates_rho=sources_rho[candidates]
candidates_theta=sources_theta[candidates]
candidates_sigma=sources_sigma[candidates]


smalls=readfits(root+klipfolder+obj+'_ifs_adiklip_YJH'+suffix+'.fits_sources.fits')
smalls=smalls[*,*,candidates]

rho_string='' & sigma_string='' & contrast_string='' & attach_string = '' & theta_string='' & candidates_string=''

cgcleanup

for mmm=0,n_elements(candidates)-1 do begin

plotcolor='Snow'
	cand=smalls[*,*,mmm]
	cand[where(finite(cand) ne 1)]=0.
	rho_string=rho_string+string(candidates_rho[mmm])
	theta_string=theta_string+string(candidates_theta[mmm])
	sigma_string=sigma_string+string(candidates_sigma[mmm])
	contrast_string=contrast_string+string(candidates_contrast[mmm])
	candidates_string=candidates_string+string(candidates[mmm]+1)
	
	cgps_open,root+klipfolder+obj+'_ifs_adiklip_YJH'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png'


   cgLoadCT, 1;, NColors=100,BOTTOM=100;,BOTTOM=100
	cgDisplay,40,40,/device

	;cgImage,cand;,stretch=stretch,xrange=[0,(size(smalls))(2)],yrange=[0,(size(smalls))(2)],/keep_aspect_ratio,minvalue=min(cand),maxvalue=max(cand),mean=mean(cand);,position=[0,0,1,1]

cgImage,cand,xrange=[0,39],yrange=[0,39],/keep_aspect_ratio;=[0,0,1,1];,minvalue=min(cand),maxvalue=max(cand),mean=mean(cand);,position=[0,0,1,1]

	cgps_close,/png

	cgcleanup

	attach_string = attach_string+ ' -A '+strcompress(root+klipfolder+obj+'_ifs_adiklip_YJH'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png',/rem)
endfor


;for cc=0,n_elements(candidates)-1 do begin

	body='ALERT!'+string(n_elements(candidates))+' Candidate exoplanets discovered around '+obj+' in adiklip IFS reduction with significance = '+sigma_string+' contrast = '+contrast_string+', rho='+rho_string+', theta='+theta_string+' source numbers ='+candidates_string+' --- Message automatically generated by ZONASPHERE on ' + String(systime())

	notify_string='echo "'+body+'" | mail -s "ALERT!!! '+string(n_elements(candidates)) +' Exoplanet Candidates Detected Around '+ obj + '!!!" kwagner@as.arizona.edu -A '+strcompress(root+klipfolder+obj+'_ifs_adiklip_YJH'+suffix+'.fits_source_map.png',/rem) +' -A '+strcompress(root+klipfolder+obj+'_ifs_adiklip_YJH'+suffix+'.fits_sources.txt',/rem) + ' -A '+strcompress(root+klipfolder+obj+'_ifs_adiklip_YJH'+suffix+'.fits_SNR_map.png',/rem)  + attach_string

	spawn, notify_string 

for mmm=0,n_elements(candidates)-1 do file_delete,strcompress(root+klipfolder+obj+'_ifs_adiklip_YJH'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png',/rem)



;endfor
endif
endif

if notify_adiklip and YJ then begin


body='IFS adiklip reduction with high pass filter width = '+string(filter)+' --- Message automatically generated on ' + String(systime())
notify_string='echo "'+body+'" | mail -s "ZONASPHERE '+ obj + ' IFS Data Reduction Results" kwagner@as.arizona.edu -A '+strcompress(root+klipfolder+obj+'_ifs_adiklip_YJ'+suffix+'.fits_source_map.png',/rem) + $
' -A '+strcompress(root+klipfolder+obj+'_ifs_adiklip_YJ'+suffix+'.fits_sources.txt',/rem) + $
' -A '+strcompress(root+klipfolder+obj+'_ifs_adiklip_YJ'+suffix+'.fits_radial_sensitivity.png',/rem) + $
' -A '+strcompress(root+klipfolder+obj+'_ifs_adiklip_YJ'+suffix+'.fits_radial_sensitivity_2.0_arcsec.png',/rem) 

spawn, notify_string


sources_rho=0. & sources_x=0. & sources_y=0. & candidates=!null ;predefining to keep things running
restore,root+klipfolder+obj+'_ifs_adiklip_YJ'+suffix+'.fits_sources.sav'

;select candidates based on contrast and separation

candidates = where(sources_rho le 2.0 and sources_rho ge 0.)

;candidates=indgen(7) 

if n_elements(candidates) eq 1 then if candidates eq -1 then candidates=!null
if sources_x[0] eq 0. or sources_y[0] eq 0. then candidates=!null


if candidates ne !NULL then begin



candidates_contrast= sources_contrast[candidates]
candidates_x=sources_x[candidates]
candidates_y=sources_y[candidates]
candidates_rho=sources_rho[candidates]
candidates_theta=sources_theta[candidates]
candidates_sigma=sources_sigma[candidates]


smalls=readfits(root+klipfolder+obj+'_ifs_adiklip_YJ'+suffix+'.fits_sources.fits')
smalls=smalls[*,*,candidates]

rho_string='' & sigma_string='' & contrast_string='' & attach_string = '' & theta_string='' & candidates_string=''

cgcleanup

for mmm=0,n_elements(candidates)-1 do begin

plotcolor='Snow'
	cand=smalls[*,*,mmm]
	cand[where(finite(cand) ne 1)]=0.
	rho_string=rho_string+string(candidates_rho[mmm])
	theta_string=theta_string+string(candidates_theta[mmm])
	sigma_string=sigma_string+string(candidates_sigma[mmm])
	contrast_string=contrast_string+string(candidates_contrast[mmm])
	candidates_string=candidates_string+string(candidates[mmm]+1)
	
	cgps_open,root+klipfolder+obj+'_ifs_adiklip_YJ'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png'


   cgLoadCT, 1;, NColors=100,BOTTOM=100;,BOTTOM=100
	cgDisplay,40,40,/device

	;cgImage,cand;,stretch=stretch,xrange=[0,(size(smalls))(2)],yrange=[0,(size(smalls))(2)],/keep_aspect_ratio,minvalue=min(cand),maxvalue=max(cand),mean=mean(cand);,position=[0,0,1,1]

cgImage,cand,xrange=[0,39],yrange=[0,39],/keep_aspect_ratio;=[0,0,1,1];,minvalue=min(cand),maxvalue=max(cand),mean=mean(cand);,position=[0,0,1,1]

	cgps_close,/png

	cgcleanup

	attach_string = attach_string+ ' -A '+strcompress(root+klipfolder+obj+'_ifs_adiklip_YJ'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png',/rem)
endfor


;for cc=0,n_elements(candidates)-1 do begin

	body='ALERT!'+string(n_elements(candidates))+' Candidate exoplanets discovered around '+obj+' in adiklip IFS reduction with significance = '+sigma_string+' contrast = '+contrast_string+', rho='+rho_string+', theta='+theta_string+' source numbers ='+candidates_string+' --- Message automatically generated by ZONASPHERE on ' + String(systime())

	notify_string='echo "'+body+'" | mail -s "ALERT!!! '+string(n_elements(candidates)) +' Exoplanet Candidates Detected Around '+ obj + '!!!" kwagner@as.arizona.edu -A '+strcompress(root+klipfolder+obj+'_ifs_adiklip_YJ'+suffix+'.fits_source_map.png',/rem) +' -A '+strcompress(root+klipfolder+obj+'_ifs_adiklip_YJ'+suffix+'.fits_sources.txt',/rem) + ' -A '+strcompress(root+klipfolder+obj+'_ifs_adiklip_YJ'+suffix+'.fits_SNR_map.png',/rem)  + attach_string

	spawn, notify_string 

for mmm=0,n_elements(candidates)-1 do file_delete,strcompress(root+klipfolder+obj+'_ifs_adiklip_YJ'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png',/rem)



;endfor
endif
endif




endif ;adi if

;----------------Starting SDI

if sdi  then begin

if fresh_sdi  then adicube=scicube else $
	adicube=readfits(root+klipfolder+obj+'_ifs_adiklip_xyln.fits')
	


if rdi_sdi  then adicube=readfits(root+klipfolder+obj+'_ifs_rdiklip_xyln.fits')

if fresh_sdi  then begin
for k=0,38 do begin
		for ii=0, cnt-1 do begin
 		adicube[*,*,k,ii]=rot(adicube[*,*,k,ii],-pa[ii]-pupoff[ii],/INTERP)
		endfor
endfor
endif


finalcube=adicube ;to be used later
for i=0,cnt-1 do begin

	tempcube=reform(adicube[*,*,*,i]) ;this is a x,y,l cube from just a single IFS frame
	tempcubesc=tempcube	;this will be the scaled cube
	
	print, size(tempcubesc)
	
	if debug  then writefits, '~/Desktop/sdi_input_prescaling.fits',tempcube, head
	
	for k=0,38 do begin ;now for each lambda slice we scale the cube then perform SDI
		for kk=0,38 do begin ;scaling the cube
			tempcubesc[*,*,kk]=rot(tempcube[*,*,kk],0., (scaling[k]/scaling[kk]), /INTERP)
		endfor
		
		tempcubesc[where(finite(tempcubesc) eq 0)]=0.
		starttime=systime(/JULIAN)

		if debug  then writefits, '~/Desktop/sdi_input_debug.fits',tempcubesc, head
	
		print, '--[ SDI-KLIPing Image ', i, ' out of ',cnt-1,' ]--'
		print, 'Object: ',obj
		print, 'On wavelength:', wavelength[k],'/',wavelength[38]
		;now use SDI klip
		if annmode  then tempcubesck = sdiklip(tempcubesc, k_sdiklip, target=k, anglemask=anglemask, distmask=distmask,posang=-angles, scaling=scaling,  wl=wavelength[k], diam=8., pixelscale=0.00746, angsep=sdisep, obj=obj, nrings=nrings_sdi,wr =wr_sdi, n_ang =n_ang,annmode_inout=annmode_inout)
		if hyp  then tempcubesck = sdiklip(tempcubesc, k_sdiklip, target=k, anglemask=anglemask, distmask=distmask,posang=-angles, scaling=scaling,  wl=wavelength[k], diam=8., pixelscale=0.00746, angsep=sdisep, obj=obj, nrings=nrings_sdi,wr =wr_sdi, n_ang =n_ang,spot_radius=spot_radius,rho=rho,phi=phi, /hyper)
		if annmode eq 0 and hyp eq 0 then tempcubesck = sdiklip(tempcubesc, k_sdiklip, target=k, anglemask=anglemask, distmask=distmask,posang=-angles, scaling=scaling,  wl=wavelength[k], diam=8., pixelscale=0.00746, angsep=sdisep, obj=obj, nrings=nrings_sdi,wr =wr_sdi, n_ang =n_ang)
		print, '------------------------[ Percent complete:', ((float(i)/float(cnt-1.))+float(k+1.)/(39.*float(cnt)))*100.
		endtime=double(systime(/JULIAN))
		print, '------------------------[ Time remaining (min):', (1.-((float(i)/float(cnt-1.))+float(k+1.)/(39.*float(cnt))))*(39.*float(cnt))*(endtime-starttime) *86400./60.
		
				if debug  then	writefits, '~/Desktop/sdi_output_debug.fits',tempcubesck, head

		finalcube[*,*,k,i]=tempcubesck
		finalcubederot=finalcube
		if k eq 38 then begin
			if rdi_sdi  then for ii=0, cnt-1 do begin
			for kk=0, 38 do begin
 			finalcubederot[*,*,kk,ii]=rot(finalcube[*,*,kk,ii],-pa[ii]-pupoff[ii],/INTERP)
 			endfor
			endfor
			medarr, reform(finalcubederot[*,*,*,i]), framei
			writefits, strcompress(root+''+klipfolder+''+obj+'_ifs_sdiklip_YJH_'+String(i)+'_prelim.fits',/rem),framei, head
		endif
	endfor

	

endfor

	finalcubederot=finalcube
	
	if rdi_sdi  then begin
for k=0,38 do begin
		for ii=0, cnt-1 do begin
 		finalcubederot[*,*,k,ii]=rot(finalcube[*,*,k,ii],-pa[ii]-pupoff[ii],/INTERP)
		endfor
endfor
endif

		;derotates the frames and median combines for each band
	for k=0,38 do begin
		if comb_type eq 'median' then medarr, reform(finalcubederot[*,*,k,*]), kframemed
		if comb_type eq 'mean' then kframemed=mean(reform(finalcubederot[*,*,k,*]),dim=3)
		if comb_type eq 'nw-mean' then kframemed=nw_ang_comb(reform(finalcubederot[*,*,k,*]),pa+pupoff)
		if k eq 0 then xylcube=kframemed else xylcube=[ [[xylcube]] , [[kframemed]] ]
		if k eq 0 then sumyjh=kframemed else sumyjh=sumyjh+kframemed
	endfor
	
	simyjh=sumyjh/39.
	
	;medarr, xylcube, medyjh
	medyjh=mean(xylcube,dim=3)
	yjhcube=[]
for kz=0, (size(finalcube))(4)-1 do begin &	kzcube=reform(finalcube[*,*,*,kz]) & medarr,kzcube,kzmed & yjhcube=[[[yjhcube]],[[kzmed]]] & endfor	


;make a convolved cube
	conv_xylcube=xylcube

for ccc=0,38 do begin
	
	sz=280.-1.
	width=(wavelength[ccc]*1E-6) / (8.2) * 206265. / 0.00746
	print, 'PSF Width: ',width
	PSF = psf_Gaussian(NPIX=sz, FWHM=[width,width])
	PSFN = PSF/MAX(PSF)
	conv_xylcube[where(finite(conv_xylcube) ne 1)]=0.
	conv_xylcube[*,*,ccc] = convolve(conv_xylcube[*,*,ccc], PSFN)

endfor
	
if yj then writefits, root+''+klipfolder+''+obj+'_ifs_sdiklip_YJ'+suffix+'.fits', medyjh, head else writefits, root+''+klipfolder+''+obj+'_ifs_sdiklip_YJH'+suffix+'.fits', medyjh, head
if yj then writefits, root+''+klipfolder+''+obj+'_ifs_sdiklip_YJ_cube'+suffix+'.fits', yjhcube, head else writefits, root+''+klipfolder+''+obj+'_ifs_sdiklip_YJH_cube'+suffix+'.fits', yjhcube, head
;writefits, root+''+klipfolder+''+obj+'_ifs_sdiklip_YJH_sum.fits', sumyjh, head
writefits, root+''+klipfolder+''+obj+'_ifs_sdiklip_xyln.fits', finalcube, head
writefits, root+''+klipfolder+''+obj+'_ifs_sdiklip_xyl'+suffix+'.fits', xylcube, head

if yj then ysdicube=xylcube[*,*,2:15] else ysdicube=xylcube[*,*,0:8]
;medarr, ysdicube, ysdi
ysdi=mean(ysdicube,dim=3)
writefits, root+''+klipfolder+''+obj+'_ifs_sdiklip_Y'+suffix+'.fits', ysdi, head

if yj then jsdicube=xylcube[*,*,17:38] else jsdicube=xylcube[*,*,10:21]
;medarr, jsdicube, jsdi
jsdi=mean(jsdicube,dim=3)
writefits, root+''+klipfolder+''+obj+'_ifs_sdiklip_J'+suffix+'.fits', jsdi, head
	
if not yj then hsdicube=xylcube[*,*,28:36]
;medarr, hsdicube, hsdi
if not yj then hsdi=mean(hsdicube,dim=3)
if not yj then writefits, root+''+klipfolder+''+obj+'_ifs_sdiklip_H'+suffix+'.fits', hsdi, head
writefits, root+''+klipfolder+''+obj+'_ifs_sdiklip_xyl_conv'+suffix+'.fits', conv_xylcube, head

if yj then ysdicube=conv_xylcube[*,*,2:15] else ysdicube=conv_xylcube[*,*,0:8]
ysdi=mean(ysdicube,dimension=3)
writefits, root+''+klipfolder+''+obj+'_ifs_sdiklip_Y_conv'+suffix+'.fits', ysdi, head

if yj then jsdicube=conv_xylcube[*,*,17:38] else jsdicube=conv_xylcube[*,*,10:21]
jsdi=mean(jsdicube,dimension=3)

writefits, root+''+klipfolder+''+obj+'_ifs_sdiklip_J_conv'+suffix+'.fits', jsdi, head
	
if not yj then hsdicube=conv_xylcube[*,*,28:36]
if not yj then hsdi=mean(hsdicube,dimension=3)

if not yj then writefits, root+''+klipfolder+''+obj+'_ifs_sdiklip_H_conv'+suffix+'.fits', hsdi, head



if not yj then ltfilt=xylcube[*,*,[14,15,16,17,18,19,20,30,31,32,33,34,35,36] ] else ltfilt=xylcube[*,*,[24:36] ]

ltfilt=mean(ltfilt,dimension=3)

yjhsdi=mean(conv_xylcube,dimension=3)

if not yj then  writefits, root+''+klipfolder+''+obj+'_ifs_sdiklip_YJH_conv'+suffix+'.fits', yjhsdi, head else writefits, root+''+klipfolder+''+obj+'_ifs_sdiklip_YJ_conv'+suffix+'.fits', yjhsdi, head

writefits, root+''+klipfolder+''+obj+'_ifs_sdiklip_LTfilt'+suffix+'.fits', ltfilt, head

path=root+'/'+klipfolder+'/'
oldfiles=FILE_SEARCH(path,'*prelim.fits',COUNT=cc)
FILE_DELETE, oldfiles

if identify and not yj then find_sources, strcompress(root+klipfolder+''+obj+'_ifs_sdiklip_YJH'+suffix+'.fits',/rem),sigma=5.0, reference=root+'products/'+'data_cube_psf_med.fits',FWHM=((1.3*1.0E-6)/(8.4))*206265./0.00746,platescale=0.00746,clip_level=0.03,stretch='linear',correction_factor=1.0,inrad=annmode_inout[0]+4, outrad=annmode_inout[1]-4
		if identify and yj then find_sources, strcompress(root+klipfolder+''+obj+'_ifs_sdiklip_YJ'+suffix+'.fits',/rem),sigma=5.0, reference=root+'products/'+'data_cube_psf_med.fits',FWHM=((1.3*1.0E-6)/(8.4))*206265./0.00746,platescale=0.00746,clip_level=0.03,stretch='linear',correction_factor=1.0,inrad=annmode_inout[0]+4, outrad=annmode_inout[1]-4

if identify then find_sources, strcompress(root+klipfolder+''+obj+'_ifs_sdiklip_LTfilt'+suffix+'.fits',/rem),sigma=5.0, reference=root+'products/'+'data_cube_psf_med.fits',FWHM=((1.3*1.0E-6)/(8.4))*206265./0.00746,platescale=0.00746,clip_level=0.03,stretch='linear',correction_factor=1.0,inrad=annmode_inout[0]+4, outrad=annmode_inout[1]-4




if notify_sdiklip and not YJ then begin


body='IFS sdiklip reduction with high pass filter width = '+string(filter)+' --- Message automatically generated on ' + String(systime())
notify_string='echo "'+body+'" | mail -s "ZONASPHERE '+ obj + ' IFS Data Reduction Results" kwagner@as.arizona.edu -A '+strcompress(root+klipfolder+obj+'_ifs_sdiklip_YJH'+suffix+'.fits_source_map.png',/rem) + $
' -A '+strcompress(root+klipfolder+obj+'_ifs_sdiklip_YJH'+suffix+'.fits_sources.txt',/rem) + $
' -A '+strcompress(root+klipfolder+obj+'_ifs_sdiklip_YJH'+suffix+'.fits_radial_sensitivity.png',/rem) + $
' -A '+strcompress(root+klipfolder+obj+'_ifs_sdiklip_YJH'+suffix+'.fits_radial_sensitivity_2.0_arcsec.png',/rem) 

spawn, notify_string


sources_rho=0. & sources_x=0. & sources_y=0. & candidates=!null ;predefining to keep things running
restore,root+klipfolder+obj+'_ifs_sdiklip_YJH'+suffix+'.fits_sources.sav'

;select candidates based on contrast and separation

candidates = where(sources_rho le 2.0 and sources_rho ge 0.)

;candidates=indgen(7) 

if n_elements(candidates) eq 1 then if candidates eq -1 then candidates=!null
if sources_x[0] eq 0. or sources_y[0] eq 0. then candidates=!null


if candidates ne !NULL then begin



candidates_contrast= sources_contrast[candidates]
candidates_x=sources_x[candidates]
candidates_y=sources_y[candidates]
candidates_rho=sources_rho[candidates]
candidates_theta=sources_theta[candidates]
candidates_sigma=sources_sigma[candidates]


smalls=readfits(root+klipfolder+obj+'_ifs_sdiklip_YJH'+suffix+'.fits_sources.fits')
smalls=smalls[*,*,candidates]

rho_string='' & sigma_string='' & contrast_string='' & attach_string = '' & theta_string='' & candidates_string=''

cgcleanup

for mmm=0,n_elements(candidates)-1 do begin

plotcolor='Snow'
	cand=smalls[*,*,mmm]
	cand[where(finite(cand) ne 1)]=0.
	rho_string=rho_string+string(candidates_rho[mmm])
	theta_string=theta_string+string(candidates_theta[mmm])
	sigma_string=sigma_string+string(candidates_sigma[mmm])
	contrast_string=contrast_string+string(candidates_contrast[mmm])
	candidates_string=candidates_string+string(candidates[mmm]+1)
	
	cgps_open,root+klipfolder+obj+'_ifs_sdiklip_YJH'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png'


   cgLoadCT, 1;, NColors=100,BOTTOM=100;,BOTTOM=100
	cgDisplay,40,40,/device

	;cgImage,cand;,stretch=stretch,xrange=[0,(size(smalls))(2)],yrange=[0,(size(smalls))(2)],/keep_aspect_ratio,minvalue=min(cand),maxvalue=max(cand),mean=mean(cand);,position=[0,0,1,1]

cgImage,cand,xrange=[0,39],yrange=[0,39],/keep_aspect_ratio;=[0,0,1,1];,minvalue=min(cand),maxvalue=max(cand),mean=mean(cand);,position=[0,0,1,1]

	cgps_close,/png

	cgcleanup

	attach_string = attach_string+ ' -A '+strcompress(root+klipfolder+obj+'_ifs_sdiklip_YJH'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png',/rem)
endfor


;for cc=0,n_elements(candidates)-1 do begin

	body='ALERT!'+string(n_elements(candidates))+' Candidate exoplanets discovered around '+obj+' in sdiklip IFS reduction with significance = '+sigma_string+' contrast = '+contrast_string+', rho='+rho_string+', theta='+theta_string+' source numbers ='+candidates_string+' --- Message automatically generated by ZONASPHERE on ' + String(systime())

	notify_string='echo "'+body+'" | mail -s "ALERT!!! '+string(n_elements(candidates)) +' Exoplanet Candidates Detected Around '+ obj + '!!!" kwagner@as.arizona.edu -A '+strcompress(root+klipfolder+obj+'_ifs_sdiklip_YJH'+suffix+'.fits_source_map.png',/rem) +' -A '+strcompress(root+klipfolder+obj+'_ifs_sdiklip_YJH'+suffix+'.fits_sources.txt',/rem) + ' -A '+strcompress(root+klipfolder+obj+'_ifs_sdiklip_YJH'+suffix+'.fits_SNR_map.png',/rem)  + attach_string

	spawn, notify_string 

for mmm=0,n_elements(candidates)-1 do file_delete,strcompress(root+klipfolder+obj+'_ifs_sdiklip_YJH'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png',/rem)



;endfor
endif
endif

if notify_sdiklip and YJ then begin


body='IFS sdiklip reduction with high pass filter width = '+string(filter)+' --- Message automatically generated on ' + String(systime())
notify_string='echo "'+body+'" | mail -s "ZONASPHERE '+ obj + ' IFS Data Reduction Results" kwagner@as.arizona.edu -A '+strcompress(root+klipfolder+obj+'_ifs_sdiklip_YJ'+suffix+'.fits_source_map.png',/rem) + $
' -A '+strcompress(root+klipfolder+obj+'_ifs_sdiklip_YJ'+suffix+'.fits_sources.txt',/rem) + $
' -A '+strcompress(root+klipfolder+obj+'_ifs_sdiklip_YJ'+suffix+'.fits_radial_sensitivity.png',/rem) + $
' -A '+strcompress(root+klipfolder+obj+'_ifs_sdiklip_YJ'+suffix+'.fits_radial_sensitivity_2.0_arcsec.png',/rem) 

spawn, notify_string


sources_rho=0. & sources_x=0. & sources_y=0. & candidates=!null ;predefining to keep things running
restore,root+klipfolder+obj+'_ifs_sdiklip_YJ'+suffix+'.fits_sources.sav'

;select candidates based on contrast and separation

candidates = where(sources_rho le 2.0 and sources_rho ge 0.)

;candidates=indgen(7) 

if n_elements(candidates) eq 1 then if candidates eq -1 then candidates=!null
if sources_x[0] eq 0. or sources_y[0] eq 0. then candidates=!null


if candidates ne !NULL then begin



candidates_contrast= sources_contrast[candidates]
candidates_x=sources_x[candidates]
candidates_y=sources_y[candidates]
candidates_rho=sources_rho[candidates]
candidates_theta=sources_theta[candidates]
candidates_sigma=sources_sigma[candidates]


smalls=readfits(root+klipfolder+obj+'_ifs_sdiklip_YJ'+suffix+'.fits_sources.fits')
smalls=smalls[*,*,candidates]

rho_string='' & sigma_string='' & contrast_string='' & attach_string = '' & theta_string='' & candidates_string=''

cgcleanup

for mmm=0,n_elements(candidates)-1 do begin

plotcolor='Snow'
	cand=smalls[*,*,mmm]
	cand[where(finite(cand) ne 1)]=0.
	rho_string=rho_string+string(candidates_rho[mmm])
	theta_string=theta_string+string(candidates_theta[mmm])
	sigma_string=sigma_string+string(candidates_sigma[mmm])
	contrast_string=contrast_string+string(candidates_contrast[mmm])
	candidates_string=candidates_string+string(candidates[mmm]+1)
	
	cgps_open,root+klipfolder+obj+'_ifs_sdiklip_YJ'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png'


   cgLoadCT, 1;, NColors=100,BOTTOM=100;,BOTTOM=100
	cgDisplay,40,40,/device

	;cgImage,cand;,stretch=stretch,xrange=[0,(size(smalls))(2)],yrange=[0,(size(smalls))(2)],/keep_aspect_ratio,minvalue=min(cand),maxvalue=max(cand),mean=mean(cand);,position=[0,0,1,1]

cgImage,cand,xrange=[0,39],yrange=[0,39],/keep_aspect_ratio;=[0,0,1,1];,minvalue=min(cand),maxvalue=max(cand),mean=mean(cand);,position=[0,0,1,1]

	cgps_close,/png

	cgcleanup

	attach_string = attach_string+ ' -A '+strcompress(root+klipfolder+obj+'_ifs_sdiklip_YJ'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png',/rem)
endfor


;for cc=0,n_elements(candidates)-1 do begin

	body='ALERT!'+string(n_elements(candidates))+' Candidate exoplanets discovered around '+obj+' in sdiklip IFS reduction with significance = '+sigma_string+' contrast = '+contrast_string+', rho='+rho_string+', theta='+theta_string+' source numbers ='+candidates_string+' --- Message automatically generated by ZONASPHERE on ' + String(systime())

	notify_string='echo "'+body+'" | mail -s "ALERT!!! '+string(n_elements(candidates)) +' Exoplanet Candidates Detected Around '+ obj + '!!!" kwagner@as.arizona.edu -A '+strcompress(root+klipfolder+obj+'_ifs_sdiklip_YJ'+suffix+'.fits_source_map.png',/rem) +' -A '+strcompress(root+klipfolder+obj+'_ifs_sdiklip_YJ'+suffix+'.fits_sources.txt',/rem) + ' -A '+strcompress(root+klipfolder+obj+'_ifs_sdiklip_YJ'+suffix+'.fits_SNR_map.png',/rem)  + attach_string

	spawn, notify_string 

for mmm=0,n_elements(candidates)-1 do file_delete,strcompress(root+klipfolder+obj+'_ifs_sdiklip_YJ'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png',/rem)



;endfor
endif
endif


if notify_sdiklip_ltfilt then begin


body='IFS sdiklip reduction with T-type filter and high pass filter width = '+string(filter)+' --- Message automatically generated on ' + String(systime())
notify_string='echo "'+body+'" | mail -s "ZONASPHERE '+ obj + ' IFS SDI Data Reduction Results with T-filter" kwagner@as.arizona.edu -A '+strcompress(root+klipfolder+obj+'_ifs_sdiklip_LTfilt'+suffix+'.fits_source_map.png',/rem) + $
' -A '+strcompress(root+klipfolder+obj+'_ifs_sdiklip_LTfilt'+suffix+'.fits_sources.txt',/rem) + $
' -A '+strcompress(root+klipfolder+obj+'_ifs_sdiklip_LTfilt'+suffix+'.fits_radial_sensitivity.png',/rem) + $
' -A '+strcompress(root+klipfolder+obj+'_ifs_sdiklip_LTfilt'+suffix+'.fits_radial_sensitivity_2.0_arcsec.png',/rem) 

spawn, notify_string


sources_rho=0. & sources_x=0. & sources_y=0. & candidates=!null ;predefining to keep things running
restore,root+klipfolder+obj+'_ifs_sdiklip_LTfilt'+suffix+'.fits_sources.sav'

;select candidates based on contrast and separation

if sources_rho ne !null then candidates = where(sources_rho le 2.0 and sources_rho ge 0.) else candidates=!null

;candidates=indgen(7) 

if n_elements(candidates) eq 1 then if candidates eq -1 then candidates=!null
if sources_x eq !null then begin sources_x=0 & sources_y=0 & endif
if sources_x[0] eq 0. or sources_y[0] eq 0. then candidates=!null


if candidates ne !NULL then begin



candidates_contrast= sources_contrast[candidates]
candidates_x=sources_x[candidates]
candidates_y=sources_y[candidates]
candidates_rho=sources_rho[candidates]
candidates_theta=sources_theta[candidates]
candidates_sigma=sources_sigma[candidates]


smalls=readfits(root+klipfolder+obj+'_ifs_sdiklip_LTfilt'+suffix+'.fits_sources.fits')
smalls=smalls[*,*,candidates]

rho_string='' & sigma_string='' & contrast_string='' & attach_string = '' & theta_string='' & candidates_string=''

cgcleanup

for mmm=0,n_elements(candidates)-1 do begin

plotcolor='Snow'
	cand=smalls[*,*,mmm]
	cand[where(finite(cand) ne 1)]=0.
	rho_string=rho_string+string(candidates_rho[mmm])
	theta_string=theta_string+string(candidates_theta[mmm])
	sigma_string=sigma_string+string(candidates_sigma[mmm])
	contrast_string=contrast_string+string(candidates_contrast[mmm])
	candidates_string=candidates_string+string(candidates[mmm]+1)
	
	cgps_open,root+klipfolder+obj+'_ifs_sdiklip_LTfilt'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png'


   cgLoadCT, 1;, NColors=100,BOTTOM=100;,BOTTOM=100
	cgDisplay,40,40,/device

	;cgImage,cand;,stretch=stretch,xrange=[0,(size(smalls))(2)],yrange=[0,(size(smalls))(2)],/keep_aspect_ratio,minvalue=min(cand),maxvalue=max(cand),mean=mean(cand);,position=[0,0,1,1]

cgImage,cand,xrange=[0,39],yrange=[0,39],/keep_aspect_ratio;=[0,0,1,1];,minvalue=min(cand),maxvalue=max(cand),mean=mean(cand);,position=[0,0,1,1]

	cgps_close,/png

	cgcleanup

	attach_string = attach_string+ ' -A '+strcompress(root+klipfolder+obj+'_ifs_sdiklip_LTfilt'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png',/rem)
endfor


;for cc=0,n_elements(candidates)-1 do begin

	body='ALERT!'+string(n_elements(candidates))+' Candidate T-type exoplanets discovered around '+obj+' in sdiklip IFS reduction with significance = '+sigma_string+' contrast = '+contrast_string+', rho='+rho_string+', theta='+theta_string+' source numbers ='+candidates_string+' --- Message automatically generated by ZONASPHERE on ' + String(systime())

	notify_string='echo "'+body+'" | mail -s "ALERT!!! '+string(n_elements(candidates)) +' T-type Exoplanet Candidates Detected Around '+ obj + '!!!" kwagner@as.arizona.edu -A '+strcompress(root+klipfolder+obj+'_ifs_sdiklip_LTfilt'+suffix+'.fits_source_map.png',/rem) +' -A '+strcompress(root+klipfolder+obj+'_ifs_sdiklip_LTfilt'+suffix+'.fits_sources.txt',/rem) + ' -A '+strcompress(root+klipfolder+obj+'_ifs_sdiklip_LTfilt'+suffix+'.fits_SNR_map.png',/rem)  + attach_string

	spawn, notify_string 

for mmm=0,n_elements(candidates)-1 do file_delete,strcompress(root+klipfolder+obj+'_ifs_sdiklip_LTfilt'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png',/rem)



;endfor
endif
endif



endif

;----------------Starting RDI

if rdi  then begin

finalcube=scicube ;to be used later

;make cubes the same size
endsize=(size(scicube))(4)-1
if endsize lt (size(refcube))(4)-1 then refcube=refcube[*,*,*,0:endsize]

print, size(refcube)
print, size(scicube)


for k=0,38 do begin

	scikcube=reform(scicube[*,*,k,*]) ;this is a x,y,l cube from just a single IFS frame
	refkcube=reform(refcube[*,*,k,*])	;this will be the scaled cube
		
	if debug  then writefits, '~/Desktop/rdi_sci_input.fits',scikcube, head
	if debug  then writefits, '~/Desktop/rdi_ref_input.fits',refkcube, refhead
	
	for i=0,cnt-1 do begin ;loop through target frames and perform RDI
		refkcube[where(finite(refkcube) eq 0)]=0.
		starttime=systime(/JULIAN)

		print, '--[ SDI-KLIPing Image ', i, ' out of ',cnt-1,' ]--'
		print, 'On wavelength:', wavelength[k]
		;now use SDI klip
		if hyp eq 0 and annmode  then rdikcube = rdiklip(scikcube, refkcube, k_rdiklip, target=i, anglemask=anglemask, distmask=distmask,posang=-angles, scaling=scaling,  wl=wavelength[k], diam=8., pixelscale=0.00746, angsep=0., obj=obj, nrings=nrings_rdi,wr =wr_rdi, n_ang =n_ang,annmode_inout=annmode_inout)
		if hyp  then rdikcube = rdiklip(scikcube, refkcube, k_rdiklip, target=i, anglemask=anglemask, distmask=distmask,posang=-angles, scaling=scaling,  wl=wavelength[k], diam=8., pixelscale=0.00746, angsep=0., obj=obj, nrings=nrings_rdi,wr =wr_rdi, n_ang =n_ang,spot_radius=spot_radius,rho=rho,phi=phi, /hyper)
		if hyp eq 0 and annmode eq 0 then rdikcube = rdiklip(scikcube, refkcube, k_rdiklip, target=i, anglemask=anglemask, distmask=distmask,posang=-angles, scaling=scaling,  wl=wavelength[k], diam=8., pixelscale=0.00746, angsep=0., obj=obj, nrings=nrings_rdi,wr =wr_rdi, n_ang =n_ang)

		if debug  then	writefits, '~/Desktop/rdi_output_debug.fits',rdikcube, head

		print, size(rdikcube)	
		
		finalcube[*,*,k,i]=rdikcube
		
		if i eq cnt-1 then begin
		finalcubederot=finalcube
			for ii = 0 , cnt-1 do finalcubederot[*,*,k,ii]=rot(finalcube[*,*,k,ii],-pa[ii]-pupoff[ii],/INTERP)
			medarr, reform(finalcubederot[*,*,k,*]), framek
			writefits, root+''+klipfolder+''+obj+'_ifs_rdiklip_k'+String(k)+'_prelim.fits',framek, head
		endif
	endfor
endfor
	finalcubederot=finalcube
	

		;derotates the frames and median combines for each band
	for k=0,38 do begin
	
	for ii=0, cnt-1 do begin
 			finalcubederot[*,*,k,ii]=rot(finalcubederot[*,*,k,ii],-pa[ii]-pupoff[ii],/INTERP)
	endfor
	
		if comb_type eq 'median' then medarr, reform(finalcubederot[*,*,k,*]), kframemed
		if comb_type eq 'mean' then kramemed=mean(reform(finalcubederot[*,*,k,*]),dim=3)
		if comb_type eq 'nw-mean' then kramemed=nw_ang_comb(reform(finalcubederot[*,*,k,*]),pa+pupoff)

		if k eq 0 then xylcube=kframemed else xylcube=[ [[xylcube]] , [[kframemed]] ]
		
		if k eq 0 then sumyjh=kframemed else sumyjh=sumyjh+kframemed
	endfor
	
	simyjh=sumyjh/39.
	
	medarr, xylcube, medyjh
	yjhcube=[]
for kz=0, (size(finalcubederot))(4)-1 do begin &	kzcube=reform(finalcubederot[*,*,*,kz]) & medarr,kzcube,kzmed & yjhcube=[[[yjhcube]],[[kzmed]]] & endfor	
	
if yj then writefits, root+''+klipfolder+''+obj+'_ifs_rdiklip_YJ'+suffix+'.fits', medyjh, head else writefits, root+''+klipfolder+''+obj+'_ifs_rdiklip_YJH'+suffix+'.fits', medyjh, head
if yj then writefits, root+''+klipfolder+''+obj+'_ifs_rdiklip_YJ_cube'+suffix+'.fits', yjhcube, head else writefits, root+''+klipfolder+''+obj+'_ifs_rdiklip_YJH_cube'+suffix+'.fits', yjhcube, head
;writefits, root+''+klipfolder+''+obj+'_ifs_rdiklip_YJH_sum.fits', sumyjh, head

writefits, root+''+klipfolder+''+obj+'_ifs_rdiklip_xyln.fits', finalcube, head
writefits, root+''+klipfolder+''+obj+'_ifs_rdiklip_xyl'+suffix+'.fits', xylcube, head

if yj then sdicube=xylcube[*,*,2:15] else ysdicube=xylcube[*,*,0:8]
;medarr, ysdicube, ysdi
ysdi=mean(ysdicube,dim=3)
writefits, root+''+klipfolder+''+obj+'_ifs_rdiklip_Y'+suffix+'.fits', ysdi, head

if yj then jsdicube=xylcube[*,*,17:38] else jsdicube=xylcube[*,*,10:21]
;medarr, jsdicube, jsdi
jsdi=mean(jsdicube,dim=3)
writefits, root+''+klipfolder+''+obj+'_ifs_rdiklip_J'+suffix+'.fits', jsdi, head
	
if not yj then hsdicube=xylcube[*,*,28:36]
;medarr, hsdicube, hsdi
if not yj then hsdi=mean(hsdicube,dim=3)
if not yj then writefits, root+''+klipfolder+''+obj+'_ifs_rdiklip_H'+suffix+'.fits', hsdi, head


path=root+'/'+klipfolder+'/'
oldfiles=FILE_SEARCH(path,'*prelim.fits',COUNT=cc)
FILE_DELETE, oldfiles




endif

end
