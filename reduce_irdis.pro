;created 2017 03 01 – updated 2017 06 26
;Author: Kevin Wagner (kwagner@as.arizona.edu)
;version: highly beta

; This is version 2 of the University of Arizona pipeline developed
; by Kevin Wagner and Daniel Apai. The first version was very robust,
; and in the end had too many modules and functions to really be 
; efficient (and was very cumbersome to debug). This new version
; is intended to keep the functionality and most used features 
; of the 'medirdis.pro' and 'klipirdis.pro' routines, but to 
; provide higher user friendliness and a cleaner script inter-
; face for future modification and debugging.


; INPUTS: 
;	OBJ: object name, and the name of the folder to look for the files in within the reduction path
;	REDUCTION_PATH: as described, the path toward the data files (expects to contain folders with
;	 		object names)
; OPTIONAL INPUTS:
;	INJECT: binary, enable to turn on fake planet injections (parameters defined below).
;	AUTOCEN: binary, turn on to enabel auto-centering of the data (only works with K12 and x-style 
;		 satellite spots at present). 
;	SKIP_START: binary, will skip initial calibrations and go straight to post-processing. Note that 
;		    the initial steps will need to be run at least on the first pass.
;	KLIP: binary, will enable KLIP processing
;	SUFFIX: extra suffix to add to the files
;	SECOND_CEN: will use second star center file for calibration.
;	
;	Other inputs defined below may be changed from target to target. 
; OUTPUTS: many files.


pro reduce_irdis,  obj, inject=inject, autocen=autocen,reduction_path=reduction_path, skip_start=skip_start, klip=klip, suffix=suffix, second_cen=second_cen,notify=notify,rho=rho,theta=theta,contrast=contrast

if not keyword_set(reduction_path) then reduction_path='/Users/kevinwagner/Data/SPHERE/'


;this loop is for performing parameter optimizations. Optimize injections with the following parameters, then see the later section test_residals to define testing parameters.
NRUNS=1.


for runs=0., NRUNS-1. do begin

;S25_2: 149.6, 0.1125, l=1.5E-3,r=1.7E-3

;currently injections are set up to test sensitivity with a uniform grid of fake planets.

if nruns gt 1 then begin 	
contrast_left=  1.0E-3  + RUNS*1E-4


theta=118.5; -  (0.1* runs) 
rhop=0.15; - (0.01*runs)	;right 0.825, left: 0.830	

contrast_right=  contrast_left;+4E-4;1.0E-3 ;-  RUNS*1E-6


planet_r=[rhop]
planet_theta=([theta]+90.)
planet_contrast_left=[-contrast_left]
planet_contrast_right=[-contrast_right]
endif
;planet_r=[rhop,rhop,rhop,rhop,rhop,rhop,rhop,rhop]
;planet_theta=[theta,theta+180.,theta+90.,theta-90.,theta+45.,theta+180.+45.,theta+90.+45.,theta-90.+45.]

;planet_r=[planet_r,planet_r-0.2,planet_r+0.2,planet_r+0.4,planet_r-0.4]
;planet_theta=[planet_theta,planet_theta-27.5,planet_theta+27.5,planet_theta+45.,planet_theta+45.]

;planet_contrast_left=[contrast_left,contrast_left,contrast_left,contrast_left,contrast_left,contrast_left,contrast_left,contrast_left]
;planet_contrast_right=[contrast_right,contrast_right,contrast_right,contrast_right,contrast_right,contrast_right,contrast_right,contrast_right]

;planet_contrast_left=[planet_contrast_left,planet_contrast_left,planet_contrast_left,planet_contrast_left,planet_contrast_left,planet_contrast_left,planet_contrast_left]
;planet_contrast_right=[planet_contrast_right,planet_contrast_right,planet_contrast_right,planet_contrast_right,planet_contrast_right,planet_contrast_right,planet_contrast_right]


if not keyword_set(rho) and 1 eq 0 then begin
nplanets=120

planet_contrast_left=fltarr(nplanets) & planet_contrast_left[*]=4.0E-5
planet_contrast_right=fltarr(nplanets) & planet_contrast_right[*]=4.0E-5



planet_theta=findgen(nplanets)*!PI/4./!DTOR

for ii=0,nplanets-1 do planet_theta[ii]=planet_theta[ii]+(float(ii)*!PI/16./!DTOR)
planet_r=findgen(nplanets)*(6.0/nplanets);/0.01225
endif

if keyword_set(rho) then begin

	planet_r=[rho]
	planet_theta=[theta]/!DTOR
	planet_contrast_left=[contrast]
	planet_contrast_right=[contrast]
	skip_right=1
endif else skip_right=0



;############################################################################
;#######################[     BEGIN USER INPUT     ]#########################
;############################################################################

;-----------------------[ Enter reduction parameters ]-----------------------


debug = 0 		;if 1 will output many diagnostic files and messages, some to User Desktop. 

identify=0		;will find sources
if keyword_set(contrast) then identify=0

pxscale=0.01225

;image combination type for final steps
comb_type='nw-mean'	;'mean', 'median', 'nw-mean' (Bottom et al. 2017)
	;median is fastest, nw-mean provides best quality combination

do_flux 	= 1		; Will reduce the  flux files. Typically
				; only needs to be on for the first run and then can be 
				; turned off. Scales the frames to DIT = 1 s and by the 
				; neutral density transmission profile in each filter. 
			
	flux_file_num=0	; Which flux file to use? Typically file 0 is fine.
	do_flux_bad_pixel_removal=0	
	IF RUNS GT 0 THEN DO_FLUX=0

	do_center = 1		; Will reduce the star center waffle frames to determine 
	if runs gt 0 then do_center=0
				; the true center. It will then align the science frames 
				; accordingly via a cross-correlation step, and finally
				; will scale to a DIT of 1s.
				; Assumes there is no ND filter.
	if keyword_set(autocen) then auto_flux_cen=1 else auto_flux_cen=0
	if obj eq 'SA1' then auto_flux_cen=0

	if keyword_set(autocen) then auto_cen=1 else auto_cen=0
				
	auto_frame_hsize=10.	;half size of centering box used to initially find the star's location (10.)

		
	cen_filter      =       1	;will filter and smooth data before centering (1)

		
	do_multi_cen	=	0	;will process extra star center frames (but won't use for science) (0)
	use_second_cen  = 	0	;will use star center frame from after the object seq. (0)
		if keyword_set(second_cen) then use_second_cen=1

	field_stab_cen  = 	0	;will align on field stabilized cube (useful for tight binaries) (0)
		

	field_stab_clean = 0 	;will check against a field stabilized median (useful for binaries/disks) (0)	
		if field_stab_cen then field_stab_clean=1

	spot_fit_hwidth	=	10.	;half-width of box around satellite spot used in mpfit (10.)


	init_centroid   =   0	;will use user click for initial centering.. then spot fits (if enabled) (0)
	no_spot_fit		=	0	; Will use user clicks instead for the waffle-fit (0)
	spot_smooth=0	;will smooth center frame before gaussian fitting (useful for low SNR spots)	(0)

	click_center	=	0	; Will use a user-click as center (with centroid). (0)

		
	disp_type='alog'		;Display stretch during star centering (alog)
					;alog typically works well, but linear may be better for bright
					;targets. Options: 'asin', 'alog', 'linear'
							
	block_rad_cen=1	& rad_cen_min=0. & rad_cen_max=50.  ;(1, 0., 50.)
			
	remove_wild=1	;if set to 1 will reject shifts that are 5 pixels or more from the previousl center (1)
		if field_stab_cen then remove_wild=0	

	use_left_for_right=0	;use left shifts for right frames (0)


	output_cubes = 0 ;will take up lots of disk space! (0)					
									
do_rotational_center=0	;rotational-based centering (0)
	censteps=41		;steps in x and y - i.e. total steps will be square of this (41.)
	censquare=2.		;test square size in pixels (please enter a float value) (2.)
	cenblockradius=5.	;size of center to block in ADI testing (5.)
	cenframewidth=60.	;(60.)

skip_rotational_center=1	;default is to skip rotational centering step (1)

do_clean=1    &    cross_thresh=0.85	;(1, 0.85)
	if not do_clean and not keyword_set(skip_start) then skip_clean=1 else skip_clean=0	;will write original cubes with '_clean' suffix for further use

if keyword_set(skip_start) then begin
	do_flux=0
	do_center=0
	;do_clean=0
endif	

if keyword_set(inject) then use_inject=1 else use_inject=0
		
if use_inject then do_inject = 1 else do_inject=0; Will inject fake planets into the data				
	use_gauss=1	;(1)	;Gaussian fit to the instrumental PSF?

; ---- next five inputs are shared by derotate and klip ----

	final_image_half_size=499;199;119;99;99 ;max=255	(119)
	high_pass_width_final=0.	; smooth after subtracting and KLIPping, only affects residuals
	smooth_width=0.;8.;0.
	destripe_frames=1	;(1)

do_rotate=1
if keyword_set(rho) then do_rotate=0
	notify_adi=0
	identify_adi=0

	;will perform ADI
	classical_rdi=0
		classical_ref_obj='S47'
	smart_adi=0 ;takes much more time (unless images are tiny), overrides RDI
	if runs gt 0 then do_rotate=0

if keyword_set(klip) then do_klip=1 else do_klip=0

;main KLIP parameters
	bin=1    	& 	  high_pass_width=13.;high pass width also defines filter for ADI and derot images	
	auto_bin=1	;bins to cube size of about fifty
	k_klip=5	&  	  k_klip_sdi=4 	& k_klip_rdi=4
	angsep=0.7  &     n_ang=6;6

;secondary parameters
	chop=0	;removes negative pixels. This never really helps anything.
	hyp = 0	;will only process the region around the planet
	if hyp then rho=[0.85]/pxscale;	spot location input in a radius in arcsec, then converted to pixels
		phi=[193.] & phi=(phi)*!DTOR;
		spot_radius=14.;	spot-width
	annmode=0 & if hyp then annmode=0	;annular mode?
		annmode_inout=[55,78]
	anglemax=360. ;set to 360 to include all frames as references
	nrings=12
	wr=fix((final_image_half_size+1.)/float(nrings))
	fill=1 & if not do_rotate then fill=0
	fill_adi=1
	fill_zero=0
	

	if keyword_set(rho) then begin
		annmode=1
		annmode_inout=[max([0.,(rho/0.01225) -13.]),(rho/0.01225)+14.]
		if annmode_inout[1] lt 20. then n_ang=1
	endif

	
	remove_bads=1 ;remove use defined bads
		bads=[]
		;e.g., bads=[15,[16:35],[45:53],[88:94],[97:104]]-1	

	notify_adiklip=0

	
do_klip_sdi=0	;will run after adi step, unless use_adi is = 0 it will use the adiklipped frames
	bin_sdi=1
	use_adi=1
	notify_sdi=0 & identify_sdi=0
	
do_klip_rdi=0	;will run after adi step, unless use_adi is = 0
	bin_rdi=3
	ref_obj='S47'
	
	
if keyword_set(suffix) then extra_suffix=suffix else extra_suffix=''
;if use_inject then suffix='_inj'+extra_suffix else suffix=''+extra_suffix
if use_inject then suffix='_inj' else suffix=''+extra_suffix

	;optimizing injection parameters?
if nruns gt 1 then test_residuals=1 else test_residuals=0

;############################################################################
;####################[     Object Specific Inputs     ]######################
;############################################################################

;f obj eq 'S24' then comb_type='median'
;if obj eq 'S34' then comb_type='median'

;if obj eq 'S59' then n_nag=1.

;if obj eq 'S12' then angsep=0.5
if obj eq 'SA23' then disp_type='linear'
if obj eq 'SA27' then disp_type='linear'
if obj eq 'SA28' then disp_type='linear'
if obj eq 'Sirius' then disp_type='linear'
if obj eq 'Sirius' then click_center=1

;if obj eq 'SA29' then disp_type='linear'

if obj eq 'S25' then begin
	
	;do_rotate=0 & identify_adi=0 & notify_adi=0 & 	notify_sdi=0 & identify_sdi=0
	;annmode=1 & if hyp then annmode=0	;annular mode?
	;	annmode_inout=[5,35]
	;n_ang=1
	;k_klip=2
	;angsep=0.


	;contrast_left=  2E-3
	;contrast_right=  2E-3


	;rhop=0.15; + 0.001*runs	;right 0.825, left: 0.830	
	;theta=148.-20.

	;planet_r=[rhop,rhop,rhop]	
	;planet_theta=[theta+180.,theta+90.,theta-90.]+90.

	;planet_contrast_left=[contrast_left,contrast_left,contrast_left]
	;planet_contrast_right=planet_contrast_left

endif



if obj eq 'S25_2' then begin
	
	;do_rotate=0 & identify_adi=0 & notify_adi=0 & 	notify_sdi=0 & identify_sdi=0
	;annmode=0 & if hyp then annmode=0	;annular mode?
		annmode_inout=[5,35]
	;n_ang=4
	;k_klip=2
	;angsep=0.

;first set of planets for b

	;contrast_left=  1E-3
	;contrast_right=  1E-3


	;rhop=0.12; + 0.001*runs	;right 0.825, left: 0.830	
	;theta=148.

	;planet_r=[rhop,rhop,rhop]	
	;planet_theta=[theta+180.,theta+90.,theta-90.]+90.

	;planet_contrast_left=[contrast_left,contrast_left,contrast_left]
	;planet_contrast_right=planet_contrast_left

;;;


endif



if obj eq '51Eri' then begin
	high_pass_width=13.
	k_klip=7 & k_klip_sdi=7	
	angsep=0.7
	n_ang=6
	do_rotate=0 & identify_adi=0 & notify_adi=0 & 	notify_sdi=1 & identify_sdi=1
	do_klip_sdi=1
	annmode=1 & annmode_inout=[0.25,0.65]/0.01225
	fill=0
endif

if obj eq '51Eri_3' then begin
	high_pass_width=13.
	k_klip=7 & k_klip_sdi=7	
	angsep=0.7
	n_ang=6
	do_rotate=0 & identify_adi=0 & notify_adi=0 & 	notify_sdi=1 & identify_sdi=1
	do_klip_sdi=1
	annmode=1 & annmode_inout=[0.25,0.65]/0.01225
	fill=0
endif
if obj eq 'SA10' then bads=[2]
if obj eq 'S61' then bads=[1:31]-1	;binary saturating detector, switched to ND1.0

if obj eq 'S12' or obj eq 'S12_2' then auto_frame_hsize=20	;intended to center on the binary		
if obj eq 'S13' then auto_frame_hsize=20  	;intended to center on the binary				
if obj eq 'S27' then auto_frame_hsize=20  	;intended to center on the binary				
if obj eq 'S24' then auto_frame_hsize= 5
if obj eq 'S32' then auto_frame_hsize= 9		
if obj eq 'S34' then auto_frame_hsize= 4

if obj eq 'HIP67036_C' then auto_frame_hsize= 12		

if obj eq 'S23_2' then cen_filter=1
if obj eq 'S32' then cen_filter=1

if obj eq 'MWC758' then begin
	auto_cen=0
	annmode=1 & annmode_inout=[10,70]
endif
;if obj eq 'SA10' then field_stab_cen=1
if obj eq 'SA20' then field_stab_cen=1
;if obj eq 'SA20' then auto_frame_hsize=5.
if obj eq 'S1' then use_second_cen=1	
if obj eq 'S2' then use_second_cen=0
if obj eq 'S11' then use_second_cen=1		
if obj eq 'S12' or obj eq 'S12_2' then use_second_cen=0	
if obj eq 'S13' then use_second_cen=0	
if obj eq 'S14_C2' then use_second_cen=1	
if obj eq 'S23_2' then use_second_cen=1
if obj eq 'S18' then use_second_cen=1
if obj eq 'S21' then use_second_cen=1
if obj eq 'S24' then use_second_cen=1
if obj eq 'S28' then use_second_cen=1
if obj eq 'S29' then use_second_cen=1
if obj eq 'S30' or obj eq 'HIP67036_C' then use_second_cen=1	
if obj eq 'S32' then use_second_cen=1	
if obj eq 'S31' then use_second_cen=1	
if obj eq 'S34' then use_second_cen=1	;tight binary
if obj eq 'HIP79366_C' then use_second_cen=1
if obj eq 'HIP80059_C' then use_second_cen=1
if obj eq 'S7_C' then use_second_cen=1
if obj eq 'S22_6' then use_second_cen=1

if obj eq 'S12'  then field_stab_cen=1
if obj eq 'S12_2' then field_stab_cen=1	
if obj eq 'S13' then field_stab_cen=1	
if obj eq 'S14' or obj eq 'S14_2' or obj eq 'S14_C2' then field_stab_cen=1	
if obj eq 'S24' then field_stab_cen=1	
if obj eq 'S27' then field_stab_cen=1	
if obj eq 'S30' or obj eq 'HIP67036_C' or obj eq 'HIP67036_X' then field_stab_cen=1	
if obj eq 'S34' then field_stab_cen=1
if obj eq 'S56' then field_stab_cen=1
if obj eq 'S57' then field_stab_cen=1
if obj eq 'S22_6' or obj eq 'S22_10' then field_stab_cen=1

if field_stab_cen then field_stab_clean=1
if obj eq 'S6' or obj eq 'S6_2' or obj eq 'S6_3' then field_stab_clean=1
if obj eq 'S31' then field_stab_clean=1
if obj eq 'S20' then field_stab_clean=1
if obj eq 'S15' then field_stab_clean=1
if obj eq 'S27' then field_stab_clean=1
if obj eq 'S33' then field_stab_clean=1
if obj eq 'S50' then field_stab_clean=1
if obj eq 'S8' or obj eq 'S8_2' then field_stab_clean=1
if obj eq 'HIP80059_C' then field_stab_clean=1	
if obj eq 'S1_2' then field_stab_clean=1

if obj eq 'S14' then spot_fit_hwidth=8
if obj eq 'S49' then spot_fit_hwidth=8
if obj eq 'S55' then spot_fit_hwidth=8
;if obj eq 'SA1' then spot_fit_hwidth=5

if obj eq 'HIP67036_C' then spot_fit_hwidth=8

if obj eq 'S14' then spot_smooth=1	

if obj eq 'S2' or obj eq 'S2_2' then click_center=1 ;mickey mouse pattern
if obj eq 'S12' then click_center=1
if obj eq 'S13' then click_center=1
if obj eq 'S27' then click_center=1
if obj eq 'SAO206462_2' then click_center=1
if obj eq 'SAO206462_3' then click_center=1
if obj eq 'PDS70_2' then click_center=1
if obj eq 'SA1' then click_center=1

if obj eq 'SA10' then auto_cen=0
if obj eq 'S34' then auto_cen=0
;if obj eq 'SA1' then auto_cen=0
;if obj eq 'SA10' then begin block_rad_cen=1 & rad_cen_min=0. & rad_cen_max=400. & endif
if obj eq 'SA20' then begin block_rad_cen=1 & rad_cen_min=10. & rad_cen_max=100. & endif
if obj eq 'S1' then begin block_rad_cen=1 & rad_cen_min=0. & rad_cen_max=7. & endif
if obj eq 'S2' or obj eq 'S2_2' then begin block_rad_cen=1 & rad_cen_min=0. & rad_cen_max=7. & endif
if obj eq 'S8' or obj eq 'S8_2' then begin block_rad_cen=1 & rad_cen_min=0. & rad_cen_max=7. & endif
if obj eq 'S11' then begin block_rad_cen=1 & rad_cen_min=0. & rad_cen_max=7. & endif	
if obj eq 'S12' then begin block_rad_cen=1 & rad_cen_min=0. & rad_cen_max=7. & endif			
if obj eq 'S12_2' then begin block_rad_cen=1 & rad_cen_min=0. & rad_cen_max=7. & endif			
if obj eq 'S10_2' then begin block_rad_cen=1 & rad_cen_min=0. & rad_cen_max=100. & endif
if obj eq 'S15' then  begin block_rad_cen=1 & rad_cen_min=0. & rad_cen_max=7. & endif
if obj eq 'S20' then  begin block_rad_cen=1 & rad_cen_min=0. & rad_cen_max=70. & endif
if obj eq 'S23'  then  begin block_rad_cen=1 & rad_cen_min=0. & rad_cen_max=20. & endif
if obj eq 'S24'  then  begin block_rad_cen=1 & rad_cen_min=0. & rad_cen_max=20. & endif
if obj eq 'S27'  then  begin block_rad_cen=1 & rad_cen_min=0. & rad_cen_max=7. & endif
if obj eq 'S31' then  rad_cen_max=7.
if obj eq 'S32' then  begin block_rad_cen=1 & rad_cen_min=0. & rad_cen_max=7. & endif
if obj eq 'S33' then  begin block_rad_cen=1 & rad_cen_min=0. & rad_cen_max=7. & endif
if obj eq 'S50' then  begin block_rad_cen=1 & rad_cen_min=0. & rad_cen_max=10. & endif
if obj eq 'S61' then  begin block_rad_cen=1 & rad_cen_min=0. & rad_cen_max=7. & endif
if obj eq 'S22_6' then  begin block_rad_cen=1 & rad_cen_min=50. & rad_cen_max=400. & endif
if obj eq 'S22_10' then  begin block_rad_cen=1 & rad_cen_min=50. & rad_cen_max=400. & endif
if obj eq 'HIP79366_C' then block_rad_cen=0		
if obj eq 'HIP80059_C' then rad_cen_max=10
if obj eq 'SA10' then cross_thresh=0.92

if obj eq 'S32' then remove_wild=1
if obj eq 'HIP79366_C' then remove_wild=0

if obj eq 'S1_2' then cross_thresh=0.5
if obj eq 'S10' then cross_thresh=0.99
if obj eq 'S32' then cross_thresh=0.75
if obj eq 'S37_2' then cross_thresh=0.8
if obj eq 'S56' then cross_thresh=0.8	

;if obj eq 'S20' then comb_type='median'

if obj eq 'HD100546_3' or obj eq 'HD100546_4' or obj eq 'HD100546_5' then field_stab_clean=1
if obj eq 'HD100546_3' or obj eq 'HD100546_4' or obj eq 'HD100546_5' then cross_thresh=0.7
if obj eq 'HD100546_3' or obj eq 'HD100546_4' or obj eq 'HD100546_5' then auto_cen=0

if obj eq 'BetaCir' then field_stab_clean = 1
if obj eq 'BetaCir' then field_stab_cen = 1
;if obj eq 'S20' then do_klip=1
;if obj eq 'S20' then do_klip_sdi=1


;Default correction factors. These must be modified for accurate photometry!!!
;A good way to do this is to inject fake planets to measure the throughput
;Then correction_factor = 1/throughput

left_derot_correction_factor=1.0
left_adi_correction_factor=1.0
left_klip_correction_factor=1.0

right_derot_correction_factor=1.0
right_adi_correction_factor=1.0
right_klip_correction_factor=1.0

total_derot_correction_factor=1.0
total_adi_correction_factor=1.0
total_klip_correction_factor=1.0



if 1 eq 0 then begin ;turning off
if obj eq 'S22_12' or obj eq 'S22' or obj eq 'S22_2' or obj eq 'S22_3' or obj eq 'S22_4' or obj eq 'S22_5' or obj eq 'S22_6' or obj eq 'S22_7' or obj eq 'S22_8' or obj eq 'S22_9' or obj eq 'S22_10' or obj eq 'S22_11'  then begin

	annmode=1 
	;if obj eq 'S22' or obj eq 'S22_2' or obj eq 'S22_3' or obj eq 'S22_4' then annmode_inout=[55,78] $
	;	else annmode_inout=[53,76]

	annmode_inout=[53,78]

	contrast_left=  0.5E-5
	contrast_right=  0.5E-5

	if obj eq 'S22_12' or obj eq 'S22' or obj eq 'S22_2' or obj eq 'S22_3' or obj eq 'S22_4' or obj eq 'S22_10' then contrast_left=1E-5
	if obj eq 'S22_12' or obj eq 'S22' or obj eq 'S22_2' or obj eq 'S22_3' or obj eq 'S22_4' or obj eq 'S22_10' then contrast_right=1E-5

	theta=193.; -  0.1* runs 
	rhop=0.8; + 0.001*runs	;right 0.825, left: 0.830	

	planet_r=[rhop,rhop,rhop,rhop,rhop,rhop,rhop,rhop]
	planet_theta=[theta,theta+180.,theta+90.,theta-90.,theta+45.,theta+180.+45.,theta+90.+45.,theta-90.+45.]

	planet_contrast_left=[0.,contrast_left,contrast_left,contrast_left,contrast_left,contrast_left,contrast_left,contrast_left]
	planet_contrast_right=planet_contrast_left

	;do_rotate=0
endif
endif


if obj eq 'S22_10' then begin
	angsep=0.25 & anglemax=10.
	high_pass_width=11.
endif

if obj eq 'S22_6' then begin
	angsep=0.25 & anglemax=10.
	high_pass_width=11.
endif

;if obj eq 'S27'  then use_gauss=0


if obj eq 'SA2_3' then begin
	block_rad_cen=1 & rad_cen_min=50. & rad_cen_max=400.
endif


if obj eq 'HD131488' then begin
	high_pass_width=20
endif

if not annmode then annmode_inout=[0,511]
if not fill and not annmode then annmode_inout=[0,final_image_half_size-1] ;these lines are helpful for find_sources
if annmode then fill=0

;default display parameters


;   0         None           No scaling whatsoever is done.
;   1         Linear         scaled = BytScl(image, MIN=minValue, MAX=maxValue)
;   2         Clip           A histogram stretch, with a percentage of pixels clipped at both the top and bottom
;   3         Gamma          scaled = GmaScl(image, MIN=minValue, MAX=maxValue, Gamma=gamma)
;   4         Log            scaled = LogScl(image, MIN=minValue, MAX=maxValue, Mean=mean, Exponent=exponent)
;   5         Asinh          scaled = AsinhScl(image, MIN=minValue, MAX=maxValue, Beta=beta)
;   6         SquareRoot     A linear stretch of the square root histogram of the image values.
;   7         Equalization   A linear stretch of the histogram equalized image histogram.
;   8         Gaussian       A Gaussian normal function is applied to the image histogram.
;   9         MODIS          Scaling done in the differential manner of the MODIS Rapid Response Team
;                            and implemented in the Coyote Library routine ScaleModis.
;   10        StdDev         Standard deviation stretch. scaled = SDevScl(image, MULTIPLIER=2).


left_adi_stretch='linear'
left_derot_stretch='linear'
left_klip_stretch='linear'


left_adi_clip_level=0.025
left_klip_clip_level=0.03
left_derot_clip_level=0.015


right_adi_stretch=left_adi_stretch
right_derot_stretch=left_derot_stretch
right_klip_stretch=left_klip_stretch


right_adi_clip_level=left_adi_clip_level
right_klip_clip_level=left_klip_clip_level
right_derot_clip_level=left_derot_clip_level

total_adi_clip_level=left_adi_clip_level
total_klip_clip_level=right_klip_clip_level
total_derot_clip_level=left_derot_clip_level


total_adi_stretch=left_adi_stretch
total_derot_stretch=left_derot_stretch
total_klip_stretch=left_klip_stretch



;############################################################################
;#########################[     END USER INPUT     ]#########################
;############################################################################
;#################[  YOU SHALL NOT EDIT BENEATH THIS LINE  ]#################
;#################(...unless you know what you are doing...)#################
;############################################################################
;############################################################################


;-------------------------[ Defining file paths ]------------------------------

path=strcompress(reduction_path+obj+'/IRDIS/products/',/remove_all)
proc_path=strcompress(reduction_path+obj+'/IRDIS/processed/',/remove_all)
ref_path=strcompress(reduction_path+ref_obj+'/IRDIS/processed/',/remove_all)
cen_path=strcompress(reduction_path+obj+'/IRDIS/center/',/remove_all)
flux_path=strcompress(reduction_path+obj+'/IRDIS/flux/',/remove_all)


leftfiles=FILE_SEARCH(path,'*left.fits',COUNT=leftfilecount)
rightfiles=FILE_SEARCH(path,'*right.fits',COUNT=rightfilecount)


leftfiles_cen=FILE_SEARCH(cen_path,'*left.fits',COUNT=leftcenfilecount)
rightfiles_cen=FILE_SEARCH(cen_path,'*right.fits',COUNT=rightcenfilecount)


leftfiles_flux=FILE_SEARCH(flux_path,'*left.fits',COUNT=leftfluxfilecount)
rightfiles_flux=FILE_SEARCH(flux_path,'*right.fits',COUNT=rightfluxfilecount)


;-------------------------[ Begin Reduction Task ]------------------------------
if do_flux eq 1 then begin
	
	left_flux=readfits(leftfiles_flux[flux_file_num],lhdr)
	right_flux=readfits(rightfiles_flux[flux_file_num],rhdr)
	
	
	
	dbf=esopar(lhdr,'ESO INS1 OPTI2 NAME ')
	if strpos(dbf,'Y23') ne -1 then dual_filter='Y23'

	if strpos(dbf,'J23') ne -1 then dual_filter='J23'
	if strpos(dbf,'H23') ne -1 then dual_filter='H23'
	if strpos(dbf,'H34') ne -1 then dual_filter='H34'
	if strpos(dbf,'K12') ne -1 then dual_filter='K12'

	if strpos(dbf,'Y23') ne -1 then dual_filter='Y23'
	
	if strpos(dbf,'J23') ne -1 then dual_filter='J23'
	if strpos(dbf,'H23') ne -1 then dual_filter='H23'
	if strpos(dbf,'H34') ne -1 then dual_filter='H34'
	if strpos(dbf,'K12') ne -1 then dual_filter='K12'


	if strpos(dbf,'CLEAR') ne -1 then begin;dual_filter='K12'
		print, 'Dual filter clear. Checking for single filter.'
		dbf2=esopar(lhdr,'ESO INS1 FILT NAME ')
		if strpos(dbf2,'Ks') ne -1 then dual_filter='Ks'
		if strpos(dbf2,'B_Y') ne -1 then dual_filter='Y'
		if strpos(dbf2,'B_J') ne -1 then dual_filter='J'
		if strpos(dbf2,'B_H') ne -1 then dual_filter='H'
	
	endif
	
	
	print, dual_filter
	if dual_filter eq 'Ks' then wavel=2.182
	if dual_filter eq 'Ks' then waver=2.182


	if dual_filter eq 'Y' then wavel=1.043
	if dual_filter eq 'Y' then waver=1.043
	
	
	if dual_filter eq 'J' then wavel=1.245
	if dual_filter eq 'J' then waver=1.245
	
	
	if dual_filter eq 'H' then wavel=1.625
	if dual_filter eq 'H' then waver=1.625

	
	if dual_filter eq 'Ks' then wavel=2.182
	if dual_filter eq 'Ks' then waver=2.182

	if dual_filter eq 'Y23' then wavel=1.022

	if dual_filter eq 'K12' then wavel=2.1025
	if dual_filter eq 'J23' then wavel=1.1895
	if dual_filter eq 'H23' then wavel=1.5888

	if dual_filter eq 'Y23' then waver=1.076

	if dual_filter eq 'K12' then waver=2.255
	if dual_filter eq 'J23' then waver=1.2698
	if dual_filter eq 'H23' then waver=1.6671

	

	;reading the neutral-density filter ID for use in determining the next scale factor 
	filt = esopar(lhdr,'ESO INS4 COMB IND')
	if isa(filt, 'string') eq 0 then filt = esopar(lhdr,'ESO INS4 FILT2 NAME')  
	readcol, reduction_path+'calib/SPHERE_CPI_ND.dat.txt', wavend, nd0, nd1, nd2, nd35
	;readcol, reduction_path+'calib/SPHERE_IRDIS_D_K12.dat.txt', wavek12, dk1, dk2

	;these are the values at the centers of the K1 and K2 filter. 
	;The ND filters are essentially flat across the bandpasses.
	;if strpos(filt, 'ND_N_3.5') ne -1 or STRPOS(filt, 'ND_3.5') ne -1 then begin 
	;	filtscale1=0.001390790 & filtscale2=0.001428780 & endif
	;if strpos(filt, 'ND_N_2.0') ne -1 or STRPOS(filt, 'ND_2.0') ne -1 then begin 
	;	filtscale1=0.022148401 & filtscale2=0.021858200 & endif
	;if strpos(filt, 'ND_N_1.0') ne -1 or STRPOS(filt, 'ND_1.0') ne -1 then begin 
	;	filtscale1=0.148992002 & filtscale2=0.139275998 & endif
	;Print,'Off-center ND filter: ',filt
	
	
	
	;new code
	
	if strpos(filt, 'ND_N_3.5') ne -1 or STRPOS(filt, 'ND_3.5') ne -1 then nd=nd35
	if strpos(filt, 'ND_N_2.0') ne -1 or STRPOS(filt, 'ND_2.0') ne -1 then nd=nd2
	if strpos(filt, 'ND_N_1.0') ne -1 or STRPOS(filt, 'ND_1.0') ne -1 then nd=nd1
	if strpos(filt, 'ND_N_0.0') ne -1 or STRPOS(filt, 'ND_0.0') ne -1 or STRPOS(filt, 'OPEN') ne -1 then nd=nd0

	Print,'Off-center ND filter: ',filt
	nearest_wavel=nearest_element(wavel*1000.,wavend,wave_posl)
	nearest_waver=nearest_element(waver*1000.,wavend,wave_posr)
	
	filtscale1=nd[wave_posl]
	filtscale2=nd[wave_posr]

	
	

	print, 'Filtscale 1:',filtscale1
	print, 'Filtscale 2:',filtscale2

if n_elements(size(left_flux)) gt 5 then begin
	 for ij=0, (size(left_flux))(3) -1 do left_flux[*,*,ij]=left_flux[*,*,ij]/float(filtscale1)
endif else  left_flux=left_flux/float(filtscale1)
if n_elements(size(right_flux)) gt 5 then begin	
	for ij=0, (size(right_flux))(3) -1 do right_flux[*,*,ij]=right_flux[*,*,ij]/float(filtscale2)	
endif else right_flux=right_flux/float(filtscale2)
	
	date=fxpar(lhdr,'DATE')
	
	dit=fxpar(lhdr,'EXPTIME') ;reading DIT and normalizing to a DIT of 1 second
	print, 'DIT: ',dit
	
if n_elements(size(left_flux)) gt 5 then begin
	 for ij=0, (size(left_flux))(3) -1 do left_flux[*,*,ij]=left_flux[*,*,ij]/float(dit)
endif else  left_flux=left_flux/float(dit)
if n_elements(size(right_flux)) gt 5 then begin	
	for ij=0, (size(right_flux))(3) -1 do right_flux[*,*,ij]=right_flux[*,*,ij]/float(dit)	
endif else right_flux=right_flux/float(dit)

;it gets averaged anyways, duplicating into cube form to make the rest of the coding easier
if n_elements(size(left_flux)) le 5 then left_flux=[[[left_flux]],[[left_flux]]]
if n_elements(size(right_flux)) le 5 then right_flux=[[[right_flux]],[[right_flux]]]
	;bad pixels in flux
	if do_flux_bad_pixel_removal eq 1 then begin 
	  			print, ' Bad pixel removal - left'

		cube=left_flux
		badmap=cube

		badmap[*]=1

		 badsig=2.5

		 smallcube = cube; We don't need all the frames to find the bad pixels
 			cs=8
  			for ix=cs, (size(cube))(1)-1-cs do begin
   			 	 for iy=cs, (size(cube))(2)-1-cs do begin
        		   box = smallcube[ix-cs:ix+cs,iy-cs:iy+cs,*]
         		   this =  box[cs,cs,*]  ; this pixel
         		   if mean(this) lt -100 then badmap[ix,iy]=0
         		   box[cs,cs,*] = !values.f_nan
          		   if mean(this) gt (mean(box,/nan)+badsig*stddev(box,/nan)) $
          		   		or mean(this) lt (mean(box,/nan)-badsig*stddev(box,/nan)) $
          		   		then bad=1 else bad=0
          		   if bad then badmap[ix,iy] = 0           
      			  endfor
    	 		  ;print, ' Left Column : ', ix
     		endfor
     

			; BAD PIXEL REMOVAL
  			print
  			for ii=0, (size(cube))(3)-1 do  begin
				;print, 'Working on frame ', ii, ' / ', (size(cube))(3)-1
  				cube[*,*,ii] = maskinterp(cube[*,*,ii],badmap,3,6,"splinterp")
  			endfor
  			  			print, ' Bad pixel removal - right'

  			left_flux=cube
  			
  				cube=right_flux
				badmap=cube

				badmap[*]=1

				 badsig=2.5

				 smallcube = cube; We don't need all the frames to find the bad pixels
 			cs=8
  			for ix=cs, (size(cube))(1)-1-cs do begin
   			 	 for iy=cs, (size(cube))(2)-1-cs do begin
        		   box = smallcube[ix-cs:ix+cs,iy-cs:iy+cs,*]
         		   this =  box[cs,cs,*]  ; this pixel
         		   if mean(this) lt -100 then badmap[ix,iy]=0
         		   box[cs,cs,*] = !values.f_nan
          		   if mean(this) gt (mean(box,/nan)+badsig*stddev(box,/nan)) $
          		   		or mean(this) lt (mean(box,/nan)-badsig*stddev(box,/nan)) $
          		   		then bad=1 else bad=0
          		   if bad then badmap[ix,iy] = 0           
      			  endfor
    	 		 ; print, ' Right Column : ', ix
     		endfor
     

			; BAD PIXEL REMOVAL
  			print
  			for ii=0, (size(cube))(3)-1 do  begin
				;print, 'Working on frame ', ii, ' / ', (size(cube))(3)-1
  				cube[*,*,ii] = maskinterp(cube[*,*,ii],badmap,3,6,"splinterp")
  			endfor
  			
  			
    		right_flux=cube
    	endif ;bad pixels if

	
	
	;corrects for distortion by 1.006 in y direction
	 print, ' Correcting for anamorphic distortion...'
	 cube=left_flux
	 print, size(cube)

	 scube = size(cube)
	 sn = size(congrid(cube[*,*,0], scube[1],scube[2]*1.006))
	 nz=scube[3]
	 print, nz
	ncube = fltarr(sn[1], sn[2], nz)
	for ii=0, nz-1 do ncube[*,*,ii] = congrid(cube[*,*,ii], scube[1],  scube[2]*1.006) 
 	left2 =ncube
 	 left_flux=ncube
 	ncube =0.
 
 	cube=right_flux
 	scube = size(cube)
	 sn = size(congrid(cube[*,*,0], scube[1],scube[2]*1.006))
	 ncube = fltarr(sn[1], sn[2], nz)
	 for ii=0, nz-1 do ncube[*,*,ii] = congrid(cube[*,*,ii], scube[1],  scube[2]*1.006) 
	 right2 =ncube
	right_flux=ncube
	 ncube =0.
	
	
		xsz=(size(right_flux))(1)
		ysz=(size(right_flux))(2) 

		window,0,xs=400,ys=400
		window,1,xs=400,ys=400
		window,2,xs=400,ys=400
		window,3,xs=400,ys=400
		
		WSET,0

		loadct,0,/silent
		if not auto_flux_cen then print, 'Click on the center of the star:'
		tvscl,-rot(left_flux[float(xsz)/2.-1.-200:float(xsz)/2.-1.+200-1.,$
						float(ysz)/2.-1.-200:float(ysz)/2.-1.+200-1.,0],0,2.0)
		if not auto_flux_cen then cursor, x, y, /DEVICE, /DOWN
		
		
		
		if auto_flux_cen then result=gauss2dfit(rot(left_flux[float(xsz)/2.-1.-200:float(xsz)/2.-1.+200-1., float(ysz)/2.-1.-200:float(ysz)/2.-1.+200-1.,0],0,2.0), A)
		if auto_flux_cen then x=A[4] & if auto_flux_cen then y=A[5]
		
		x=float(x)/2.
		y=float(y)/2.
		print, 'Center at ',x,y
		Print, 'Determining fine center with centroid'
				print, x,y

		cntrd, left_flux[*,*,0],float(xsz)/2.-1.-100.+x,float(ysz)/2.-1.-100.+y,x,y, 4.35
		print, 'Real image coords: ',float(xsz)/2.-1.-100.+x,float(ysz)/2.-1.-100.+y
		print, 'Centering by ',x-512.,y-512.
		for nn=0,(size(left_flux))(3) -1 do left_flux[*,*,nn]=fshift(left_flux[*,*,nn],-(x-512.),-(y-512.))

		tvscl,-rot(left_flux[float(xsz)/2.-1.-200:float(xsz)/2.-1.+200-1.,$
						float(ysz)/2.-1.-200:float(ysz)/2.-1.+200-1.,0],0,2.0)
		!p.color=0 & !p.charsize=1.5
		xyouts, [30], [360], 'IRDIS Flux Calibration - '+obj,/device
		
		      plots,[200,200]+2.5,[0,401],color=250,/device
		      plots,[0,401],[200,200]-2.5,color=250,/device

		
				WSET,1

		
		if not auto_flux_cen then print, 'Click on the center of the star:'
		tvscl,-rot(right_flux[float(xsz)/2.-1.-200:float(xsz)/2.-1.+200-1.,$
						float(ysz)/2.-1.-200:float(ysz)/2.-1.+200-1.,0],0,2.0)
		if not auto_flux_cen then cursor, x, y, /DEVICE, /DOWN
		
		
		if auto_flux_cen then result=gauss2dfit(rot(right_flux[float(xsz)/2.-1.-200:float(xsz)/2.-1.+200-1.,$
						float(ysz)/2.-1.-200:float(ysz)/2.-1.+200-1.,0],0,2.0), A)
		if auto_flux_cen then 	x=A[4] & if auto_flux_cen then y=A[5]
		
		x=float(x)/2.
		y=float(y)/2.
		print, 'Center at ',x,y
		Print, 'Determining fine center with centroid'
		cntrd, right_flux[*,*,0],float(xsz)/2.-1.-100.+x,float(ysz)/2.-1.-100.+y,x,y, 4.35
		print, x,y
		print, 'Centering by ',x-512.,y-512.
		for nn=0,(size(right_flux))(3)-1 do right_flux[*,*,nn]=fshift(right_flux[*,*,nn],-(x-512.),-(y-512.))

		tvscl,-rot(right_flux[float(xsz)/2.-1.-200:float(xsz)/2.-1.+200-1.,$
						float(ysz)/2.-1.-200:float(ysz)/2.-1.+200-1.,0],0,2.0)
			!p.color=0 & !p.charsize=1.5
		xyouts, [30], [360], 'IRDIS Flux Calibration - '+obj,/device
		
		      plots,[200,200]+2.5,[0,401],color=250,/device
		      plots,[0,401],[200,200]-2.5,color=250,/device
		
		;hak
		;cgcleanup
		
		print, 'Reducing flux calibrations....'
		
		for jj=0, (size(right_flux))(3)-1 do begin
			; rcorr=crosscorr(right_flux[*,*,0],right_flux[*,*,jj],pmax, /cut)
			;if pmax(0) ge 0. then dx = (size(right_flux))(1)/2.-pmax(0) else dx = -((size(right_flux))(1)/2.)+abs(pmax(0))
			;if pmax(1) ge 0. then dy = (size(right_flux))(2)/2.-pmax(1) else dy = -((size(right_flux))(2)/2.)+abs(pmax(1)) ;si was ww everywhere before subimages
			;if jj eq 0 then dx=0
			;if jj eq 0 then dy=0
			
			cntrd,right_flux[*,*,jj],512.,512.,x,y,4.2
			
			dx=511.5-x
			dy=511.5-y
			
			
			right_flux[*,*,jj]=fshift(right_flux[*,*,jj],dx,dy)

		endfor
		
		for jj=0, (size(left_flux))(3)-1 do begin
			 ;rcorr=crosscorr(left_flux[*,*,0],left_flux[*,*,jj],pmax, /cut)
			;if pmax(0) ge 0. then dx = (size(left_flux))(1)/2.-pmax(0) else dx = -((size(left_flux))(1)/2.)+abs(pmax(0))
			;if pmax(1) ge 0. then dy = (size(left_flux))(2)/2.-pmax(1) else dy = -((size(left_flux))(2)/2.)+abs(pmax(1)) ;si was ww everywhere before subimages
			;if jj eq 0 then dx=0
			;if jj eq 0 then dy=0
			
			
			cntrd,left_flux[*,*,jj],512.,512.,x,y,4.2
			
			dx=511.5-x
			dy=511.5-y
			
			left_flux[*,*,jj]=fshift(left_flux[*,*,jj],dx,dy)

		endfor
		
		sz=50.
		subleft=left_flux[512.-sz:512.+sz-1.,512.-sz:512.+sz-1.,*]
		subright=right_flux[512.-sz:512.+sz-1.,512.-sz:512.+sz-1.,*]
		
		

		
		medarr,subleft,left_flux
		medarr,subright,right_Flux

		
		
		
		;left_flux=fshift(left_flux,-0.5,-0.5)
		;right_flux=fshift(right_flux,-0.5,-0.5)
		
		writefits,reduction_path+obj+'/IRDIS/flux/'+obj+'_left_median.fits',left_flux,lhdr
		writefits,reduction_path+obj+'/IRDIS/flux/'+obj+'_right_median.fits',right_flux,lhdr

		writefits,reduction_path+obj+'/IRDIS/flux/'+obj+'_left+right_median.fits',(left_flux+right_flux)/2.,lhdr


endif   ; do_flux if


;-------------------------[ Begin Reduction Task ]------------------------------
if do_center eq 1 then begin

first_left=readfits(leftfiles_cen[0],lhdr0)

if n_elements(size(first_left)) gt 5 then medarr,first_left,first_left

if leftcenfilecount le 1 then use_second_cen=0

if do_multi_cen ne 1 then begin leftcenfilecount=1 & rightcenfilecount=1 & endif

if leftcenfilecount gt 1 or use_second_cen then last_left=readfits(leftfiles_cen[1],lhdr1)
if leftcenfilecount gt 1 or use_second_cen then medarr, last_left,last_left

first_right=readfits(rightfiles_cen[0],rhdr0)
if n_elements(size(first_right)) gt 5 then medarr, first_right,first_right
if rightcenfilecount gt 1 or use_second_cen then last_right=readfits(rightfiles_cen[1],rhdr1)
if rightcenfilecount gt 1 or use_second_cen then medarr, last_right, last_right

if use_second_cen then first_left=last_left
if use_second_cen then first_right=last_right

if use_second_cen then begin lhdr0=lhdr1 & rhdr0=rhdr1 & endif

if rightcenfilecount gt 1 and leftcenfilecount gt 1 then cen_steps_max = 3 $
	else cen_steps_max=1
	
	
	
	;corrects for distortion by 1.006 in y direction
	 print, ' Correcting for anamorphic distortion...'
	 cube=first_left
	 print, size(cube)
	 scube = size(cube)
	 sn = size(congrid(cube, scube[1],scube[2]*1.006))
	 nz=scube[3]
	 print, nz
	ncube = fltarr(sn[1], sn[2])
	ncube = congrid(cube, scube[1],  scube[2]*1.006) 
 	left2 =ncube
 	ncube =0.
 	first_left=cube
 	
 	 cube=first_right
	 print, size(cube)
	 scube = size(cube)
	 sn = size(congrid(cube, scube[1],scube[2]*1.006))
	 nz=scube[3]
	 print, nz
	ncube = fltarr(sn[1], sn[2])
	ncube = congrid(cube, scube[1],  scube[2]*1.006) 
 	left2 =ncube
 	ncube =0.
 	first_right=cube
 	
 	if rightcenfilecount gt 1 and leftcenfilecount gt 1 then begin
 		;corrects for distortion by 1.006 in y direction
	 print, ' Correcting for anamorphic distortion...'
	 cube=last_left
	 print, size(cube)
	 scube = size(cube)
	 sn = size(congrid(cube, scube[1],scube[2]*1.006))
	 nz=scube[3]
	 print, nz
	ncube = fltarr(sn[1], sn[2])
	ncube = congrid(cube, scube[1],  scube[2]*1.006) 
 	left2 =ncube
 	ncube =0.
	 last_left=cube
 	
 	 cube=last_right
	 print, size(cube)
	 scube = size(cube)
	 sn = size(congrid(cube, scube[1],scube[2]*1.006))
	 nz=scube[3]
	 print, nz
	ncube = fltarr(sn[1], sn[2])
	ncube = congrid(cube, scube[1],  scube[2]*1.006) 
 	left2 =ncube
 	ncube =0.
 	last_right=cube
 	
 	first_left=first_left[0:1023,0:1023,*]
 	first_right=first_right[0:1023,0:1023,*]
 	last_left=last_left[0:1023,0:1023,*]
 	last_right=last_right[0:1023,0:1023,*]
 	endif
	
for cen_steps=0,cen_steps_max do begin


if cen_steps eq 0 then image=first_left
if cen_steps eq 1 then image=first_right

if cen_steps eq 2 then image=last_left
if cen_steps eq 3 then image=last_right




		sz=size(image)
		nframes=sz[3]
		xsz=sz[1] & ysz=sz[2]
		
		
		if cen_steps eq 0 then	WSET,2
		if cen_steps eq 1 then	WSET,3

		if dual_filter ne 'K12' then auto_cen=0
		
		new_cen=512.	;should stay as 512. as click_cen uses cntrd
				;spot fit center later uses 511.5


		if not auto_cen then print, 'Click on the center of the star:' else print, 'Determine star center automatically....'
		if disp_type eq 'alog' then tvscl,-alog(rot(image[float(xsz)/2.-1.-200:float(xsz)/2.-1.+200-1.,$
						float(ysz)/2.-1.-200:float(ysz)/2.-1.+200-1.],0,2.0))
		if disp_type eq 'asin' then tvscl,-asin(rot(image[float(xsz)/2.-1.-200:float(xsz)/2.-1.+200-1.,$
						float(ysz)/2.-1.-200:float(ysz)/2.-1.+200-1.],0,2.0))
		if disp_type eq 'linear' then tvscl,-(rot(image[float(xsz)/2.-1.-200:float(xsz)/2.-1.+200-1.,$
						float(ysz)/2.-1.-200:float(ysz)/2.-1.+200-1.],0,2.0))
		!p.color=0 & !p.charsize=1.5
		xyouts, [30], [360], 'IRDIS Star Centering - '+obj,/device
		if auto_cen eq 1 then begin
			if cen_steps eq 0 or cen_steps eq 2 then x=480. else x=483.
			if cen_steps eq 0 or cen_steps eq 2 then y=525. else y=512.


		
			frame=image[x-auto_frame_hsize:x+auto_frame_hsize-1,y-auto_frame_hsize:y+auto_frame_hsize-1,0]
			
			
			result=gauss2dfit(frame, A)
			x=A[4]+x-auto_frame_hsize &  y=A[5]+y-auto_frame_hsize

			;cntrd, image[*,*,0],x,y,x,y, 4.35

		endif else begin
			cursor, x, y, /DEVICE, /DOWN
			x=float(x)/2. &	y=float(y)/2.
			
			if init_centroid then Print, 'Determining fine center with centroid'
			if init_centroid then cntrd, image[*,*,0],float(xsz)/2.-1.-100.+x,float(ysz)/2.-1.-100.+y,x,y, 4.35
			if not init_centroid then begin x=float(xsz)/2.-1.-100.+x & y=float(ysz)/2.-1.-100.+y & endif
		
			
		endelse
		print, 'Centroid: ',x,y
		if auto_cen ne 1 then print, 'Real coords: ',float(xsz)/2.-1.-100.+x,float(ysz)/2.-1.-100.+y
			;x=float(xsz)/2.-1.-100.+x
			;y=float(ysz)/2.-1.-100.+y
				print, 'Center at ',x,y


		print, 'Centering by ',x-new_cen,y-new_cen
		xinit=x-new_cen &	yinit=y-new_cen
		image=fshift(image,-(x-new_cen),-(y-new_cen))

		if disp_type eq 'alog' then tvscl,-alog(rot(image[float(xsz)/2.-1.-200:float(xsz)/2.-1.+200-1.,$
						float(ysz)/2.-1.-200:float(ysz)/2.-1.+200-1.,0],0,2.0))
		if disp_type eq 'asin' then tvscl,-asin(rot(image[float(xsz)/2.-1.-200:float(xsz)/2.-1.+200-1.,$
						float(ysz)/2.-1.-200:float(ysz)/2.-1.+200-1.,0],0,2.0))
		if disp_type eq 'linear' then tvscl,-(rot(image[float(xsz)/2.-1.-200:float(xsz)/2.-1.+200-1.,$
						float(ysz)/2.-1.-200:float(ysz)/2.-1.+200-1.,0],0,2.0))
														
		if click_center ne 1 then begin			
		xyouts, [30], [360], 'IRDIS Star Centering - '+obj,/device

		cut=float(xsz)/2.-1.-100.
		
		;hak		;hit any key				
		Print,'Click on left spot'
		if auto_cen eq 1 then begin
			if cen_steps eq 0 or cen_steps eq 2 then x1=467. else x1=464.
			if cen_steps eq 0 or cen_steps eq 2 then y1=556. else y1=559.
		endif else begin
			cursor, x1,y1, /DEVICE, /DOWN
			x1=x1/2.+cut & y1=y1/2.+cut
		endelse
		      Print, 'First spot approx. at',x1,y1


		Print,'Click on right spot'
		if auto_cen eq 1 then begin
			if cen_steps eq 0 or cen_steps eq 2 then x2=557. else x2=559.
			if cen_steps eq 0 or cen_steps eq 2 then y2=468. else y2=465.
		endif else begin
			cursor, x2,y2, /DEVICE, /DOWN
			x2=x2/2.+cut & y2=y2/2.+cut
		endelse
		     Print, 'Second spot approx. at',x2,y2
		
      plots,[(x1-cut)*2.,(x2-cut)*2.],[(y1-cut)*2.,(y2-cut)*2.],color=250,/device
		
		;hak
		Print,'Click on top spot'
		if auto_cen eq 1 then begin
			if cen_steps eq 0 or cen_steps eq 2 then x3=467. else x3=467.
			if cen_steps eq 0 or cen_steps eq 2 then y3=464. else y3=464.
		endif else begin
			cursor, x3,y3, /DEVICE, /DOWN
			x3=x3/2.+cut & y3=y3/2.+cut
		endelse
		            Print, 'Third spot approx. at',x3,y3

		;hak
		Print,'Click on bottom spot'
		if auto_cen eq 1 then begin
			if cen_steps eq 0 or cen_steps eq 2 then x4=557. else x4=556.
			if cen_steps eq 0 or cen_steps eq 2 then y4=560. else y4=559.
		endif else begin
			cursor, x4,y4, /DEVICE, /DOWN
			x4=x4/2.+cut & y4=y4/2.+cut
		endelse
		            Print, 'Fourth spot approx. at',x4,y4

	 plots,[(x3-cut)*2.,(x4-cut)*2.],[(y3-cut)*2.,(y4-cut)*2.],color=250,/device
	 
	 
	 	 	 lint,[x1,y1],[x2,y2],[x3,y3],[x4,y4],center1
	 	 	 
	 	 	 if center1 eq !NULL then center1=[x1,y1]
	 	endif else center1=[512.,512.] 	 	 	 
	Print, 'Approximate center: ',center1
	
	if no_spot_fit eq 1 then center=center1
	if click_center ne 1 and no_spot_fit ne 1 then begin
	 Print, 'Fitting 2D Gaussian... Adjusting fit...'


		image2=image-smooth(image,15.)
		if spot_smooth then image2=smooth(image2,3.)
		;image2=image2>0.
		image2[where(image2 lt 0)]=image2[where(image2 lt 0)]/1000000000.
	 
	 ext=spot_fit_hwidth	;half size of box to fit the spot peaks (may need to modify depending
	 		;on bad-pixels and binaries).
	;sub1=image2[fix(x1)-ext:fix(x1)+ext-1,fix(y1)-ext:fix(y1)+ext-1]


	;try making new sub by filling in zero for outside of x1,y1 region
	sub=image2
	delvar,a
	for xx=0,1023 do begin
		for yy=0,1023 do begin
			if sqrt((xx-x1)^2. + (yy-y1)^2.) gt spot_fit_hwidth then sub[xx,yy]=0.
		endfor	
	endfor
	fit1 = mpfit2dpeak(sub,a, /tilt)

      x1 = a[4];+float(fix(x1))-float(ext)
      y1 = a[5];+float(fix(y1))-float(ext)
      
      Print, 'First spot at',x1,y1

	sub=fshift(sub,511.5-x1,511.5-y1)

	if cen_steps eq 0 then writefits,reduction_path+obj+'/IRDIS/center/'+obj+'_left-sub1c.fits',sub
      	if cen_steps eq 1 then writefits,reduction_path+obj+'/IRDIS/center/'+obj+'_right-sub1c.fits',sub
     
	 sub=image2
	delvar,a
	for xx=0,1023 do begin
		for yy=0,1023 do begin
			if sqrt((xx-x2)^2. + (yy-y2)^2.) gt spot_fit_hwidth then sub[xx,yy]=0.
		endfor	
	endfor
	fit1 = mpfit2dpeak(sub,a, /tilt)

      x2 = a[4];+float(fix(x2))-float(ext)
      y2 = a[5];+float(fix(y2))-float(ext)


	sub=fshift(sub,511.5-x2,511.5-y2)
	if cen_steps eq 0 then writefits,reduction_path+obj+'/IRDIS/center/'+obj+'_left-sub2c.fits',sub
      	if cen_steps eq 1 then writefits,reduction_path+obj+'/IRDIS/center/'+obj+'_right-sub2c.fits',sub


     Print, 'Second spot at',x2,y2
      
      
     sub=image2
	delvar,a
	for xx=0,1023 do begin
		for yy=0,1023 do begin
			if sqrt((xx-x3)^2. + (yy-y3)^2.) gt spot_fit_hwidth then sub[xx,yy]=0.
		endfor	
	endfor
	fit1 = mpfit2dpeak(sub,a, /tilt)

      
      x3 = a[4];+float(fix(x3))-float(ext)
      y3 = a[5];+float(fix(y3))-float(ext)


	sub=fshift(sub,511.5-x3,511.5-y3)
	if cen_steps eq 0 then writefits,reduction_path+obj+'/IRDIS/center/'+obj+'_left-sub3c.fits',sub
      	if cen_steps eq 1 then writefits,reduction_path+obj+'/IRDIS/center/'+obj+'_right-sub3c.fits',sub
      
            Print, 'Third spot at',x3,y3

	 sub=image2
	delvar,a
	for xx=0,1023 do begin
		for yy=0,1023 do begin
			if sqrt((xx-x4)^2. + (yy-y4)^2.) gt spot_fit_hwidth then sub[xx,yy]=0.
		endfor	
	endfor
	fit1 = mpfit2dpeak(sub,a, /tilt)

      
      x4 = a[4];+float(fix(x4))-float(ext)
      y4 = a[5];+float(fix(y4))-float(ext)


	sub=fshift(sub,511.5-x4,511.5-y4)
	if cen_steps eq 0 then writefits,reduction_path+obj+'/IRDIS/center/'+obj+'_left-sub4c.fits',sub
      	if cen_steps eq 1 then writefits,reduction_path+obj+'/IRDIS/center/'+obj+'_right-sub4c.fits',sub
      
            Print, 'Fourth spot at',x4,y4
    

	cut=float(cut)
      plots,[(x1+1.-cut)*2.,(x2+1.-cut)*2.],[(y1+1.-cut)*2.,(y2+1.-cut)*2.],color=0,/device
	 plots,[(x3+1.-cut)*2.,(x4+1.-cut)*2.],[(y3+1.-cut)*2.,(y4+1.-cut)*2.],color=0,/device
	 		
	 	 lint,[x1,y1],[x2,y2],[x3,y3],[x4,y4],center

	 	 
	 	  	 	 	 if center eq !NULL then center=[x1,y1]	;if no solution 
	endif else center=center1

	Print, 'Center after centroid at:',center
	print, 'Initial center at:', center[0]+xinit,' ',center[1]+yinit
;	Print, 'Final shifting by', -(center[0]-new_cen),-(center[1]-new_cen)
	
		;hak
		;cgcleanup	;ends the GUI

if cen_steps eq 0 then first_left_center = center
if cen_steps eq 1 then first_right_center = center

if cen_steps eq 2 then last_left_center = center
if cen_steps eq 3 then last_right_center = center


if cen_steps eq 0 then first_left=image
if cen_steps eq 1 then first_right=image

if cen_steps eq 2 then last_left=image
if cen_steps eq 3 then last_right=image

if cen_steps eq 0 then xxlc=center[0]+xinit-1.	;used later in cross-correlation window setup 
if cen_steps eq 0 then yylc=center[1]+yinit+2.
if cen_steps eq 1 then xxrc=center[0]+xinit-1.
if cen_steps eq 1 then yyrc=center[1]+yinit+2.



endfor
if not keyword_set(autocen) then cgcleanup


; now shift star center frames and write them out for diagnositcs
; using 511.5 so the center is in the middle of the cube

first_left_centered=fshift(first_left,511.5-first_left_center[0] , 511.5-first_left_center[1])
												
first_right_centered=fshift(first_right,511.5-first_right_center[0] , 511.5-first_right_center[1])

date=fxpar(lhdr0,'DATE')


	yyyys=strmid(date,0,4)
	mms=strmid(date,5,2)
	dds=strmid(date,8,2)
print, date

	
	if float(yyyys) lt 2017 then begin
		if float(mms) le 7 then fix_time_sync=1 else fix_time_sync=0
		if float(yyyys) lt 2016 then fix_time_sync=1
	endif else fix_time_sync=0

	if field_stab_cen then begin 
		
 	
 	parangstart = esopar(lhdr0,'ESO TEL PARANG START ')
 	parangend = esopar(lhdr0,'ESO TEL PARANG END ')
		
 		parangle=(parangstart+parangend)/2.
	altj=esopar(lhdr0, 'ESO TEL ALT')
 	
 	
 	
 	DROT2BEGIN=esopar(lhdr0,'ESO INS4 DROT2 BEGIN')
 	DROT2END=esopar(lhdr0,'ESO INS4 DROT2 END')
	
 		drangle=(DROT2BEGIN+DROT2END)/2.
 		

 		if fix_time_sync then print, 'Altitude: ', altj
 		if fix_time_sync then print, 'DROT2 : ',drangle

 		if fix_time_sync then 	deroterror=parangle+atan( tan(  (!PI/180.) * (altj - parangle- (2. $
 			   											*drangle ) ))  ) /(!PI/180.)
 			   											
 			if fix_time_sync then if abs(deroterror) gt 100 and deroterror lt 0 then deroterror=deroterror+180.  							
 			 if fix_time_sync then if abs(deroterror) gt 100 and deroterror gt 0 then deroterror=deroterror-180.  							
 	
 			
 	if fix_time_sync eq 1 then print,'----------[ Calcualted derotation error: ',deroterror, ' ]----------------'

 	if fix_time_sync eq 1 then true_north=-1.75 - 1.1105*deroterror else true_north=-1.75
 	
 	 	print,'----------[ True North : ', true_north, ' ]----------------'

 	offset=135.87 
	offset=offset+true_north
 		
 		cen_angle=parangle+offset

		first_left_centered=rot(first_left_centered,-cen_angle,/interp)
		first_right_centered=rot(first_right_centered,-cen_angle,/interp)
	endif
												
writefits,reduction_path+obj+'/IRDIS/center/'+obj+'_first_left_median_centered.fits',first_left_centered
writefits,reduction_path+obj+'/IRDIS/center/'+obj+'_first_right_median_centered.fits',first_right_centered

												
	if rightcenfilecount gt 1 and leftcenfilecount gt 1 then begin
			last_left_centered=fshift(last_left,511.5-last_left_center[0] $
														, 511.5-last_left_center[1])
														
			last_right_centered=fshift(last_right,511.5-last_right_center[0] $
														, 511.5-last_right_center[1])			
		
		
			writefits,reduction_path+obj+'/IRDIS/center/'+obj+'_last_left_median_centered.fits',last_left_centered
			writefits,reduction_path+obj+'/IRDIS/center/'+obj+'_last_right_median_centered.fits',last_right_centered				
	endif			
	
	
	
	;now stacking science frames, correcting for anamorphism, and generating angle arrays
	k=0	;counting frames
	
	
	
times=[]	

for i=0,leftfilecount-1, 1 do begin
	
	;setting up for left and right stacks (after first frame, frames go into left2 and 
	;right2, to be stacked into left and right in a later step).
	

	if i eq 0 then begin 
		left=readfits(leftfiles[i],lefthdr)
		right=readfits(rightfiles[i],righthdr)
		
		nextleft=readfits(leftfiles[i+1],nextlefthdr)
		nextright=readfits(rightfiles[i+1],nextrighthdr)
		print, 'Size of left array',size(left)
	endif
	if i ne 0 then begin 
		left2=readfits(leftfiles[i],lefthdr)
		right2=readfits(rightfiles[i],righthdr)
		
		if i ne leftfilecount-1 then begin
		nextleft=readfits(leftfiles[i+1],nextlefthdr)
		nextright=readfits(rightfiles[i+1],nextrighthdr)
		endif
		
		lastleft=readfits(leftfiles[i-1],lastlefthdr)
		lastright=readfits(rightfiles[i-1],lastrighthdr)
	endif

	;print, leftfiles[ii]
	
	
 	
 	derot = strcompress(String(esopar(lefthdr,'ESO INS4 COMB ROT')),/remove_all)
 	 	print, 'Derotator mode:',derot

 	if derot eq 'FIELD' then begin
 	
 	;FIELD orientation images are NOT true north calibrated. They are usually used for 
 	;calibrating true North, so this is actually desired.
 	
 	parangstart=0.
 	parangend=0.
 	offset=0.
 	
 	endif else begin
 	
 	
 	
 	parangstart = esopar(lefthdr,'ESO TEL PARANG START ')
 	parangend = esopar(lefthdr,'ESO TEL PARANG END ')
 	
 	ALTSTART=esopar(lefthdr, 'ESO TEL ALT')
 	
 	if fix_time_sync and i lt leftfilecount-1 then ALTEND=esopar(nextlefthdr,'ESO TEL ALT')
 	
 	if fix_time_sync and i eq leftfilecount - 1 then begin
 		ALTLAST=esopar(lastlefthdr,'ESO TEL ALT')
 		ALTSTEP=ALTSTART-ALTLAST
 		ALTEND=ALTSTART+ALTSTEP
 	endif
 	
 	DROT2BEGIN=esopar(lefthdr,'ESO INS4 DROT2 BEGIN')
 	DROT2END=esopar(lefthdr,'ESO INS4 DROT2 END')
 	
 	endelse

	;what size are the images?
 	if i eq 0 then szl2=size(left) else szl2=size(left2)
	nstack=szl2[3]
	ww=szl2[2]/2.
	
	;corrects for distortion by 1.006 in y direction
	 print, ' Correcting for anamorphic distortion...'
	 if i eq 0 then cube=left else cube=left2 
	 print, size(cube)

	 scube = size(cube)
	 	if n_elements(size(cube)) gt 5 then sn = size(congrid(cube[*,*,0], scube[1],scube[2]*1.006)) else sn = size(congrid(cube, scube[1],scube[2]*1.006)) 
	 nz=scube[3]
	 print, nz
	if n_elements(size(cube)) gt 5 then ncube = fltarr(sn[1], sn[2], nz) else ncube = fltarr(sn[1], sn[2])
	if n_elements(size(cube)) gt 5 then for iij=0, nz-1 do ncube[*,*,iij] = congrid(cube[*,*,iij], scube[1],  scube[2]*1.006) else ncube = congrid(cube, scube[1],  scube[2]*1.006)
 	left2 =ncube
 	if i eq 0 then left=ncube
 	ncube =0.
 
 	 if i eq 0 then cube=right else cube=right2
 	scube = size(cube)
	if n_elements(size(cube)) gt 5 then  sn = size(congrid(cube[*,*,0], scube[1],scube[2]*1.006)) else sn = size(congrid(cube, scube[1],scube[2]*1.006))
	if n_elements(size(cube)) gt 5 then  ncube = fltarr(sn[1], sn[2], nz) else  ncube = fltarr(sn[1], sn[2]) 
	 	if n_elements(size(cube)) gt 5 then for ii=0, nz-1 do ncube[*,*,ii] = congrid(cube[*,*,ii], scube[1],  scube[2]*1.006)  else ncube = congrid(cube, scube[1],  scube[2]*1.006)
	 right2 =ncube
	 if i eq 0 then right=ncube
	 ncube =0.
	 
	 left=left[0:1023,0:1023,*]
	 right=right[0:1023,0:1023,*]
	 
	 left2=left2[0:1023,0:1023,*]
	 right2=right2[0:1023,0:1023,*]
	
 	;this loop takes the original input cubes, de-rotates them to N-up E-left, to be 
 	;stacked into the single cube in the next loop
 	parangle=fltarr((size(left2))(3))
 	 drangle=fltarr((size(left2))(3))
 	 altj=fltarr((size(left2))(3))
	if n_elements(size(left2)) gt 5 then begin
 	for j = 0, (size(left2))(3)-1 do begin

		
		;store times
		time=fxpar(lefthdr,'UTC') 

		times=[times,time]
	
 		parangle=parangstart+(parangend-parangstart)*float(j)/float(nstack)
	
 		drangle=DROT2BEGIN+(DROT2END-DROT2BEGIN)*float(j)/float(nstack)
 		
 		if fix_time_sync then altj=ALTSTART+(ALTEND-ALTSTART)*float(j)/float(nstack)
 		;^^^ interpolation between start of cube and end of cube

 			print, 'Parang: ', parangle
 		if fix_time_sync then print, 'Altitude: ', altj
 		if fix_time_sync then print, 'DROT2 : ',drangle

 		if fix_time_sync then 	deroterror=parangle+atan( tan(  (!PI/180.) * (altj - parangle- (2. $
 			   											*drangle ) ))  ) /(!PI/180.)
 			   											
 			if fix_time_sync then if abs(deroterror) gt 100 and deroterror lt 0 then deroterror=deroterror+180.  							
 			 if fix_time_sync then if abs(deroterror) gt 100 and deroterror gt 0 then deroterror=deroterror-180.  							

 	
 			
 	if fix_time_sync eq 1 then print,'----------[ Calcualted derotation error: ',deroterror, ' ]----------------'

 	if fix_time_sync eq 1 then true_north=-1.75 - 1.1105*deroterror else true_north=-1.75
 	
 	 	print,'----------[ True North : ', true_north, ' ]----------------'

 	offset=135.87 
	offset=offset+true_north
 		
 		;if j eq 0 then angle=fltarr(nstack)
 		if j eq 0 then angle=parangle+offset else angle=[angle, parangle+offset]

 
 		k=k+1
 		print, 'Frame # ',k,' angle = ',angle[j],' degrees'
 		
 	endfor
	endif else begin
		
parangle=(parangstart+parangstart)/2.
	
 		drangle=(DROT2BEGIN+DROT2END)/2.
 		
 		if fix_time_sync then altj=(ALTSTART+ALTEND)/2.
 		;^^^ interpolation between start of cube and end of cube

 			print, 'Parang: ', parangle
 		if fix_time_sync then print, 'Altitude: ', altj
 		if fix_time_sync then print, 'DROT2 : ',drangle

 		if fix_time_sync then 	deroterror=parangle+atan( tan(  (!PI/180.) * (altj - parangle- (2. $
 			   											*drangle ) ))  ) /(!PI/180.)
 			   											
 			if fix_time_sync then if abs(deroterror) gt 100 and deroterror lt 0 then deroterror=deroterror+180.  							
 			 if fix_time_sync then if abs(deroterror) gt 100 and deroterror gt 0 then deroterror=deroterror-180.  						
 	
 			
 	if fix_time_sync eq 1 then print,'----------[ Calcualted derotation error: ',deroterror, ' ]----------------'

 	if fix_time_sync eq 1 then true_north=-1.75 - 1.1105*deroterror else true_north=-1.75
 	
 	 	print,'----------[ True North : ', true_north, ' ]----------------'

 	offset=135.87 
	offset=offset+true_north
 		
 		;if j eq 0 then angle=fltarr(nstack)
 		angle=parangle+offset; else; angle=[angle, parangle+offset]


 		k=k+1
 		print, 'Frame # ',k,' angle = ',angle,' degrees'
 		


	endelse
 	
 	;now we put the current stacks into the bigger cubes with the next two if statements
	if i ne 0 then begin 
		left=[[[left]],[[left2]]]
		right=[[[right]],[[right2]]]	
		angles=[angles,angle]
	endif 
	if i eq 0 then begin  
		left=left2
		right=right2
		angles=angle
	endif


endfor


	left2=left
	right2=right

	;angles are known and stored above
	
	if cen_filter then begin
		first_left_centered=first_left_centered-smooth(first_left_centered,5.)
		first_left_centered=smooth(first_left_centered,3.)
		
		first_right_centered=first_right_centered-smooth(first_right_centered,5.)
		first_right_centered=smooth(first_right_centered,3.)
		
		
	for iii=0,(size(left))(3)-1 do begin
		;print, size(first_left_centered)
		;print, size(left)
		left2[*,*,iii]=left2[*,*,iii]-smooth(left2[*,*,iii],5.)
		right2[*,*,iii]=right2[*,*,iii]-smooth(right2[*,*,iii],5.)
		
		right2[*,*,iii]=smooth(right2[*,*,iii],3.)
		left2[*,*,iii]=smooth(left2[*,*,iii],3.)
	endfor
	endif
	
	if block_rad_cen eq 1 then begin
	
	
		
		for xx=0,1023 do begin
			for yy=0,1023 do begin
				xd=xx-xxlc
				yd=yy-yylc
				
				xdc=xx-511.
				ydc=yy-511.
				if sqrt(xd*xd+yd*yd) lt rad_cen_min or sqrt(xd*xd+yd*yd) gt rad_cen_max then left2[xx,yy,*]=left2[xx,yy,*]/1000000.
				if sqrt(xdc*xdc+ydc*ydc) lt rad_cen_min or sqrt(xdc*xdc+ydc*ydc) gt rad_cen_max then first_left_centered[xx,yy]=first_left_centered[xx,yy]/1000000.

			endfor
		endfor
		
		for xx=0,1023 do begin
			for yy=0,1023 do begin
				xd=xx-xxrc
				yd=yy-yyrc
				
				
				xdc=xx-511.
				ydc=yy-511.
				if sqrt(xd*xd+yd*yd) lt rad_cen_min or sqrt(xd*xd+yd*yd) gt rad_cen_max then right2[xx,yy,*]=right2[xx,yy,*]/1000000.
				if sqrt(xdc*xdc+ydc*ydc) lt rad_cen_min or sqrt(xdc*xdc+ydc*ydc) gt rad_cen_max then first_right_centered[xx,yy]=first_right_centered[xx,yy]/1000000.

			endfor
		endfor
	
	endif
	
	

;field stabilized centering
	if field_stab_cen then begin

		;first derotate
		for ii=0, (size(left2))(3)-1 do begin
			left2[*,*,ii]=rot(left2[*,*,ii],-angles[ii],/interp)
			right2[*,*,ii]=rot(right2[*,*,ii],-angles[ii],/interp)
			left[*,*,ii]=rot(left[*,*,ii],-angles[ii],/interp)
			right[*,*,ii]=rot(right[*,*,ii],-angles[ii],/interp)
		endfor

	endif

	left2[where(left2 lt 0)]=left2[where(left2 lt 0)]/10000000.
	right2[where(right2 lt 0)]=right2[where(right2 lt 0)]/10000000.
	first_left_centered[where(first_left_centered lt 0)]=first_left_centered[where(first_left_centered lt 0)]/10000000.
	first_right_centered[where(first_right_centered lt 0)]=first_right_centered[where(first_right_centered lt 0)]/10000000.
	
	; now center the science frames
	
	
	
		diff=512.

	
		
	;ldxlast=0. & ldylast=0. & rdxlast=0. & rdylast=0.	;set the wild shift tests to zero at start
	wildflagl=0 & wildflagr=0
	for iii=0,(size(left))(3)-1 do begin
		;print, size(first_left_centered)
		;print, size(left)
		lcorr=crosscorr(first_left_centered,left2[*,*,iii],pmax, /cut)
		if pmax(0) ge 0. then dx = diff-pmax(0) else dx = -diff+abs(pmax(0))
		if pmax(1) ge 0. then dy = diff-pmax(1) else dy = -diff+abs(pmax(1))
			
		if iii gt 0 and remove_wild and wildflagl eq 0 then begin
			if abs(dx-ldxlast) gt 5 or abs(dy-ldylast) gt 5 then begin dx=ldxlast & dy=ldylast & print, '!!! Wild shift detected! Using previous frame center instead!' &  wildflagl=1 & endif else wildflagl=0
		endif else wildflagl=0
		print, 'Shifting left x,y by: ',String(-dx),String(-dy)
		left[*,*,iii]=fshift(left[*,*,iii],-dx,-dy)
		left2[*,*,iii]=fshift(left2[*,*,iii],-dx,-dy)		
		
		ldxlast=dx
		ldylast=dy
	
		rcorr=crosscorr(first_right_centered,right2[*,*,iii],pmax, /cut)
		if pmax(0) ge 0. then dx = diff-pmax(0) else dx = -diff+abs(pmax(0))
		if pmax(1) ge 0. then dy = diff-pmax(1) else dy = -diff+abs(pmax(1))
		
		if iii gt 0 and remove_wild and wildflagr eq 0 and not use_left_for_right then begin
			if abs(dx-rdxlast) gt 5 or abs(dy-rdylast) gt 5 then begin dx=rdxlast & dy=rdylast & print, '!!! Wild shift detected! Using previous frame center instead!' &  wildflagr=1 & endif else wildflagr=0
		endif  else wildflagr=0

		

		if use_left_for_right then begin

			if iii eq 0 then ldx=ldxlast
			if iii eq 0 then ldy=ldylast

			if iii eq 0 then rdx=dx
			if iii eq 0 then rdy=dy

			if iii eq 0 then rldiffx=rdx-ldx
			if iii eq 0 then rldiffy=rdy-ldy

			dx=rldiffx+ldxlast
			dy=rldiffy+ldylast

		endif

		print, 'Shifting right x,y by: ',String(-dx),String(-dy)

		right[*,*,iii]=fshift(right[*,*,iii],-dx,-dy)
		right2[*,*,iii]=fshift(right2[*,*,iii],-dx,-dy)
		
		rdxlast=dx
		rdylast=dy
		
	endfor

	

	if field_stab_cen then begin

		;now go back to pupil stabilization
		for ii=0, (size(left))(3)-1 do begin
			left[*,*,ii]=rot(left[*,*,ii],angles[ii],/interp)
			right[*,*,ii]=rot(right[*,*,ii],angles[ii],/interp)
		endfor

	endif
	
	; now scale the science frames
	
		
	dit=fxpar(lefthdr,'EXPTIME') ;reading DIT and normalizing to a DIT of 1 second
	print, 'Science DIT: ',dit
	for ij=0, (size(left))(3) -1 do left[*,*,ij]=left[*,*,ij]/float(dit)
	for ij=0, (size(right))(3) -1 do right[*,*,ij]=right[*,*,ij]/float(dit)

	;ND code
	;reading the neutral-density filter ID for use in determining the next scale factor 
	filt = esopar(lefthdr,'ESO INS4 COMB IND')
	if isa(filt, 'string') eq 0 then filt = esopar(lefthdr,'ESO INS4 FILT2 NAME')  

	print, 'Science ND filter: ',filt
	;hak
	readcol, reduction_path+'calib/SPHERE_CPI_ND.dat.txt', wavend, nd0, nd1, nd2, nd35
	;readcol, reduction_path+'calib/SPHERE_IRDIS_D_K12.dat.txt', wavek12, dk1, dk2

	
	
	if strpos(filt, 'ND_N_3.5') ne -1 or STRPOS(filt, 'ND_3.5') ne -1 then nd=nd35
	if strpos(filt, 'ND_N_2.0') ne -1 or STRPOS(filt, 'ND_2.0') ne -1 then nd=nd2
	if strpos(filt, 'ND_N_1.0') ne -1 or STRPOS(filt, 'ND_1.0') ne -1 then nd=nd1
	if strpos(filt, 'ND_N_0.0') ne -1 or STRPOS(filt, 'ND_0.0') ne -1 or STRPOS(filt, 'OPEN') ne -1 then nd=nd0

	;Print,'Science ND filter: ',filt
	nearest_wavel=nearest_element(wavel*1000.,wavend,wave_posl)
	nearest_waver=nearest_element(waver*1000.,wavend,wave_posr)
	
	filtscale1=nd[wave_posl]
	filtscale2=nd[wave_posr]

	
	

	print, 'Filtscale 1:',filtscale1 
	print, 'Filtscale 2:',filtscale2




	for ij=0, (size(left))(3) -1 do left[*,*,ij]=left[*,*,ij]/float(filtscale1)
	for ij=0, (size(right))(3) -1 do right[*,*,ij]=right[*,*,ij]/float(filtscale2)

	
	; now output the science frames									

	writefits,proc_path+obj+'_left_cube.fits',left,lefthdr
	writefits,proc_path+obj+'_right_cube.fits',right,righthdr

	save, FILENAME=proc_path+obj+'_angles.sav',angles ;this stores the parallactic angle of the images, needed to rotate the images back to original field orientation when performing ADI
	save, FILENAME=proc_path+obj+'_times.sav',times ;this stores the parallactic angle of the images, needed to

endif ; do_center if 

cgcleanup

;------------------------------[ Begin Thought Process  ]---------------------------------


if do_rotational_center eq 1 then begin
print, 'Finding precise cneter...'


left_cube=readfits(proc_path+obj+'_left_cube.fits',hdr)
restore, filename=proc_path+obj+'_angles.sav'
left_angles=angles
right_angles=angles
;left_cube2=left_cube

right_cube=readfits(proc_path+obj+'_right_cube.fits',hdr)

;right_cube2=right_cube

bestresiduals=99e99
for kk=0, censteps-1 do begin
for jj=0, censteps-1 do begin

angles=left_angles

print, 'On position ', (kk)*censteps + jj+1 , ' out of ', censteps*censteps

left_cube2=left_cube[512.-cenframewidth/2.:512.+cenframewidth/2.-1.,512.-cenframewidth/2.:512.+cenframewidth/2.-1.,*]


writefits,reduction_path+'/subleft.fits',left_cube2

;shift cube frames
shiftx=float(kk)*censquare/float(censteps-1) - censquare/2. 
shifty=float(jj)*censquare/float(censteps-1) - censquare/2. 

print, 'Testing shift x', shiftx, ' shift y', shifty

for ll=0, (size(left_cube))(3)-1 do left_cube2[*,*,ll]=fshift(left_cube2[*,*,ll],shiftx,shifty)


;take median
medarr, left_cube2, medframe

;derotate images

for ii=0, (size(left_cube2))(3)-1 do begin
	;print, 'Rotating by ', angles[ii]
	;left_cube2[*,*,ii]=left_cube2[*,*,ii]
	left_cube2[*,*,ii]=rot(left_cube2[*,*,ii],-angles[ii],/interp)-medframe
endfor

left_cube2[where(finite(left_cube2) eq 0)]=0.

;block center
for xx=0.,cenframewidth-1. do begin
	for yy=0.,cenframewidth-1. do begin
		xc=xx-cenframewidth/2.
		yc=yy-cenframewidth/2.
		
		if sqrt(xc*xc + yc*yc) lt cenblockradius then left_cube2[xx,yy,*]=0.
	endfor
endfor

residuals=total(abs(left_cube2))
print, residuals

;if kk eq 0 and jj eq 0 then bestresiduals=residuals
if residuals le bestresiduals then begin
	
	print, 'Found better residual!'

	bestresiduals=residuals
	
	bestshiftx=shiftx
	bestshifty=shifty
endif



endfor	;kk censteps for
endfor	;jj censteps for

print, 'best x shift: ', bestshiftx
print, 'best y shift: ', bestshifty

for ll=0, (size(left_cube))(3)-1 do left_cube[*,*,ll]=fshift(left_cube[*,*,ll],bestshiftx,bestshifty)


		writefits,proc_path+obj+'_left_cube_cen.fits', left_cube, hdr

		medarr, left_cube, medframe
		writefits,proc_path+obj+'_left_pupil.fits', medframe
		
		
angles=right_angles

		
bestresiduals=99e99
for kk=0, censteps-1 do begin
for jj=0, censteps-1 do begin

print, 'On position ', (kk)*censteps + jj+1 , ' out of ', censteps*censteps


right_cube2=right_cube[512.-cenframewidth/2.:512.+cenframewidth/2.-1.,512.-cenframewidth/2.:512.+cenframewidth/2.-1.,*]



;shift cube frames
shiftx=float(kk)*censquare/float(censteps-1) - censquare/2. 
shifty=float(jj)*censquare/float(censteps-1) - censquare/2. 

print, 'Testing shift x', shiftx, ' shift y', shifty

for ll=0, (size(right_cube))(3)-1 do right_cube2[*,*,ll]=fshift(right_cube2[*,*,ll],shiftx,shifty)


;take median
medarr, right_cube2, medframe

;derotate images

for ii=0, (size(right_cube2))(3)-1 do begin
	;print, 'Rotating by ', angles[ii]
	;right_cube2[*,*,ii]=right_cube2[*,*,ii]

	right_cube2[*,*,ii]=rot(right_cube2[*,*,ii],-angles[ii],/interp)-medframe
endfor

right_cube2[where(finite(right_cube2) eq 0)]=0.


;block center
for xx=0.,cenframewidth-1. do begin
	for yy=0.,cenframewidth-1. do begin
		xc=xx-cenframewidth/2.
		yc=yy-cenframewidth/2.
		
		if sqrt(xc*xc + yc*yc) lt cenblockradius then right_cube2[xx,yy,*]=0.
	endfor
endfor


residuals=total(abs(right_cube2))
print, residuals


if residuals le bestresiduals then begin
	
	print, 'Found better residual!'

	bestresiduals=residuals
	
	bestshiftx=shiftx
	bestshifty=shifty
	;medarr,obj_cube2, bestresframe
	;writefits,reduction_path+'/test-adi-residuals.fits', bestresframe
endif



endfor	;kk censteps for
endfor	;jj censteps for

print, 'best x shift: ', bestshiftx
print, 'best y shift: ', bestshifty

for ll=0, (size(right_cube))(3)-1 do right_cube[*,*,ll]=fshift(right_cube[*,*,ll],bestshiftx,bestshifty)


		writefits,proc_path+obj+'_right_cube_cen.fits', right_cube, hdr

		medarr, right_cube, medframe
		writefits,proc_path+obj+'_right_pupil.fits', medframe
		
endif ;precise recenter if

;------------------------------[ Begin Thought Process  ]---------------------------------



if do_clean then begin



	print, 'Loading cubes...'

	files=file_search(proc_path+obj+'_left_cube.fits',count=cnt)
if cnt ge 1 then begin

	 lcube=readfits(proc_path+obj+'_left_cube.fits',lhd) 
	 rcube=readfits(proc_path+obj+'_right_cube.fits',rhd) 

	restore, filename=proc_path+obj+'_angles.sav'
	;restore, filename=proc_path+obj+'_times.sav'

;block outer regions
	clean_hsize=50.

	lcube2=lcube
	rcube2=rcube

	lcube2=lcube2[512.-clean_hsize:512.+clean_hsize-1.,512.-clean_hsize:512.+clean_hsize-1.,*]
	rcube2=rcube2[512.-clean_hsize:512.+clean_hsize-1.,512.-clean_hsize:512.+clean_hsize-1.,*]

;derotate if set to field stabilized
	if field_stab_clean then begin

		for ii=0,(size(lcube))(3)-1 do lcube2[*,*,ii]=rot(lcube2[*,*,ii],-angles[ii],/interp)
		for ii=0,(size(rcube))(3)-1 do rcube2[*,*,ii]=rot(rcube2[*,*,ii],-angles[ii],/interp)
	endif

;calculate median

	medarr,lcube2,lmed
	medarr,rcube2,rmed

;loop through images cross-correlating with the median, only working on one side since pupils are linked

	for ii=0,(size(lcube))(3)-1 do begin

		maxcor=max(crosscorr(lcube2[*,*,ii],lmed))
			print, 'Frame ',ii,' has max correlation value of ', maxcor, ' with pupil median.'
		if ii eq 0 then corrs=maxcor else corrs=[corrs,maxcor]
	endfor


	

;keep those above cross_thresh percentile (defined out of unity), e.g. 0.8=80%

	lcube=lcube[*,*,where(corrs ge cross_thresh)]
	rcube=rcube[*,*,where(corrs ge cross_thresh)]

	angles=angles[where(corrs ge cross_thresh)]
	;times=times[where(corrs ge cross_thresh)]

	print, 'Cleaned ', n_elements(where(corrs lt cross_thresh)) ,' frames.'
	print, 'Bad frames:', where(corrs lt cross_thresh)

	if  n_elements(where(corrs lt cross_thresh)) ge 0.5*n_elements(corrs) then begin
		print, '!!! Warning! Autoclean has rejected over 50% of frames !!! 
		print, 'You will porbably want to go back to the reduce_irdis.pro routine and change the cross_thresh parameter to a lower value (0 will accept all frames).'
		print, 'If this is a binary star, this may have happened because of comparing single frames (with a bright star) to a median-combined pupil containing a smeared image of the star. Try enabling pup_stab_clean in this case.'
		print, 'It is recommended to apply the /skip keyword next run to reduce tiem spet on initial steps.'
		print, 'Hit any key to continue anyways, or ESC to stop the reduction.'
		hak
	endif

	writefits,proc_path+obj+'_left_cube_clean.fits',lcube,lhd
	writefits,proc_path+obj+'_right_cube_clean.fits',rcube,rhd
	save, filename=proc_path+obj+'_angles_clean.sav',angles
	;save, filename=proc_path+obj+'_times_clean.sav',times


	file_delete, proc_path+obj+'_left_cube.fits'
	file_delete, proc_path+obj+'_right_cube.fits'
	file_delete, proc_path+obj+'_angles.sav'
	print, 'Done cleaning.'

endif 

endif

if skip_clean then begin 
	;restore, filename=proc_path+obj+'_times.sav'
	restore, filename=proc_path+obj+'_angles.sav'
	 lcube=readfits(proc_path+obj+'_left_cube.fits',lhd) 
	 rcube=readfits(proc_path+obj+'_right_cube.fits',rhd) 

	writefits,proc_path+obj+'_left_cube_clean.fits',lcube,lhd
	writefits,proc_path+obj+'_right_cube_clean.fits',rcube,rhd
	save, filename=proc_path+obj+'_angles_clean.sav',angles
	;save, filename=proc_path+obj+'_times_clean.sav',times

	file_delete, proc_path+obj+'_left_cube.fits'
	file_delete, proc_path+obj+'_right_cube.fits'
	file_delete, proc_path+obj+'_angles.sav'
endif;clean if 





;-------------------------[ Begin Reduction Task ]------------------------------
if do_inject eq 1 then begin


ref_r=readfits(flux_path+obj+'_right_median.fits',rhdrpsf)
ref_l=readfits(flux_path+obj+'_left_median.fits',lhdrpsf)

;first subtract sky
mmm,ref_r,skyr
mmm,ref_l,skyl

ref_l=ref_l-skyl
ref_r=ref_r-skyr


if use_gauss eq 1 then begin

;Result = GAUSS2DFIT( Z, A [, X, Y] [, FITA=vector] [, MASK=array] [, /NEGATIVE] [ ] )


gauss_ref_r=gauss2dfit(ref_r, A)

writefits, flux_path+obj+'_original_refr.fits', ref_r



ref_r=gauss_ref_r

;A=[0.,1,1,1,1,1,1]

gauss_ref_l=gauss2dfit(ref_l, A)

 writefits,flux_path+obj+'_original_refl.fits', ref_l



ref_l=gauss_ref_l

;remove sky again

mmm, ref_l,skyl
mmm, ref_r,skyr

ref_l=ref_l-skyl
ref_r=ref_r-skyr

writefits, flux_path+obj+'_gauss_refr.fits', ref_r
writefits, flux_path+obj+'_gauss_refl.fits', ref_l



endif

if skip_rotational_center then begin
	left_cube=readfits(proc_path+obj+'_left_cube_clean.fits',lefthdr)
	right_cube=readfits(proc_path+obj+'_right_cube_clean.fits',righthdr)
endif else begin
	left_cube=readfits(proc_path+obj+'_left_cube_cen_clean.fits',lefthdr)
	right_cube=readfits(proc_path+obj+'_right_cube_cen_clean.fits',righthdr)
endelse

restore, filename=proc_path+obj+'_angles_clean.sav'
print, proc_path+obj+'_angles_clean.sav'
print, size(angles)
;hak

if not skip_right then llmax=1 else llmax=0
for ll=0,llmax do begin

if ll eq 0 then obj_cube=left_cube
if ll eq 1 then obj_cube=right_cube

if ll eq 0 then ref=ref_l else ref=ref_r


refsz=15
for xx=0,(size(ref))(1)-1 do begin
	for yy=0,(size(ref))(2)-1 do begin
		dx=xx-(size(ref))(1)/2.+1.
		dy=yy-(size(ref))(2)/2.+1.
		if sqrt(dx*dx+dy*dy) gt refsz then ref[xx,yy]=0.
	endfor
endfor

;airmasses=left_airmasses

;detrotate
print, (size(obj_cube))(3)-1
for ii=0, (size(obj_cube))(3)-1 do begin
	print, 'Derotating ',ii,' by ', -angles[ii]
	obj_cube[*,*,ii]=rot(obj_cube[*,*,ii],-angles[ii],/interp)
endfor

big_ref=obj_cube[*,*,0]
big_ref[*]=0

;build bigger reference image
for xx=0, (size(ref))(1) - 1 do begin
	for yy=0, (size(ref))(1) - 1 do begin

		big_ref[512.-(size(ref))(1)/2.+xx,512.-(size(ref))(1)/2.+yy]=ref[xx,yy]
		

	endfor
endfor
		;big_ref=fshift(big_ref,+0.5,+0.5)

if debug then writefits,reduction_path+'/bigref.fits',big_ref

if ll eq 0 then planet_theta=(planet_theta)*!DTOR

if ll eq 0 then planet_contrast=planet_contrast_left
if ll eq 1 then planet_contrast=planet_contrast_right


;inject planets
for ii=0, (size(planet_r))(1)-1 do begin
	Print, 'Injecting planet', ii
	big_ref_ii=big_ref*planet_contrast[ii]
	
		xshift= planet_r[ii] * (1./pxscale) * Cos(planet_theta[ii])
		yshift= planet_r[ii] * (1./pxscale) * Sin(planet_theta[ii])
		big_ref_ii=fshift(big_ref_ii, xshift,yshift)
	
	for jj=0,(size(obj_cube))(3) - 1 do begin
	
		obj_cube[*,*,jj]=obj_cube[*,*,jj]+big_ref_ii
	
	endfor
	
endfor

;rotate back to pupil stabilized orientation
for ii=0, (size(obj_cube))(3)-1 do begin
	print, 'Rotating by ', angles[ii]
	;print, 'Airmass:', airmasses[ii]
	;print, 'Airmass scaling:',( 10.^( (0.114/2.5) * (1.195-airmasses[ii] ) ))

	obj_cube[*,*,ii]=rot(obj_cube[*,*,ii],angles[ii],/interp)
endfor


	if ll eq 0 then writefits,proc_path+obj+'_left_cube_inject.fits', obj_cube,lefthdr
	if ll eq 1 then writefits,proc_path+obj+'_right_cube_inject.fits', obj_cube,righthdr

	
endfor	
endif ;inject if



;------------------------------[ Begin Thought Process  ]---------------------------------


if do_rotate eq 1 then begin



ww=final_image_half_size

if skip_rotational_center then begin
	left_cube=readfits(proc_path+obj+'_left_cube_clean.fits',lefthdr)
	right_cube=readfits(proc_path+obj+'_right_cube_clean.fits',righthdr)
endif else begin
	left_cube=readfits(proc_path+obj+'_left_cube_cen.fits',lefthdr)
	right_cube=readfits(proc_path+obj+'_right_cube_cen.fits',righthdr)
endelse



if classical_rdi  then begin
	ref_left_cube=readfits(ref_path+classical_ref_obj+'_left_cube_clean.fits',lefthdr)
	ref_right_cube=readfits(ref_path+classical_ref_obj+'_right_cube_clean.fits',righthdr)
endif  else begin
	ref_left_cube=left_cube
	ref_right_cube=right_cube
endelse

restore, filename=proc_path+obj+'_angles_clean.sav'




;for ii=0, (size(left_cube))(3)-1 do left_cube[*,*,ii]=left_cube[*,*,ii]-smooth(left_cube[*,*,ii],10.)
if use_inject eq 1 then left_cube=readfits(proc_path+obj+'_left_cube_inject.fits',lefthdr)
if use_inject eq 1 then ref_left_cube=left_cube





goods=angles ;just check how many frames there are

if auto_bin then begin
	if n_elements(goods) gt 60 then bin=2
	if n_elements(goods) gt 100 then bin=3
	if n_elements(goods) gt 120 then bin=4
	if n_elements(goods) gt 150 then bin=5
endif


pre_binned_angles=angles

if bin gt 1 then begin
	st=1
	for ii=0.,(size(left_cube))(3)-1. do begin
		print, fix(ii+1.) mod fix(bin)
		if fix(ii+1.) mod fix(bin) eq 1 then binned=left_cube[*,*,ii] else binned=[ [[binned]],[[ left_cube[*,*,ii] ]] ]
		if fix(ii+1.) mod fix(bin) eq 0 then begin
		
			print,'Binning left frames...'
			
			if st eq 1 then begin 
			
			medarr,binned,binned
			binned_cube = binned 
			st=0 
			
			 endif else begin
			 medarr,binned,binned

			 binned_cube=[[[binned_cube]],[[binned]]]
			 
			 endelse
			print, size(binned_cube)
		endif
	endfor	
	print, size(binned_cube)
	st=1
	for ii=0,(size(left_cube))(3)-1. do begin
		if fix(ii+1.) mod fix(bin) eq 1 then binned_angle=angles[ii] else binned_angle=binned_angle+angles[ii]
		if fix(ii+1.) mod fix(bin) eq 0 then begin
			if st eq 1 then begin binned_angles = binned_angle/bin & st=0 & endif else binned_angles=[[[binned_angles]],[[binned_angle/bin]]]
		endif
	endfor	
	
	left_cube=binned_cube
	angles=binned_angles
endif


;lmask=readfits(reduction_path+'/calib/irdis_left_mask.fits')
;rmask=readfits(reduction_path+'/calib/irdis_right_mask.fits')

;for ii=0,(size(left_cube))(3)-1 do left_cube[*,*,ii]=left_cube[*,*,ii]*lmask
;for ii=0,(size(right_cube))(3)-1 do right_cube[*,*,ii]=right_cube[*,*,ii]*rmask


;block bad regions
;left_cube[0:50,*,*]=0.
;left_cube[1023-50:1023,*,*]=0.
;left_cube[*,1023-50:1023,*]=0.
;left_cube[*,0:50,*]=0.

dbf=esopar(lefthdr,'ESO INS1 OPTI2 NAME ')
	if strpos(dbf,'Y23') ne -1 then dual_filter='Y23'

	if strpos(dbf,'J23') ne -1 then dual_filter='J23'
	if strpos(dbf,'H23') ne -1 then dual_filter='H23'
	if strpos(dbf,'H34') ne -1 then dual_filter='H34'
	if strpos(dbf,'K12') ne -1 then dual_filter='K12'

	if strpos(dbf,'Y23') ne -1 then dual_filter='Y23'
	
	if strpos(dbf,'J23') ne -1 then dual_filter='J23'
	if strpos(dbf,'H23') ne -1 then dual_filter='H23'
	if strpos(dbf,'H34') ne -1 then dual_filter='H34'
	if strpos(dbf,'K12') ne -1 then dual_filter='K12'


	if strpos(dbf,'CLEAR') ne -1 then begin;dual_filter='K12'
		print, 'Dual filter clear. Checking for single filter.'
		dbf2=esopar(lhdr,'ESO INS1 FILT NAME ')
		print, dbf,dbf2
		if strpos(dbf2,'Ks') ne -1 then dual_filter='Ks'
		if strpos(dbf2,'B_Y') ne -1 then dual_filter='Y'
		if strpos(dbf2,'B_J') ne -1 then dual_filter='J'
		if strpos(dbf2,'B_H') ne -1 then dual_filter='H'
	
	endif

	
	print, 'IRDIS Filter = ',dual_filter
	if dual_filter eq 'Ks' then wavel=2.182
	if dual_filter eq 'Ks' then waver=2.182


	if dual_filter eq 'Y' then wavel=1.043
	if dual_filter eq 'Y' then waver=1.043
	
	
	if dual_filter eq 'J' then wavel=1.245
	if dual_filter eq 'J' then waver=1.245
	
	
	if dual_filter eq 'H' then wavel=1.625
	if dual_filter eq 'H' then waver=1.625

	
	if dual_filter eq 'Ks' then wavel=2.182
	if dual_filter eq 'Ks' then waver=2.182

	if dual_filter eq 'Y23' then wavel=1.022

	if dual_filter eq 'K12' then wavel=2.1025
	if dual_filter eq 'J23' then wavel=1.1895
	if dual_filter eq 'H23' then wavel=1.5888

	if dual_filter eq 'Y23' then waver=1.076

	if dual_filter eq 'K12' then waver=2.255
	if dual_filter eq 'J23' then waver=1.2698
	if dual_filter eq 'H23' then waver=1.6671

if destripe_frames ne 0 then begin
	for ii=0, (size(ref_left_cube))(3)-1 do ref_left_cube[*,*,ii]=destripe(ref_left_cube[*,*,ii],90.,clip_level=0.0,/nodisp)

	print, 'destriping left frames 90 degrees...'
	for ii=0, (size(left_cube))(3)-1 do left_cube[*,*,ii]=destripe(left_cube[*,*,ii],90.,clip_level=0.0,/nodisp)
	;print, 'destriping 0 degrees...'
	;for ii=0, (size(left_cube))(3)-1 do left_cube[*,*,ii]=destripe(left_cube[*,*,ii],0.,clip_level=0.0,/nodisp)

endif

restore, filename=proc_path+obj+'_angles_clean.sav'

if high_pass_width gt 0 then for ii= 0, (size(left_cube))(3)-1 do $
		left_cube[*,*,ii]=left_cube[*,*,ii]-smooth(left_cube[*,*,ii],high_pass_width)
if smooth_width gt 0 then for ii= 0, (size(left_cube))(3)-1 do $
		left_cube[*,*,ii]=smooth(left_cube[*,*,ii],smooth_width)
left_cube[where(left_cube eq 0)]=!values.f_nan				
left_cube_pupil=left_cube
left_adi_cube=left_cube


if high_pass_width gt 0 then for ii= 0, (size(ref_left_cube))(3)-1 do $
		ref_left_cube[*,*,ii]=ref_left_cube[*,*,ii]-smooth(ref_left_cube[*,*,ii],high_pass_width)
if smooth_width gt 0 then for ii= 0, (size(ref_left_cube))(3)-1 do $
		ref_left_cube[*,*,ii]=smooth(ref_left_cube[*,*,ii],smooth_width)
ref_left_cube[where(ref_left_cube eq 0)]=!values.f_nan				
;ref_left_cube_pupil=ref_left_cube


;for ii=0, (size(right_cube))(3)-1 do right_cube[*,*,ii]=right_cube[*,*,ii]-smooth(right_cube[*,*,ii],10.)
if use_inject eq 1 then right_cube=readfits(proc_path+obj+'_right_cube_inject.fits',righthdr)
if use_inject eq 1 then ref_right_cube=right_cube

angles=pre_binned_angles

if bin gt 1 then begin
	st=1
	for ii=0.,(size(right_cube))(3)-1. do begin
		print, fix(ii+1.) mod fix(bin)
		if fix(ii+1.) mod fix(bin) eq 1 then binned=right_cube[*,*,ii] else binned=[ [[binned]],[[ right_cube[*,*,ii] ]] ]
		if fix(ii+1.) mod fix(bin) eq 0 then begin
		
			print,'Binning right frames...'
			if st eq 1 then begin 
			medarr,binned,binned
			binned_cube = binned 
			 st=0 
			 endif else begin
			 medarr,binned,binned

			 binned_cube=[[[binned_cube]],[[binned]]]
			 endelse
			print, size(binned_cube)
		endif
	endfor	
	print, size(binned_cube)
	st=1
	for ii=0,(size(right_cube))(3)-1. do begin
		if fix(ii+1.) mod fix(bin) eq 1 then binned_angle=angles[ii] else binned_angle=binned_angle+angles[ii]
		if fix(ii+1.) mod fix(bin) eq 0 then begin
			if st eq 1 then begin binned_angles = binned_angle/bin & st=0 & endif else binned_angles=[[[binned_angles]],[[binned_angle/bin]]]
		endif
	endfor	
	
	right_cube=binned_cube
	angles=binned_angles
endif


if destripe_frames ne 0 then begin

	print, 'destriping right frames 90 degrees...'
		for ii=0, (size(ref_right_cube))(3)-1 do ref_right_cube[*,*,ii]=destripe(ref_right_cube[*,*,ii],90.,clip_level=0.0,/nodisp)

	for ii=0, (size(right_cube))(3)-1 do right_cube[*,*,ii]=destripe(right_cube[*,*,ii],90.,clip_level=0.0,/nodisp)
	;print, 'destriping 0 degrees...'
	;for ii=0, (size(right_cube))(3)-1 do right_cube[*,*,ii]=destripe(right_cube[*,*,ii],0.,clip_level=0.0,/nodisp)

endif


bigright=fltarr(512,512,(size(right_cube))(3))

if high_pass_width gt 0 then for ii= 0, (size(right_cube))(3)-1 do $
		right_cube[*,*,ii]=right_cube[*,*,ii]-smooth(right_cube[*,*,ii],high_pass_width)
if smooth_width gt 0 then for ii= 0, (size(right_cube))(3)-1 do $
		right_cube[*,*,ii]=smooth(right_cube[*,*,ii],smooth_width)
right_cube[where(right_cube eq 0)]=!values.f_nan				
right_adi_cube=right_cube
right_cube_pupil=right_cube

if high_pass_width gt 0 then for ii= 0, (size(ref_right_cube))(3)-1 do $
		ref_right_cube[*,*,ii]=ref_right_cube[*,*,ii]-smooth(ref_right_cube[*,*,ii],high_pass_width)
if smooth_width gt 0 then for ii= 0, (size(ref_right_cube))(3)-1 do $
		ref_right_cube[*,*,ii]=smooth(ref_right_cube[*,*,ii],smooth_width)
ref_right_cube[where(ref_right_cube eq 0)]=!values.f_nan				


if smart_adi eq 1 then begin
for ii=0, (size(left_cube))(3)-1 do begin
	print, 'Performing ADI on left frame', ii
	medarr, ref_left_cube[*,*,where(abs(angles-angles[ii]) gt 5.)], medframel
	left_adi_cube[*,*,ii]=left_adi_cube[*,*,ii]-medframel

endfor

for ii=0, (size(right_cube))(3)-1 do begin
	print, 'Performing ADI on right frame', ii

	medarr, ref_right_cube[*,*,where(abs(angles-angles[ii]) gt 5.)], rightmed
	right_adi_cube[*,*,ii]=right_adi_cube[*,*,ii]-rightmed
	
endfor
endif else begin	;classical ADI
	medarr,ref_left_cube,medframel
	medarr,ref_right_cube,rightmed
	
	
	medarr,left_cube,medframel_orig
	medarr,right_cube,rightmed_orig
	

	;this is an ADI scale factor which is only non-unity when using a reference cube
	medframel=medframel*max(medframel_orig[500:519,500:519])/max(medframel[500:519,500:519])
	rightmed=rightmed*max(rightmed_orig[500:519,500:519])/max(rightmed[500:519,500:519])

	
for ii=0, (size(left_cube))(3)-1 do begin
	print, 'Performing ADI on left frame', ii
	left_adi_cube[*,*,ii]=left_adi_cube[*,*,ii]-medframel
	
endfor

for ii=0, (size(right_cube))(3)-1 do begin
	print, 'Performing ADI on right frame', ii
	right_adi_cube[*,*,ii]=right_adi_cube[*,*,ii]-rightmed
	
endfor
endelse


;block bad regions
left_cube[0:100,*,*]=0.
left_cube[1023-100:1023,*,*]=0.
left_cube[*,1023-100:1023,*]=0.
left_cube[*,0:100,*]=0.

left_adi_cube[0:100,*,*]=0.
left_adi_cube[1023-100:1023,*,*]=0.
left_adi_cube[*,1023-100:1023,*]=0.
left_adi_cube[*,0:100,*]=0.


right_cube[0:100,*,*]=0.
right_cube[1023-100:1023,*,*]=0.
right_cube[*,1023-100:1023,*]=0.
right_cube[*,0:100,*]=0.

right_adi_cube[0:100,*,*]=0.
right_adi_cube[1023-100:1023,*,*]=0.
right_adi_cube[*,1023-100:1023,*]=0.
right_adi_cube[*,0:100,*]=0.


	
for ii=0, (size(left_cube))(3)-1 do begin

	left_cube[*,*,ii]=rot(left_cube[*,*,ii],-angles[ii],/interp)
	left_adi_cube[*,*,ii]=rot(left_adi_cube[*,*,ii],-angles[ii],/interp)

endfor

for ii=0, (size(right_cube))(3)-1 do begin

	right_cube[*,*,ii]=rot(right_cube[*,*,ii],-angles[ii],/interp)
	right_adi_cube[*,*,ii]=rot(right_adi_cube[*,*,ii],-angles[ii],/interp)

endfor


left_cube[0:100,*,*]=0.
left_cube[1023-100:1023,*,*]=0.
left_cube[*,1023-100:1023,*]=0.
left_cube[*,0:100,*]=0.

left_adi_cube[0:100,*,*]=0.
left_adi_cube[1023-100:1023,*,*]=0.
left_adi_cube[*,1023-100:1023,*]=0.
left_adi_cube[*,0:100,*]=0.


right_cube[0:100,*,*]=0.
right_cube[1023-100:1023,*,*]=0.
right_cube[*,1023-100:1023,*]=0.
right_cube[*,0:100,*]=0.

right_adi_cube[0:100,*,*]=0.
right_adi_cube[1023-100:1023,*,*]=0.
right_adi_cube[*,1023-100:1023,*]=0.
right_adi_cube[*,0:100,*]=0.




if output_cubes then writefits,proc_path+obj+'_left_cube_derot'+suffix+'.fits', left_cube, hdr
if output_cubes then writefits,proc_path+obj+'_right_cube_derot'+suffix+'.fits', right_cube, hdr


writefits,proc_path+obj+'_left_pupil.fits', medframel, hdr
writefits,proc_path+obj+'_right_pupil.fits', rightmed, hdr




if comb_type eq 'mean' then leftderotmed=mean(left_cube,dim=3)
if comb_type eq 'mean' then rightderotmed=mean(right_cube,dim=3)



if comb_type eq 'median' then medarr,left_cube,leftderotmed
if comb_type eq 'median' then medarr,right_cube,rightderotmed


if comb_type eq 'nw-mean' then leftderotmed=nw_ang_comb(left_cube,angles)
if comb_type eq 'nw-mean' then rightderotmed=nw_ang_comb(right_cube,angles)


writefits,proc_path+obj+'_left_derot'+suffix+'.fits', leftderotmed, lefthdr
writefits,proc_path+obj+'_right_derot'+suffix+'.fits', rightderotmed, righthdr

writefits,proc_path+obj+'_left+right_derot'+suffix+'.fits', (leftderotmed+rightderotmed)/2., hdr


if identify then find_sources, proc_path+obj+'_left_derot'+suffix+'.fits',sigma=5.0,inrad=annmode_inout[0]+2, reference=reduction_path+obj+'/IRDIS/flux/'+obj+'_left_median.fits',FWHM=((wavel*1.0E-6)/(8.4))*206265./0.01225,platescale=0.01225,clip_level=left_derot_clip_level,stretch=left_derot_stretch,correction_factor=left_derot_correction_factor,outrad=annmode_inout[1]-2

;if identify then find_sources, proc_path+obj+'_right_derot'+suffix+'.fits',sigma=5.0,inrad=annmode_inout[0]+2, reference=reduction_path+obj+'/IRDIS/flux/'+obj+'_right_median.fits',FWHM=((waver*1.0E-6)/(8.4))*206265./0.01225,platescale=0.01225,clip_level=right_derot_clip_level,stretch=right_derot_stretch,correction_factor=right_derot_correction_factor,outrad=annmode_inout[1]-2

;if identify then find_sources, proc_path+obj+'_left+right_derot'+suffix+'.fits',sigma=5.0,inrad=annmode_inout[0]+2, reference=reduction_path+obj+'/IRDIS/flux/'+obj+'_left+right_median.fits',FWHM=((((wavel+waver)/2.)*1.0E-6)/(8.4))*206265./0.01225,platescale=0.01225,clip_level=total_derot_clip_level,stretch=total_derot_stretch,correction_factor=total_derot_correction_factor,outrad=annmode_inout[1]-2

;if keyword_set(notify) then notify=1 else notify=1
notify=0
if notify then begin


body='Raw data reduction with high pass fitler width = '+string(high_pass_width)+' --- Message automatically generated on ' + String(systime())

notify_string='echo "'+body+'" | mail -s "ZONASPHERE '+ obj + ' Data Reduction Results" kwagner@as.arizona.edu -A '+strcompress(proc_path+obj+'_left+right_derot'+suffix+'.fits_source_map.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_left_derot'+suffix+'.fits_source_map.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_right_derot'+suffix+'.fits_source_map.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_left+right_derot'+suffix+'.fits_sources.txt',/rem) + $
' -A '+strcompress(proc_path+obj+'_left_derot'+suffix+'.fits_sources.txt',/rem) + $
' -A '+strcompress(proc_path+obj+'_right_derot'+suffix+'.fits_sources.txt',/rem) + $
' -A '+strcompress(proc_path+obj+'_right_derot'+suffix+'.fits_radial_sensitivity.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_right_derot'+suffix+'.fits_radial_sensitivity_2.0_arcsec.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_left_derot'+suffix+'.fits_radial_sensitivity.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_left_derot'+suffix+'.fits_radial_sensitivity_2.0_arcsec.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_left+right_derot'+suffix+'.fits_radial_sensitivity.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_left+right_derot'+suffix+'.fits_radial_sensitivity_2.0_arcsec.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_right_derot'+suffix+'.fits_radial_sensitivity.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_right_derot'+suffix+'.fits_radial_sensitivity_2.0_arcsec.png',/rem) 

spawn, notify_string





endif


if comb_type eq 'median' then medarr,left_adi_cube,leftadimed
if comb_type eq 'median' then medarr,right_adi_cube,rightadimed

if comb_type eq 'mean' then leftadimed=mean(left_adi_cube,dim=3)
if comb_type eq 'mean' then rightadimed=mean(right_adi_cube,dim=3)


if comb_type eq 'nw-mean' then leftadimed=nw_ang_comb(left_adi_cube,angles)
if comb_type eq 'nw-mean' then rightadimed=nw_ang_comb(right_adi_cube,angles)





if output_cubes then writefits,proc_path+obj+'_left_cADI_cube'+suffix+'.fits', left_adi_cube, hdr
if output_cubes then writefits,proc_path+obj+'_right_cADI_cube'+suffix+'.fits', right_adi_cube, hdr



writefits,proc_path+obj+'_left_cADI'+suffix+'.fits', leftadimed, lefthdr
writefits,proc_path+obj+'_right_cADI'+suffix+'.fits', rightadimed, righthdr

if identify_adi then find_sources, proc_path+obj+'_left_cADI'+suffix+'.fits',sigma=5.0,inrad=annmode_inout[0]+2, reference=reduction_path+obj+'/IRDIS/flux/'+obj+'_left_median.fits',FWHM=((wavel*1.0E-6)/(8.4))*206265./0.01225,platescale=0.01225,clip_level=left_adi_clip_level,stretch=left_adi_stretch,correction_factor=left_adi_correction_factor,outrad=annmode_inout[1]-2

if identify_adi then find_sources, proc_path+obj+'_right_cADI'+suffix+'.fits',sigma=5.0,inrad=annmode_inout[0]+2, reference=reduction_path+obj+'/IRDIS/flux/'+obj+'_right_median.fits',FWHM=((waver*1.0E-6)/(8.4))*206265./0.01225,platescale=0.01225,clip_level=right_adi_clip_level,stretch=right_adi_stretch,correction_factor=right_adi_correction_factor,outrad=annmode_inout[1]-2


scale=12.25/7.46
leftadimed_ifsc=rot(leftadimed,0.,scale,/interp)
rightadimed_ifsc=rot(rightadimed,0.,scale,/interp)

if not fill then for xx=0.,1023. do for yy=0.,1023. do if sqrt( (xx-512.)^2. + (yy-512.)^2 ) gt scale*wr*nrings then leftadimed_ifsc[xx,yy]=!values.f_nan
if not fill then for xx=0.,1023. do for yy=0.,1023. do if sqrt( (xx-512.)^2. + (yy-512.)^2 ) gt scale*wr*nrings then rightadimed_ifsc[xx,yy]=!values.f_nan


writefits,strcompress(proc_path+obj+'_left+right_cADI_IFS-scale'+suffix+'.fits'), leftadimed_ifsc+rightadimed_ifsc,hdr
writefits,strcompress(proc_path+obj+'_left_cADI_IFS-scale'+suffix+'.fits'), leftadimed_ifsc,hdr
writefits,strcompress(proc_path+obj+'_right_cADI_IFS-scale'+suffix+'.fits'), rightadimed_ifsc,hdr

sdi=((leftadimed-rightadimed)+(rightadimed-leftadimed))


writefits,proc_path+obj+'_left+right_cADI'+suffix+'.fits', (leftadimed+rightadimed)/2., hdr

if identify_adi then find_sources, proc_path+obj+'_left+right_cADI'+suffix+'.fits',sigma=5.0,inrad=annmode_inout[0]+2, reference=reduction_path+obj+'/IRDIS/flux/'+obj+'_left+right_median.fits',FWHM=((((wavel+waver)/2.)*1.0E-6)/(8.4))*206265./0.01225,platescale=0.01225,clip_level=total_adi_clip_level,stretch=total_adi_stretch,correction_factor=adi_correction_factor,outrad=annmode_inout[1]-2

if notify_adi then begin


body='Classical ADI reduction with high pass filter width = '+string(high_pass_width)+' --- Message automatically generated on ' + String(systime())
notify_string='echo "'+body+'" | mail -s "ZONASPHERE '+ obj + ' Data Reduction Results" kwagner@as.arizona.edu -A '+strcompress(proc_path+obj+'_left+right_cADI'+suffix+'.fits_source_map.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_left_cADI'+suffix+'.fits_source_map.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_right_cADI'+suffix+'.fits_source_map.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_left+right_cADI'+suffix+'.fits_sources.txt',/rem) + $
' -A '+strcompress(proc_path+obj+'_left_cADI'+suffix+'.fits_sources.txt',/rem) + $
' -A '+strcompress(proc_path+obj+'_right_cADI'+suffix+'.fits_sources.txt',/rem) + $
' -A '+strcompress(proc_path+obj+'_right_cADI'+suffix+'.fits_radial_sensitivity.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_right_cADI'+suffix+'.fits_radial_sensitivity_2.0_arcsec.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_left_cADI'+suffix+'.fits_radial_sensitivity.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_left_cADI'+suffix+'.fits_radial_sensitivity_2.0_arcsec.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_left+right_cADI'+suffix+'.fits_radial_sensitivity.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_left+right_cADI'+suffix+'.fits_radial_sensitivity_2.0_arcsec.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_right_cADI'+suffix+'.fits_radial_sensitivity.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_right_cADI'+suffix+'.fits_radial_sensitivity_2.0_arcsec.png',/rem) 

spawn, notify_string



restore,proc_path+obj+'_left+right_cADI'+suffix+'.fits_sources.sav'

;select candidates based on contrast and separation

candidates = where(sources_rho le 2.0 and sources_rho ge 0.18)

;candidates=indgen(7) 

if n_elements(candidates) eq 1 then if candidates eq -1 then candidates=!null


if n_elements(sources_x) eq 1 then if sources_x eq 0. and sources_y eq 0. then candidates=!null
if candidates ne !NULL then begin



candidates_contrast= sources_contrast[candidates]
candidates_x=sources_x[candidates]
candidates_y=sources_y[candidates]
candidates_rho=sources_rho[candidates]
candidates_theta=sources_theta[candidates]
candidates_sigma=sources_sigma[candidates]


smalls=readfits(proc_path+obj+'_left+right_cADI'+suffix+'.fits_sources.fits')
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
	
	cgps_open,proc_path+obj+'_left+right_cADI'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png'


   cgLoadCT, 1;, NColors=100,BOTTOM=100;,BOTTOM=100
	cgDisplay,40,40,/device

	;cgImage,cand;,stretch=stretch,xrange=[0,(size(smalls))(2)],yrange=[0,(size(smalls))(2)],/keep_aspect_ratio,minvalue=min(cand),maxvalue=max(cand),mean=mean(cand);,position=[0,0,1,1]

cgImage,cand,xrange=[0,39],yrange=[0,39],/keep_aspect_ratio;=[0,0,1,1];,minvalue=min(cand),maxvalue=max(cand),mean=mean(cand);,position=[0,0,1,1]

	cgps_close,/png

	cgcleanup

	attach_string = attach_string+ ' -A '+strcompress(proc_path+obj+'_left+right_cADI'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png',/rem)
endfor


;for cc=0,n_elements(candidates)-1 do begin

	body='ALERT!'+string(n_elements(candidates))+' Candidate exoplanets discovered around '+obj+' in cADI reduction with significance = '+sigma_string+' contrast = '+contrast_string+', rho='+rho_string+', theta='+theta_string+' source numbers ='+candidates_string+' --- Message automatically generated by ZONASPHERE on ' + String(systime())

	notify_string='echo "'+body+'" | mail -s "ALERT!!! '+string(n_elements(candidates)) +' Exoplanet Candidates Detected Around '+ obj + '!!!" kwagner@as.arizona.edu -A '+strcompress(proc_path+obj+'_left+right_cADI'+suffix+'.fits_source_map.png',/rem) +' -A '+strcompress(proc_path+obj+'_left+right_cADI'+suffix+'.fits_sources.txt',/rem) + ' -A '+strcompress(proc_path+obj+'_left+right_cADI'+suffix+'.fits_SNR_map.png',/rem)  + attach_string

	spawn, notify_string 

for mmm=0,n_elements(candidates)-1 do file_delete,strcompress(proc_path+obj+'_left+right_cADI'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png',/rem)



;endfor
endif



endif

;classical SDI

dbf=esopar(lefthdr,'ESO INS1 OPTI2 NAME ')
if strpos(dbf,'Y23') ne -1 then dual_filter='Y23'

if strpos(dbf,'J23') ne -1 then dual_filter='J23'
if strpos(dbf,'H23') ne -1 then dual_filter='H23'
if strpos(dbf,'H34') ne -1 then dual_filter='H34'
if strpos(dbf,'K12') ne -1 then dual_filter='K12'


	if strpos(dbf,'CLEAR') ne -1 then begin;dual_filter='K12'
		print, 'Dual filter clear. Checking for single filter.'
		dbf2=esopar(lhdr,'ESO INS1 FILT NAME ')
		print, dbf,dbf2
		if strpos(dbf2,'Ks') ne -1 then dual_filter='Ks'
		if strpos(dbf2,'B_Y') ne -1 then dual_filter='Y'
		if strpos(dbf2,'B_J') ne -1 then dual_filter='J'
		if strpos(dbf2,'B_H') ne -1 then dual_filter='H'
	
	endif


	print, dbf
	if dual_filter eq 'Ks' then wavel=2.182
	if dual_filter eq 'Ks' then waver=2.182


if dual_filter eq 'Y23' then wavel=1.022

if dual_filter eq 'K12' then wavel=2.1025
if dual_filter eq 'J23' then wavel=1.1895
if dual_filter eq 'H23' then wavel=1.5888

if dual_filter eq 'Y23' then waver=1.076

if dual_filter eq 'K12' then waver=2.255
if dual_filter eq 'J23' then waver=1.2698
if dual_filter eq 'H23' then waver=1.6671



	if dual_filter eq 'H' then wavel=1.625
	if dual_filter eq 'H' then waver=1.625


sz=1024.

widthl=(wavel*1E-6) / (8.) * 206265. / 0.01225
widthr=(waver*1E-6) / (8.) * 206265. / 0.01225

print, 'Width left: ',widthl
print, 'Width right: ', widthr

PSFl = psf_Gaussian(NPIX=sz, FWHM=[widthl,widthl])
PSFr = psf_Gaussian(NPIX=sz, FWHM=[widthr,widthr])

; Normalize the PSF

PSFNl = PSFl/MAX(PSFl)
PSFNr = PSFr/MAX(PSFr)

scale12=waver/wavel
scale21=wavel/waver

left_scaled=rot(leftadimed,0.,scale12 )
right_scaled=rot(rightadimed,0.,scale21 )


leftdr_scaled=rot(leftderotmed,0.,scale12 )
rightdr_scaled=rot(rightderotmed,0.,scale21 )

left_sdi=leftadimed-right_scaled

right_sdi=rightadimed-left_scaled

left_csdi=leftderotmed-rightdr_scaled
right_csdi=rightderotmed-leftdr_scaled

writefits,  proc_path+obj+'_left_cADI_SDI'+suffix+'.fits', left_sdi, lefthdr
writefits,  proc_path+obj+'_right_cADI_SDI'+suffix+'.fits', right_sdi, righthdr


writefits,  proc_path+obj+'_left_cSDI'+suffix+'.fits', left_csdi, lefthdr
writefits,  proc_path+obj+'_right_cSDI'+suffix+'.fits', right_csdi, righthdr
writefits,  proc_path+obj+'_left+right_cSDI'+suffix+'.fits', (left_csdi+right_csdi)/2.0, righthdr


writefits,  proc_path+obj+'_left+right_cADI_SDI'+suffix+'.fits', (left_sdi+right_sdi)/2.0, righthdr

;sdi_c=convolve((left_sdi+right_sdi)/2.0, PSFr)
;writefits,  proc_path+obj+'_left+right_cADI_SDI_convolved'+suffix+'.fits', sdi_c, righthdr

leftadimed[where(finite(leftadimed) eq 0 )]=0.
rightadimed[where(finite(leftadimed) eq 0 )]=0.

;convolve images


wave=(wavel+waver)/2. ;µm


sz=final_image_half_size*2.+1.;511

width=(wave*1.0E-6) / (8.2) * 206265. / 0.01225

print, 'PSF Width: ',width

PSF = psf_Gaussian(NPIX=sz, FWHM=[width,width])

; Normalize the PSF

PSFN = PSF/MAX(PSF)

leftadimed[where(finite(leftadimed) eq 0)]=0.
rightadimed[where(finite(rightadimed) eq 0)]=0.


leftadi_c = convolve(leftadimed, PSFN)
rightadi_c = convolve(rightadimed, PSFN)

hdr=lefthdr ;headers for left and right are the same so it doesn't matter which is used

writefits,proc_path+obj+'_left_cADI_conv'+suffix+'.fits', leftadi_c, lefthdr
writefits,proc_path+obj+'_right_cADI_conv'+suffix+'.fits', rightadi_c, righthdr

total_conv_image=(leftadi_c+rightadi_c)/2.
writefits,proc_path+obj+'_left+right_cADI_conv'+suffix+'.fits', total_conv_image, hdr



left_sdi[where(finite(left_sdi) eq 0)]=0.
right_sdi[where(finite(right_sdi) eq 0)]=0.

left_csdi[where(finite(left_csdi) eq 0)]=0.
right_csdi[where(finite(right_csdi) eq 0)]=0.



lcsdi_c=convolve(left_csdi,PSFN)
rcsdi_c=convolve(right_csdi,PSFN)


writefits,  proc_path+obj+'_left_cSDI_conv'+suffix+'.fits', lcsdi_c, lefthdr
writefits,  proc_path+obj+'_right_cSDI_conv'+suffix+'.fits', rcsdi_c, righthdr

lcsdi_c=convolve(left_sdi,PSFN)
rcsdi_c=convolve(right_sdi,PSFN)



csdi_c=convolve((left_sdi+right_sdi)/2.0,PSFN)


writefits,  proc_path+obj+'_left_cADI_SDI_conv'+suffix+'.fits', lcsdi_c, lefthdr
writefits,  proc_path+obj+'_right_cADI_SDI_conv'+suffix+'.fits', rcsdi_c, righthdr

writefits,  proc_path+obj+'_left+right_cADI_SDI_conv'+suffix+'.fits', csdi_c, righthdr




endif ;rotate if



;------------------------------[ Begin Thought Process  ]---------------------------------


if do_klip eq 1 then begin

	;suffix=strcompress( suffix+ '_'+string(sigfig(k_klip,2))+'k_'+string(sigfig(angsep,2))+'as_'+string(sigfig(high_pass_width,2))+'sf_'+string(sigfig(nrings,2))+'nrings_'+string(sigfig(n_ang,1))+'nang_'+string(sigfig(bin,2))+'bin',/rem)


ww=final_image_half_size

if skip_rotational_center then begin
	left_cube=readfits(proc_path+obj+'_left_cube_clean.fits',lefthdr)
	right_cube=readfits(proc_path+obj+'_right_cube_clean.fits',righthdr)
endif else begin
	left_cube=readfits(proc_path+obj+'_left_cube_cen.fits',lefthdr)
	right_cube=readfits(proc_path+obj+'_right_cube_cen.fits',righthdr)
endelse
if use_inject eq 1 then left_cube=readfits(proc_path+obj+'_left_cube_inject.fits',lefthdr)


right_cube[where(finite(right_cube) ne 1)]=0.
left_cube[where(finite(left_cube) ne 1)]=0.
restore, filename=proc_path+obj+'_angles_clean.sav'


nimages=(size(left_cube))(3)

goods=indgen(nimages)

if remove_bads then remove,bads,goods,angles

if remove_bads then left_cube=left_cube[*,*,goods]


if auto_bin then begin
	if n_elements(goods) gt 60 then bin=2
	if n_elements(goods) gt 100 then bin=3
	if n_elements(goods) gt 120 then bin=4
	if n_elements(goods) gt 150 then bin=5
endif


if bin gt 1 then begin
	st=1
	for ii=0.,(size(left_cube))(3)-1. do begin
		print, fix(ii+1.) mod fix(bin)
		if fix(ii+1.) mod fix(bin) eq 1 then binned=left_cube[*,*,ii] else binned=[ [[binned]],[[ left_cube[*,*,ii] ]] ]
		if fix(ii+1.) mod fix(bin) eq 0 then begin
		
			print,'Binning left frames...'
			
			if st eq 1 then begin 
			
			medarr,binned,binned
			binned_cube = binned 
			st=0 
			
			 endif else begin
			 medarr,binned,binned

			 binned_cube=[[[binned_cube]],[[binned]]]
			 
			 endelse
			print, size(binned_cube)
		endif
	endfor	
	print, size(binned_cube)
	st=1
	for ii=0,(size(left_cube))(3)-1. do begin
		if fix(ii+1.) mod fix(bin) eq 1 then binned_angle=angles[ii] else binned_angle=binned_angle+angles[ii]
		if fix(ii+1.) mod fix(bin) eq 0 then begin
			if st eq 1 then begin binned_angles = binned_angle/bin & st=0 & endif else binned_angles=[[[binned_angles]],[[binned_angle/bin]]]
		endif
	endfor	
	
	left_cube=binned_cube
	angles=binned_angles
endif



if destripe_frames ne 0 then begin

	print, 'destriping 90 degrees...'
	for ii=0, (size(left_cube))(3)-1 do left_cube[*,*,ii]=destripe(left_cube[*,*,ii],90.,clip_level=0.0,/nodisp)
	;print, 'destriping 0 degrees...'
	;for ii=0, (size(left_cube))(3)-1 do left_cube[*,*,ii]=destripe(left_cube[*,*,ii],0.,clip_level=0.0,/nodisp)

endif

if high_pass_width gt 0 then for ii= 0, (size(left_cube))(3)-1 do $
		left_cube[*,*,ii]=left_cube[*,*,ii]-smooth(left_cube[*,*,ii],high_pass_width)
if smooth_width gt 0 then for ii= 0, (size(left_cube))(3)-1 do $
		left_cube[*,*,ii]=smooth(left_cube[*,*,ii],smooth_width)				
;left_cube_pupil=left_cube
left_klip_cube=left_cube


;right_cube=readfits(proc_path+obj+'_right_cube.fits',righthdr)
if use_inject eq 1 and not skip_right then right_cube=readfits(proc_path+obj+'_right_cube_inject.fits',righthdr)


dbf=esopar(lefthdr,'ESO INS1 OPTI2 NAME ')
if strpos(dbf,'Y23') ne -1 then dual_filter='Y23'

if strpos(dbf,'J23') ne -1 then dual_filter='J23'
if strpos(dbf,'H23') ne -1 then dual_filter='H23'
if strpos(dbf,'H34') ne -1 then dual_filter='H34'
if strpos(dbf,'K12') ne -1 then dual_filter='K12'

	if strpos(dbf,'CLEAR') ne -1 then begin;dual_filter='K12'
		print, 'Dual filter clear. Checking for single filter.'
		dbf2=esopar(lefthdr,'ESO INS1 FILT NAME ')
		print, dbf,dbf2
		if strpos(dbf2,'Ks') ne -1 then dual_filter='Ks'
		if strpos(dbf2,'B_Y') ne -1 then dual_filter='Y'
		if strpos(dbf2,'B_J') ne -1 then dual_filter='J'
		if strpos(dbf2,'B_H') ne -1 then dual_filter='H'
	
	endif

	print, dbf
	if dual_filter eq 'Ks' then wavel=2.182
	if dual_filter eq 'Ks' then waver=2.182


	if dual_filter eq 'H' then wavel=1.625
	if dual_filter eq 'H' then waver=1.625

if dual_filter eq 'Y23' then wavel=1.022

if dual_filter eq 'K12' then wavel=2.1025
if dual_filter eq 'J23' then wavel=1.1895
if dual_filter eq 'H23' then wavel=1.5888

if dual_filter eq 'Y23' then waver=1.076

if dual_filter eq 'K12' then waver=2.255
if dual_filter eq 'J23' then waver=1.2698
if dual_filter eq 'H23' then waver=1.6671


restore, filename=proc_path+obj+'_angles_clean.sav'


nimages=(size(right_cube))(3)

goods=indgen(nimages)


if remove_bads then remove,bads,goods,angles

if remove_bads and not skip_right then right_cube=right_cube[*,*,goods]

if bin gt 1 then begin
	st=1
	for ii=0.,(size(right_cube))(3)-1. do begin
		print, fix(ii+1.) mod fix(bin)
		if fix(ii+1.) mod fix(bin) eq 1 then binned=right_cube[*,*,ii] else binned=[ [[binned]],[[ right_cube[*,*,ii] ]] ]
		if fix(ii+1.) mod fix(bin) eq 0 then begin
		
			print,'Binning right frames...'
			if st eq 1 then begin 
			medarr,binned,binned
			binned_cube = binned 
			 st=0 
			 endif else begin
			 medarr,binned,binned

			 binned_cube=[[[binned_cube]],[[binned]]]
			 endelse
			print, size(binned_cube)
		endif
	endfor	
	print, size(binned_cube)

	right_cube=binned_cube
	angles=binned_angles
endif

save, filename=proc_path+obj+'_adiklip_angles.sav',angles



if destripe_frames ne 0 and not skip_right then begin

	print, 'destriping 90 degrees...'
	for ii=0, (size(right_cube))(3)-1 do right_cube[*,*,ii]=destripe(right_cube[*,*,ii],90.,clip_level=0.0,/nodisp)
	;print, 'destriping 0 degrees...'
	;for ii=0, (size(right_cube))(3)-1 do right_cube[*,*,ii]=destripe(right_cube[*,*,ii],0.,clip_level=0.0,/nodisp)

endif

if not skip_right then begin
if high_pass_width gt 0 then for ii= 0, (size(right_cube))(3)-1 do $
		right_cube[*,*,ii]=right_cube[*,*,ii]-smooth(right_cube[*,*,ii],high_pass_width)
if smooth_width gt 0 then for ii= 0, (size(right_cube))(3)-1 do $
		right_cube[*,*,ii]=smooth(right_cube[*,*,ii],smooth_width)
endif
right_klip_cube=right_cube

;right_cube_pupil=right_cube

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



posang=angles
for ii=0, (size(left_cube))(3)-1 do begin

Print, '-------[Object: ',obj,' -- KLIPing left image ', ii, ' out of ', (size(left_cube))(3)-1, ' ]----------'

if not hyp and not annmode then left_klip_cube[*,*,ii] = adiklip(left_cube, k_klip, target=ii, anglemask=anglemask, distmask=distmask, posang=posang, wl=2.2,diam=8.2, pixelscale=0.01225, angsep=angsep, anglemax=anglemax,obj=obj,nrings=nrings, wr =wr, n_ang =n_ang) 

if hyp and not annmode then left_klip_cube[*,*,ii] = adiklip(left_cube, k_klip, target=ii, anglemask=anglemask, distmask=distmask, posang=posang, wl=2.2,diam=8.2, pixelscale=0.01225, angsep=angsep, anglemax=anglemax,obj=obj,nrings=nrings, wr =wr, n_ang =n_ang, spot_radius=spot_radius,rho=rho,phi=phi, /hyper) 

if annmode then left_klip_cube[*,*,ii] = adiklip(left_cube, k_klip, target=ii, anglemask=anglemask, distmask=distmask, posang=posang, wl=2.2,diam=8.2, pixelscale=0.01225, angsep=angsep, anglemax=anglemax,obj=obj,nrings=nrings, wr =wr, n_ang =n_ang, annmode_inout=annmode_inout) 

endfor

left_klip_cube[where(left_klip_cube eq 0)]=!values.f_nan
right_klip_cube[where(right_klip_cube eq 0)]=!values.f_nan

	;no output_cubes step because this is used in sdi.
	writefits,proc_path+obj+'_left_cube_klip.fits', left_klip_cube,lefthdr



;block bad regions
left_klip_cube[0:100,*,*]=0.
left_klip_cube[1023-100:1023,*,*]=0.
left_klip_cube[*,1023-100:1023,*]=0.
left_klip_cube[*,0:100,*]=0.




for ii=0, (size(left_cube))(3)-1 do begin
	print, 'Rotating by ', angles[ii]
	left_klip_cube[*,*,ii]=rot(left_klip_cube[*,*,ii],-angles[ii],/interp)

endfor


;block bad regions
left_klip_cube[0:100,*,*]=0.
left_klip_cube[1023-100:1023,*,*]=0.
left_klip_cube[*,1023-100:1023,*]=0.
left_klip_cube[*,0:100,*]=0.

if chop eq 1 then left_klip_cube=left_klip_cube>0


;if not do_rotate then fill=0


	
	if comb_type eq 'median' then medarr, left_klip_cube, medframel
	if comb_type eq 'mean' then medframel=mean(left_klip_cube,dim=3)
	if comb_type eq 'nw-mean' then medframel=nw_ang_comb(left_klip_cube,angles)

	IF fill and not do_rotate then begin
		leftderotmed = readfits(proc_path+obj+'_left_derot.fits')
		leftadimed = readfits(proc_path+obj+'_left_cADI.fits')
	endif
	
	if fill eq 1 and fill_adi eq 0  then medframel[where(finite(medframel) eq 0 )]=leftderotmed[where(finite(medframel) eq 0 )]
	if fill eq 1  and fill_adi eq 1 then	medframel[where(finite(medframel) eq 0 )]=leftadimed[where(finite(medframel) eq 0 )]


if not fill then for xx=0.,1023. do for yy=0.,1023. do if sqrt( (xx-512.)^2. + (yy-512.)^2 ) gt wr*nrings then medframel[xx,yy]=!values.f_nan
	
	writefits,proc_path+obj+'_left_klip'+suffix+'.fits', medframel,lefthdr


	
	;medframel[where(finite(medframel) eq 0 )]=0.;rightadimed[where(finite(medframel) eq 0 )]

	;identical code for right




if not skip_right then begin	

for ii=0, (size(right_cube))(3)-1 do begin

Print, '-------[Object: ',obj,' -- KLIPing right image ', ii, ' out of ', (size(right_cube))(3)-1, ' ]----------'

if not hyp and not annmode then right_klip_cube[*,*,ii] = adiklip(right_cube, k_klip, target=ii, anglemask=anglemask, distmask=distmask, posang=posang, wl=2.2,diam=8.2, pixelscale=0.01225, angsep=angsep, anglemax=anglemax,obj=obj,nrings=nrings, wr =wr, n_ang =n_ang) 

if hyp and not annmode then right_klip_cube[*,*,ii] = adiklip(right_cube, k_klip, target=ii, anglemask=anglemask, distmask=distmask, posang=posang, wl=2.2,diam=8.2, pixelscale=0.01225, angsep=angsep, anglemax=anglemax,obj=obj,nrings=nrings, wr =wr, n_ang =n_ang, spot_radius=spot_radius,rho=rho,phi=phi, /hyper) 

if annmode then right_klip_cube[*,*,ii] = adiklip(right_cube, k_klip, target=ii, anglemask=anglemask, distmask=distmask, posang=posang, wl=2.2,diam=8.2, pixelscale=0.01225, angsep=angsep, anglemax=anglemax,obj=obj,nrings=nrings, wr =wr, n_ang =n_ang, annmode_inout=annmode_inout) 


endfor
	;no output_cubes switch since this is used in sdi step
		 writefits,proc_path+obj+'_right_cube_klip.fits', right_klip_cube,righthdr


;block bad regions
right_klip_cube[0:100,*,*]=0.
right_klip_cube[1023-100:1023,*,*]=0.
right_klip_cube[*,1023-100:1023,*]=0.
right_klip_cube[*,0:100,*]=0.



for ii=0, (size(right_cube))(3)-1 do begin
	print, 'Rotating by ', angles[ii]
	right_klip_cube[*,*,ii]=rot(right_klip_cube[*,*,ii],-angles[ii],/interp)

endfor


;block bad regions
right_klip_cube[0:100,*,*]=0.
right_klip_cube[1023-100:1023,*,*]=0.
right_klip_cube[*,1023-100:1023,*]=0.
right_klip_cube[*,0:100,*]=0.


if chop eq 1 then right_klip_cube=right_klip_cube>0


	



	if comb_type eq 'median' then medarr, right_klip_cube, medframer
	if comb_type eq 'mean' then medframer=mean(right_klip_cube,dim=3)
	if comb_type eq 'nw-mean' then medframer=nw_ang_comb(right_klip_cube,angles)

IF fill and not do_rotate then begin
		rightderotmed = readfits(proc_path+obj+'_right_derot.fits')
		rightadimed = readfits(proc_path+obj+'_right_cADI.fits')
	endif
	
		if fill eq 1  and fill_adi eq 0 then	medframer[where(finite(medframer) eq 0 )]=rightderotmed[where(finite(medframer) eq 0 )]
		if fill eq 1 and fill_adi eq 1 then	medframer[where(finite(medframer) eq 0 )]=rightadimed[where(finite(medframer) eq 0 )]



if not fill then for xx=0.,1023. do for yy=0.,1023. do if sqrt( (xx-512.)^2. + (yy-512.)^2 ) gt wr*nrings then medframer[xx,yy]=!values.f_nan
	
	
	writefits,proc_path+obj+'_right_klip'+suffix+'.fits', medframer,righthdr
hdr=righthdr




if identify then find_sources, proc_path+obj+'_left_klip'+suffix+'.fits',sigma=5.0,inrad=annmode_inout[0]+2, reference=reduction_path+obj+'/IRDIS/flux/'+obj+'_left_median.fits',FWHM=((wavel*1.0E-6)/(8.4))*206265./0.01225,platescale=0.01225,clip_level=left_klip_clip_level,stretch=left_klip_stretch,correction_factor=left_klip_correction_factor,outrad=annmode_inout[1]-2


if identify then find_sources, proc_path+obj+'_right_klip'+suffix+'.fits',sigma=5.0,inrad=annmode_inout[0]+2, reference=reduction_path+obj+'/IRDIS/flux/'+obj+'_right_median.fits',FWHM=((waver*1.0E-6)/(8.4))*206265./0.01225,platescale=0.01225,clip_level=right_klip_clip_level,stretch=right_klip_derot_stretch,correction_factor=right_klip_correction_factor,outrad=annmode_inout[1]-2


	dbf=esopar(righthdr,'ESO INS1 OPTI2 NAME ')
	
	print, dbf
	
if strpos(dbf,'Y23') ne -1 then dual_filter='Y23'

if strpos(dbf,'J23') ne -1 then dual_filter='J23'
if strpos(dbf,'H23') ne -1 then dual_filter='H23'
if strpos(dbf,'H34') ne -1 then dual_filter='H34'
if strpos(dbf,'K12') ne -1 then dual_filter='K12'


	if strpos(dbf,'CLEAR') ne -1 then begin;dual_filter='K12'
		print, 'Dual filter clear. Checking for single filter.'
		dbf2=esopar(lefthdr,'ESO INS1 FILT NAME ')
		if strpos(dbf2,'Ks') ne -1 then dual_filter='Ks'
		if strpos(dbf2,'B_Y') ne -1 then dual_filter='Y'
		if strpos(dbf2,'B_J') ne -1 then dual_filter='J'
		if strpos(dbf2,'B_H') ne -1 then dual_filter='H'
		
		

	endif
	
	
	print, dual_filter
	if dual_filter eq 'Ks' then wavel=2.182
	if dual_filter eq 'Ks' then waver=2.182


	if dual_filter eq 'Y' then wavel=1.043
	if dual_filter eq 'Y' then waver=1.043
	
	
	if dual_filter eq 'J' then wavel=1.245
	if dual_filter eq 'J' then waver=1.245
	
	
	if dual_filter eq 'H' then wavel=1.625
	if dual_filter eq 'H' then waver=1.625


if dual_filter eq 'Y23' then wavel=1.022

if dual_filter eq 'K12' then wavel=2.1025
if dual_filter eq 'J23' then wavel=1.1895
if dual_filter eq 'H23' then wavel=1.5888

if dual_filter eq 'Y23' then waver=1.076

if dual_filter eq 'K12' then waver=2.255
if dual_filter eq 'J23' then waver=1.2698
if dual_filter eq 'H23' then waver=1.6671
	
	medframer_sc=rot(medframer,0.,wavel/waver,/interp)
	medframel_sc=rot(medframel,0.,waver/wavel,/interp)
	
	;medframel_sc[where(finite(medframel_sc) eq 0)]=0.
	;medframer_sc[where(finite(medframer_sc) eq 0)]=0.

	
		writefits,proc_path+obj+'_right_klip_cSDI'+suffix+'.fits', medframer-medframel_sc,hdr
		writefits,proc_path+obj+'_left_klip_cSDI'+suffix+'.fits', medframel-medframer_sc,hdr
		
		writefits,proc_path+obj+'_left+right_klip_cSDI'+suffix+'.fits', ((medframel-medframer_sc)+(medframer-medframel_sc))/2.


	
	
	
	;medframer[where(finite(medframer) eq 0 )]=0.;rightadimed[where(finite(medframer) eq 0 )]

	
	total_image=(medframel+medframer)/2.

	writefits,strcompress(proc_path+obj+'_left+right_klip'+suffix+'.fits'), total_image,hdr


if identify then find_sources, proc_path+obj+'_left+right_klip'+suffix+'.fits',sigma=5.0,inrad=annmode_inout[0]+2, reference=reduction_path+obj+'/IRDIS/flux/'+obj+'_left+right_median.fits',FWHM=((((wavel+waver)/2.)*1.0E-6)/(8.4))*206265./0.01225,platescale=0.01225,clip_level=total_klip_clip_level,stretch=total_klip_stretch,correction_factor=total_klip_correction_factor,outrad=annmode_inout[1]-2

endif ; skip right if	

if notify_adiklip then begin


body='KLIP ADI reduction with high pass filter width = '+string(high_pass_width)+' --- Message automatically generated on ' + String(systime())
notify_string='echo "'+body+'" | mail -s "ZONASPHERE '+ obj + ' Data Reduction Results" kwagner@as.arizona.edu -A '+strcompress(proc_path+obj+'_left_klip'+suffix+'.fits_source_map.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_left+right_klip'+suffix+'.fits_source_map.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_right_klip'+suffix+'.fits_source_map.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_left_klip'+suffix+'.fits_sources.txt',/rem) + $
' -A '+strcompress(proc_path+obj+'_left+right_klip'+suffix+'.fits_sources.txt',/rem) + $
' -A '+strcompress(proc_path+obj+'_right_klip'+suffix+'.fits_sources.txt',/rem) + $
' -A '+strcompress(proc_path+obj+'_right_klip'+suffix+'.fits_radial_sensitivity.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_right_klip'+suffix+'.fits_radial_sensitivity_2.0_arcsec.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_left_klip'+suffix+'.fits_radial_sensitivity.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_left_klip'+suffix+'.fits_radial_sensitivity_2.0_arcsec.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_left+right_klip'+suffix+'.fits_radial_sensitivity.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_left+right_klip'+suffix+'.fits_radial_sensitivity_2.0_arcsec.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_right_klip'+suffix+'.fits_radial_sensitivity.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_right_klip'+suffix+'.fits_radial_sensitivity_2.0_arcsec.png',/rem) 

spawn, notify_string



restore,proc_path+obj+'_left_klip'+suffix+'.fits_sources.sav'

;select candidates based on contrast and separation

candidates = where(sources_rho le 2.0 and sources_rho ge 0.18)

;candidates=indgen(7) 

if n_elements(candidates) eq 1 then if candidates eq -1 then candidates=!null

print, sources_rho
print, sources_x
print, sources_y
print, candidates

if n_elements(sources_x) eq 1 then if sources_x eq 0. and sources_y eq 0. then candidates=!null
if candidates ne !NULL then begin



candidates_contrast= sources_contrast[candidates]
candidates_x=sources_x[candidates]
candidates_y=sources_y[candidates]
candidates_rho=sources_rho[candidates]
candidates_theta=sources_theta[candidates]
candidates_sigma=sources_sigma[candidates]


smalls=readfits(proc_path+obj+'_left_klip'+suffix+'.fits_sources.fits')
smalls=smalls[*,*,candidates]

rho_string='' & sigma_string='' & contrast_string='' & attach_string = '' & theta_string='' & candidates_string=''

cgcleanup

for mmm=0,n_elements(candidates)-1 do begin

plotcolor='Snow'
	cand=smalls[*,*,mmm]
	rho_string=rho_string+string(candidates_rho[mmm])
	theta_string=theta_string+string(candidates_theta[mmm])
	sigma_string=sigma_string+string(candidates_sigma[mmm])
	contrast_string=contrast_string+string(candidates_contrast[mmm])
	candidates_string=candidates_string+string(candidates[mmm]+1)
	
	cgps_open,proc_path+obj+'_left_klip'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png'


   cgLoadCT, 1;, NColors=100,BOTTOM=100;,BOTTOM=100
	cgDisplay,40,40,/device

	;cgImage,cand;,stretch=stretch;,xrange=[0,(size(smalls))(2)],yrange=[0,(size(smalls))(2)],/keep_aspect_ratio,minvalue=min(cand),maxvalue=max(cand),mean=mean(cand);,position=[0,0,1,1]

cgImage,cand,xrange=[0,39],yrange=[0,39],/keep_aspect_ratio


	cgps_close,/png

	cgcleanup

	attach_string = attach_string+ ' -A '+strcompress(proc_path+obj+'_left_klip'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png',/rem)
endfor


;for cc=0,n_elements(candidates)-1 do begin

	body='ALERT!'+string(n_elements(candidates))+' Candidate exoplanets discovered around '+obj+' in KLIP ADI reduction with significance = '+sigma_string+' contrast = '+contrast_string+', rho='+rho_string+', theta='+theta_string+' source numbers ='+candidates_string+' --- Message automatically generated by ZONASPHERE on ' + String(systime())

	notify_string='echo "'+body+'" | mail -s "ALERT!!! '+string(n_elements(candidates)) +' Exoplanet Candidates Detected Around '+ obj + '!!!" kwagner@as.arizona.edu -A '+strcompress(proc_path+obj+'_left_klip'+suffix+'.fits_source_map.png',/rem) +' -A '+strcompress(proc_path+obj+'_left_klip'+suffix+'.fits_sources.txt',/rem) + ' -A '+strcompress(proc_path+obj+'_left_klip'+suffix+'.fits_SNR_map.png',/rem)  + attach_string

	spawn, notify_string 

for mmm=0,n_elements(candidates)-1 do file_delete,strcompress(proc_path+obj+'_left_klip'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png',/rem)



;endfor
endif



endif



;scale to IFS platescale
scale=12.25/7.46
if not skip_right then total_image_ifsc=rot(total_image,0.,scale,/interp)
medframel_ifsc=rot(medframel,0.,scale,/interp)
if not skip_right then medframer_ifsc=rot(medframer,0.,scale,/interp)


if not fill then for xx=0.,1023. do for yy=0.,1023. do if sqrt( (xx-512.)^2. + (yy-512.)^2 ) gt scale*wr*nrings then medframel_ifsc[xx,yy]=!values.f_nan
if not skip_right then if not fill then for xx=0.,1023. do for yy=0.,1023. do if sqrt( (xx-512.)^2. + (yy-512.)^2 ) gt scale*wr*nrings then medframer_ifsc[xx,yy]=!values.f_nan
if not skip_right then if not fill then for xx=0.,1023. do for yy=0.,1023. do if sqrt( (xx-512.)^2. + (yy-512.)^2 ) gt scale*wr*nrings then total_image_ifsc[xx,yy]=!values.f_nan

if not skip_right then writefits,strcompress(proc_path+obj+'_left+right_klip_IFS-scale'+suffix+'.fits'), total_image_ifsc,hdr
writefits,strcompress(proc_path+obj+'_left_klip_IFS-scale'+suffix+'.fits'), medframel_ifsc,hdr
if not skip_right then writefits,strcompress(proc_path+obj+'_right_klip_IFS-scale'+suffix+'.fits'), medframer_ifsc,hdr


;convolve images


wave=(wavel+waver)/2 ;µm


sz=final_image_half_size*2.+1.;511

width=(wave*1E-6) / (8.2) * 206265. / 0.01225

print, 'PSF Width: ',width

PSF = psf_Gaussian(NPIX=sz, FWHM=[width,width])

; Normalize the PSF

PSFN = PSF/MAX(PSF)


medframel_c = convolve(medframel, PSFN)
if not skip_right then rightmed_c = convolve(medframer, PSFN)


	if not fill then for xx=0.,1023. do for yy=0.,1023. do if sqrt( (xx-512.)^2. + (yy-512.)^2 ) gt wr*nrings then medframel_c[xx,yy]=!values.f_nan

if not skip_right then if not fill then for xx=0.,1023. do for yy=0.,1023. do if sqrt( (xx-512.)^2. + (yy-512.)^2 ) gt wr*nrings then rightmed_c[xx,yy]=!values.f_nan

writefits,proc_path+obj+'_left_klip_conv'+suffix+'.fits', medframel_c, hdr
if not skip_right then writefits,proc_path+obj+'_right_klip_conv'+suffix+'.fits', rightmed_c, hdr

if not skip_right then total_conv_image=(medframel_c+rightmed_c)/2.
if not skip_right then writefits,proc_path+obj+'_left+right_klip_conv'+suffix+'.fits', total_conv_image, hdr



if not skip_right then 		writefits,proc_path+obj+'_left+right_klip_cSDI_conv'+suffix+'.fits', convolve(((medframel-medframer_sc)+(medframer-medframel_sc))/2.,PSFN),hdr




endif ;klip if

;--------------------------------------------------------------------------------------------

if do_klip_sdi then begin


ww=final_image_half_size



if skip_rotational_center then begin
	left_cube=readfits(proc_path+obj+'_left_cube_clean.fits',lefthdr)
	right_cube=readfits(proc_path+obj+'_right_cube_clean.fits',righthdr)
endif else begin
	left_cube=readfits(proc_path+obj+'_left_cube_cen.fits',lefthdr)
	right_cube=readfits(proc_path+obj+'_right_cube_cen.fits',righthdr)
endelse


if use_inject then left_cube=readfits(proc_path+obj+'_left_cube_inject.fits',lefthdr)
if use_adi then left_cube=readfits(proc_path+obj+'_left_cube_klip.fits',lefthdr)
if use_adi then restore, filename=proc_path+obj+'_adiklip_angles.sav' else restore, filename=proc_path+obj+'_angles_clean.sav'


if use_inject then right_cube=readfits(proc_path+obj+'_right_cube_inject.fits',righthdr)
if use_adi then right_cube=readfits(proc_path+obj+'_right_cube_klip.fits',righthdr)

if use_adi then restore, filename=proc_path+obj+'_adiklip_angles.sav' else restore, filename=proc_path+obj+'_angles_clean.sav'


angles=reform(angles)

;print, size(left_cube)
;print, size(angles)
;print, 'Check cube sizes!!!'
;hak


dbf=esopar(lefthdr,'ESO INS1 OPTI2 NAME ')
if strpos(dbf,'Y23') ne -1 then dual_filter='Y23'

if strpos(dbf,'J23') ne -1 then dual_filter='J23'
if strpos(dbf,'H23') ne -1 then dual_filter='H23'
if strpos(dbf,'H34') ne -1 then dual_filter='H34'
if strpos(dbf,'K12') ne -1 then dual_filter='K12'



if dual_filter eq 'Y23' then wavel=1.022

if dual_filter eq 'K12' then wavel=2.1025
if dual_filter eq 'J23' then wavel=1.1895
if dual_filter eq 'H23' then wavel=1.5888

if dual_filter eq 'Y23' then waver=1.076

if dual_filter eq 'K12' then waver=2.255
if dual_filter eq 'J23' then waver=1.2698
if dual_filter eq 'H23' then waver=1.6671


if not use_adi and auto_bin then begin
	goods=(size(left_cube))(3)
	if n_elements(goods) gt 100 then bin_sdi=2
	if n_elements(goods) gt 150 then bin_sdi=3
	if n_elements(goods) gt 200 then bin_sdi=4
endif

if bin_sdi gt 1 then begin
	st=1
	for ii=0.,(size(left_cube))(3)-1. do begin
		print, fix(ii+1.) mod fix(bin_sdi)
		if fix(ii+1.) mod fix(bin_sdi) eq 1 then binned=left_cube[*,*,ii] else binned=[ [[binned]],[[ left_cube[*,*,ii] ]] ]
		if fix(ii+1.) mod fix(bin_sdi) eq 0 then begin
		
			print,'Binning left frames...'
			
			if st eq 1 then begin 
			
			medarr,binned,binned
			binned_cube = binned 
			st=0 
			
			 endif else begin
			 medarr,binned,binned

			 binned_cube=[[[binned_cube]],[[binned]]]
			 
			 endelse
			print, size(binned_cube)
		endif
	endfor	
	print, size(binned_cube)
	st=1
	for ii=0,(size(left_cube))(3)-1. do begin
		if fix(ii+1.) mod fix(bin_sdi) eq 1 then binned_angle=angles[ii] else binned_angle=binned_angle+angles[ii]
		if fix(ii+1.) mod fix(bin_sdi) eq 0 then begin
			if st eq 1 then begin binned_angles = binned_angle/bin_sdi & st=0 & endif else binned_angles=[[[binned_angles]],[[binned_angle/bin_sdi]]]
		endif
	endfor	
	
	left_cube=binned_cube
	angles=binned_angles
endif



if destripe_frames ne 0 then begin

	print, 'destriping 90 degrees...'
	for ii=0, (size(left_cube))(3)-1 do left_cube[*,*,ii]=destripe(left_cube[*,*,ii],90.,clip_level=0.0,/nodisp)
	;print, 'destriping 0 degrees...'
	;for ii=0, (size(left_cube))(3)-1 do left_cube[*,*,ii]=destripe(left_cube[*,*,ii],0.,clip_level=0.0,/nodisp)

endif

if not use_adi then begin

if high_pass_width gt 0 then for ii= 0, (size(left_cube))(3)-1 do $
		left_cube[*,*,ii]=left_cube[*,*,ii]-smooth(left_cube[*,*,ii],high_pass_width)
if smooth_width gt 0 then for ii= 0, (size(left_cube))(3)-1 do $
		left_cube[*,*,ii]=smooth(left_cube[*,*,ii],smooth_width)	
		
endif			
;left_cube_pupil=left_cube
left_klip_cube=left_cube




if bin_sdi gt 1 then begin
	st=1
	for ii=0.,(size(right_cube))(3)-1. do begin
		print, fix(ii+1.) mod fix(bin_sdi)
		if fix(ii+1.) mod fix(bin_sdi) eq 1 then binned=right_cube[*,*,ii] else binned=[ [[binned]],[[ right_cube[*,*,ii] ]] ]
		if fix(ii+1.) mod fix(bin_sdi) eq 0 then begin
		
			print,'Binning right frames...'
			if st eq 1 then begin 
			medarr,binned,binned
			binned_cube = binned 
			 st=0 
			 endif else begin
			 medarr,binned,binned

			 binned_cube=[[[binned_cube]],[[binned]]]
			 endelse
			print, size(binned_cube)
		endif
	endfor	
	print, size(binned_cube)
	st=1
	for ii=0,(size(right_cube))(3)-1. do begin
		if fix(ii+1.) mod fix(bin_sdi) eq 1 then binned_angle=angles[ii] else binned_angle=binned_angle+angles[ii]
		if fix(ii+1.) mod fix(bin_sdi) eq 0 then begin
			if st eq 1 then begin binned_angles = binned_angle/bin_sdi & st=0 & endif else binned_angles=[[[binned_angles]],[[binned_angle/bin_sdi]]]
		endif
	endfor	
	
	right_cube=binned_cube
	angles=binned_angles
endif


if destripe_frames ne 0 then begin

	print, 'destriping 90 degrees...'
	for ii=0, (size(right_cube))(3)-1 do right_cube[*,*,ii]=destripe(right_cube[*,*,ii],90.,clip_level=0.0,/nodisp)
	;print, 'destriping 0 degrees...'
	;for ii=0, (size(right_cube))(3)-1 do right_cube[*,*,ii]=destripe(right_cube[*,*,ii],0.,clip_level=0.0,/nodisp)

endif

if not use_adi then begin

if high_pass_width gt 0 then for ii= 0, (size(right_cube))(3)-1 do $
		right_cube[*,*,ii]=right_cube[*,*,ii]-smooth(right_cube[*,*,ii],high_pass_width)
if smooth_width gt 0 then for ii= 0, (size(right_cube))(3)-1 do $
		right_cube[*,*,ii]=smooth(right_cube[*,*,ii],smooth_width)
		
endif
right_klip_cube=right_cube
;right_cube_pupil=right_cube

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;function RDIKLIP, cube, refcube, k_klip, target=target, checkpoint = checkpoint, anglemask=anglemask, $
;distmask=distmask, posang=posang,scaling=scaling, diam=diam,pixelscale=pixelscale, wl=wl,angsep=angsep,hyper=hyper,$
;obj=obj,nrings=nrings,wr =wr,n_ang =n_ang

right_cube_sc=right_cube
left_cube_sc=left_cube

right_cube_sc[where(finite(right_cube_sc) eq 0)]=0.
left_cube_sc[where(finite(left_cube_sc) eq 0)]=0.


right_cube[where(finite(right_cube) eq 0)]=0.
left_cube[where(finite(left_cube) eq 0)]=0.

print,size(left_cube)
print, size(right_cube_sc)

for ii=0, (size(right_cube))(3)-1 do begin
	right_cube_sc[*,*,ii]=rot(right_cube[*,*,ii],0.,wavel/waver)
endfor

for ii=0, (size(left_cube))(3)-1 do begin
	left_cube_sc[*,*,ii]=rot(left_cube[*,*,ii],0.,waver/wavel)
endfor


posang=angles
for ii=0, (size(left_cube))(3)-1 do begin

Print, '-------[ SDI KLIPing left image ', ii, ' out of ', (size(left_cube))(3)-1, ' ]----------'

if not hyp and not annmode then left_klip_cube[*,*,ii] = rdiklip(left_cube, right_cube_sc,k_klip_sdi, target=ii, anglemask=anglemask, distmask=distmask, posang=posang, wl=wavel,diam=8.2, pixelscale=0.01225, angsep=angsep, obj=obj,nrings=nrings, wr =wr, n_ang =n_ang) 

if hyp and not annmode then left_klip_cube[*,*,ii] = rdiklip(left_cube, right_cube_sc, k_klip_sdi, target=ii, anglemask=anglemask, distmask=distmask, posang=posang, wl=wavel,diam=8.2, pixelscale=0.01225, angsep=angsep, obj=obj,nrings=nrings, wr =wr, n_ang =n_ang, spot_radius=spot_radius,rho=rho,phi=phi, /hyper) 

if annmode then left_klip_cube[*,*,ii] = rdiklip(left_cube, right_cube_sc,k_klip_sdi, target=ii, anglemask=anglemask, distmask=distmask, posang=posang, wl=wavel,diam=8.2, pixelscale=0.01225, angsep=angsep, obj=obj,nrings=nrings, wr =wr, n_ang =n_ang, annmode_inout=annmode_inout) 

endfor

left_klip_cube[where(left_klip_cube eq 0)]=!values.f_nan
right_klip_cube[where(right_klip_cube eq 0)]=!values.f_nan


for ii=0, (size(left_cube))(3)-1 do begin
	print, 'Rotating by ', angles[ii]
	left_klip_cube[*,*,ii]=rot(left_klip_cube[*,*,ii],-angles[ii],/interp)

endfor

if chop eq 1 then left_klip_cube=left_klip_cube>0


	if output_cubes then writefits,proc_path+obj+'_left_cube_klip_sdi'+suffix+'.fits', left_klip_cube
	
	;medarr, left_klip_cube, medframel


	if comb_type eq 'median' then medarr, left_klip_cube, medframel
	if comb_type eq 'mean' then medframel=mean(left_klip_cube,dim=3)
	if comb_type eq 'nw-mean' then medframel=nw_ang_comb(left_klip_cube,angles)

IF fill and not do_rotate then begin
		leftderotmed = readfits(proc_path+obj+'_left_derot.fits')
		leftadimed = readfits(proc_path+obj+'_left_cADI.fits')
	endif
	
	if fill eq 1  and fill_adi eq 0  then	medframel[where(finite(medframel) eq 0 )]=leftderotmed[where(finite(medframel) eq 0 )]
	if fill eq 1  and fill_adi eq 1  then	medframel[where(finite(medframel) eq 0 )]=leftadimed[where(finite(medframel) eq 0 )]


	if not fill then for xx=0.,1023. do for yy=0.,1023. do if sqrt( (xx-512.)^2. + (yy-512.)^2 ) gt wr*nrings then medframel[xx,yy]=!values.f_nan
	
	
	writefits,proc_path+obj+'_left_klip_sdi'+suffix+'.fits', medframel



if identify_sdi then find_sources, proc_path+obj+'_left_klip_sdi'+suffix+'.fits',sigma=5.0,inrad=annmode_inout[0]+2, reference=reduction_path+obj+'/IRDIS/flux/'+obj+'_left_median.fits',FWHM=((wavel*1.0E-6)/(8.4))*206265./0.01225,platescale=0.01225,clip_level=left_klip_sdi_clip_level,stretch=left_klip_sdi_stretch,correction_factor=left_klip_sdi_correction_factor,outrad=annmode_inout[1]-2
	
	medframel[where(finite(medframel) eq 0 )]=0.;rightadimed[where(finite(medframel) eq 0 )]

	;identical code for right
	

for ii=0, (size(right_cube))(3)-1 do begin

Print, '-------[ SDIKLIPing right image ', ii, ' out of ', (size(right_cube))(3)-1, ' ]----------'

if not hyp and not annmode then right_klip_cube[*,*,ii] = rdiklip(right_cube,left_cube_sc, k_klip_sdi, target=ii, anglemask=anglemask, distmask=distmask, posang=posang, wl=waver,diam=8.2, pixelscale=0.01225, angsep=angsep,obj=obj,nrings=nrings, wr =wr, n_ang =n_ang) 

if hyp and not annmode then right_klip_cube[*,*,ii] = rdiklip(right_cube,left_cube_sc, k_klip_sdi, target=ii, anglemask=anglemask, distmask=distmask, posang=posang, wl=waver,diam=8.2, pixelscale=0.01225, angsep=angsep,obj=obj,nrings=nrings, wr =wr, n_ang =n_ang, spot_radius=spot_radius,rho=rho,phi=phi, /hyper) 

if annmode then right_klip_cube[*,*,ii] = rdiklip(right_cube,left_cube_sc, k_klip_sdi, target=ii, anglemask=anglemask, distmask=distmask, posang=posang, wl=waver,diam=8.2, pixelscale=0.01225, angsep=angsep,obj=obj,nrings=nrings, wr =wr, n_ang =n_ang,  annmode_inout=annmode_inout) 


endfor

for ii=0, (size(right_cube))(3)-1 do begin

	print, 'Rotating by ', angles[ii]
	right_klip_cube[*,*,ii]=rot(right_klip_cube[*,*,ii],-angles[ii],/interp)

endfor

if chop eq 1 then right_klip_cube=right_klip_cube>0


	if output_cubes then writefits,proc_path+obj+'_right_cube_klip_sdi'+suffix+'.fits', right_klip_cube
	
	;medarr, right_klip_cube, medframer


	if comb_type eq 'median' then medarr, right_klip_cube, medframer
	if comb_type eq 'mean' then medframer=mean(right_klip_cube,dim=3)
	if comb_type eq 'nw-mean' then medframer=nw_ang_comb(right_klip_cube,angles)

IF fill and not do_rotate then begin
		rightderotmed = readfits(proc_path+obj+'_right_derot.fits')
		rightadimed = readfits(proc_path+obj+'_right_cADI.fits')
	endif
	
		if fill eq 1   and fill_adi eq 0 then	medframer[where(finite(medframer) eq 0 )]=rightderotmed[where(finite(medframer) eq 0 )]
		if fill eq 1 and fill_adi eq 1then	medframer[where(finite(medframer) eq 0 )]=rightadimed[where(finite(medframer) eq 0 )]

	
	writefits,proc_path+obj+'_right_klip_sdi'+suffix+'.fits', medframer



if identify_sdi then find_sources, proc_path+obj+'_right_klip_sdi'+suffix+'.fits',sigma=5.0,inrad=annmode_inout[0]+2, reference=reduction_path+obj+'/IRDIS/flux/'+obj+'_right_median.fits',FWHM=((waver*1.0E-6)/(8.4))*206265./0.01225,platescale=0.01225,clip_level=right_klip_sdi_clip_level,stretch=right_klip_sdi_derot_stretch,correction_factor=right_klip_sdi_correction_factor,outrad=annmode_inout[1]-2

	
	
	medframer_sc=rot(medframer,0.,wavel/waver,/interp)
	medframel_sc=rot(medframel,0.,waver/wavel,/interp)
	
	medframel_sc[where(finite(medframel_sc) eq 0)]=0.
	medframer_sc[where(finite(medframer_sc) eq 0)]=0.

	
	
	
	medframer[where(finite(medframer) eq 0 )]=0.;rightadimed[where(finite(medframer) eq 0 )]

	
	total_image=(medframel+medframer)/2.

	writefits,strcompress(proc_path+obj+'_left+right_klip_sdi'+suffix+'.fits'), total_image



if identify_sdi then find_sources, proc_path+obj+'_left+right_klip_sdi'+suffix+'.fits',sigma=5.0,inrad=annmode_inout[0]+2, reference=reduction_path+obj+'/IRDIS/flux/'+obj+'_left+right_median.fits',FWHM=((((wavel+waver)/2.)*1.0E-6)/(8.4))*206265./0.01225,platescale=0.01225,clip_level=total_klip_sdi_clip_level,stretch=total_klip_sdi_stretch,correction_factor=total_klip_sdi_correction_factor,outrad=annmode_inout[1]-2


;convolve images


wave=(wavel+waver)/2 ;µm


sz=final_image_half_size*2.+1.;511

width=(wave*1E-6) / (8.2) * 206265. / 0.01225

print, 'PSF Width: ',width

PSF = psf_Gaussian(NPIX=sz, FWHM=[width,width])

; Normalize the PSF

PSFN = PSF/MAX(PSF)


medframel_c = convolve(medframel, PSFN)
rightmed_c = convolve(medframer, PSFN)

writefits,proc_path+obj+'_left_klip_sdi_conv'+suffix+'.fits', medframel_c, hdr
writefits,proc_path+obj+'_right_klip_sdi_conv'+suffix+'.fits', rightmed_c, hdr

total_conv_image=(medframel_c+rightmed_c)/2.
writefits,proc_path+obj+'_left+right_klip_sdi_conv'+suffix+'.fits', total_conv_image, hdr



if notify_sdi then begin


body='KLIP SDI reduction with high pass filter width = '+string(high_pass_width)+' --- Message automatically generated on ' + String(systime())
notify_string='echo "'+body+'" | mail -s "ZONASPHERE '+ obj + ' Data Reduction Results" kwagner@as.arizona.edu -A '+strcompress(proc_path+obj+'_left_klip_sdi'+suffix+'.fits_source_map.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_left+right_klip_sdi'+suffix+'.fits_source_map.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_right_klip_sdi'+suffix+'.fits_source_map.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_left_klip_sdi'+suffix+'.fits_sources.txt',/rem) + $
' -A '+strcompress(proc_path+obj+'_left+right_klip_sdi'+suffix+'.fits_sources.txt',/rem) + $
' -A '+strcompress(proc_path+obj+'_right_klip_sdi'+suffix+'.fits_sources.txt',/rem) + $
' -A '+strcompress(proc_path+obj+'_right_klip_sdi'+suffix+'.fits_radial_sensitivity.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_right_klip_sdi'+suffix+'.fits_radial_sensitivity_2.0_arcsec.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_left_klip_sdi'+suffix+'.fits_radial_sensitivity.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_left_klip_sdi'+suffix+'.fits_radial_sensitivity_2.0_arcsec.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_left+right_klip_sdi'+suffix+'.fits_radial_sensitivity.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_left+right_klip_sdi'+suffix+'.fits_radial_sensitivity_2.0_arcsec.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_right_klip_sdi'+suffix+'.fits_radial_sensitivity.png',/rem) + $
' -A '+strcompress(proc_path+obj+'_right_klip_sdi'+suffix+'.fits_radial_sensitivity_2.0_arcsec.png',/rem) 

spawn, notify_string



restore,proc_path+obj+'_left_klip_sdi'+suffix+'.fits_sources.sav'

;select candidates based on contrast and separation

candidates = where(sources_rho le 2.0 and sources_rho ge 0.18)

;candidates=indgen(7) 

if n_elements(candidates) eq 1 then if candidates eq -1 then candidates=!null


print, candidates
if sources_x eq 0. and sources_y eq 0. then candidates=!null
if candidates ne !NULL then begin



candidates_contrast= sources_contrast[candidates]
candidates_x=sources_x[candidates]
candidates_y=sources_y[candidates]
candidates_rho=sources_rho[candidates]
candidates_theta=sources_theta[candidates]
candidates_sigma=sources_sigma[candidates]


smalls=readfits(proc_path+obj+'_left_klip_sdi'+suffix+'.fits_sources.fits')
smalls=smalls[*,*,candidates]

rho_string='' & sigma_string='' & contrast_string='' & attach_string = '' & theta_string='' & candidates_string=''

cgcleanup

for mmm=0,n_elements(candidates)-1 do begin

plotcolor='Snow'
	cand=smalls[*,*,mmm]
	rho_string=rho_string+string(candidates_rho[mmm])
	theta_string=theta_string+string(candidates_theta[mmm])
	sigma_string=sigma_string+string(candidates_sigma[mmm])
	contrast_string=contrast_string+string(candidates_contrast[mmm])
	candidates_string=candidates_string+string(candidates[mmm]+1)
	
	cgps_open,proc_path+obj+'_left_klip_sdi'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png'


   cgLoadCT, 1;, NColors=100,BOTTOM=100;,BOTTOM=100
	cgDisplay,40,40,/device

	;cgImage,cand;,stretch=stretch;,xrange=[0,(size(smalls))(2)],yrange=[0,(size(smalls))(2)],/keep_aspect_ratio,minvalue=min(cand),maxvalue=max(cand),mean=mean(cand);,position=[0,0,1,1]

cgImage,cand,xrange=[0,39],yrange=[0,39],/keep_aspect_ratio


	cgps_close,/png

	cgcleanup

	attach_string = attach_string+ ' -A '+strcompress(proc_path+obj+'_left_klip_sdi'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png',/rem)
endfor


;for cc=0,n_elements(candidates)-1 do begin

	body='ALERT!'+string(n_elements(candidates))+' Candidate exoplanets discovered around '+obj+' in KLIP SDI reduction with significance = '+sigma_string+' contrast = '+contrast_string+', rho='+rho_string+', theta='+theta_string+' source numbers ='+candidates_string+' --- Message automatically generated by ZONASPHERE on ' + String(systime())

	notify_string='echo "'+body+'" | mail -s "ALERT!!! '+string(n_elements(candidates)) +' Exoplanet Candidates Detected Around '+ obj + '!!!" kwagner@as.arizona.edu -A '+strcompress(proc_path+obj+'_left_klip_sdi'+suffix+'.fits_source_map.png',/rem) +' -A '+strcompress(proc_path+obj+'_left_klip_sdi'+suffix+'.fits_sources.txt',/rem) + ' -A '+strcompress(proc_path+obj+'_left_klip_sdi'+suffix+'.fits_SNR_map.png',/rem)  + attach_string

	spawn, notify_string 

for mmm=0,n_elements(candidates)-1 do file_delete,strcompress(proc_path+obj+'_left_klip_sdi'+suffix+'.fits_sources_'+strcompress((string(mmm+1)),/rem)+'.png',/rem)



;endfor
endif



endif



endif ;end sdi

;--------------------------------------------------------------------------------------------

if do_klip_rdi then begin


ww=final_image_half_size



if skip_rotational_center then begin
	left_cube=readfits(proc_path+obj+'_left_cube_clean.fits',lefthdr)
	right_cube=readfits(proc_path+obj+'_right_cube_clean.fits',righthdr)
endif else begin
	left_cube=readfits(proc_path+obj+'_left_cube_cen.fits',lefthdr)
	right_cube=readfits(proc_path+obj+'_right_cube_cen.fits',righthdr)
endelse


if use_inject then left_cube=readfits(proc_path+obj+'_left_cube_inject.fits',lefthdr)

if use_inject then right_cube=readfits(proc_path+obj+'_right_cube_inject.fits',righthdr)

ref_left_cube=readfits(ref_path+ref_obj+'_left_cube_clean.fits',lefthdr)
ref_right_cube=readfits(ref_path+ref_obj+'_right_cube_clean.fits',righthdr)

restore, filename=proc_path+obj+'_angles_clean.sav'


angles=reform(angles)

;print, size(left_cube)
;print, size(angles)
;print, 'Check cube sizes!!!'
;hak


dbf=esopar(lefthdr,'ESO INS1 OPTI2 NAME ')
if strpos(dbf,'Y23') ne -1 then dual_filter='Y23'

if strpos(dbf,'J23') ne -1 then dual_filter='J23'
if strpos(dbf,'H23') ne -1 then dual_filter='H23'
if strpos(dbf,'H34') ne -1 then dual_filter='H34'
if strpos(dbf,'K12') ne -1 then dual_filter='K12'



if dual_filter eq 'Y23' then wavel=1.022

if dual_filter eq 'K12' then wavel=2.1025
if dual_filter eq 'J23' then wavel=1.1895
if dual_filter eq 'H23' then wavel=1.5888

if dual_filter eq 'Y23' then waver=1.076

if dual_filter eq 'K12' then waver=2.255
if dual_filter eq 'J23' then waver=1.2698
if dual_filter eq 'H23' then waver=1.6671

if bin_rdi gt 1 then begin
	st=1
	for ii=0.,(size(left_cube))(3)-1. do begin
		print, fix(ii+1.) mod fix(bin_rdi)
		if fix(ii+1.) mod fix(bin_rdi) eq 1 then binned=left_cube[*,*,ii] else binned=[ [[binned]],[[ left_cube[*,*,ii] ]] ]
		if fix(ii+1.) mod fix(bin_rdi) eq 0 then begin
		
			print,'Binning left frames...'
			
			if st eq 1 then begin 
			
			medarr,binned,binned
			binned_cube = binned 
			st=0 
			
			 endif else begin
			 medarr,binned,binned

			 binned_cube=[[[binned_cube]],[[binned]]]
			 
			 endelse
			print, size(binned_cube)
		endif
	endfor	
	print, size(binned_cube)
	st=1
	for ii=0,(size(left_cube))(3)-1. do begin
		if fix(ii+1.) mod fix(bin_rdi) eq 1 then binned_angle=angles[ii] else binned_angle=binned_angle+angles[ii]
		if fix(ii+1.) mod fix(bin_rdi) eq 0 then begin
			if st eq 1 then begin binned_angles = binned_angle/bin_sdi & st=0 & endif else binned_angles=[[[binned_angles]],[[binned_angle/bin_sdi]]]
		endif
	endfor	
	
	left_cube=binned_cube
	angles=binned_angles
endif

if bin_rdi gt 1 then begin
	st=1
	for ii=0.,(size(ref_left_cube))(3)-1. do begin
		print, fix(ii+1.) mod fix(bin_rdi)
		if fix(ii+1.) mod fix(bin_rdi) eq 1 then binned=ref_left_cube[*,*,ii] else binned=[ [[binned]],[[ ref_left_cube[*,*,ii] ]] ]
		if fix(ii+1.) mod fix(bin_rdi) eq 0 then begin
		
			print,'Binning left frames...'
			
			if st eq 1 then begin 
			
			medarr,binned,binned
			binned_cube = binned 
			st=0 
			
			 endif else begin
			 medarr,binned,binned

			 binned_cube=[[[binned_cube]],[[binned]]]
			 
			 endelse
			print, size(binned_cube)
		endif
	endfor	
	print, size(binned_cube)
	
	
	ref_left_cube=binned_cube
endif



if destripe_frames ne 0 then begin

	print, 'destriping 90 degrees...'
		for ii=0, (size(ref_left_cube))(3)-1 do ref_left_cube[*,*,ii]=destripe(ref_left_cube[*,*,ii],90.,clip_level=0.0,/nodisp)

	for ii=0, (size(left_cube))(3)-1 do left_cube[*,*,ii]=destripe(left_cube[*,*,ii],90.,clip_level=0.0,/nodisp)
	;print, 'destriping 0 degrees...'
	;for ii=0, (size(left_cube))(3)-1 do left_cube[*,*,ii]=destripe(left_cube[*,*,ii],0.,clip_level=0.0,/nodisp)

endif


if high_pass_width gt 0 then for ii= 0, (size(left_cube))(3)-1 do $
		left_cube[*,*,ii]=left_cube[*,*,ii]-smooth(left_cube[*,*,ii],high_pass_width)
if smooth_width gt 0 then for ii= 0, (size(left_cube))(3)-1 do $
		left_cube[*,*,ii]=smooth(left_cube[*,*,ii],smooth_width)	
		
		
if high_pass_width gt 0 then for ii= 0, (size(ref_left_cube))(3)-1 do $
		ref_left_cube[*,*,ii]=ref_left_cube[*,*,ii]-smooth(ref_left_cube[*,*,ii],high_pass_width)
if smooth_width gt 0 then for ii= 0, (size(ref_left_cube))(3)-1 do $
		ref_left_cube[*,*,ii]=smooth(ref_left_cube[*,*,ii],smooth_width)	
		
;left_cube_pupil=left_cube
left_klip_cube=left_cube




if bin_rdi gt 1 then begin
	st=1
	for ii=0.,(size(right_cube))(3)-1. do begin
		print, fix(ii+1.) mod fix(bin_rdi)
		if fix(ii+1.) mod fix(bin_rdi) eq 1 then binned=right_cube[*,*,ii] else binned=[ [[binned]],[[ right_cube[*,*,ii] ]] ]
		if fix(ii+1.) mod fix(bin_rdi) eq 0 then begin
		
			print,'Binning right frames...'
			if st eq 1 then begin 
			medarr,binned,binned
			binned_cube = binned 
			 st=0 
			 endif else begin
			 medarr,binned,binned

			 binned_cube=[[[binned_cube]],[[binned]]]
			 endelse
			print, size(binned_cube)
		endif
	endfor	
	print, size(binned_cube)
	
	
	right_cube=binned_cube
endif

if bin_rdi gt 1 then begin
	st=1
	for ii=0.,(size(ref_right_cube))(3)-1. do begin
		print, fix(ii+1.) mod fix(bin_rdi)
		if fix(ii+1.) mod fix(bin_rdi) eq 1 then binned=ref_right_cube[*,*,ii] else binned=[ [[binned]],[[ ref_right_cube[*,*,ii] ]] ]
		if fix(ii+1.) mod fix(bin_rdi) eq 0 then begin
		
			print,'Binning left frames...'
			
			if st eq 1 then begin 
			
			medarr,binned,binned
			binned_cube = binned 
			st=0 
			
			 endif else begin
			 medarr,binned,binned

			 binned_cube=[[[binned_cube]],[[binned]]]
			 
			 endelse
			print, size(binned_cube)
		endif
	endfor	
	print, size(binned_cube)
	
	
	ref_right_cube=binned_cube
endif


if destripe_frames ne 0 then begin

	print, 'destriping 90 degrees...'
		for ii=0, (size(ref_right_cube))(3)-1 do ref_right_cube[*,*,ii]=destripe(ref_right_cube[*,*,ii],90.,clip_level=0.0,/nodisp)

	for ii=0, (size(right_cube))(3)-1 do right_cube[*,*,ii]=destripe(right_cube[*,*,ii],90.,clip_level=0.0,/nodisp)
	;print, 'destriping 0 degrees...'
	;for ii=0, (size(right_cube))(3)-1 do right_cube[*,*,ii]=destripe(right_cube[*,*,ii],0.,clip_level=0.0,/nodisp)

endif


if high_pass_width gt 0 then for ii= 0, (size(ref_right_cube))(3)-1 do $
		ref_right_cube[*,*,ii]=ref_right_cube[*,*,ii]-smooth(ref_right_cube[*,*,ii],high_pass_width)
if smooth_width gt 0 then for ii= 0, (size(ref_right_cube))(3)-1 do $
		ref_right_cube[*,*,ii]=smooth(ref_right_cube[*,*,ii],smooth_width)
		

if high_pass_width gt 0 then for ii= 0, (size(right_cube))(3)-1 do $
		right_cube[*,*,ii]=right_cube[*,*,ii]-smooth(right_cube[*,*,ii],high_pass_width)
if smooth_width gt 0 then for ii= 0, (size(right_cube))(3)-1 do $
		right_cube[*,*,ii]=smooth(right_cube[*,*,ii],smooth_width)
		
right_klip_cube=right_cube
;right_cube_pupil=right_cube

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;function RDIKLIP, cube, refcube, k_klip, target=target, checkpoint = checkpoint, anglemask=anglemask, $
;distmask=distmask, posang=posang,scaling=scaling, diam=diam,pixelscale=pixelscale, wl=wl,angsep=angsep,hyper=hyper,$
;obj=obj,nrings=nrings,wr =wr,n_ang =n_ang




posang=angles
for ii=0, (size(left_cube))(3)-1 do begin

Print, '-------[ RDI KLIPing left image ', ii, ' out of ', (size(left_cube))(3)-1, ' ]----------'

if not hyp and not annmode then left_klip_cube[*,*,ii] = rdiklip(left_cube, ref_left_cube,k_klip_rdi, target=ii, anglemask=anglemask, distmask=distmask, posang=posang, wl=wavel,diam=8.2, pixelscale=0.01225, angsep=angsep, obj=obj,nrings=nrings, wr =wr, n_ang =n_ang) 

if hyp and not annmode then left_klip_cube[*,*,ii] = rdiklip(left_cube, ref_left_cube, k_klip_rdi, target=ii, anglemask=anglemask, distmask=distmask, posang=posang, wl=wavel,diam=8.2, pixelscale=0.01225, angsep=angsep, obj=obj,nrings=nrings, wr =wr, n_ang =n_ang, spot_radius=spot_radius,rho=rho,phi=phi, /hyper) 

if annmode then left_klip_cube[*,*,ii] = rdiklip(left_cube, ref_left_cube,k_klip_rdi, target=ii, anglemask=anglemask, distmask=distmask, posang=posang, wl=wavel,diam=8.2, pixelscale=0.01225, angsep=angsep, obj=obj,nrings=nrings, wr =wr, n_ang =n_ang, annmode_inout=annmode_inout) 

endfor

left_klip_cube[where(left_klip_cube eq 0)]=!values.f_nan
right_klip_cube[where(right_klip_cube eq 0)]=!values.f_nan


for ii=0, (size(left_cube))(3)-1 do begin
	print, 'Rotating by ', angles[ii]
	left_klip_cube[*,*,ii]=rot(left_klip_cube[*,*,ii],-angles[ii],/interp)

endfor

if chop eq 1 then left_klip_cube=left_klip_cube>0


	if output_cubes then writefits,proc_path+obj+'_left_cube_klip_rdi.fits', left_klip_cube
	
	;medarr, left_klip_cube, medframel


	if comb_type eq 'median' then medarr, left_klip_cube, medframel
	if comb_type eq 'mean' then medframel=mean(left_klip_cube,dim=3)
	if comb_type eq 'nw-mean' then medframel=nw_ang_comb(left_klip_cube,angles)
	
	if fill eq 1 and hyp eq 1 and do_rotate eq 1 and fill_adi eq 0  then	medframel[where(finite(medframel) eq 0 )]=leftderotmed[where(finite(medframel) eq 0 )]
	if fill eq 1 and hyp eq 1 and fill_adi eq 1 and do_rotate eq 1 then	medframel[where(finite(medframel) eq 0 )]=leftadimed[where(finite(medframel) eq 0 )]


	
	writefits,proc_path+obj+'_left_klip_rdi.fits', medframel
	
	medframel[where(finite(medframel) eq 0 )]=0.;rightadimed[where(finite(medframel) eq 0 )]

	;identical code for right
	

for ii=0, (size(right_cube))(3)-1 do begin

Print, '-------[ RDIKLIPing right image ', ii, ' out of ', (size(right_cube))(3)-1, ' ]----------'

if not hyp and not annmode then right_klip_cube[*,*,ii] = rdiklip(right_cube,ref_right_cube, k_klip_rdi, target=ii, anglemask=anglemask, distmask=distmask, posang=posang, wl=waver,diam=8.2, pixelscale=0.01225, angsep=angsep,obj=obj,nrings=nrings, wr =wr, n_ang =n_ang) 

if hyp and not annmode then right_klip_cube[*,*,ii] = rdiklip(right_cube,ref_right_cube, k_klip_rdi, target=ii, anglemask=anglemask, distmask=distmask, posang=posang, wl=waver,diam=8.2, pixelscale=0.01225, angsep=angsep,obj=obj,nrings=nrings, wr =wr, n_ang =n_ang, spot_radius=spot_radius,rho=rho,phi=phi, /hyper) 

if annmode then right_klip_cube[*,*,ii] = rdiklip(right_cube,ref_right_cube, k_klip_rdi, target=ii, anglemask=anglemask, distmask=distmask, posang=posang, wl=waver,diam=8.2, pixelscale=0.01225, angsep=angsep,obj=obj,nrings=nrings, wr =wr, n_ang =n_ang,  annmode_inout=annmode_inout) 


endfor

for ii=0, (size(right_cube))(3)-1 do begin
	print, 'Rotating by ', angles[ii]
	right_klip_cube[*,*,ii]=rot(right_klip_cube[*,*,ii],-angles[ii],/interp)

endfor

if chop eq 1 then right_klip_cube=right_klip_cube>0


	if output_cubes then writefits,proc_path+obj+'_right_cube_klip_rdi.fits', right_klip_cube
	
	;medarr, right_klip_cube, medframer


	if comb_type eq 'median' then medarr, right_klip_cube, medframer
	if comb_type eq 'mean' then medframer=mean(right_klip_cube,dim=3)
	if comb_type eq 'nw-mean' then medframer=nw_ang_comb(right_klip_cube,angles)
	
		if fill eq 1 and hyp eq 1 and do_rotate eq 1  and fill_adi eq 0 then	medframer[where(finite(medframer) eq 0 )]=rightderotmed[where(finite(medframer) eq 0 )]
		if fill eq 1 and hyp eq 1 and fill_adi eq 1 and do_rotate eq 1 then	medframer[where(finite(medframer) eq 0 )]=rightadimed[where(finite(medframer) eq 0 )]

	
	writefits,proc_path+obj+'_right_klip_rdi.fits', medframer
	
	
	medframer_sc=rot(medframer,0.,wavel/waver,/interp)
	medframel_sc=rot(medframel,0.,waver/wavel,/interp)
	
	medframel_sc[where(finite(medframel_sc) eq 0)]=0.
	medframer_sc[where(finite(medframer_sc) eq 0)]=0.

	
	
	
	medframer[where(finite(medframer) eq 0 )]=0.;rightadimed[where(finite(medframer) eq 0 )]

	
	total_image=(medframel+medframer)/2.

	writefits,strcompress(proc_path+obj+'_left+right_klip_rdi'+suffix+'.fits'), total_image

;convolve images


wave=(wavel+waver)/2 ;µm


sz=final_image_half_size*2.+1.;511

width=(wave*1E-6) / (8.2) * 206265. / 0.01225

print, 'PSF Width: ',width

PSF = psf_Gaussian(NPIX=sz, FWHM=[width,width])

; Normalize the PSF

PSFN = PSF/MAX(PSF)


medframel_c = convolve(medframel, PSFN)
rightmed_c = convolve(medframer, PSFN)

writefits,proc_path+obj+'_left_klip_rdi_conv.fits', medframel_c, hdr
writefits,proc_path+obj+'_right_klip_rdi_conv.fits', rightmed_c, hdr

total_conv_image=(medframel_c+rightmed_c)/2.
writefits,proc_path+obj+'_left+right_klip_rdi_conv.fits', total_conv_image, hdr






endif


;------------------------------[ Begin Thought Process  ]---------------------------------

if test_residuals eq 1 then begin


targx=526.-1.
targy=449.-1.
resrad=20000.

residual_mask=medframel
residual_mask[*]=0
cenx=float((size(medframel))(1)/2.-1.)
ceny=float((size(medframel))(1)/2.-1.)
for xx=0., (size(medframel))(1)-1. do begin
	for yy=0., (size(medframel))(2)-1. do begin
		dx=targx-xx
		dy=targy-yy
		if sqrt(dx*dx+dy*dy) lt resrad then residual_mask[xx,yy]=1.		
	endfor
endfor


leftres=residual_mask*medframel
if not skip_right then rightres=residual_mask*medframer

combres=leftres+rightres


writefits,proc_path+obj+'_left_klip_res.fits', leftres, hdr
writefits,proc_path+obj+'_right_klip_res.fits', rightres, hdr

if high_pass_width_final gt 0. then leftres=leftres-smooth(leftres,high_pass_width_final)
if high_pass_width_final gt 0. then rightres=rightres-smooth(rightres,high_pass_width_final)
if high_pass_width_final gt 0. then combres=combres-smooth(combres,high_pass_width_final)

;resl=total(abs(leftres))
;resr=total(abs(rightres))
;resc=total(abs(combres))

resl=stddev(leftres,/nan)
resr=stddev(rightres,/nan)
resc=stddev(combres,/nan)

if runs eq 0 then bestl=resl
if runs eq 0 then bestr=resr
if runs eq 0 then lresar=resl else lresar=[lresar,resl]
if runs eq 0 then rresar=resr else rresar=[rresar,resr]
if runs eq 0 then cresar=resc else cresar=[cresar,resc]



if runs eq 0 then rhoar=rhop else rhoar=[rhoar,rhop]
if runs eq 0 then thetaar=theta else thetaar=[thetaar,theta]
if runs eq 0 then contarl=contrast_left else contarl=[contarl,contrast_left]
if runs eq 0 then contarr=contrast_right else contarr=[contarr,contrast_right]

if runs eq 0 then lefts=leftres else lefts=[[[lefts]],[[leftres]]]
if runs eq 0 then rights=rightres else rights=[[[rights]],[[rightres]]]



if resl le bestl then begin
	
	bestleft=medframel
	bestl=resl


endif

if resr le bestr then begin
	
	bestright=medframer
	bestr=resr


endif



endif
if test_residuals then !p.linestyle=0
;left=solid

if runs gt 0 then plot, (lresar-mean(lresar))/max(lresar-mean(lresar)), TITLE='Left Residuals (Solid), Right Residuals (Dashed), Combined Residuals (Dot-dashed)',YRANGE=[-1,1]

if test_residuals eq 1 then !p.linestyle=2
;right=dashed
if runs gt 0 then oplot, (rresar-mean(rresar))/max(rresar-mean(rresar)) ;right is red

if test_residuals eq 1 then !p.linestyle=3

if runs gt 0 then oplot, (cresar-mean(cresar))/max(cresar-mean(cresar)) ;right is red



endfor; runs for

if test_residuals eq 1 then begin

;write the best files
;print the residuals

;show the best residual
writefits,  proc_path+obj+'_best_left+right_klip_median'+suffix+'.fits', (bestleft+bestright)/2.   

if nruns gt 1 then writefits,  proc_path+obj+'_tests_left'+suffix+'.fits', lefts
if nruns gt 1 then writefits,  proc_path+obj+'_tests_right'+suffix+'.fits', rights


print, 'Residuals from all left tests:', lresar
print, 'Residuals from all right tests:', rresar

min_lres=min(lresar,min_lel)
min_rres=min(rresar,min_rel)
min_cres=min(cresar,min_cel)



print, 'Min residual Left = ',min_lres, ' @ position ', min_lel
print, '--[ Best Left parameters:'
print, '----Contrast: ', contarl[min_lel]
print, '----Radius: ',rhoar[min_lel]
print, '----Theta:',thetaar[min_lel]


print, 'Min residual Right = ',min_rres, ' @ position ', min_rel
print, '--[ Best Right parameters:'
print, '----Contrast: ', contarr[min_rel]
print, '----Radius: ',rhoar[min_rel]
print, '----Theta:',thetaar[min_rel]



print, 'Min residual combined = ',min_cres, ' @ position ', min_cel
print, '--[ Best Right parameters:'
print, '----Contrast: ', contarr[min_cel]
print, '----Radius: ',rhoar[min_cel]
print, '----Theta:',thetaar[min_cel]



endif



end
