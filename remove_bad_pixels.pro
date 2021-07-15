pro remove_bad_pixels, badmap, path, fullmode=fullmode

;A bad pixel removal script for use with the Scorpion data reduction pipeline
;Author: Kevin Wagner

print, 'Bad pixel remover called upon.'

files=FILE_SEARCH(path,'*_.fits',COUNT=cnt)
if cnt eq 0 then files=FILE_SEARCH(path,COUNT=cnt)
print, cnt

for ff=0,cnt-1 do begin
cube=readfits(files[ff],hdr)

print, 'Removing bad pixels from file ',ff+1,'/',cnt,'.'

if  keyword_set(fullmode) then begin 

badmap=fltarr((size(cube))(1),(size(cube))(2))

badmap[*]=1

if n_elements(size(cube)) gt 5 then smallcube = mean(cube,dim=3) else smallcube=cube
 print, 'Finding bad pixels... may take a while....'
  for ix=cs, (size(cube))(1)-1-cs do begin
     for iy=cs, (size(cube))(2)-1-cs do begin
           box = smallcube[ix-cs:ix+cs,iy-cs:iy+cs]
           this =  box[cs,cs]  ; this pixel
           box[cs,cs] = !values.f_nan
           if mean(this) gt (mean(box,/nan)+badsig*stddev(box,/nan)) or mean(this) lt (mean(box,/nan)-badsig*stddev(box,/nan)) then bad=1 else bad=0
           if bad then badmap[ix,iy] = 0           
        endfor
     endfor
endif

; BAD PIXEL REMOVAL
  print
  print, ' Removing bad pixels...'
  
  
if n_elements(size(cube)) gt 5 then begin
  for ii=0, (size(cube))(3)-1 do  begin
	print, 'Working on frame ', ii, ' / ', (size(cube))(3)-1
	
  		fixpix,cube[*,*,ii],badmap,outimg
  		;outimg = maskinterp(cube[*,*,ii],badmap,3,6,"splinterp",gpix=10,gpoints=5,cdis=2)
  		cube[*,*,ii]=outimg


  	endfor
endif else begin
		print, 'Working on single frame ...'
	
  		fixpix,cube[*,*],badmap,outimg
  		;outimg = maskinterp(cube[*,*],badmap,3,6,"splinterp",gpix=10,gpoints=5,cdis=2)
  		cube[*,*]=outimg

endelse
 writefits,files[ff],cube,hdr
  
  endfor
  
end
