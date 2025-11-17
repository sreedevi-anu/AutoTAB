FUNCTION HDR2COORD, header, lat=lat, lon=lon
	
	wcs   =  FITSHEAD2WCS(header)
	coord =  WCS_GET_COORD(wcs)	;spherical cordinate
	WCS_CONVERT_FROM_COORD, wcs, coord,'hg', lon, lat
	RETURN,0
end

;Step 1 : Extract the binary files of BMR's detected
data_dir='1cycle_save/bmr_save/'
sav = file_search(data_dir,'*.sav',count=nn)

;Create a new text file for purpose of tracking
;openw,1,'BMR_id.txt' & close,1

kk=3613

for i=41594ul, nn-1 do begin

	;Recreate the binary file
	restore, sav[i]
	data = fltarr(hdr_los.naxis1, hdr_los.naxis2)
	data[bmr_ind] = 1
	ll = label_region(data,/all)
	header_day1 =hdr_los
	time1 = anytim2tai(hdr_los.date_obs)
	
	;Create an array to contain tracked information
	mask = fltarr(hdr_los.naxis1, hdr_los.naxis2)

	
	for j=1, max(ll) do begin
		;isolate each BMR 
		mask_los = fltarr(hdr_los.naxis1, hdr_los.naxis2)
		index = where(ll eq j,cnt)
		mask_los[index] = 1
		area = cnt
		
		;Check if the BMR has been already tracked, if yes we remove that
		readcol,'BMR_id.txt', spot_id,id, begin_id,/silent,/NAN,format='(a,ul,a)' 
		if N_elements(spot_id) ne 0 then begin
			new_id = header_day1.date_obs+'_'+strtrim(j,2)
			void = where(spot_id eq new_id, cnt)
			;BMR has been already tracked
			if cnt ne 0 then begin
				print,'Already tracked with id',id[void[0]]
				print,i,cnt, id[void[0]]
				mask_los[index]=0
				continue
			endif
		endif
		
		;making map of binary file
		index2map, header_day1, mask_los, mask_map
		
		!p.multi=[0,4,5]
		window, 0, xs=500, ys=500
		plot_map, mask_map,  dmax=1,  /grid, grid_spacing=30,/notitle
		
		update = hdr_los.date_obs+'_'+strtrim(j,2)
		
		void=hdr2coord(hdr_los, lat=lat, lon=lon)
		lat = mean(lat[index]) & lon = mean(lon[index])
		rate = 13.45 - 3.0 * sin(lat*!DTOR)^2                           ; diff_rot rate
		max_mag = fix(((90.0 - lon)/rate*24.0)/1.6)                    ; maximum possible number of magneotagram bmr can be tracked 
		if max_mag eq 0 then continue
		tmax = ((90.0 - lon)/rate*24.0)                                 ; t_max

		track = fltarr(max_mag)
	
		missed_cnt = 0                                                  ; Number of magnetograms where BMR's are not identified in between
		for k=i+1, i+max_mag do begin
			restore, sav[k]
			track_los = fltarr(hdr_los.naxis1, hdr_los.naxis2)      ; make copy of the days magnetogram
			track_los[bmr_ind] = 1
			time_track = anytim2tai(hdr_los.date_obs)
			
			dt=(time_track - time1)/3600.0                          ; check time      
			if dt gt tmax then break
			
			index2map, hdr_los, track_los, track_map
			plot_map, track_map,  dmax=1,  /grid, grid_spacing=30,/notitle
			
			lbl = label_region(track_los,/all)
			drotmap = drot_map(mask_map, ref_map=track_map)
			
		
			sum = drotmap.data + track_map.data
			ind = where(sum eq 2,cnt1)                               ;Check overlap
			
			
			if cnt1 lt 150 then begin                             
				missed_cnt += 1
				track[k-(i+1)] = 0
			endif
			
			if cnt1 ge 150 then begin
				if missed_cnt lt 30 then begin
					track[k-(i+1)]=1
					over = max(lbl[ind])
					openw, 1,'BMR_id.txt',/append
					printf, 1, hdr_los.date_obs+'_'+strtrim(over,2), 10000+kk, '    '+header_day1.date_obs+'_'+strtrim(j,2)
					close,1
					missed_cnt = 0
				endif
				
				if missed_cnt ge 30 then begin
					if cnt1 lt 250 then begin
						missed_cnt += 1
						track[k-(i+1)] = 0
					endif
					
					if cnt1 ge 250 then begin
						track[k-(i+1)]=1
						over = max(lbl[ind])
						openw, 1,'BMR_id.txt',/append
						printf, 1, hdr_los.date_obs+'_'+strtrim(over,2), 10000+kk, '    '+header_day1.date_obs+'_'+strtrim(j,2)
						close,1
						missed_cnt = 0
					endif
			        endif
			endif
			
		endfor
		
		b = where(track eq 1, cnt2)
		if cnt2 gt 0 then  begin
		  	mask[index] = 10000+kk
		  	kk += 1
		endif
		
		if cnt2 lt 0 then continue
		
		
	endfor

	name = 'mdi_bmr_ind_'+time2file(header_day1.date_obs)+'.sav'
	save, mask, header_day1, filename='1cycle_save/track_save/' + name
endfor

			
end

			
