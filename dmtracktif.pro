; dmtracktif
; written by David Meer 6/2026
; Tracks droplets moving in two dimensions and finds many of their dynamic properties
;
; data is stored in 89 columns.
; List of variables, stored as data(n,*):
; 0-> x center of mass
; 1-> y center of mass
; 2-> frame number (duplicates out)
; 3-> index of blob within frame (vaguely sorted by y values)
; 4-> area
; 5-> perimeter
; 6-> Symmetry about neck axis. Used to be second perimeter method (not used)
; 7-> shape parameter, p^2/4*pi*A
; 8-> distance to nearest obstacle
; 9-> symmetry about the x axis<-typical one
; 10-> symmetry about the y axis
; 11-> thinnest portion of neck IF neck exists. Assuming perfectly downstream neck
; 12-> is the blob near the edge? 0=no (and therefore good), 1=yes
; 13-> number of neighbors
; 14-> index in frame of neighbor 1
; 15-> index in frame of neighbor 2
; 16-> index in frame of neighbor 3
; 17-> index in frame of neighbor 4
; 18-> index in frame of neighbor 5
; 19-> index in frame of neighbor 6
; 20-> interface length with neighbor 1
; 21-> interface length with neighbor 2
; 22-> interface length with neighbor 3
; 23-> interface length with neighbor 4
; 24-> interface length with neighbor 5
; 25-> interface length with neighbor 6
; 26-> neck distance
; 27-> neck angle
; 28-> position of neck 1x
; 29-> position of neck 1y
; 30-> position of neck 2x
; 31-> position of neck 2y
; 32-> Are they in contact with an obstacle? 0=no, 1=yes
; 33-> number of obstacles touched
; 34-> index of obstacle 1
; 35-> index of obstacle 2
; 36-> index of obstacle 3
; 37-> index of obstacle 4
; 38-> index of obstacle 5
; 39-> index of obstacle 6
; 40-> Actual time (at 30fps)
; ***41-> blob index through time***
; 42-> x velocity between this frame and previous frame
; 43-> y velocity
; 44-> magnitude of velocity
; 45-> fate parameter. 0=1 frame, 1=offscreen 2=in frame, 3=combin, 4=brkup, 5=too big
; 46-> time at begin
; 47-> time at death
; 48-> time since birth
; 49-> distance to endpoint
; 50-> index in time of neighbor 1
; 51-> index in time of neighbor 2
; 52-> index in time of neighbor 3
; 53-> index in time of neighbor 4
; 54-> index in time of neighbor 5
; 55-> index in time of neighbor 6
; 56-> total change in area
; 57-> normalized change in area
; 58-> life parameter. 0=1 frame, 1=offscreen, 2=in frame, 3=combin, 4=brkup
; 59-> area change velocity
; 60-> time in seconds
; 61-> area in meters^2
; 62-> mag. velocity in m/s
; 63-> distance perpendicular since inception in m
; 64-> distance parallel since inception in m
; 65-> wrapped angle
; 66-> are they thin? 0->no, 1->yes
; 67-> are they near something gross? 0->no, 1->yes
; 68-> effective radius
; 69-> velocity of neck
; 70-> index of closest obstacle
; 71-> velocity angle
; 72-> smoothed velocity angle
; 73-> distance to final obstacle
; 74-> smoothed total velocity
; 75-> TODO - SYMMETRY AROUND VELOCITY ANGLE
; 76-> TODO - SYMMETRY AROUND VELOCITY ANGLE OF FIRST CONTACT
; 77-> Maximum extend in X
; 78-> Smoothed symmetry around velocity
; 79-> Smoothed symmetry around angle of first contact
; 80-> symmetry about axis of first contact
; 81-> first x-coordinate in droplet
; 82-> first y-coordinate in droplet
; 83-> offset from obstacle
; 84-> offset at collision
; 85-> index of nearest upstream neighbor
; 86-> edge-to-edge distance to nearest upstream neighbor
; 87-> moment of area symmetry
; 88-> smoothed moment of area
; 89-> moment of area symmetry at first contact
; [xcom,ycom,t,blobindx(i),mass,perim,perimc,calA,sqrt(obsdist),xsymm,ysymm,neckl,edgeflag,nneigh,neighs,t]

;within deth and brth
; 0->index through time
; 1->time
; 2->x
; 3->y
; 4->area
; 5->fate
; 6->life
; 7-12->6 neighbors given by time index
function dmtracktif,fname,directory,collect=collect,xdim=xdim,ydim=ydim,brek=brek,combin=combin,conec=conec,ddrops=ddrops,contacttim=contacttim,ddropcheck=ddropcheck,antoine=antoine,maxdisp=maxdisp,mov=mov,nochug=nochug,bglst=bglst,interface=interface,coal=coal,newptedge=newptedge,nnn=nnn
	if keyword_set(antoine) then direc='/data/antoine/'
	if keyword_set(coal) then direc='/data/david/coal/'+directory+'/'
	if (keyword_set(collect) EQ 0) then begin &$
		direc='/data/david/drops/'+directory+'/' &$
		if keyword_set(coal) then direc='/data/david/coal/'+directory+'/' &$
		if keyword_set(antoine) then direc='/data/antoine/' &$
		ptdata=read_gdf(direc+fname+'ptcollect') &$
	endif else begin &$
		ptdata=collect &$
	endelse
	ptdata1=ptdata([0:40],*)
	ptdata1(2,*)=ptdata(4,*)/4.0
	ptdata1(4,*)=ptdata(40,*)
	ptdata1(40,*)=ptdata(2,*)
	
	pt2=ptdata1(*,where(ptdata1(12,*) EQ 1))

	maxt=max(pt2(-1,*),min=mint)
	if maxt EQ mint then maxt=maxt+1
	nx=n_elements(pt2(*,0))
	tmp=fltarr(nx,maxt-mint)
	tmp[0,*] = max(pt2[0,*])+100.0
	tmp[1,*] = max(pt2[1,*])+100.0
	tmp[-1,*] = findgen(maxt-mint) + mint
 
	pt3=[[pt2],[tmp]]
	s=sort(pt3(-1,*))
	pt3=pt3(*,s)
	if maxt EQ mint+1 then pt3(-1)=pt3(-1)+1
	if NOT keyword_set(maxdisp) then begin &$
		maxdisp=50 &$
		maxdisp1=35 &$
	endif else begin &$
		maxdisp1=maxdisp/2 &$
	endelse &$
	datae=track(pt3,maxdisp1,mem=1,good=1,/quiet)

	data1=track(ptdata1,maxdisp,mem=1,good=0,dim=3,/quiet)
	data2=data1
	data3=datae
	data2(2,*)=data1(40,*)
	data2(4,*)=data1(2,*)*4.0
	data2(40,*)=data1(4,*)
	data3(2,*)=datae(40,*)
	data3(4,*)=datae(2,*)*4.0
	data3(40,*)=datae(4,*)
	leeway=12
	data2=data2(*,where(data2(12,*) EQ 0))
	for i=0,max(data3(41,*)) do begin &$
		locs=where(data3(41,*) EQ i) &$
		if locs(0) NE -1 then locs=locs(sort(data3(2,locs))) &$
		if data3(0,locs(0)) LT 10 then begin &$
			tmppos=[data3(0,locs(-1)),data3(1,locs(-1))] &$
			chckt=data3(2,locs(-1))+1 &$
		endif else begin &$
			tmppos=[data3(0,locs(0)),data3(1,locs(0))] &$
			chckt=data3(2,locs(0))-1 &$
		endelse &$
		mtch=where((abs(tmppos(0)-data2(0,*)) LT leeway) AND (abs(tmppos(1)-data2(1,*)) LT leeway) AND (chckt EQ data2(2,*))) &$
		if length(mtch) gt 1 then mtch=where((abs(tmppos(0)-data2(0,*)) LT leeway/2.0)AND(abs(tmppos(1)-data2(1,*)) LT leeway/2.0)AND(chckt EQ data2(2,*))) &$
		spcr=length(locs) &$
		data4=replicate(!values.F_NAN,length(data2),length(transpose(data2))+spcr) &$
		if mtch(0) NE -1 then begin &$
			data4(*,[0:mtch-1])=data2(*,[0:mtch-1]) &$
			data4(*,[mtch:mtch+spcr-1])=data3(*,[locs]) &$
			data4(41,[mtch:mtch+spcr-1])=data2(41,mtch) &$
			data4(*,[mtch+spcr:length(transpose(data4))-1])=data2(*,[mtch:length(transpose(data2))-1]) &$
		endif else begin &$
			data4(*,[0:length(transpose(data2))-1])=data2 &$
			data4(*,[length(transpose(data2)):length(transpose(data2))+spcr-1])=data3(*,[locs]) &$
			data4(41,[length(transpose(data2)):length(transpose(data2))+spcr-1])=max(data2(41,*))+1 &$
		endelse &$
		data2=data4 &$
		if keyword_set(antoine) then print,[(i+0.0)/max(data3(41,*)),1] &$
	endfor
	mov=read_gdf(direc+'pt'+fname+'gdf')
	b=size(mov)
	xdim=b(1)
	ydim=b(2)
	tdim=b(3)
	data2=data2(*,where((data2(0,*) GE 0) AND (data2(0,*) LE xdim) AND (data2(1,*) GE 0) AND (data2(1,*) LE ydim)))
	sr=FILE_SEARCH(direc+fname+'obs',count=countobs)
	if (countobs GT 0) then begin &$
		obs=float(readtext(direc+fname+'obs')) &$
		xobs1=obs(*,0) &$
;		yobs1=max(obs(*,1))-obs(*,1) &$
		yobs1=obs(*,1) &$
		stup=read_gdf(direc+fname+'obssetup') &$
		scl=stup(0) &$
		xshift=stup(1) &$
		yshift=stup(2) &$
		ang=stup(3) &$
		obsr=stup(4) &$
		pix2um=stup(5) &$
		FPS=stup(6) &$
		xobs2=(xobs1*scl+xshift)*cos(ang*!pi/180.0)-(yobs1*scl+yshift)*sin(ang*!pi/180.0) &$
		yobs2=(xobs1*scl+xshift)*sin(ang*!pi/180.0)+(yobs1*scl+yshift)*cos(ang*!pi/180.0) &$
		xobs=xobs2(where((xobs2 GE 0) and (xobs2 LT xdim) AND (yobs2 GE 0) and (yobs2 LT ydim))) &$
		yobs=yobs2(where((xobs2 GE 0) and (xobs2 LT xdim) AND (yobs2 GE 0) and (yobs2 LT ydim))) &$
	endif else begin &$
		sr=FILE_SEARCH(direc+fname+'cobs',count=countcobs) &$
		if countcobs GT 0 then begin &$
			stup=read_gdf(direc+fname+'cobs') &$
			lstup=n_elements(stup)-3 &$
			xobs=stup([0:(lstup/2)-1]) &$
			yobs=stup([lstup/2:(lstup-1)]) &$
			obsr=stup(-3) &$
			pix2um=stup(-2) &$
			FPS=stup(-1) &$
		endif else begin &$
			if keyword_set(antoine) then print,'obstacle error' &$
			pix2um=5.28 &$
			FPS=60 &$
			xobs=[-100] &$
			yobs=[-100] &$
			obsr=[0] &$
		endelse &$
	endelse
	if NOT keyword_set(gross) then gross=[520,145]
	data2=data2(*,where(data2(2,*) LT tdim))
	lng=length(transpose(data2(0,*)))
	numblob=max(data2(41,*))
	data=replicate(!Values.F_NAN,90,lng)
	data([0:41],*)=data2
	hasneigh=where(data(13,*) GT 0)
	data(9,*)=1-abs(2*(data(9,*)-0.5))
	data(10,*)=1-abs(2*(data(10,*)-0.5))
	data(68,*)=sqrt(data(4,*)/!PI)
	data(67,where(sqrt((data(0,*)-gross(0))^2+(data(1,*)-gross(1))^2) LT data(68,*)+5))=1
	data(67,where(sqrt((data(0,*)-gross(0))^2+(data(1,*)-gross(1))^2) GE data(68,*)+5))=0
	data(66,where(data(26,*)/data(68,*) LT 1))=1
	data(66,where(data(66,*) NE 1))=0
;	tdim=max(data(2,*))+1
	pix2m=pix2um/1000000.0
	data(60,*)=data(40,*)*FPS/tdim
	data(61,*)=data(4,*)*pix2m^2
	visc=50.0*10.0^(-6)
	surftens=11.2*10.0^(-3)
	dens=0.96*10.0^(3)
	dvisc=visc*dens
;t=0
;tv,bytscl(mov(*,*,t))
;cgtext,data(0,where(data(2,*) EQ t))-87,data(1,where(data(2,*) EQ t)),string(round(data(41,where(data(2,*) EQ t)))),/device,charthick=2,color=[500000]
;cgtext,data(0,where(data(2,*) EQ t)),data(1,where(data(2,*) EQ t)),'.',/device,charthick=2,color=[500000]
;dmmov2,mov,data(0,*),data(1,*),data(2,*),data(41,*)
	for i=0,lng-1 do begin &$
		dstemp=sqrt((xobs-data(0,i))^2+(yobs-data(1,i))^2) &$
		data(70,i)=where(dstemp EQ min(dstemp)) &$
	end
;mova=dmreadmp4('/home/dameer/micro/videos/'+fname+'.mp4')
;t=0
;tv,mova(0,*,*,t)
;cgtext,xobs-95,yobs,string([0:length(xobs)-1]),/device,charthick=2,color=[500000]
	for i=0,length(hasneigh)-1 do begin &$
		tim=data(2,hasneigh(i)) &$
		if data(13,hasneigh(i)) EQ 1 then begin &$
			tmpn=where( (data(14,hasneigh(i)) EQ data(3,*)) AND (data(2,*) EQ tim) ) &$
			data(50,hasneigh(i))=data(41,tmpn(0)) &$
		endif else begin &$
			for j=0,data(13,hasneigh(i))-1 do begin &$
				tmph=data(41,where( (data(14+j,hasneigh(i)) EQ data(3,*)) AND (data(2,*) EQ tim) )) &$
				tmph1=tmph(uniq(tmph)) &$
				data(50+j,hasneigh(i))=tmph1 &$
			end &$
		endelse &$
		if keyword_set(antoine) then print,[(i+0.0)/length(hasneigh),2] &$
	endfor
	brth=[]
	deth=[]
	for i=0,numblob do begin &$
		locs=where(data(41,*) EQ i) &$
		if locs(0) NE -1 then locs=locs(sort(data(2,locs))) &$
		if length(locs) GT 0 then begin &$
			data(46,locs)=min(data(2,locs)) &$
			data(47,locs)=max(data(2,locs)) &$
			data(48,locs)=data(47,locs)-data(46,locs) &$
			finalpos=[data(0,locs(where(data(2,locs) EQ max(data(2,locs))))),data(1,locs(where(data(2,locs) EQ max(data(2,locs)))))] &$
			finalspot=locs(where(data(2,locs) EQ max(data(2,locs)))) &$
			firstpos=[data(0,locs(where(data(2,locs) EQ min(data(2,locs))))),data(1,locs(where(data(2,locs) EQ min(data(2,locs)))))] &$
			firstspot=locs(where(data(2,locs) EQ min(data(2,locs)))) &$
			data(49,locs)=sqrt((data(0,locs)-finalpos(0))^2+(data(1,locs)-finalpos(1))^2) &$
			data(63,locs)=pix2m*(data(0,locs)-firstpos(0)) &$
			data(64,locs)=pix2m*(data(1,locs)-firstpos(1)) &$
			if length(locs) GT 1 then begin &$
				srch=where(mov(*,*,data(2,firstspot(0))) EQ data(3,firstspot(0))) &$
				srchx=srch mod xdim &$
				if min(srchx) LT 10 then data(58,firstspot)=1 &$
				srch=where(mov(*,*,data(2,finalspot(0))) EQ data(3,finalspot(0))) &$
				srchx=srch mod xdim &$
				if max(srchx) GT xdim-10 then data(45,finalspot)=1 &$
				if data(58,firstspot(0)) NE 1 then data(58,firstspot)=2 &$
				if data(45,finalspot(0)) NE 1 then data(45,finalspot)=2 &$
				if max(data(4,locs)) GT 70000.000 then data(45,locs)=5 &$
				brth=[[brth],[data(41,locs(0)),data(2,locs(0)),data(0,locs(0)),data(1,locs(0)),data(4,locs(0)),data(45,locs(0)),data(58,locs(0)),$
					data(50,locs(0)),data(51,locs(0)),data(52,locs(0)),data(53,locs(0)),data(54,locs(0)),data(55,locs(0))]] &$
				deth=[[deth],[data(41,locs(-1)),data(2,locs(-1)),data(0,locs(-1)),data(1,locs(-1)),data(4,locs(-1)),data(45,locs(-1)),data(58,locs(-1)),$
					data(50,locs(-1)),data(51,locs(-1)),data(52,locs(-1)),data(53,locs(-1)),data(54,locs(-1)),data(55,locs(-1))]] &$
				for j=1,length(locs)-1 do begin &$
					data(42,locs(j))=(data(0,locs(j))-data(0,locs(j-1)))/(data(2,locs(j))-data(2,locs(j-1))) &$
					data(43,locs(j))=(data(1,locs(j))-data(1,locs(j-1)))/(data(2,locs(j))-data(2,locs(j-1))) &$
					if (data(40,locs(j-1)) mod 1 GT 0) OR (data(40,locs(j)) mod 1 GT 0) then begin &$
						if (data(42,locs(j))-data(42,locs(j-1)) GT data(42,locs(j-1))/1.5) then data(42,locs(j))=data(42,locs(j))/2 &$
						if (data(43,locs(j))-data(43,locs(j-1)) GT data(43,locs(j-1))/1.5) then data(43,locs(j))=data(43,locs(j))/2 &$
					endif &$
					data(44,locs(j))=sqrt(data(42,locs(j))^2+data(43,locs(j-1))^2) &$
					data(56,locs(j))=data(4,locs(j))-data(4,locs(j-1)) &$
					data(57,locs(j))=data(56,locs(j))/data(4,locs(0)) &$
					data(59,locs(j))=(data(4,locs(j))-data(4,locs(j-1)))/(data(40,locs(j))-data(40,locs(j-1))) &$
					data(62,locs(j))=(sqrt(data(42,locs(j))^2+data(43,locs(j-1))^2))*pix2m/(data(60,locs(j))-data(60,locs(j-1))) &$
					data(69,locs(j))=data(26,locs(j))-data(26,locs(j-1)) &$
				end &$
				angl=data(27,locs([0:length(locs)-2]))-data(27,locs([1:length(locs)-1])) &$
				fx=where(abs(angl) GT !pi/2) &$
				nangl=data(27,locs) &$
				if fx(0) GT -1 then begin &$
					for k=0,length(fx)-1 do begin &$
						if angl(fx(k)) GT 0 then begin &$
							nangl([fx(k)+1:(length(nangl)-1)])=nangl([fx(k)+1:(length(nangl)-1)])+!pi &$
						endif else begin &$
							nangl([fx(k)+1:(length(nangl)-1)])=nangl([fx(k)+1:(length(nangl)-1)])-!pi &$
						endelse &$
					end &$
				endif &$
				data(65,locs)=nangl &$
			endif &$
		endif else begin &$
			if locs GE 0 then begin &$
				data(45,locs)=0 &$
				data(58,locs)=0 &$
			endif &$
		endelse &$
		if keyword_set(antoine) then print,[(i+0.0)/numblob,3] &$
	end
	wdirec=where(data(43,*) LT 0)
	data(71,*)=acos(data(42,*)/data(44,*))
	data(71,wdirec)=2*!PI-acos(data(42,wdirec)/data(44,wdirec))
	data(71,where(data(71,*) GT !pi))=data(71,where(data(71,*) GT !pi))-2*!PI
	for i=0,numblob do begin &$
		locs=where(data(41,*) EQ i) &$
		if locs(0) NE -1 then locs=locs(sort(data(2,locs))) &$
		chck=where(data(71,locs) GT -10) &$
		if length(chck) GT 2 then data(72,locs)=smooth(data(71,locs),3,/edge_zero,/NaN) &$
	end
;all the values i need to fix because of dropped frames
	crrct=[42,43,44,62,69]
	tlst=[]
	for i=1,tdim-2 do begin &$
		m1=mean(data(44,where((data(2,*) EQ i) AND (data(44,*) GE 0)))) &$
		m2=mean(data(44,where((data(2,*) EQ i+1) AND (data(44,*) GE 0)))) &$
		m0=mean(data(44,where((data(2,*) EQ i-1) AND (data(44,*) GE 0)))) &$
		if m1/[(m0+m2)/2.0] GT 1.6 then tlst=[tlst,i] &$
	end
	for i=0,numblob do begin &$
		locs=where(data(41,*) EQ i) &$
		if locs(0) NE -1 then locs=locs(sort(data(2,locs))) &$
		for j=0,length(tlst)-1 do begin &$
			ttemp=where(data(2,locs) EQ tlst(j)) &$
			ttemp1=where(data(2,locs) EQ tlst(j)-1) &$
			ttemp2=where(data(2,locs) EQ tlst(j)+1) &$
			if (ttemp(0) NE -1) AND (ttemp1(0) NE -1) AND (ttemp2(0) NE -1) then begin &$
				for k=0,length(crrct)-1 do begin &$
					data(crrct(k),locs(ttemp))=mean([transpose(data(crrct(k),locs(ttemp1))),transpose(data(crrct(k),locs(ttemp2)))]) &$
				endfor &$
			endif &$
		end &$
		if keyword_set(antoine) then print,[(i+0.0)/numblob,4] &$
	end
;smooth out the final velocities
	for i=0,numblob do begin &$
		locs=where(data(41,*) EQ i) &$
		if locs(0) NE -1 then locs=locs(sort(data(2,locs))) &$
		if locs(0) NE -1 then begin &$
			if length(locs) gt 5 then begin &$
				smth=smooth(transpose([data(44,locs)]),5,/EDGE_TRUNCATE) &$
				smth(where(~finite(smth), /null))=data(44,locs(where(~finite(smth), /null))) &$
				data(74,locs)=smth &$
			endif &$
		endif &$
		if keyword_set(antoine) then print,[(i+0.0)/numblob,5] &$
	end
;symemtry parameter based on velocity angle
;while I'm digging around in there, the nearest upstream neighbor too
	if NOT keyword_set(nochug) then begin &$
		for i=0,numblob do begin &$
			locs=where(data(41,*) EQ i) &$
			if locs(0) NE -1 then locs=locs(sort(data(2,locs))) &$
			if length(locs) GT 1 then begin &$
				for j=0,length(locs)-1 do begin &$
					mov1=reform(mov(*,*,data(2,locs(j)))) &$
					strtx=where(mov1 EQ data(3,locs(j))) mod xdim &$
					strty=floor(where(mov1 EQ data(3,locs(j)))/xdim) &$
					if (strtx(0) GE 0) AND (strty(0) GE 0) then begin &$
						data(81,i)=strtx(0) &$
						data(82,i)=strty(0) &$
						srch=search2d(mov1,strtx(0),strty(0),data(3,locs(j)),data(3,locs(j))) &$
						cent=[xobs(data(70,locs(j))),yobs(data(70,locs(j)))] &$
						data(77,locs(j))=max(strtx) &$
						if data(71,locs(j)) GT -10 then begin &$
							slop=tan(data(71,locs(j))) &$
							theta2=atan((cent(1)-data(1,locs(j)))/(cent(0)-data(0,locs(j))))-data(71,locs(j)) &$
							off=sqrt((cent(1)-data(1,locs(j)))^2+(cent(0)-data(0,locs(j)))^2)*sin(theta2) &$
							data(83,locs(j))=off &$
							blw=where(cent(1)+(strtx-cent(0))*slop LT strty,complement=abv) &$
							data(75,locs(j))=1.0-abs((data(4,locs(j))-2.0*length(blw))/data(4,locs(j))) &$
;for a line of the form y+bx+c=0, and the point (x1,y1),
;the minimum distance from that line to the point is
;abs(bx1+y1+c)/sqrt(1+b^2)
							coeffc=slop*cent(0)-cent(1) &$
							coeffb=-slop &$
							Iarea=abs(coeffb*strtx+strty+coeffc)/sqrt(1+coeffb^2) &$
							Iarea1=total(Iarea(blw)) &$
							Iarea2=total(Iarea(abv)) &$
							data(87,locs(j))=1-abs(Iarea1-Iarea2)/(Iarea1+Iarea2) &$
							if (data(70,locs(j)) GE -1) AND (data(33,locs(j)) GT 0) then begin &$
								locs2=where(data(70,locs) EQ data(70,locs(j)) AND (data(33,locs) GT 0)) &$
								locs3=where(data(75,locs2) GE 0) &$
								data(76,locs(j))=data(75,locs(locs2(locs3(0)))) &$
								data(84,locs(j))=data(83,locs(locs2(locs3(0)))) &$
							endif &$
						endif &$
						if (data(70,locs(j)) GE -1) AND (data(33,locs(j)) GT 0) then begin &$
							locs2=where(data(70,locs) EQ data(70,locs(j)) AND (data(33,locs) GT 0)) &$
							locs3=where(data(9,locs2) GE 0) &$
							data(80,locs(j))=data(9,locs(locs2(locs3(0)))) &$
						endif &$
					endif &$
				endfor &$
				if length(where(data(75,locs) GT -1)) GT 7 then begin &$
					data(78,locs)=smooth(data(75,locs),7,/edge_truncate,/nan) &$
					data(88,locs)=smooth(data(87,locs),7,/edge_truncate,/nan) &$
					for j=0,length(xobs)-1 do begin &$
						tmploc=where((data(34,locs) EQ j) OR (data(35,locs) EQ j) OR (data(36,locs) EQ j) OR (data(37,locs) EQ j) OR (data(38,locs) EQ j) OR (data(39,locs) EQ j)) &$
						data(79,locs(tmploc))=data(78,locs(tmploc(0))) &$
						data(89,locs(tmploc))=data(88,locs(tmploc(0))) &$
					endfor &$
;					for j=0,length(locs)-1 do begin &$
;						tmploc=where(data())
;						if data(76,locs(j)) GT -1 then begin &$
;							tmploc=where(data(75,locs) EQ data(76,locs(j))) &$
;							data(79,locs(j))=data(78,locs(tmploc(0))) &$
;						endif &$
;					endfor &$
				endif &$
			endif &$
		endfor &$
		for t=0,tdim-1 do begin &$
;first define the perimeter blobs
			section=mov(*,*,t) &$
			indxcount=histogram(section,binsize=1,locations=indx) &$
			blobindx=round(indx(where((indxcount GT 0) AND (indx GE 0)))) &$
			maxp=1 &$
			c=-1 &$
			mastp=[] &$
			for i=0,length(blobindx)-1 do begin &$
				region=where(section EQ blobindx(i)) &$
				xr=region mod xdim &$
				yr=floor(region/xdim) &$
				if (max(yr) LT ydim-1) AND (max(xr) LT xdim-1) AND (min(yr) GT 1) AND (min(xr) GT 1) then begin &$
					attack=dmedgewalk(xr,yr,xdim,ydim,plist=plist,ls=1) &$
					if length(transpose(plist)) GT 1 then begin &$
						c=c+1 &$
						if length(transpose(plist)) GE maxp then begin &$
							mastp1=bytarr(2,length(transpose(plist)),c+1) &$
							mastp1(*,*,c)=plist &$ &$
							if c NE 0 then begin &$
								for k=0,c-1 do begin &$
									mastp1(*,[0:maxp-1],k)=mastp(*,*,k) &$
								endfor &$
							endif &$
							maxp=length(transpose(plist)) &$
								mastp=mastp1 &$
						endif else begin &$
							plist1=[[plist],[bytarr(2,maxp-length(transpose(plist)))]] &$
							mastp1=bytarr(2,maxp,c+1) &$
								mastp1(*,*,[0:c-1])=mastp &$
							mastp1(*,*,c)=plist1 &$
							mastp=mastp1 &$
						endelse &$
					endif &$
				endif &$
			endfor &$
;			if keyword_set(antoine) then print,'finish step 1 for time'+string(t) &$
;now define nearest neighbor
			szp=size(mastp) &$
			if isa(nnn) EQ 0 then begin &$
				for i=0,length(blobindx)-1 do begin &$
					defperim=where((mastp(0,*,i) EQ 0) AND (mastp(1,*,i) EQ 0)) &$
					if defperim(0) EQ -1 then begin &$
						lenp=maxp &$
					endif else begin &$
						lenp=defperim(0) &$
					endelse &$
					minp=replicate(100d,2,1) &$
					for j=0,lenp-1 do begin &$
						for k=0,length(blobindx)-1 do begin &$
							if (k NE i) and (isa(mastp) EQ 1) then begin &$
								dstp=sqrt((mastp(0,j,i)-mastp(0,*,k))^2+(mastp(1,j,i)-mastp(1,*,k))^2) &$
								if minp(0) GT min(dstp) then minp=[min(dstp),blobindx(k)] &$
							endif else begin &$
								dstp=-1 &$
								minp=[-1,-1] &$
							endelse &$
						endfor &$
					endfor &$
					place1=where((data(2,*) EQ t) AND (data(3,*) EQ blobindx(i))) &$
					place2=where((data(2,*) EQ t) AND (data(3,*) EQ minp(1))) &$
					data(86,place1)=minp(0) &$
					data(85,place1)=data(41,place2(0)) &$
				endfor &$
			endif &$
		endfor &$
	endif
;begin breakup and combination detection. index 6 is 0 if clean, and 1 if disturbed
	combin=dmcombtif(brth,deth,tdim,xdim,obsr,data)
;write in code for if combin=empty set
	if length(combin) EQ 0 then begin &$
		combin=replicate(!values.f_nan,8,1) &$
		tempb=-1 &$
		tempc=-1 &$
	endif else begin &$
		tempb=where(combin(6,*) EQ 5) &$
		tempc=where(combin(6,*) EQ 4) &$
	endelse
	ddrops=[]
	if tempb(0) NE -1 then begin &$
		brek=combin(*,where(combin(6,*) EQ 5)) &$
		for i=0,length(tempb)-1 do begin &$
			loc=where((data(41,*) EQ brek(0,i)) AND (data(2,*) EQ brek(3,i))) &$
			firsttouch1=data(34,where(data(41,*) EQ brek(0,i))) &$
			firsttouch2=where((data(41,*) EQ brek(0,i)) AND (data(34,*) EQ firsttouch1(-1))) &$
			firsttouch=firsttouch2(0) &$
			if length(loc) EQ 1 then begin &$
				if keyword_set(ddropcheck) EQ 0 then begin &$
					brek(6,i)=0 &$
				endif else begin &$
					if data(13,loc) EQ 0 then brek(6,i)=0 &$
					if data(13,loc) EQ 1 then begin &$
						brek(6,i)=1 &$
					endif &$
					if data(12,loc) EQ 1 then begin &$
						brek(6,i)=1 &$
					endif &$
					if data(33,loc) NE 1 then begin &$
						brek(6,i)=1 &$
					endif &$
					if firsttouch2(0) NE -1 then begin &$
						if data(13,firsttouch) EQ 1 then begin &$
							brek(6,i)=1 &$
							endif &$
						if data(12,firsttouch) EQ 1 then begin &$
							brek(6,i)=1 &$
						endif &$
						if brth(1,where(brth(0,*) EQ brek(0,i))) EQ data(40,firsttouch) then begin &$
							brek(6,i)=1 &$
						endif &$
					endif &$
				endelse &$
			endif else begin &$
				if keyword_set(antoine) then print,'breakup '+string(i)+' has '+string(length(loc))+' many valid points.' &$
			endelse &$
		endfor &$
		brek1=replicate(!values.F_nan,13,length(transpose(brek))) &$
		brek1[0:7,*]=brek &$
		brek=brek1 &$
		for i=0,length(transpose(brek))-1 do begin &$
			tempp=where((data(41,*) EQ brek(0,i)) AND (data(2,*) EQ brek(3,i))) &$
			if data(12,tempp(0)) EQ 0 then begin &$
				brek(7,i)=data(76,tempp(0)) &$
				brek(8,i)=data(4,tempp(0)) &$
				brek(9,i)=data(74,tempp(0)) &$
				brek(10,i)=data(7,tempp(0)) &$
				brek(11,i)=data(85,tempp(0)) &$
				brek(12,i)=length(where((data(41,*) EQ brek(0,i)) AND ((data(34,*) EQ data(70,*)) $
OR (data(35,*) EQ data(70,*)) OR (data(36,*) EQ data(70,*)) OR (data(37,*) EQ data(70,*)) OR (data(38,*) EQ data(70,*)) $
OR (data(39,*) EQ data(70,*))))) &$
			endif &$
			if brek(6,i) EQ 0 then begin &$
				ddrop1=data(4,where((data(41,*) EQ brek(1,i)) AND (data(2,*) EQ brek(3,i)+1) AND (brek(6,i) EQ 0))) &$
				ddrop2=data(4,where((data(41,*) EQ brek(2,i)) AND (data(2,*) EQ brek(3,i)+1) AND (brek(6,i) EQ 0))) &$
				ddrop0=data(76,where((data(41,*) EQ brek(0,i)) AND (data(2,*) EQ brek(3,i)) AND (brek(6,i) EQ 0))) &$
				ddrop3=data(79,where((data(41,*) EQ brek(0,i)) AND (data(2,*) EQ brek(3,i)) AND (brek(6,i) EQ 0))) &$
				ddrop4=data(70,where((data(41,*) EQ brek(0,i)) AND (data(2,*) EQ brek(3,i)) AND (brek(6,i) EQ 0))) &$
				ddropr=(ddrop1+0.0)/ddrop2 &$
				pddr=ddrop3/(2-ddrop3) &$
				if ddropr GT 1 then ddropr=1/ddropr &$
				alpha=1-2*ddropr/(ddrop0+ddrop0*ddropr) &$
				ddrops=[[ddrops],[ddropr,ddrop0,i, brek(8,i), brek(9,i), brek(10,i),ddrop3,ddrop4,brek(12,i),alpha,min([ddrop1,ddrop2]),max([ddrop1,ddrop2]),brek(11,i)]] &$
;DDR, symmetry, index, area, velocity, shape parameter, smooth symmetry, index of closest obstacle, contact time,smaller drop, larger drop,distance to nearest neighbor
			endif &$
		endfor &$
;		brek=brek(*,where(brek(6,*) EQ 0)) &$
	endif else begin &$
		brek=replicate(!values.f_nan,12,1) &$
	endelse
	if length(ddrops) EQ 0 then ddrops=replicate(!values.f_nan,8,1)
	if tempc(0) NE -1 then begin &$
		combin=combin(*,where(combin(6,*) EQ 4)) &$
		combin=combin([0:5],*) &$
	endif else begin &$
		combin=replicate(!values.f_nan,7,1) &$
	endelse
	for i=0,numblob do begin &$
		if tempc(0) NE -1 then begin &$
			tmpcc=where(((combin(1,*) EQ i) OR (combin(2,*) EQ i))) &$
			if tmpcc(0) GT -1 then data(45,where(data(41,*) EQ i))=3 &$
		endif &$
		if tempb(0) NE -1 then begin &$
			tmpbb=where(brek(0,*) EQ i) &$
			if tmpbb(0) GT -1 then data(45,where(data(41,*) EQ i))=4 &$
		endif &$
	endfor
;create the ancestor network
	ancestry=replicate(-10,33,length(transpose(combin)))
	if combin(0) GT -1 then begin &$
		for i=0,length(transpose(combin))-1 do begin &$
			anc=[combin(1,i),combin(2,i)] &$
			satis=0 &$
			while satis EQ 0 do begin &$
				anc1=[] &$
				for j=0,length(anc)-1 do begin &$
					tmpa=where(combin(0,*) EQ anc(j)) &$
					if tmpa(0) NE -1 then anc1=[anc1,combin(1,tmpa),combin(2,tmpa)] &$
				end &$
				anc2=[] &$
				for j=0,length(anc1)-1 do begin &$
					tmpb=where(anc1(j) EQ anc) &$
					if tmpb(0) EQ -1 then anc2=[anc2,anc1(j)] &$
				end &$
				anc=[anc,anc2] &$
				if length(anc2) EQ 0 then satis=1 &$
			end &$
			if length(anc) GT 16 then begin &$
				if keyword_set(antoine) then print,'ancestor history more than 32' &$
			endif else begin &$
				ancestry(0:length(anc),i)=[combin(0,i),anc] &$
			endelse &$
		endfor &$
		if brek(0) GT -1 then begin &$
			for i=0,length(transpose(brek))-1 do begin &$
				for j=0,1 do begin &$
					tmpb=where(ancestry(0,*) EQ brek(1+j,i)) &$
					tmpc=where(ancestry(0,*) EQ brek(0,i)) &$
					tmpd=where(ancestry EQ brek(1+j,i)) &$
					if tmpb(0) EQ -1 then begin &$
						if tmpc(0) EQ -1 then begin &$
							ancestry=[[ancestry],[brek(1+j,i),brek(0,i),replicate(-10,31,1)]] &$
						endif else begin &$
							tmpc1=ancestry(where(ancestry(*,tmpc) GT -10),tmpc(0)) &$
							if (length(tmpc1) LE 31) then begin &$
								ancestry=[[ancestry],[brek(1+j,i),tmpc1,replicate(-10,32-length(tmpc1),1)]] &$
							endif else begin &$
								ancestry=[[ancestry],[brek(1+j,i),tmpc1([0:31])]] &$
							endelse &$
						endelse &$
					endif else begin &$
						tmpb1=where(ancestry(*,tmpb) EQ -1) &$
						ancestry(tmpb1(0),tmpb(0))=brek(0,i) &$
					endelse &$
					if tmpd(0) NE -1 then begin &$
						pullfrom=where(ancestry(0,*) EQ brek(1+j,i)) &$
						pullfrom1=ancestry([1:length(where(ancestry(*,pullfrom) GT -10))-1],pullfrom(0)) &$
						for k=0,length(tmpd)-1 do begin &$
							xd=tmpd(k) mod 33 &$
							yd=tmpd(k)/33 &$
							if xd NE 0 then begin &$
								plc=where(ancestry(*,yd) EQ -10) &$
								ancestry(plc([0:length(pullfrom1)-1]),yd)=pullfrom1 &$
							endif &$
						endfor &$
					endif &$
				endfor &$
			endfor &$
		endif &$
		ancestry=fix(ancestry) &$
	endif
;blob connection network
;conec=[index 1, index 2, time start, time end, what happens to 1, what happens to 2, max interface length,total memory time,edgeflag or not at final time]
;positive index ->coalesces into that droplet. -1=nothing happens, -2=wanders offscreen, -3=breaks, -4=video ends, -5=sub-two frame connection
	if not keyword_set(nochug) then begin &$
		if combin(0) GT -1 then begin &$
			conec=dmconecnetwork1(data,numblob,combin,brek,tdim,xdim,ancestry) &$
		endif else begin &$
			conec=replicate(!values.f_nan,8,1) &$
		endelse &$
		contacttim=[] &$
		for i=0,numblob do begin &$
			locs=where(data(41,*) EQ i) &$
			if locs(0) NE -1 then locs=locs(sort(data(2,locs))) &$
			for j=0,length(xobs)-1 do begin &$
				tmp1=where((data(34,locs) EQ j) OR (data(35,locs) EQ j) OR (data(36,locs) EQ j) OR (data(37,locs) EQ j) OR (data(38,locs) EQ j) OR (data(39,locs) EQ j)) &$
				if tmp1(0) NE -1 then begin &$
					contacttim=[[contacttim],[i,j,data(2,locs(tmp1(-1)))-data(2,locs(tmp1(0)))+1,length(tmp1)]] &$
					if contacttim(2,-1) NE contacttim(3,-1) then print,'droplet '+string(i)+' reattaches to obstacle '+string(j) &$
				endif &$
			endfor &$
		endfor &$
		sr=FILE_SEARCH(direc+fname+'ptedge',count=count) &$
;0 -> x-com position
;1 -> y-com position
;2 -> time
;3 -> index in frame of drop 1
;4 -> index in frame of drop 2
;5 -> mass
;6 -> slope
;7 -> area
;8 -> edgeflag
;9-> time (again)
;10-> index through time
;11-> index in time of drop 1
;12-> index in time of drop 2
;13-> does the droplet coalesce?
;14-> what index breakup was most recent
;15-> time since breakup
;16-> mass of drop 1
;17-> mass of drop 2
;18-> velocity of drop 1
;19-> velocity angle of drop 1
;20-> velocity of drop 2
;21-> velocity angle of drop 2
;22-> mass change rate
;23-> mass angle
;24-> mass angle velcity
		if (count NE 0) AND keyword_set(interface) then begin &$
			ptedge=read_gdf(direc+fname+'ptedge') &$
			if (isa(newptedge) EQ 1) then begin &$
				if ((ptedge(0) NE -1) EQ 1) then begin &$
					sizeptdata=size(ptedge) &$
					extredge=fltarr(sizeptdata(1),tdim) &$
					extredge(0,*)=1 &$
					extredge(1,*)=1 &$
					extredge(5,*)=-1 &$
					extredge(9,*)=[0:tdim-1] &$
					ptedge=[[ptedge],[extredge]] &$
					ptedge=ptedge(*,sort(ptedge(9,*))) &$
					dataedge=track(ptedge,15,mem=1,good=3,/quiet) &$
					dataedge=dataedge(*,where(dataedge(5,*) NE -1,dataedgel)) &$
					if dataedgel GT 0 then begin &$
						edgesz=size(dataedge) &$
						datae=fltarr(edgesz(1)+15,edgesz(2)) &$
						datae([0:(edgesz(1)-1)],*)=dataedge &$
						datae(8,where(datae(2,*) EQ 0))=1 &$
						datae(8,where(datae(2,*) EQ tdim))=1 &$
						for t=0,tdim-1 do begin &$
							fxme=where(datae(2,*) EQ t,fxmel) &$
							if fxmel GT 0 then begin &$
								for i=0,fxmel-1 do begin &$
									bleh=datae(*,fxme(i)) &$
									b1=where((data(3,*) EQ bleh(3)) and data(2,*) EQ bleh(2),b1l) &$
									if b1l EQ 1 then begin &$
										datae(11,fxme(i))=data(41,b1) &$
										datae(16,fxme(i))=data(4,b1) &$
										datae(18,fxme(i))=data(44,b1) &$
										datae(19,fxme(i))=data(71,b1) &$
									endif &$
									b2=where((data(3,*) EQ bleh(4)) and data(2,*) EQ bleh(2),b2l) &$
									if b2l EQ 1 then begin &$
										datae(12,fxme(i))=data(41,b2) &$
										datae(17,fxme(i))=data(4,b2) &$
										datae(20,fxme(i))=data(44,b2) &$
										datae(21,fxme(i))=data(71,b2) &$
									endif &$
									fxmetmp=where(((data(41,*) EQ datae(11,fxme(i))) OR (data(41,*) EQ datae(12,fxme(i)))) AND (data(2,*) EQ datae(2,fxme(i)))) &$
									if length(fxmetmp) GT 1 then begin &$
										for j=0,length(fxmetmp)-1 do begin &$
											if data(12,fxmetmp(j)) EQ 1 then datae(8,fxme(i))=1 &$
										endfor &$
									endif else begin &$
										if length(fxmetemp) EQ 1 then begin &$
											if data(12,fxmetmp) EQ 1 then datae(8,fxme(i))=1 &$
										endif &$
									endelse &$
									aatmp=( (max(data(2,where(data(41,*) EQ datae(11,fxme(i))))) EQ t) OR (max(data(2,where(data(41,*) EQ datae(12,fxme(i))))) EQ t) OR ((data(12,where((data(41,*) EQ datae(11,fxme(i))) AND (data(2,*) EQ t+1)))) EQ 1) OR ((data(12,where((data(41,*) EQ datae(12,fxme(i))) AND (data(2,*) EQ t+1)))) EQ 1) ) &$
									if ((aatmp(0) EQ 1) and (length(aatmp) EQ 1)) then begin &$
										fxmetmp=where(((data(41,*) EQ datae(11,fxme(i))) OR (data(41,*) EQ datae(12,fxme(i)))) AND (data(2,*) EQ datae(2,fxme(i))+1)) &$
										if length(fxmetmp) GE 1 then datae(8,fxme(i))=1 &$
									endif &$
									datae(23,fxme(i))=atan(datae(6,fxme(i))) &$
									if t NE 0 then begin &$
										fxmet=where(((datae(2,*) EQ t-1) AND (datae(10,*) EQ datae(10,fxme(i)))),fxmelt) &$
										if fxmelt EQ 1 then begin &$
											datae(22,fxme(i))=datae(5,fxme(i))-datae(5,fxmet) &$
											datae(24,fxme(i))=datae(23,fxme(i))-datae(23,fxmet) &$
										endif &$
									endif &$
								endfor &$
							endif &$
						endfor &$
						numedge=max(datae(10,*)) &$
						if length(transpose(combin)) GT 1 then begin &$
							for i=0,length(transpose(combin))-1 do begin &$
								tmp1=datae(10,where( ((datae(11,*) EQ combin(1,i)) AND (datae(12,*) EQ combin(2,i)) AND (datae(2,*) EQ combin(3,i)-1)) OR ((datae(11,*) EQ combin(2,i)) AND (datae(12,*) EQ combin(1,i)) AND ((datae(2,*) EQ combin(3,i)-1))) )) &$
								tmp1=tmp1(sort(tmp1)) &$
								tmp1=tmp1(uniq(tmp1)) &$
								tmp1l=length([tmp1]) &$
								if tmp1l GT 0 then begin &$
									if tmp1l EQ 1 then begin &$
										datae(13,where(datae(10,*) EQ tmp1))=1 &$
									endif else begin &$
										if keyword_set(antoine) then print,'strange error with'+string(i) &$
									endelse &$
								endif &$
							endfor &$
						endif &$
;						plot,[0,1],[0,1],/NODATA,xrange=[0,600],yrange=[-0.2,1.2],charsize=3 &$
						bglst=[] &$
						for i=0,numedge do begin &$
							locs=where(datae(10,*) EQ i,locsl) &$
							velf=datae(22,locs(-1)) &$
							if length(locs) GT 2 then begin &$
								velfs=datae(22,locs) &$
								if velf(-1) GT mean(abs(velfs([0:length(locs)-2])))+2*stdev(abs(velfs([0:length(locs)-2]))) then begin &$
									velf=velfs(-2) &$
									locs=locs(sort(datae(2,locs))) &$
									locs=locs[0:length(locs)-2] &$
									locsl=locsl-1 &$
								endif &$
							endif &$
							if locsl GE 2 then begin &$
								locs=locs(sort(datae(2,locs))) &$
								maxedge=max(datae(5,locs),met) &$
								maxedge=maxedge &$
								ttime=length(locs) &$
								edgefl=where(datae(8,locs) EQ 1,edgefll) &$
								for j=0,length(transpose(brek))-1 do begin &$
									chckhist=where(((datae(11,locs) EQ brek(0,j)) OR (datae(12,locs) EQ brek(0,j))) AND (datae(2,locs) EQ brek(3,j)),chckhistl) &$
									if chckhistl GT 0 then begin &$
										if chckhistl EQ 1 then begin &$
											datae(12,locs)=brek(0,j) &$
											datae(14,locs)=brek(3,j)-datae(2,locs) &$
										endif else begin &$
											if keyword_set(antoine) then print,'multiple breakup events for i,j='+string([i,j]) &$
										endelse &$
									endif &$
								endfor &$
								edgebr=where(datae(15,locs) NE 0,edgebrl) &$
								if ((edgebrl EQ 0) and (edgefll EQ 0)) then begin &$
									massr=datae(16,locs(0))/datae(17,locs(0)) &$
									if massr GT 1 then massr=1/massr &$
									comv=[0] &$
									comvp=[0] &$
								totl=[0] &$
									totn=[0] &$
									for k=0,locsl-1 do begin &$
										lv1=where((data(2,*) EQ datae(2,locs(k))) AND (data(41,*) EQ datae(11,locs(k)) AND (data(44,*) GT 0)),lvl1) &$
										lv2=where((data(2,*) EQ datae(2,locs(k))) AND (data(41,*) EQ datae(12,locs(k)) AND (data(44,*) GT 0)),lvl2) &$
										if (lvl1 EQ 1) and (lvl2 EQ 1) then begin &$
											avv=[(data(4,lv1)*data(42,lv1)+data(4,lv2)*data(42,lv2))/(data(4,lv1)+data(4,lv2)),$
												(data(4,lv1)*data(43,lv1)+data(4,lv2)*data(43,lv2))/(data(4,lv1)+data(4,lv2))] &$
											v2x=data(42,lv1)-avv(0) &$
											v2y=data(43,lv1)-avv(1) &$
											v1x=data(42,lv2)-avv(0) &$
											v1y=data(43,lv2)-avv(1) &$
											cv1x=data(0,lv1)-data(0,lv2) &$
											cv1y=data(1,lv1)-data(1,lv2) &$
											cv2x=data(0,lv2)-data(0,lv1) &$
											cv2y=data(1,lv2)-data(1,lv1) &$
											project1=(v1x*cv1x+v1y*cv1y)/sqrt(cv1y^2+cv1x^2) &$
											project2=(v2x*cv2x+v2y*cv2y)/sqrt(cv2y^2+cv2x^2) &$
											pproject1=abs(v1x*cv1y-v1y*cv1x)/sqrt(cv1x^2+cv1y^2) &$
											pproject2=abs(v2x*cv2y-v2y*cv2x)/sqrt(cv2x^2+cv2y^2) &$
											comv=[comv,comv(-1)+project1+project2] &$
											comvp=[comvp,comvp(-1)+pproject1+pproject2] &$
											totl=[totl,datae(5,locs(k))+totl(-1)] &$
											totn=[totn,datae(5,locs(k))/data(5,lv1)+datae(5,locs(k))/data(5,lv2)+totn(-1)] &$
										endif &$
									endfor &$
									rij=1/((1/sqrt(datae(17,locs(0))/!PI))+(1/sqrt(datae(16,locs(0))/!PI))) &$
									if stdev(datae(22,locs)) LT 10 then begin &$
										bglst=[[bglst],[datae(13,locs(0)),maxedge,ttime,massr,datae(16,locs(0)),datae(17,locs(0)),$
											datae(2,locs(-1))-datae(2,locs(met)),1/rij+(1/rij)^2,datae(18,locs(-1)),datae(19,locs(-1)),$
											datae(20,locs(-1)),datae(21,locs(-1)),comv(-1),comvp(-1),totl(-1),totn(-1),velf,$
											stdev(datae(22,locs)),datae(10,locs(0))]] &$
									endif &$
								endif &$
							endif &$
						endfor &$
					endif &$
				endif else begin &$
					bglst=replicate(-1d,6,1) &$
				endelse &$
			endif &$
		endif &$
	endif
;dmmov3,mov,data(0,*),data(1,*),data(2,*),data(41,*),brek,combin
;dmmov5,mov,data(0,*),data(1,*),data(2,*),data(28,*),data(29,*),data(30,*),data(31,*),data(41,*)
	return,data
end

