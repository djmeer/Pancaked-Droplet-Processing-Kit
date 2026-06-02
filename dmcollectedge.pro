; dmcollectedge
; written by David Meer 6/2026
; collects static information for droplet-droplet interfaces
;
function dmcollectedge,fname,directory,first=first,chngobs=chngobs,antoine=antoine,coal=coal
;fname='10_24_70um_21'
;directory='10_24_24'
direc='/data/david/drops/'+directory+'/'
        if keyword_set(antoine) then direc='/data/antoine/'
	 if keyword_set(coal) then direc='/data/david/coal/'+directory+'/'
	if (keyword_set(first) EQ 0) then begin
		vid=read_gdf(direc+'pt'+fname+'gdf')
	endif else begin
		sr=FILE_SEARCH(direc+'setup'+fname+'gdf*',count=count)
		if count GT 0 then begin
			vid=dmptpink4(fname,directory)
		endif else begin
			vid=dmptpink4(fname,directory,firsttime=1)
		endelse
	endelse
	vids=size(vid)
	xdim=vids(1)
	ydim=vids(2)
	tdim=vids(3)
	tim=[0:tdim-1]
	sr=FILE_SEARCH(direc+fname+'obs',count=countobs)
	if (countobs GT 0) then begin
		obs=float(readtext(direc+fname+'obs'))
		xobs1=obs(*,0)
;		yobs1=max(obs(*,1))-obs(*,1)
		yobs1=obs(*,1)
		stup=read_gdf(direc+fname+'obssetup')
		scl=stup(0)
		xshift=stup(1)
		yshift=stup(2)
		ang=stup(3)
		obsr=stup(4)
		pix2um=stup(5)
		fps=stup(6)
		xobs2=(xobs1*scl+xshift)*cos(ang*!pi/180.0)-(yobs1*scl+yshift)*sin(ang*!pi/180.0)
		yobs2=(xobs1*scl+xshift)*sin(ang*!pi/180.0)+(yobs1*scl+yshift)*cos(ang*!pi/180.0)
		xobs=xobs2(where((xobs2 GE 0) and (xobs2 LT xdim) AND (yobs2 GE 0) and (yobs2 LT ydim)))
		yobs=yobs2(where((xobs2 GE 0) and (xobs2 LT xdim) AND (yobs2 GE 0) and (yobs2 LT ydim)))
	endif else begin
		sr=FILE_SEARCH(direc+fname+'cobs',count=countcobs)
		if countcobs GT 0 then begin
			stup=read_gdf(direc+fname+'cobs')
			lstup=n_elements(stup)-3
			xobs=stup([0:(lstup/2)-1])
			yobs=stup([lstup/2:(lstup-1)])
			obsr=stup(-3)
			pix2um=stup(-2)
			fps=stup(-1)
		endif else begin
			print,'obstacle error'
		endelse
	endelse
	data=[]
	for t=0,tdim-1 do begin &$
		data1=[] &$
		section=vid(*,*,t) &$
		indxcount=histogram(section,binsize=1,locations=indx) &$
		blobindx=round(indx(where((indxcount GT 0) AND (indx GE 0)))) &$
		print,string((t+0.0)/tdim)+'collect step' &$
		for i=0,length(blobindx)-1 do begin &$
			region=where(section EQ blobindx(i)) &$
			xr=region mod xdim &$
			yr=floor(region/xdim) &$
			xmax=max(xr) &$
			ymax=max(yr) &$
			xmin=min(xr) &$
			ymin=min(yr) &$
			edgeflag=0 &$
			if ((ymax GE ydim-2) OR (xmax GE xdim-2) OR (xmin LE 1) OR (ymin LE 1)) then edgeflag=1 &$
			if edgeflag EQ 0 then begin &$
				neighs=[section(xr+1,yr),section(xr-1,yr),section(xr,yr-1),section(xr,yr+1)] &$
				neighs=neighs(sort(neighs)) &$
				neighs=neighs(uniq(neighs)) &$
				neighs=neighs(where((neighs GT 0) AND (neighs NE blobindx(i)),neighl)) &$
				if neighl GT 0 then begin &$
					nk=dmedgewalk2(xr,yr,xdim,ydim,plist=edges,ls=1) &$
					perim=length(transpose(edges)) &$
					for j=0,neighl-1 do begin &$
						onedge=[] &$
						for k=0,perim-1 do begin &$
							chck=[section(edges(0,k)+1,edges(1,k)),section(edges(0,k)-1,edges(1,k)),section(edges(0,k),edges(1,k)-1),section(edges(0,k),edges(1,k)+1)] &$
							tmpchck=where(chck EQ neighs(j)) &$
							if tmpchck(0) GT -1 then onedge=[[onedge],[edges(0,k),edges(1,k)]] &$
						endfor &$
						if isa(onedge) eq 1 then begin &$
							x=transpose(onedge(0,*)) &$
							y=transpose(onedge(1,*)) &$
							xcom=mean(x) &$
							ycom=mean(y) &$
							mass=length(x) &$
							lfit=linfit(x,y) &$
							slop=lfit(1) &$
							ang=atan(slop) &$
							data1=[[data1],[xcom,ycom,t,blobindx(i),neighs(j),mass,slop,ang,edgeflag,t]] &$
						endif &$
					endfor &$
				endif &$
			endif &$
		endfor &$
		if length(data1) GT 0 then begin &$
			for i=0,length(blobindx)-1 do begin &$
				tmp1=where(((data1(3,*) EQ blobindx(i)) AND (data1(3,*) LT data1(4,*))),tmpl) &$
				if tmpl GT 0 then begin &$
					for j=0,tmpl-1 do begin &$
						tmp2=where((data1(4,*) EQ blobindx(i)) AND (data1(3,*) EQ data1(4,tmp1(j))),tmpl) &$
						if tmpl GT 0 then begin &$
							pt1=data1(*,tmp1(j)) &$
							pt2=data1(*,tmp2) &$
							ptef=0 &$
							if (pt1(8) EQ 1) OR (pt2(8) EQ 1) then ptef=1 &$
							data=[[data],[(pt1(0)+pt2(0))/2,(pt1(1)+pt2(1))/2,pt1(2),pt1(3),pt1(4),(pt1(5)+pt2(5))/2,(pt1(6)+pt2(6))/2,(pt1(7)+pt2(7))/2,ptef,pt1(9)]] &$
						endif &$
					endfor &$
				endif &$
			endfor &$
		endif &$
	endfor
	if isa(data) EQ 0 then data=replicate(-1d,9,1)
	write_gdf,data,direc+fname+'ptedge'
	return,data
end

