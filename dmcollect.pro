; dmcollect
; written by David Meer 6/2026
; collects static information for droplets
;
function dmcollect,fname,directory,first=first,chngobs=chngobs
;fname='10_24_70um_21'
;directory='10_24_24'
direc='/data/david/drops/'+directory+'/'
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
	for t=0,tdim-1 do begin
		section=vid(*,*,t)
		indxcount=histogram(section,binsize=1,locations=indx)
		blobindx=round(indx(where((indxcount GT 0) AND (indx GE 0))))
		print,string((t+0.0)/tdim)+'collect step'
		for i=0,length(blobindx)-1 do begin
			region=where(section EQ blobindx(i))
			xr=region mod xdim
			yr=floor(region/xdim)
			mass=length(region)
			xcom=mean(xr)
			ycom=mean(yr)
			perim=0
			obsdist=min((xcom-xobs)^2+(ycom-yobs)^2)
			obsdistloc=where(((xcom-xobs)^2+(ycom-yobs)^2) EQ obsdist)
			obsdistfront=where(((xcom-xobs)^2+(ycom-yobs)^2) EQ min(obsdist(where(xobs GT xcom))))
			xobsnear=xobs(obsdistloc)
			yobsnear=yobs(obsdistloc)
			xobsnearf=xobs(obsdistfront)
			yobsnearf=yobs(obsdistfront)
			touchobs=0
			nobs=0
			obstrack=replicate(!Values.F_NAN,6,1)
			for j=0,length(xobs)-1 do begin &$
				tempobs=where(sqrt((xr-xobs(j))^2+(yr-yobs(j))^2) LE obsr+5.0) &$
				if tempobs(0) GT -1 then begin &$
					touchobs=1 &$
					nobs=nobs+1 &$
					if nobs LT 6 then obstrack(nobs-1)=j &$
				endif &$
			endfor
			xsymm=length(where(yobsnearf(0) LT yr))/(mass+0.0)
			ysymm=length(where(xobsnearf(0) LT xr))/(mass+0.0)
			xsymm=length(where(yobsnear(0) LT yr))/(mass+0.0)
			ysymm=length(where(xobsnear(0) LT xr))/(mass+0.0)
			edgeflag=0
			if max(xr) GE xdim-2 then edgeflag=1
			if min(xr) LE 2 then edgeflag=1
			if max(yr) GE ydim-2 then edgeflag=1
			if min(yr) LE 2 then edgeflag=1
			neighlst=[]
			neighs=replicate(!Values.F_NAN,6,1)
			nneigh=0
			calA=!Values.F_NAN
;NOT USED		perimc=!Values.F_NAN
			if (mass GT 20) AND (edgeflag EQ 0) then begin &$
				attack=dmedgewalk3(xr,yr,xdim,ydim) &$
			endif else begin &$
				attack=[!Values.F_NAN,!Values.F_NAN,!Values.F_NAN,!Values.F_NAN,!Values.F_NAN,!Values.F_NAN] &$
			endelse &$
			symslop=!Values.F_NAN
			if attack(1) GT -100 then begin &$
				obsx=xobs(where(((xcom-xobs)^2+(ycom-yobs)^2) EQ obsdist)) &$
				obsy=yobs(where(((xcom-xobs)^2+(ycom-yobs)^2) EQ obsdist)) &$
				slop=tan(attack(1)) &$
				symup=length(where((obsy(0)+(xr-obsx(0))*slop) GT yr)) &$
				symdn=length(where((obsy(0)+(xr-obsx(0))*slop) LE yr)) &$
				symslop=1.0-abs((symup-symdn+0.0)/mass) &$
			end
			if edgeflag EQ 0 then begin &$
				ylist1=yr(sort(yr)) &$
				ylist=ylist1(uniq(ylist1)) &$
				neckls=[] &$
				for j=0,length(ylist)-1 do begin &$
					xsuby=xr(where(yr EQ ylist(j))) &$
					neckls=[neckls,max(xsuby)-min(xsuby)] &$
				endfor &$
				neckmax=[] &$
				for j=1,length(ylist)-2 do begin &$
					if (neckls(j) GT neckls(j+1)) AND (neckls(j) GT neckls(j+1)) then neckmax=[neckmax,j] &$
				endfor &$
				if length(neckmax) GT 1 then begin &$
					neckl=min(neckls([neckmax(1):neckmax(-1)])) &$
				endif else begin &$
					neckl=!Values.F_NAN &$
				endelse &$
				for j=0,mass-1 do begin &$
					testedge=[section(xr(j)+1,yr(j)+1),section(xr(j),yr(j)+1),section(xr(j)-1,yr(j)+1)$
,section(xr(j)+1,yr(j)),section(xr(j)-1,yr(j)),section(xr(j)+1,yr(j)-1),section(xr(j),yr(j)-1),section(xr(j)-1,yr(j)-1)] &$
					perimchk=where(testedge NE blobindx(i)) &$
					if perimchk(0) NE -1 then begin &$
						perim=perim+1 &$
						neighchk=testedge(perimchk) &$
						boundchk=where(neighchk NE -1) &$
						if boundchk(0) GE 0 then begin &$
							if length(bounchk) EQ 1 then begin &$
								neighlst=[neighlst,neighchk(boundchk)] &$
							endif else begin &$
								for k=0,length(perimchk)-1 do begin &$
									neighlst=[neighlst,neighchk(boundchk)] &$
								endfor &$
							endelse &$
						endif &$
					endif &$
				endfor &$
			endif else begin &$
				neckl=!Values.F_NAN &$
				perim=!Values.F_NAN &$
			endelse
			perimc=(!pi/4)*(perim+8)
			calA=perimc^2/(mass*4*!pi)
			if length([neighlst]) GE 2 then begin &$
				tmplst1=neighlst(sort(neighlst)) &$
				tmplst2=tmplst1(uniq(tmplst1)) &$
				neighlst=tmplst2 &$
			endif
			nneigh=length([neighlst])
			if nneigh GT 6 then begin &$
				print,'more than 6 neighbors on blob'+string(blobindx(i))+' at time'+string(t) &$
				print,'there are '+string(nneigh)+' neighbors' &$
				masslist=replicate(!VALUES.F_NAN,length(neighlst),1) &$
				for k=0,length(neighlst)-1 do begin &$
					if neighlst(k) LT blobindx(i) then begin &$
						masslist(k)=data(4,where((data(3,*) EQ neighlst(k)) AND (data(2,*) EQ t))) &$
					endif else begin &$
						masslist(k)=length(where(section EQ neighlst(k))) &$
					endelse &$
				endfor &$
				neighmass=neighlst(sort(masslist)) &$
				neighlst=neighmass([0:5]) &$
				nneigh=6 &$
			endif
			if nneigh GT 0 then begin
				if nneigh GT 1 then begin
					neighs([0:(nneigh-1)])=neighlst
				endif else begin
					neighs(0)=neighlst
				endelse
			endif
			touching=replicate(!values.F_NAN,6)
			if nneigh gt 0 then begin &$
				excess=[[xr+1,xr,xr+1,xr-1,xr+1,xr,xr-1,xr-1],[yr+1,yr+1,yr,yr+1,yr-1,yr-1,yr,yr-1]] &$
				thng=1000.0*excess(*,0)+excess(*,1) &$
				thng1=thng(sort(thng)) &$
				thng2=thng1(uniq(thng1)) &$
				chckme=[[floor(thng2/1000.0)],[thng2-floor(thng2/1000.0)*1000.0]] &$
				if nneigh LE 6 then begin &$
                                        for j=0,nneigh-1 do begin &$
                                                tmp1=where(section(excess(*,0),excess(*,1)) EQ neighlst(j)) &$
                                                if tmp1(0) GT -1 then touching(j)=length(tmp1) &$
                                        endfor &$
                                endif else begin &$
                                        for j=0,5 do begin &$
                                                tmp1=where(section(excess(*,0),excess(*,1)) EQ neighlst(j)) &$
                                                if tmp1(0) GT -1 then touching(j)=length(tmp1) &$
                                        endfor &$
                                endelse &$
			endif
			blobdata=[xcom,ycom,t,blobindx(i),mass,perim,symslop,calA,sqrt(obsdist),xsymm,ysymm,neckl,edgeflag,nneigh,neighs,touching,$
				attack(0),attack(1),attack(2),attack(3),attack(4),attack(5),touchobs,nobs,obstrack(0),obstrack(1),$
				obstrack(2),obstrack(3),obstrack(4),obstrack(5),tim(t),obsdistfront(0)]
			data=[[data],[blobdata]]
		endfor
	endfor
	write_gdf,data,direc+fname+'ptcollect'
	return,data
end

