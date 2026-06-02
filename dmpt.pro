; dmpt
; written by David Meer 6/2026
; takes tif stacks of droplets moving in 2D and turns them into a data array set up for pretracking
; Each frame takes a value of 0 everywhere except for droplets, which take a value according to their assigned index.
; Indexes are assigned vertically, with droplets that have pixels at y=0 taking the lowest.
;
;keywords
;If keyword format is not specified, it takes the value 1 to turn on.
; firsttime - runs every subroutine instead of pulling from saved data
; chngobs - modifies the obstacle placement instead of pulling from a saved data
; colorbalance - applies an additional color balancing to the image
; invert - applies a reflection to the obstacle placements

function dmpt,fname,directory,firsttime=firsttime,chngobs=chngobs,colorbalance=colorbalance,invert=invert
;fname='10_24_70um_21'
;directory='10_24_24'
	direc='/data/David/drops/'+directory+'/'
;step 1 is putting the image together
	sr=FILE_SEARCH(direc+fname+'*'+'tif',count=count)
	count=count-1
	aaa=readtiffstack(direc+fname+'.tif')
	b=size(aaa)
	aaa=replicate(!values.F_NAN,3,b(1),b(2),b(3))
	aaa(0,*,*,*)=readtiffstack(direc+fname+'.tif',color=0)
	aaa(1,*,*,*)=readtiffstack(direc+fname+'.tif',color=1)
	aaa(2,*,*,*)=readtiffstack(direc+fname+'.tif',color=2)
	a=replicate(!values.F_NAN,3,b(1),b(2),b(3)*count)
	xdim=b(1)
	ydim=b(2)
	tdim=b(3)
	a(*,*,*,[0:tdim-1])=aaa
	for i=1,count-1 do begin &$
		a(0,*,*,[(tdim*i):tdim*(i+1)-1])=readtiffstack((direc+fname+'_'+string(floor(i))+'.tif').compress(),color=0) &$
		a(1,*,*,[(tdim*i):tdim*(i+1)-1])=readtiffstack((direc+fname+'_'+string(floor(i))+'.tif').compress(),color=1) &$
		a(2,*,*,[(tdim*i):tdim*(i+1)-1])=readtiffstack((direc+fname+'_'+string(floor(i))+'.tif').compress(),color=2) &$
	endfor
	b=size(a)
	tdim=b(4)
	if isa(colorbalance) then begin &$
		a(0,*,*,*)=a(0,*,*,*)*200.0/median(a(0,*,*,*)) &$
		a(1,*,*,*)=a(1,*,*,*)*200.0/median(a(1,*,*,*)) &$
		a(2,*,*,*)=a(2,*,*,*)*200.0/median(a(2,*,*,*)) &$
		a(where(a GT 255.0))=255.0 &$
	endif &$


;then go through the process of median subtracting
	medi=replicate(-1d,3,xdim,ydim)
	if keyword_set(firsttime) then begin
		satis=0
		print,'red median and stdev'
		print,[median(a(0,*,*,*)),stdev(a(0,*,*,*))]
		print,'green median and stdev'
		print,[median(a(1,*,*,*)),stdev(a(1,*,*,*))]
		print,'blue median and stdev'
		print,[median(a(2,*,*,*)),stdev(a(2,*,*,*))]
		while satis EQ 0 do begin &$
			READ, grndiff, PROMPT='Enter green difference: ' &$
			for i=0,xdim-1 do begin &$
				for j=0,ydim-1 do begin &$
					tmp1=[] &$
					m1=[] &$
					m2=[] &$
					k=0 &$
					tmp1=where((max(a(1,i,j,*),dimension=1) GT (a(0,i,j,*)-grndiff)),kk) &$
					if kk NE 0 then begin &$
						if kk LT 11 then begin &$
							tmp2=tmp1 &$
						endif else begin &$
							m1=reform(mean(a(*,i,j,tmp1),dimension=1)) &$
							tmpsort=sort(m1) &$
							tmp2=tmp1(tmpsort(length(tmpsort)-[1:10])) &$
						endelse &$
					endif else begin &$
						m2=reform(mean(a(*,i,j,*),dimension=1)) &$
						tmpsort=sort(m2) &$
						tmp2=tmp1(tmpsort(length(tmpsort)-[1:10])) &$
					endelse &$
					medi(0,i,j)=median(a(0,i,j,tmp2)) &$
					medi(1,i,j)=median(a(1,i,j,tmp2)) &$
					medi(2,i,j)=median(a(2,i,j,tmp2)) &$
				end &$
				print, (i+0.0)/xdim &$
			end &$
			tv,medi,true=1 &$
			READ, satis, PROMPT='satisfied?: ' &$
			write_gdf,medi,direc+fname+'medi' &$
		end
	endif else begin
		stup1=read_gdf(direc+'setup'+fname+'gdf')
		bckgrnd=stup1(0)
		grndiff=stup1(1)
		sr=FILE_SEARCH(direc+fname+'medi',count=count)
		if count EQ 0 then begin &$
			for i=0,xdim-1 do begin &$
				for j=0,ydim-1 do begin &$
					tmp1=[] &$
					m1=[] &$
					m2=[] &$
					k=0 &$
					tmp1=where((max(a(1,i,j,*),dimension=1) GT (a(0,i,j,*)-grndiff)),kk) &$
					if kk NE 0 then begin &$
							if kk LT 11 then begin &$
						tmp2=tmp1 &$
							endif else begin &$
							m1=reform(mean(a(*,i,j,tmp1),dimension=1)) &$
							tmpsort=sort(m1) &$
							tmp2=tmp1(tmpsort(length(tmpsort)-[1:10])) &$
						endelse &$
					endif else begin &$
						m2=reform(mean(a(*,i,j,*),dimension=1)) &$
						tmpsort=sort(m2) &$
						tmp2=tmp1(tmpsort(length(tmpsort)-[1:10])) &$
					endelse &$
					medi(0,i,j)=median(a(0,i,j,tmp2)) &$
					medi(1,i,j)=median(a(1,i,j,tmp2)) &$
					medi(2,i,j)=median(a(2,i,j,tmp2)) &$
				end &$
			end &$
		endif else begin &$
			medi=read_gdf(direc+fname+'medi') &$
		endelse &$
	endelse

; Median subtraction and thresholding
	grad=a
        medi2=medi
        for i=0,ydim-1 do begin &$
                grad(*,*,i,*)=192.0*a(*,*,i,*)/mean(medi(*,*,i)) &$
                medi2(*,*,i)=192.0*medi(*,*,i)/mean(medi(*,*,i)) &$
        endfor
        grad(where(grad GT 255.0))=255.0
        medi2(where(medi2 GT 255.0))=255.0
        medsub=grad
        for i=0,tdim-1 do begin &$
                medsub(0,*,*,i)=-(grad(0,*,*,i)-medi2(0,*,*)) &$
                medsub(1,*,*,i)=-(grad(1,*,*,i)-medi2(1,*,*)) &$
                medsub(2,*,*,i)=-(grad(2,*,*,i)-medi2(2,*,*)) &$
        end
        medsub(where(medsub GT 255))=255
        medsub(where(medsub LT 0))=0
        wtthrs=240
        minm=min(medi2)
        maxm=max(medi2)
        for i=0,xdim-1 do begin &$
                for j=0,ydim-1 do begin &$
                        for k=0,2 do begin &$
                                medsub(k,i,j,*)=((medi2(k,i,j)-minm)/(maxm-minm))*(255.0/medi2(k,i,j))*medsub(k,i,j,*) &$
                        end &$
                end &$
        end
        medsub(where(medsub LT 0))=0
        medsub=255.0-medsub

; Placing all the obstacles
	obs=float(readtext(direc+fname+'obs'))
	xobs1=obs(*,0)
	yobs1=obs(*,1)
	if isa(invert) then yobs1=max(obs(*,1))-obs(*,1)
	if (keyword_set(chngobs)) then begin
		ogvid=a
		satist=0 &$
		while satist NE 1 do begin &$
			READ, obst, PROMPT='Pick a Time: ' &$
			if obst LT tdim then begin &$
				tv,ogvid(*,*,*,obst),true=1 &$
				READ, satist, PROMPT='Good?: ' &$
			endif else begin &$
				print,'ERROR: tdim is'+string(tdim) &$
			endelse &$
		end &$
		satis=0
		rngtheta=[0:2*!PI:0.001]
		while satis NE 1 do begin &$
			READ, scl, PROMPT='Enter Scale: ' &$
			READ, xshift, PROMPT='Enter Xshift: ' &$
			READ, yshift, PROMPT='Enter Yshift: ' &$
			READ, ang, PROMPT='Enter Angle in Degrees: ' &$
			READ, obsr, PROMPT='Enter Obstacle Radius: ' &$
			xobs2=(xobs1*scl+xshift)*cos(ang*!pi/180.0)-(yobs1*scl+yshift)*sin(ang*!pi/180.0) &$
			yobs2=(xobs1*scl+xshift)*sin(ang*!pi/180.0)+(yobs1*scl+yshift)*cos(ang*!pi/180.0) &$
			xobs=xobs2(where((xobs2 GE 0) and (xobs2 LT xdim) AND (yobs2 GE 0) and (yobs2 LT ydim))) &$
			yobs=yobs2(where((xobs2 GE 0) and (xobs2 LT xdim) AND (yobs2 GE 0) and (yobs2 LT ydim))) &$
			ogvid1=reform(ogvid(*,*,*,obst)) &$
			rng=[[obsr*sin(rngtheta),(obsr+0.75)*sin(rngtheta),(obsr-0.75)*sin(rngtheta)],$
				[obsr*cos(rngtheta),(obsr+0.75)*cos(rngtheta),(obsr-0.75)*cos(rngtheta)]] &$
			for i=0,length(xobs)-1 do begin &$
				rngtmpx=rng(where((rng(*,0)+xobs(i) GE 0) AND (rng(*,0)+xobs(i) LE xdim-1) AND $
					(rng(*,1)+yobs(i) GE 0) AND (rng(*,1)+yobs(i) LE ydim-1)),0) &$
				rngtmpy=rng(where((rng(*,0)+xobs(i) GE 0) AND (rng(*,0)+xobs(i) LE xdim-1) AND $
					(rng(*,1)+yobs(i) GE 0) AND (rng(*,1)+yobs(i) LE ydim-1)),1) &$
				for j=0,length(rngtmpx)-1 do begin &$
					ogvid1(1,round(rngtmpx(j)+xobs(i)),round(rngtmpy(j)+yobs(i)))=0 &$
				endfor &$
			endfor &$
			tv,ogvid1,true=1 &$
			READ, satis, PROMPT='satisfied?: ' &$
		end 
		READ, obsum, PROMPT='How many microns across is an obstacle? ' &$
		pix2um=obsr/(obsum/2.0)
		READ, FPS, PROMPT='What was the framerate? ' &$
		print,string([scl,xshift,yshift,ang,obsr,pix2um])
		write_gdf,[scl,xshift,yshift,ang,obsr,pix2um,fps],direc+fname+'obssetup'
	endif else begin
		stup=read_gdf(direc+fname+'obssetup')
		scl=stup(0)
		xshift=stup(1)
		yshift=stup(2)
		ang=stup(3)
		obsr=stup(4)
		pix2um=stup(5)
		fps=stup(6)
	endelse
	xobs2=(xobs1*scl+xshift)*cos(ang*!pi/180.0)-(yobs1*scl+yshift)*sin(ang*!pi/180.0)
	yobs2=(xobs1*scl+xshift)*sin(ang*!pi/180.0)+(yobs1*scl+yshift)*cos(ang*!pi/180.0)
	xobs=xobs2(where((xobs2 GE 0) and (xobs2 LT xdim) AND (yobs2 GE 0) and (yobs2 LT ydim)))
	yobs=yobs2(where((xobs2 GE 0) and (xobs2 LT xdim) AND (yobs2 GE 0) and (yobs2 LT ydim)))

;fill in medsub with rings where the obstacles are
	xobs=xobs(where((xobs GT 0) AND (xobs LT xdim)))
	xobs=xobs(where((yobs GT 0) AND (yobs LT ydim)))
	rngtheta=[0:2*!PI:0.001]
	rng=[[obsr*sin(rngtheta),(obsr+0.75)*sin(rngtheta),(obsr-0.75)*sin(rngtheta)],[obsr*cos(rngtheta),(obsr+0.75)*cos(rngtheta),(obsr-0.75)*cos(rngtheta)]]
	for i=0,length(xobs)-1 do begin &$
		rngtmpx=rng(where((rng(*,0)+xobs(i) GE 0) AND (rng(*,0)+xobs(i) LE xdim-1) AND (rng(*,1)+yobs(i) GE 0) AND (rng(*,1)+yobs(i) LE ydim-1)),0) &$
		rngtmpy=rng(where((rng(*,0)+xobs(i) GE 0) AND (rng(*,0)+xobs(i) LE xdim-1) AND (rng(*,1)+yobs(i) GE 0) AND (rng(*,1)+yobs(i) LE ydim-1)),1) &$
		bleh=where((round(rngtmpy+yobs(i)) LT ydim) AND (round(rngtmpy+yobs(i)) GT 0) AND (round(rngtmpx+xobs(i)) LT xdim) AND (round(rngtmpx+xobs(i)) GT 0),rl) &$
		rngtmpx=rngtmpx(bleh) &$
                rngtmpy=rngtmpy(bleh) &$
		if rl GT 0 then begin &$
			for t=0,tdim-1 do begin &$
				for j=0,length(rngtmpx)-1 do begin &$
					medsub(1,round(rngtmpx(j)+xobs(i)),round(rngtmpy(j)+yobs(i)),t)=0 &$
				endfor &$
			endfor &$
		endif &$
		print,'placing obstacles'+string((i+0.0)/(length(xobs)-1)) &$
	endfor

;begin thresholding
	if keyword_set(firsttime) then begin
		satis=0
		thrsr=100
		thrsg=100
		thrsb=100
		thrsg1=100
		while satis EQ 0 do begin &$
			READ, thrsr, PROMPT='Enter red threshold. Previous value='+string(thrsr)+':' &$
			tv,bytscl(reform(medsub(0,*,*,*) GT thrsr)) &$
			READ, satis, PROMPT='satisfied?: ' &$
		end &$
		satis=0 &$
		while satis EQ 0 do begin &$
			READ, thrsg, PROMPT='Enter green threshold. Previous value='+string(thrsg)+':' &$
			tv,bytscl(reform(medsub(1,*,*,*) GT thrsg)) &$
			READ, satis, PROMPT='satisfied?: ' &$
		end &$
		satis=0 &$
		while satis EQ 0 do begin &$
			READ, thrsg1, PROMPT='Enter thicker green threshold. Previous value='+string(thrsg1)+':' &$
			tv,bytscl(reform(medsub(1,*,*,*) GT thrsg1)) &$
			READ, satis, PROMPT='satisfied?: ' &$
		end &$
		write_gdf,[-1,grndiff,thrsr,thrsg,thrsb,thrsg1],direc+'setup'+fname+'gdf'
	endif else begin
		thrsr=stup1(2)
		thrsg=stup1(3)
		thrsb=stup1(4)
		thrsg1=stup1(5)
	endelse
	niceg=reform(medsub(1,*,*,*) GT thrsg)
	nicer=reform(medsub(0,*,*,*) GT thrsr)

;the main processing of the thresholdded images are stored in the dmfill and dmassign programs
	output4=dmfill(nicer,niceg,xdim,ydim,tdim,tol=-0.6)
	output3=dmassign(-1.0*output4,reform(medsub(1,*,*,*)),thrsr,xdim,ydim,tdim)

;take out the rings we added in
	for i=0,length(xobs)-1 do begin &$
		rngtmpx=rng(where((rng(*,0)+xobs(i) GE 0) AND (rng(*,0)+xobs(i) LE xdim-1) AND (rng(*,1)+yobs(i) GE 0) AND (rng(*,1)+yobs(i) LE ydim-1)),0) &$
		rngtmpy=rng(where((rng(*,0)+xobs(i) GE 0) AND (rng(*,0)+xobs(i) LE xdim-1) AND (rng(*,1)+yobs(i) GE 0) AND (rng(*,1)+yobs(i) LE ydim-1)),1) &$
		bleh=where((round(rngtmpy+yobs(i)) LT ydim) AND (round(rngtmpy+yobs(i)) GT 0) AND (round(rngtmpx+xobs(i)) LT xdim-1) AND (round(rngtmpx+xobs(i)) GT 0),rl) &$
                rngtmpx=rngtmpx(bleh) &$
                rngtmpy=rngtmpy(bleh) &$
                if rl GT 0 then begin &$
			for t=0,tdim-1 do begin &$
				for j=0,length(rngtmpx)-1 do begin &$
					output3(round(rngtmpx(j)+xobs(i)),round(rngtmpy(j)+yobs(i)),t)=-1 &$
					output3(round(rngtmpx(j)+xobs(i)+1),round(rngtmpy(j)+yobs(i)),t)=-1 &$
				endfor &$
			endfor &$
		endif &$
		print,'removing obstacles'+string((i+0.0)/(length(xobs)-1)) &$
	endfor
	
;one last finishing check to clean up ugly droplets
	for t=0,tdim-1 do begin &$
		section=output3(*,*,t) &$
		indxcount=histogram(section,binsize=1,locations=indx) &$
		blobindx=indx(where((indxcount GT 0) AND (indx GE 0))) &$
		for i=0,length(blobindx)-1 do begin &$
			region=where(section EQ blobindx(i)) &$
			regionx=region mod xdim &$
			regiony=floor(region/xdim) &$
			region1=region &$
			while length(region1) GT 0 do begin &$
				tmpr=search2d(section,regionx(0),regiony(0),blobindx(i),blobindx(i)) &$
				tmprx=tmpr mod xdim &$
				tmpry=floor(tmpr/xdim) &$
				if length(tmpr) GT length(region)*0.5 then begin &$
					for j=0,length(regionx)-1 do begin &$
						output3(regionx(j),regiony(j),t)=-1 &$
					endfor &$
					for j=0,length(tmprx)-1 do begin &$
						output3(tmprx(j),tmpry(j),t)=blobindx(i) &$
					endfor &$
					region1=[] &$
					if (length(tmpr)+0.0)/length(region) LT 0.8 then print,'blob'+string(blobindx(i))+'at time'+string(t)+'is weird' &$ 
				endif else begin &$
					for j=0,length(tmprx)-1 do begin &$
						output3(tmprx(j),tmpry(j),t)=-1 &$
					endfor &$
					section=output3(*,*,t) &$
					region=where(section EQ blobindx(i)) &$
					regionx=region mod xdim &$
					regiony=floor(region/xdim) &$
					region1=region &$
				endelse &$
			end &$
		end &$
		print,'finishing step'+string((t+0.0)/tdim) &$
	end

write_gdf,output3,direc+'pt'+fname+'gdf'
	return,output3
end
