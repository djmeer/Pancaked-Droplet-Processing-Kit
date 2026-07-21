;	David Meer 8/19/2024

function dmedgewalk3,x,y,xmax,ymax,plist=plist,ls=ls
	b=replicate(0,xmax,ymax)
	b(x,y)=1
	xcom=mean(x)
	ycom=mean(y)
	strtx=x(where(y EQ max(y)))
	lim=length(x)
	strt=reform([strtx(0),max(y)])
	sv=strt
	i=double(0)
	satis=0
	k=0
	p=0
	while ((satis EQ 0) AND (i LT lim)) do begin &$
		j=0 &$
		satis1=0 &$
		while (satis1 EQ 0) AND (j LT 8) do begin &$
			case (j+p) mod 8 of &$
				0: nxt=[strt(0),strt(1)+1] &$
				1: nxt=[strt(0)+1,strt(1)+1] &$
				2: nxt=[strt(0)+1,strt(1)] &$
				3: nxt=[strt(0)+1,strt(1)-1] &$
				4: nxt=[strt(0),strt(1)-1] &$
				5: nxt=[strt(0)-1,strt(1)-1] &$
				6: nxt=[strt(0)-1,strt(1)] &$
				7: nxt=[strt(0)-1,strt(1)+1] &$
			endcase &$
			if ((nxt(0) GT -1) AND (nxt(0) LT xmax) AND (nxt(1) GT -1) AND (nxt(1) LT ymax)) then begin &$
				if (b(nxt(0),nxt(1)) EQ 1) then begin &$
					if (nxt(0)+1 LT xmax) then begin &$
						stp1=b(nxt(0)+1,nxt(1)) &$
					endif else begin &$
						stp1=-2 &$
					endelse &$
					if (nxt(1)+1 LT ymax) then begin &$
						stp2=b(nxt(0),nxt(1)+1) &$
					endif else begin &$
						stp2=-2 &$
					endelse &$
					if (nxt(0)-1 GT -1) then begin &$
						stp3=b(nxt(0)-1,nxt(1)) &$
					endif else begin &$
						stp3=-2 &$
					endelse &$
					if (nxt(1)-1 GT -1) then begin &$
						stp4=b(nxt(0),nxt(1)-1) &$
					endif else begin &$
						stp4=-2 &$
					endelse &$
					if (stp1 EQ 0) OR (stp2 EQ 0) OR (stp3 EQ 0) OR (stp4 EQ 0) then begin &$
						chck=where((sv(0,*) EQ nxt(0)) AND (sv(1,*) EQ nxt(1))) &$
						if chck(0) EQ -1 then begin &$
							strt=nxt &$
							sv=[[sv],[strt]] &$
							satis1=1 &$
							k=0 &$
							p=(((j+p) mod 8)-2+8) mod 8 &$
						end &$
						if (strt(0) EQ strtx(0)) AND (strt(1) EQ max(y)) then begin &$
							satis=1 &$
							satis1=1 &$
						endif &$
					endif &$
				endif &$
			endif &$
			j=j+1 &$
		end &$
		if (satis1 EQ 0) and (j EQ 8) then begin &$
			k=k+1.0 &$
			strt=[sv(0,-k),sv(1,-k)] &$
		end &$
		i=i+1.0 &$
	end
	plist=sv
	lon=length(transpose(sv))
	xx = REBIN(sv[0,*], lon, lon)
	yy = REBIN(sv[1,*], lon, lon)
	darr = SQRT((xx - TRANSPOSE(xx))^2 + (yy - TRANSPOSE(yy))^2)
	sarr=smooth(darr,10.0,/edge_wrap)
	msarr=255.0/max(sarr)
;	plot,[0],[0],/NODATA,xrange=[0,lon-1],yrange=[0,lon-1],charsize=3,/isotropic
;	for i=0,lon-1 do begin &$
;		for j=0,lon-1 do begin &$
;			oplot,[i],[j],psym=symcat(16),color=msarr*sarr(i,j) &$
;		endfor &$
;	endfor
	mns=where((sarr LT SHIFT(sarr, 1, 0)) AND (sarr LT SHIFT(sarr, -1, 0)) AND (sarr LT SHIFT(sarr, 0, 1)) AND (sarr LT SHIFT(sarr, 0, -1)) AND $
		(sarr LT SHIFT(sarr, 1, 1)) AND (sarr LT SHIFT(sarr, -1, -1)) AND (sarr LT SHIFT(sarr, 1, -1)) AND (sarr LT SHIFT(sarr, -1, 1)))
	mnsx=mns mod lon
	mnsy=mns/lon
;	oplot,[mnsx],[mnsy],color=[255.0*256.0],psym=symcat(16)
	tempm=where((abs(mnsx-mnsy) GT 2) AND (mnsx GT mnsy))
	if tempm(0) GT -1 then begin &$
		mnsx1=mnsx(tempm) &$
		mnsy1=mnsy(tempm) &$
;		oplot,[mnsx1],[mnsy1],color=[255.0*256.0^2],psym=symcat(16) &$
;now i check for moving within the system rather than outside it
		midx=round((sv(0,mnsx1)+sv(0,mnsy1))/2) &$
		midy=round((sv(1,mnsx1)+sv(1,mnsy1))/2) &$
		wlst=[] &$
		for i=0,length(mnsx1)-1 do begin &$
			tmpw=where((x EQ midx(i)) and (y EQ midy(i))) &$
			if tmpw(0) GT -1 then wlst=[wlst,i] &$
		endfor &$
		if isa(wlst) EQ 1 then begin &$
			mnsx2=mnsx1(wlst) &$
			mnsy2=mnsy1(wlst) &$
			ans=where(sarr(mnsx2,mnsy2) EQ min(sarr(mnsx2,mnsy2))) &$
			ans=ans(0) &$
			ang=atan( (sv(1,(mnsx2(ans)))-sv(1,(mnsy2(ans))))/(sv(0,(mnsx2(ans)))-sv(0,(mnsy2(ans)))) ) &$
			neck=[sarr(mnsx2(ans),mnsy2(ans)),mnsx2(ans),mnsy2(ans),ang] &$
			output=[neck(0),neck(3),sv(0,neck(1)),sv(1,neck(1)),sv(0,neck(2)),sv(1,neck(2))] &$
		endif else begin &$
			output=[!Values.F_NAN,!Values.F_NAN,!Values.F_NAN,!Values.F_NAN,!Values.F_NAN,!Values.F_NAN] &$
		endelse &$
	endif else begin &$
		output=[!Values.F_NAN,!Values.F_NAN,!Values.F_NAN,!Values.F_NAN,!Values.F_NAN,!Values.F_NAN] &$
	endelse &$
;oplot,[sv(0,neck(1)),sv(0,neck(2))],[sv(1,neck(1)),sv(1,neck(2))],psym=5
	return,output
end
