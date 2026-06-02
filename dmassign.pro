; dmpt
; written by David Meer 6/2026
; fills in droplet centers with an index, then slowly fills in bordering regions of the black ring with those indexes until all pixels are assigned.

function dmassign,output1,medd,thrsr,xdim,ydim,tdim
	output3=replicate(!values.F_nan,xdim,ydim,tdim)
	for t=0,tdim-1 do begin
		output=output1(*,*,t)
		grn1=medd(*,*,t) LT thrsr
		gtemp=reform(grn1)
		len=where(gtemp EQ 0)
		c=0
		while length(len) GT 0 do begin &$
			srch=search2d(grn1,len(0) mod xdim,floor(len(0)/xdim),0,0) &$
			srchx=srch mod xdim &$
			srchy=floor(srch/xdim) &$
			if ((max(srchx) GE xdim-1) AND (min(srchx) LE 1) AND (max(srchy) GE ydim-1) AND (min(srchy) LE 1)) then begin &$
				xdim=xdim &$
			endif else begin &$
				c=c+1 &$
				if (mean(output1(srchx,srchy,t)) GT -0.4) then begin &$
					output(srchx,srchy)=c &$
				endif &$
			endelse &$
			gtemp(srchx,srchy)=1 &$
			len=where(gtemp EQ 0) &$
;			tv,bytscl(output) &$
;			print,c &$
		end
		for i=1,c do begin &$
			tmpa=where(output EQ i) &$
			if length(tmpa) LE 10 then begin &$
				tmpax=tmpa mod xdim &$
				tmpay=floor(tmpa/xdim) &$
				for j=0,length(tmpa)-1 do begin &$
					output(tmpax(j),tmpay(j))=-1 &$
				end &$
			endif &$
		end
		areas=replicate(-1,c+1,1)
		for i=1,c do begin &$
			areas(i)=length(where(output EQ i)) &$
		end
		fillz=where(output EQ 0)
		fillx=fillz mod xdim
		filly=floor(fillz/xdim)
		output2=output
		while (fillz(0) NE -1) do begin &$ 
			bord=search2d(output,fillx(0),filly(0),0,0) &$
			bordx=bord mod xdim &$
			bordy=floor(bord/xdim) &$
			sav1=[output(bordx,bordy+1),output(bordx+1,bordy),output(bordx-1,bordy),output(bordx,bordy-1)] &$
			kep1=uniq(sav1(sort(sav1))) &$
			kep2=sav1(sort(sav1)) &$
			kep=kep2(kep1) &$
			kep=kep(where(kep GT 0)) &$
;			print,kep &$
			if length(kep) GT 1 then begin &$
				bord1=bord+0.0 &$
				bordx1=bordx &$
				bordy1=bordy &$
				jj=0 &$
				output5=output2 &$
;				output5t=output5 &$
				while (length(bord1) GT 0) AND (jj LT 20) do begin &$
					for i=0,length(bord1)-1 do begin &$
						xtemps=[bordx(i),bordx(i)+1,bordx(i)-1,bordx(i)] &$
						ytemps=[bordy(i)+1,bordy(i),bordy(i),bordy(i)-1] &$
						ctemps=where((xtemps GE -1) AND (ytemps GE -1) AND (xtemps LT xdim) AND (ytemps LT ydim)) &$
						chck=[output2(xtemps(ctemps),ytemps(ctemps))] &$
						wchck1=where(chck GT 0) &$
						wchck2=wchck1(sort(wchck1)) &$
						wchck=chck(wchck2(unique(wchck2))) &$
						if wchck1(0) NE -1 then begin &$
							if length(wchck) EQ 1 then begin &$
								if output2(bordx(i),bordy(i)) EQ 0 then output5(bordx(i),bordy(i))=wchck(0) &$
							endif else begin &$
								wchck3=wchck(where(areas(wchck) EQ max(areas(wchck)))) &$
								if output2(bordx(i),bordy(i)) EQ 0 then output5(bordx(i),bordy(i))=wchck3(0) &$
							endelse &$
						endif &$
					endfor &$
;					tv,bytscl(output5t) &$
					btemp=where(output5(bordx1,bordy1) EQ 0) &$
					bord1=bord(btemp) &$
					bordx=bordx1(btemp) &$
					bordy=bordy1(btemp) &$
					output2=output5 &$
					jj=jj+1 &$
					if jj GT 20 then print, 'jj GT 20 at'+string([bordx(0),bordy(0)]) &$
				end &$
			endif else begin &$
				if kep(0) GT 0 then begin &$
					output2(bordx,bordy)=kep(0) &$
				endif else begin &$
					output2(bordx,bordy)=-1 &$
				endelse &$
			endelse &$
			fillz=where(output2 EQ 0) &$
			fillx=fillz mod xdim &$
			filly=floor(fillz/xdim) &$
;			print,length(fillz) &$
;			tv,bytscl(output2) &$
		end
		output3(*,*,t)=output2
		print,string((t+0.0)/tdim)+' index assign'
	end
	return,output3
end
