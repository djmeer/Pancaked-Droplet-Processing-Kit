; dmpt
; written by David Meer 6/2026
; checks every region of every image and fills it in if there are colored droplets filling in that region
;
;keywords
; If keyword format is not specified, it takes the value 1 to turn on.
; tol - the tolerance of how much space within a region needs to be filled to automatically fill in the whole region. Use negative numbers between 0 and -1
;
function dmfill,nicer,niceg,xdim,ydim,tdim,tol=tol
	;output=[]
	if isa(tol) EQ 0 then tol=-0.4
	output=bytarr(xdim,ydim,tdim)
	for t=0,tdim-1 do begin &$
		nicegg=niceg(*,*,t) &$
		nicerr=nicer(*,*,t)
;		print,'step 1' &$
		rtemp=nicer(*,*,t) &$
		full=rtemp &$
		len=where(rtemp EQ 1) &$
		while length(len) GT 0 do begin &$
;			print,'step 2' &$
			srch=search2d(rtemp,len(0) mod xdim,floor(len(0)/xdim),1,1) &$
			srchx=srch mod xdim &$
			srchy=floor(srch/xdim) &$
;			print,'step 3' &$
			mxx=max(srchx,min=mnx) &$
			mxy=max(srchy,min=mny) &$
			if ((mxx GE xdim-1) AND (mnx LE 1) AND ((mxy GE ydim-1) OR (mny LE 1))) then begin &$
				xdim=xdim &$
			endif else begin &$
;				print,'step 4' &$
;				print,n_elements(srchx) &$
;				print,total(nicegg(srchx,srchy)+0.0-nicerr(srchx,srchy))/n_elements(srchx) &$
				if total(nicegg(srchx,srchy)+0.0-nicerr(srchx,srchy))/n_elements(srchx) LT tol then begin &$
;					print,'step 5' &$
					full(srchx,srchy)=0 &$
;					print,mean(niceg(srchx,srchy,t)+0.0-nicer(srchx,srchy,t)) &$
				endif &$
			endelse &$
;			print,'step 6' &$
			rtemp(srchx,srchy)=0 &$
;			print,'step 7' &$
			len=where(rtemp EQ 1) &$
;			tv,rtemp*255 &$
		end &$
;		print,'step 8' &$
		output(*,*,t)=full &$
		print,string((t+0.0)/tdim)+' fillred' &$
	end
	return,output
end
