; dmmov3
; written by David Meer 6/2026
; displays tracked videos of droplets moving with some extra bit of data about them "trackdat" being displayed, as well as droplet breakup and coalescence
pro dmmov3,imag,xdat,ydat,tdat,trackdat,brk,combin
	tmax=max(tdat)
	delay = 0.1
	i=0
	c=0
	sz=size(imag)
	while 1 do case strupcase(get_kbrd(0)) of
		'': begin
			tv,bytscl(imag(*,*,i))
			cgtext,xdat(where(tdat eq i)),ydat(where(tdat eq i)),'.',/device,charthick=2,color=[500000]
			if (trackdat(0) mod 1) EQ 0 then begin
				cgtext,xdat(where(tdat eq i))-87,ydat(where(tdat eq i)),string(round(trackdat(where(tdat eq i)))),/device,charthick=2,color=[500000]
			endif else begin
				cgtext,xdat(where(tdat eq i))-50,ydat(where(tdat eq i)),string(trackdat(where(tdat eq i))),/device,charthick=2,color=[500000]
			endelse
			cgtext,-40,10,string(i),/device,charthick=2,color=[500000]
			for j=0,19 do begin
				lstb=where(brk(3,*) EQ i+j)
				lstc=where(combin(3,*) EQ i+j)
				if lstb(0) NE -1 then cgtext,brk(4,lstb),brk(5,lstb),'X',/device,charthick=5,color=[600000]
				if lstc(0) NE -1 then cgtext,combin(4,lstc),combin(5,lstc),'O',/device,charthick=5,color=[8000200]
			endfor
			i=i+1
			if i EQ tmax then i=0
			wait,delay
		endcase
		'S':delay = delay * 1.5
		'F':delay = delay / 1.5
		'Q':begin
			print,'time'+string(i)
			return
		endcase
		else:
	endcase
end

