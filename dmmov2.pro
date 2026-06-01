; dmmov2
; written by David Meer 6/2026
; displays tracked videos of droplets moving with some extra bit of data about them "trackdat" being displayed
pro dmmov2,imag,xdat,ydat,tdat,trackdat
	tmax=max(tdat)
	trackdat=round(trackdat)
	delay = 0.03
	i=0
	c=0
	sz=size(imag)
	while 1 do case strupcase(get_kbrd(0)) of
		'': begin
			tv,imag(*,*,i)*255
			cgtext,xdat(where(tdat eq i)),ydat(where(tdat eq i)),'.',/device,charthick=2,color=[500000]
			cgtext,xdat(where(tdat eq i))-90,ydat(where(tdat eq i)),string(trackdat(where(tdat eq i))),/device,charthick=2,color=[500000]
			cgtext,10,10,string(i),/device,charthick=2,color=[500000]
			i=i+1
			if i EQ tmax+1 then i=0
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

