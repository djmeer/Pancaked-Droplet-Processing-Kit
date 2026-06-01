; dmcombtif.pro
; written by David Meer 6/2026
; detects droplet mergers and breakups
;
; within deth and brth
; 0->index through time
; 1->time
; 2->x
; 3->y
; 4->area
; 5->fate
; 6->life
; 7-12->6 neighbors given by time index
; 13->symmetry at time
function dmcombtif,brth,deth,tdim,xdim,obsr,data
	combin=[]
	for t=0,tdim-1 do begin &$
		chckb=where((brth(1,*) EQ t+1) AND brth(6,*) EQ 2) &$
		chckd=where((deth(1,*) EQ t) AND ((deth(5,*) EQ 2) OR (deth(5,*) EQ 1))) &$
		distmapb=[] &$
		distmapd=[] &$
;breakup
		if length(chckd) GT 0 and length(chckb) GT 1 then begin &$
			for j=0,length(chckb)-2 do begin &$
				for k=j+1,length(chckb)-1 do begin &$
					distmapb=[[distmapb],[j,k,$
						(brth(2,chckb(j))*brth(4,chckb(j))+brth(2,chckb(k))*brth(4,chckb(k)))/(brth(4,chckb(j))+brth(4,chckb(k))),$
						(brth(3,chckb(j))*brth(4,chckb(j))+brth(3,chckb(k))*brth(4,chckb(k)))/(brth(4,chckb(j))+brth(4,chckb(k))),-1,$
						sqrt((brth(2,chckb(j))-brth(2,chckb(k)))^2+(brth(3,chckb(j))-brth(3,chckb(k)))^2)]] &$
				endfor &$
			endfor &$
			for i=0,length(chckd)-1 do begin &$
				distmapb1=distmapb &$
				distmapb1(4,*)=sqrt((distmapb1(2,*)-deth(2,chckd(i)))^2+(distmapb1(3,*)-deth(3,chckd(i)))^2) &$
				trgt=where((distmapb1(4,*) EQ min(distmapb1(4,*))) AND (min(distmapb1(4,*)) LT sqrt(deth(4,chckd(i))/!pi))$
					AND (distmapb(5,*) LT 2*obsr+sqrt(brth(4,chckb(distmapb1(0,*)))/!pi)+sqrt(brth(4,chckb(distmapb1(1,*)))/!pi))$
					AND (brth(4,chckb(distmapb1(0,*))) LT deth(4,chckd(i))) AND (brth(4,chckb(distmapb1(1,*))) LT deth(4,chckd(i)))) &$
				if trgt(0) GE 0 then begin &$
					combin=[[combin],$
					[deth(0,chckd(i)),brth(0,chckb(distmapb(0,trgt))),brth(0,chckb(distmapb(1,trgt))),t,deth(2,chckd(i)),deth(3,chckd(i)),$
5,data(76,where((data(2,*) EQ deth(0,chckd(i))) AND (data(41,*) EQ brth(0,chckb(distmapb(0,trgt))))))]] &$
				endif &$
			endfor &$
		end &$
;combine
		if length(chckb) GT 0 and length(chckd) GT 1 then begin &$
			for j=0,length(chckd)-2 do begin &$
				for k=j+1,length(chckd)-1 do begin &$
					distmapd=[[distmapd],[j,k,$
						(deth(2,chckd(j))*deth(4,chckd(j))+deth(2,chckd(k))*deth(4,chckd(k)))/(deth(4,chckd(j))+deth(4,chckd(k))),$
						(deth(3,chckd(j))*deth(4,chckd(j))+deth(3,chckd(k))*deth(4,chckd(k)))/(deth(4,chckd(j))+deth(4,chckd(k))),-1,$
						sqrt((deth(2,chckd(j))-deth(2,chckd(k)))^2+(deth(3,chckd(j))-deth(3,chckd(k)))^2)]] &$
				endfor &$
			endfor &$
			for i=0,length(chckb)-1 do begin &$
				distmapd1=distmapd &$
				distmapd1(4,*)=sqrt((distmapd1(2,*)-brth(2,chckb(i)))^2+(distmapd1(3,*)-brth(3,chckb(i)))^2) &$
				trgt=where((distmapd1(4,*) EQ min(distmapd1(4,*))) AND (min(distmapd1(4,*)) LT sqrt(brth(4,chckb(i))/!pi))$
					AND (distmapd(5,*) LT 4*obsr+4*sqrt(deth(4,chckb(distmapd1(0,*)))/!pi)+4*sqrt(deth(4,chckd(distmapd1(1,*)))/!pi))$
					AND (deth(4,chckd(distmapd1(0,*))) LT brth(4,chckb(i))) AND (deth(4,chckd(distmapd1(1,*))) LT brth(4,chckb(i)))) &$
				if trgt(0) GE 0 then begin &$
					combin=[[combin],$
						[brth(0,chckb(i)),deth(0,chckd(distmapd(0,trgt))),deth(0,chckd(distmapd(1,trgt))),t,brth(2,chckb(i))$
,brth(3,chckb(i)),4,-1]] &$
				endif &$
			endfor &$
		end &$
	end
	return,combin
end
