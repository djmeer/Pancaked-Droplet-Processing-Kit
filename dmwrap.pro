; dmtracktif
; written by David Meer 6/2026
; Runs dmtracktif on many different files and compiles all the relevant information into a large set of arrays
function dmwrap,bupall=bupall,plst=plst
	directories=['3_24_26','4_17_26','4_23_26','5_5_26']
	fnames=['3_24_26','4_17_26b','4_23_26a','5_5_26c']
	properties=[[85,85,50],[85,85,50],[85,85,50],[85,85,50]]
	pix2um=2500/473.0
	pix2m=pix2um*1e-6
	fps=60.0 &$
	nexps=length(fnames)
	dlst=[2,14,15,16,17,18,19,34,35,36,37,38,39,40,41,50,51,52,53,54,55] &$
	clst=[0,1,2,3] & blst=[0,1,2,3] & nlst=[0,1,2,3] & ctlst=[0,1] & nnlst=[4,5] & ddlst=[2] & ilst=[2] &$
	if length(fnames) NE length(directories) then print,'ERROR: Mismatch in number of experiments and directories'
	surftens=20.1*10.0^(-3) & dens=0.96*10.0^(3)
	xdim=640 & ydim=512
	breakall=[] & bupall=[] & plst=[] & ddrop=[] & f=[] & bgg=[]
	for q=0,3 do begin &$
		fname=fnames(q) &$
		directory=directories(q) &$
		direc='/data/david/coal/'+directory+'/' &$
		sra=FILE_SEARCH(direc+fname+'*bgg',count=count) &$
		if count LT 1 then begin &$
			print,'begin processing for '+direc+fname &$
			visc=properties(2,q)*10.0^(-6) &$
			dvisc=visc*dens &$
			viscr=0.001/dvisc &$
			zheight=properties(0,q)*10.0^(-6) & ro=properties(1,q)*10.0^(-6) &$
			dlst=[2,14,15,16,17,18,19,34,35,36,37,38,39,40,41,50,51,52,53,54,55] &$
			clst=[0,1,2,3] & blst=[0,1,2,3] & nlst=[0,1,2,3] & ctlst=[0,1] & nnlst=[4,5] & ddlst=[2] & ilst=[-1] &$
			data=[] &$
			c=[] & b=[] & n=[] & dd=[] & ct=[] & bg=[] &$
			srd=FILE_SEARCH(direc+fname+'*trackbglst',count=count) &$
			if count LT 1 then begin &$
				srb=FILE_SEARCH(direc+fname+'*ptcollect',count=count) &$
				npull=strlen(direc)-1 &$
				for i=0,count-1 do begin &$
					if i NE 22 then begin &$
						stemp=srb(i) &$
						pull=stemp.Remove(0,npull) &$
						fnam=pull.Remove(-9) &$
						dat=dmtracktif(fnam,directory,combin=combin,brek=brek,conec=conec,ddrops=ddrops,contacttim=contacttim,bglst=bglst,interface=1,coal=1,newptedge=1,nnn=1) &$
						dat(dlst,*)=dat(dlst,*)+10000.0*i &$
						if combin(0) GT -1 then	begin &$
							combin(clst,*)=combin(clst,*)+10000.0*i &$
							c=[[c],[combin]] &$
						endif &$
						szct=size(contacttim) &$
						if szct(0) GT 0 then begin &$
							contacttim(ctlst,*)=contacttim(ctlst,*)+10000.0*i &$
							ct=[[ct],[contacttim]] &$
						endif &$
						if brek(0) GT -1 then begin &$
							brek=brek(*,where(brek(6,*) EQ 0)) &$
							brek(blst,*)=brek(blst,*)+10000.0*i &$
							b=[[b],[brek]] &$
						endif &$
						if ddrops(0) GT -1 then begin &$
							dropl=replicate(!values.F_NAN,length(ddrops)+2,length(transpose(ddrops))) &$
							dropl([0:length(ddrops)-1],*)=ddrops &$
							dropl(length(ddrops),*)=i &$
							dropl(length(ddrops)+1,*)=q &$
							dd=[[dd],[dropl]] &$
						endif &$
						if conec(0) GT -1 then begin &$
							conec(nlst,*)=conec(nlst,*)+10000.0*i &$
							tmp1=where(conec(nnlst,*) GE 0) &$
							if tmp1(0) NE -1 then begin &$
								conec(4,where(conec(4,*) GE 0))=conec(4,where(conec(4,*) GE 0))+10000.0*i &$
								conec(5,where(conec(5,*) GE 0))=conec(5,where(conec(5,*) GE 0))+10000.0*i &$
							endif &$
							n=[[n],[conec]] &$
						endif &$
						print,i &$
						if isa(bglst) EQ 1 then begin &$
							if bglst(0) GT -1 then begin &$
								bglst(-1,*)=bglst(-1,*)+10000.0*i &$
								bglst(1,*)=bglst(1,*)*pix2um &$
;								print,size(bglst) &$
								bgll=length(transpose(bglst)) &$
;								print,bgll &$
								bg=[[bg],[bglst,replicate(q,1,bgll)]] &$
							endif &$
						endif &$
						data=[[data],[dat]] &$
						print,[(i+0.0)/count,(q+0.0)/nexps] &$
					endif &$
				endfor &$
				combin=c &$
				if isa(combin) EQ 1 then write_gdf,combin,direc+fname+'trackc' &$
				brek=b &$
				write_gdf,brek,direc+fname+'trackb' &$
				conec=n &$
				if isa(combin) EQ 1 then write_gdf,conec,direc+fname+'trackn' &$
				ddrops=dd &$
				if isa(ddrops) EQ 1 then write_gdf,ddrops,direc+fname+'trackd' &$
				contacttim=ct &$
				if isa(contacttim) EQ 1 then write_gdf,contacttim,direc+fname+'trackt' &$
				bglst=bg &$
				if isa(bglst) EQ 1 then write_gdf,bglst,direc+fname+'trackbglst' &$
				write_gdf,data,direc+fname+'trackdata' &$
			endif else begin &$
				sre=FILE_SEARCH(direc+fname+'trackc',count=count) &$
				if count EQ 1 then begin &$
					combin=read_gdf(direc+fname+'trackc') &$
				endif else begin &$
					combin=replicate(-1,6,1) &$
				endelse &$
				brek=read_gdf(direc+fname+'trackb') &$
				sre=FILE_SEARCH(direc+fname+'trackn',count=count) &$
				if count EQ 1 then begin &$
					conec=read_gdf(direc+fname+'trackn') &$
				endif else begin &$
					conec=replicate(-1,10,1) &$
				endelse &$
				sre=FILE_SEARCH(direc+fname+'trackd',count=count) &$
				if count EQ 1 then begin &$
					ddrops=read_gdf(direc+fname+'trackd') &$
				endif else begin &$
					ddrops=replicate(-1,length(ddrop),1) &$
				endelse &$
				sre=FILE_SEARCH(direc+fname+'trackbglst',count=count) &$
				if count EQ 1 then begin &$
					bglst=read_gdf(direc+fname+'trackbglst') &$
				endif else begin &$
					bglst=replicate(-1,length(13),1) &$
				endelse &$
				contacttim=read_gdf(direc+fname+'trackt') &$
				data=read_gdf(direc+fname+'trackdata') &$
			endelse &$
			bgg=[[bgg],[bglst]] &$
			write_gdf,bglst,direc+fname+'trackbgg' &$
		endif else begin &$
			bgl=read_gdf(direc+fname+'trackbgg') &$
			bgg=[[bgg],[bgl]] &$
		endelse
	end
	return,bgg
end
