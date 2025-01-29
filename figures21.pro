PRO figures21, ifig

!P.font = 0
!P.charsize = 1.5
thk = 5

if (fix(ifig) eq 0) then begin
   if (ifig eq 0.1) then Nt = 153
   if (ifig eq 0.2) then Nt = 174
   if (ifig eq 0.3) then Nt = 169
   if (ifig eq 0.4) then Nt = 162
   rdata = fltarr(3,Nt)
   tt = fltarr(Nt)
   mdot = fltarr(Nt)
   if (ifig eq 0.1) then openr,1,'mdot_ThinHR.dat'
   if (ifig eq 0.2) then openr,1,'mdot_ThinHR0.5new.dat'
   if (ifig eq 0.3) then openr,1,'mdot_ThinHR0.9.dat'
   if (ifig eq 0.4) then openr,1,'mdot_ThinHR0.99new.dat'
   readf,1,rdata
   close,1
   tt(*) = rdata(0,*)
   mdot(*) = rdata(1,*)
   plot_io,tt,mdot,xrange=[0,1.8d4],yrange=[1d-5,2d-4],xstyle=1,ystyle=1,$
           xtitle='!3t (M)', ytitle = 'Accretion Rate',thick=thk
   if (ifig eq 0.1) then xyouts,1000,1.5d-4,'a/M=0'
   if (ifig eq 0.2) then xyouts,1000,1.5d-4,'a/M=0.5'
   if (ifig eq 0.3) then xyouts,1000,1.5d-4,'a/M=0.9'
   if (ifig eq 0.4) then xyouts,1000,1.5d-4,'a/M=0.99'
endif

if (fix(ifig) eq 10) then begin
   if (ifig eq 10.1) then Nr = 912
   if (ifig eq 10.2) then Nr = 960
   if (ifig eq 10.3) then Nr = 1020
   if (ifig eq 10.4) then Nr = 1056
   rdata = fltarr(4,Nr)
   rr = fltarr(Nr)
   mdot = fltarr(Nr)
   if (ifig eq 10.1) then openr,1,'mdot_tavg_00.dat'
   if (ifig eq 10.2) then openr,1,'mdot_tavg_05.dat'
   if (ifig eq 10.3) then openr,1,'mdot_tavg_09.dat'
   if (ifig eq 10.4) then openr,1,'mdot_tavg_99.dat'
   if (ifig eq 10.1) then aa = 0.0
   if (ifig eq 10.2) then aa = 0.5
   if (ifig eq 10.3) then aa = 0.9
   if (ifig eq 10.4) then aa = 0.99
   a2 = aa*aa
   Z1 = 1.+(1-aa*aa)^(1./3.)*((1.+aa)^(1./3.)+(1-aa)^(1./3.))
   Z2 = sqrt(3.*aa*aa+Z1^2.)
   Risco = (3.+Z2-((3.-Z1)*(3+Z1+2.*Z2))^0.5)
   readf,1,rdata
   close,1
   rr(*) = rdata(0,*)
   mdot(*) = -rdata(1,*)
   inner = where(rr lt 1.1*Risco)
   mdot_in = mean(mdot(inner))
   mdot = mdot/mdot_in
   plot_oo,rr/Risco,mdot,xrange=[0.3,4],yrange=[0.3,3],xstyle=1,ystyle=1,$
           xticks = 2,xtickname = ['0.3','1.0','3.0'],$
           xtitle='!3r/r!LISCO!N', ytitle = 'Local <Mdot>/Inner <Mdot>',thick=thk
   if (ifig eq 10.1) then xyouts,0.35,2.5,'a/M=0'
   if (ifig eq 10.2) then xyouts,0.35,2.5,'a/M=0.5'
   if (ifig eq 10.3) then xyouts,0.35,2.5,'a/M=0.9'
   if (ifig eq 10.4) then xyouts,0.35,2.5,'a/M=0.99'
endif

if (fix(ifig) eq 1) then begin
    runids = [8151,8251,8352,8452]
    spec_runids = [8199,8299,8398,8498]
    aa_ = [0.0,0.5,0.9,0.99]
    sigma_SB = 5.6705d-5
    N1 = 0
    N2 = 0
    N3 = 0
    Nth = 41
    Nspec = 101
    for iaa = 0,3 do begin
        aa = aa_(iaa)
        a2 = aa*aa
        Z1 = 1.+(1-aa*aa)^(1./3.)*((1.+aa)^(1./3.)+(1-aa)^(1./3.))
        Z2 = sqrt(3.*aa*aa+Z1^2.)
        Risco = (3.+Z2-((3.-Z1)*(3+Z1+2.*Z2))^0.5)
        Eisco = (Risco*Risco-2*Risco+aa*sqrt(Risco))/ $
          (Risco*sqrt(Risco*Risco-3*Risco+2*aa*sqrt(Risco))) 
        eta = 1.-Eisco
        Rhor = 1.+sqrt(1.-a2)

        R_attach = 15.
        if (aa eq 0.0) then begin
            mdotfile = 'mdot_tavg_00.dat'
            Nr2 = 912
            ntrunid = 6000
        endif
        if (aa eq 0.5) then begin
            mdotfile = 'mdot_tavg_05.dat'
            Nr2 = 960
            ntrunid = 6010
        endif
        if (aa eq 0.9) then begin
            mdotfile = 'new_mdot_tavg_09.dat'
;            Nr2 = 1020
            Nr2 = 204
            ntrunid = 6018
        endif
        if (aa eq 0.99) then begin
            mdotfile = 'new_mdot_tavg_99.dat'
;            Nr2 = 1056
            Nr2 = 210
            ntrunid = 6027
            R_attach = 10.
        endif
        if ((aa eq 0.0)or(aa eq 0.5)) then rdata2 = fltarr(4,Nr2)
        if ((aa eq 0.9)or(aa eq 0.99)) then rdata2 = fltarr(2,Nr2)
        openr,1,mdotfile
        readf,1,rdata2
        close,1

        runid = runids(iaa)
        spec_runid = spec_runids(iaa)
        rdatafile = 'data/gr_0000.dat'
        dumpstr = string(runid,format='(I4.4)')
        strput,rdatafile,dumpstr,8
        openr,1,rdatafile
        readf,1,N1,N2,N3
        r = fltarr(N1)
        t = fltarr(N2)
        p = fltarr(N3)
        readf,1,r,t,p
        close,1
        Nr = N1
        drr = deriv(r)
        Tdisk = fltarr(N1,N3)
        blnk = fltarr(N1,N3)

        Rspec = fltarr(Nspec,Nr)
        Ispec = fltarr(Nspec,Nth)
        Ispecr = fltarr(Nspec,Nth,Nr)
        Qspec = fltarr(Nspec,Nth)
        Qspecr = fltarr(Nspec,Nth,Nr)
        Uspec = fltarr(Nspec,Nth)
        Uspecr = fltarr(Nspec,Nth,Nr)
        Rspecr = fltarr(Nspec,Nth,Nr)
        Cspecr = fltarr(Nspec,Nr)
        spec = fltarr(Nspec)
        rdata = fltarr(3,Nr)
        spec_s = dblarr(6,Nspec,Nth)
        spec2 = dblarr(Nspec,Nth)
        Qspec2 = dblarr(Nspec,Nth)
        Uspec2 = dblarr(Nspec,Nth)

        Mdotr = fltarr(N1)
        for ir=0,Nr2-1 do begin
            if (rdata2(0,ir) gt r(0)) then begin
                inrr = where(r le rdata2(0,ir))
                rdex = max(inrr)
                Mdotr(rdex) = rdata2(1,ir)
            endif
        endfor
        if ((aa eq 0.9)or(aa eq 0.99)) then Mdotr(0:Nr/2-1)=rdata2(1,*)
        Mdot0 = mean(Mdotr(where((r gt 2)and(r lt 6))))
        Mdot0 = mean(Mdotr(where((r gt Rhor)and(r lt 2*Rhor))))

        Mdotr(where(r gt R_attach)) = Mdot0
        if ((aa eq 0.9)or(aa eq 0.99)) then mdotr = smooth(mdotr,21)
;stop
        rdatafile = 'data/ph_0000.dat'
        dumpstr = string(runid,format='(I4.4)')
        strput,rdatafile,dumpstr,8
        openr,1,rdatafile
        readf,1,blnk,blnk,Tdisk,blnk,blnk,blnk,blnk
        close,1

        llr = fltarr(N1)
        llr(*) = 2.*sigma_SB*Tdisk(*,0)^4.
        llr = llr*(mdot0/mdotr)

        rdatafile = 'data/scat_spec.0000.dat'
        dumpstr = string(spec_runid,format='(I4.4)')
        strput,rdatafile,dumpstr,15
        openr,1,rdatafile
        readf,1,spec,Qspec,Uspec
        for isort=0,5 do begin
            readf,1,spec2,Qspec2,Uspec2
            spec_s(isort,*,*)=spec2
        endfor
        close,1
        rdatafile = 'data/scat_spcr.0000.dat'
        dumpstr = string(spec_runid,format='(I4.4)')
        strput,rdatafile,dumpstr,15
        openr,1,rdatafile
        readf,1,Ispecr,rdata,Qspecr,Uspecr,Rspecr,Cspecr
        close,1
        
        nu = fltarr(Nspec)
        dnu = nu
        e_min = 0.001
        e_max = 1000.0
        nu_min = e_min*2.
        nu_max = e_max/2.
        nu(0) = e_min
        for j=1.,Nspec-1 do begin
            nu(j) = e_min*10.0^(j/(Nspec-1.)*alog10(e_max/e_min))
            dnu(j-1) = nu(j)-nu(j-1)
        endfor
        dnu = deriv(nu)
        dLdR = dblarr(Nr)
        dLdR_scale = dblarr(Nr)
        dLCdR = dblarr(Nr)
        for ir=0,Nr-1 do begin
            for i_dex = 0,Nth-1 do begin
                dLdR(ir)=dLdr(ir)+total(Ispecr(*,i_dex,ir)*dnu)/drr(ir)
;                dLdR_scale(ir)=dLdr_scale(ir)+total(Ispecr(*,i_dex,ir)*dnu)/drr(ir)*(mdot0/mdotr(ir))
            endfor
            dLCdR(ir)=total(Cspecr(*,ir)*dnu)/drr(ir)
        endfor
        ;scale by mdot and Ledd 0.3 -> 0.1
;        if (aa gt 0.5) then mdotr(*)=mdot0 ;no r-scaling
        dLdR = dLdR*(mdot0/mdotr)/3.
        dLCdR = dLCdR*(mdot0/mdotr)/3.
        llr = llr/3.
        llr = (dldr+dlcdr)*2.4d17
        ll_100 = llr(where((r-100.)^2 eq min((r-100.)^2.)))
        ll_100 = ll_100(0)
stop
        x_title = 'r/M'
        y_title = 'dL/dr'
        if ((aa eq 0)or(aa eq 0.5)) then x_title=''
        if ((aa eq 0.5)or(aa eq 0.99)) then y_title=''
        plot_oo,r,dldr*2.4d17,xrange=[3,100],yrange=[1d35,3d37],$
          xtitle=x_title,ytitle = y_title,ystyle=1,linestyle=2,thick=thk
        oplot,r,dldr*2.4d17,thick=thk,linestyle=2,color=120
        oplot,r,dlcdr*2.4d17,linestyle=1,thick=thk,color=120
        ;oplot,r,llr*r*2.*!PI*(1.45d6)^2.,thick=thk
        oplot,r,llr,thick=thk,color=120

        runid = runids(iaa)
        spec_runid = spec_runids(iaa)
        rdatafile = 'data/gr_0000.dat'
        dumpstr = string(ntrunid,format='(I4.4)')
        strput,rdatafile,dumpstr,8
        openr,1,rdatafile
        readf,1,N1,N2,N3
        r = fltarr(N1)
        t = fltarr(N2)
        p = fltarr(N3)
        readf,1,r,t,p
        close,1
        Nr = N1
        drr = deriv(r)
        Tdisk = fltarr(N1,N3)
        blnk = fltarr(N1,N3)

        rdata = fltarr(3,Nr)
        Rspec = fltarr(Nspec,Nr)
        Ispecr = fltarr(Nspec,Nth,Nr)
        Qspecr = fltarr(Nspec,Nth,Nr)
        Uspecr = fltarr(Nspec,Nth,Nr)
        Rspecr = fltarr(Nspec,Nth,Nr)
        Cspecr = fltarr(Nspec,Nr)

        rdatafile = 'data/scat_spcr.0000.dat'
        dumpstr = string(ntrunid,format='(I4.4)')
        strput,rdatafile,dumpstr,15
        openr,1,rdatafile
        readf,1,Ispecr,rdata,Qspecr,Uspecr,Rspecr,Cspecr
        close,1
        
        dLdR = dblarr(Nr)
        dLCdR = dblarr(Nr)
        for ir=0,Nr-1 do begin
            for i_dex = 0,Nth-1 do begin
                dLdR(ir)=dLdr(ir)+total(Ispecr(*,i_dex,ir)*dnu)/drr(ir)
            endfor
            dLCdR(ir)=total(Cspecr(*,ir)*dnu)/drr(ir)
        endfor

;        rdatafile = 'data/ph_0000.dat'
;        dumpstr = string(ntrunid,format='(I4.4)')
;        strput,rdatafile,dumpstr,8
;        openr,1,rdatafile
;        readf,1,blnk,blnk,Tdisk,blnk,blnk,blnk,blnk
;        close,1

        llr = (dldr+dlcdr)*2.4d17
        ll_100nt = llr(where((r-100.)^2 eq min((r-100.)^2.)))
        ll_100nt = ll_100nt(0)
        llr = llr/ll_100nt*ll_100
        llr(where(r lt Risco)) = 1.

;        oplot,r,llr*r*2.*!PI*(1.45d6)^2.,thick=thk,color=120
        oplot,r,llr,thick=thk
        if (aa eq 0) then xyouts,30,10^36.5,'a/M=0'
        if (aa eq 0.5) then xyouts,30,10^36.5,'a/M=0.5'
        if (aa eq 0.9) then xyouts,30,10^36.5,'a/M=0.9'
        if (aa eq 0.99) then xyouts,30,10^36.5,'a/M=0.99'

        stop
        
    endfor
endif

if (fix(ifig) eq 2) then begin
    runids = [8151,8251,8352,8452]
    spec_runids = [8199,8299,8398,8498]
;    spec_runids = [8151,8251,8351,8451]
    aa_ = [0.0,0.5,0.9,0.99]
    sigma_SB = 5.6705d-5
    N1 = 0
    N2 = 0
    N3 = 0
    Nth = 41
    Nspec = 101
    for iaa = 0,3 do begin
        aa = aa_(iaa)
        a2 = aa*aa
        Z1 = 1.+(1-aa*aa)^(1./3.)*((1.+aa)^(1./3.)+(1-aa)^(1./3.))
        Z2 = sqrt(3.*aa*aa+Z1^2.)
        Risco = (3.+Z2-((3.-Z1)*(3+Z1+2.*Z2))^0.5)
        Eisco = (Risco*Risco-2*Risco+aa*sqrt(Risco))/ $
          (Risco*sqrt(Risco*Risco-3*Risco+2*aa*sqrt(Risco))) 
        eta = 1.-Eisco
        Rhor = 1.+sqrt(1.-a2)

        R_attach = 15.
        if (aa eq 0.0) then begin
            mdotfile = 'mdot_tavg_00.dat'
            Nr2 = 912
            ntrunid = 6000
        endif
        if (aa eq 0.5) then begin
            mdotfile = 'mdot_tavg_05.dat'
            Nr2 = 960
            ntrunid = 6010
        endif
        if (aa eq 0.9) then begin
            mdotfile = 'new_mdot_tavg_09.dat'
;            Nr2 = 1020
            Nr2 = 204
            ntrunid = 6018
        endif
        if (aa eq 0.99) then begin
            mdotfile = 'new_mdot_tavg_99.dat'
;            Nr2 = 1056
            Nr2 = 210
            ntrunid = 6027
            R_attach = 10.
        endif
        if ((aa eq 0.0)or(aa eq 0.5)) then rdata2 = fltarr(4,Nr2)
        if ((aa eq 0.9)or(aa eq 0.99)) then rdata2 = fltarr(2,Nr2)
        openr,1,mdotfile
        readf,1,rdata2
        close,1

        runid = runids(iaa)
        spec_runid = spec_runids(iaa)
        rdatafile = 'data/gr_0000.dat'
        dumpstr = string(runid,format='(I4.4)')
        strput,rdatafile,dumpstr,8
        openr,1,rdatafile
        readf,1,N1,N2,N3
        r = fltarr(N1)
        t = fltarr(N2)
        p = fltarr(N3)
        readf,1,r,t,p
        close,1
        Nr = N1
        drr = deriv(r)

        Rspec = fltarr(Nspec,Nr)
        Ispec = fltarr(Nspec,Nth)
        Ispecr = fltarr(Nspec,Nth,Nr)
        Qspec = fltarr(Nspec,Nth)
        Qspecr = fltarr(Nspec,Nth,Nr)
        Uspec = fltarr(Nspec,Nth)
        Uspecr = fltarr(Nspec,Nth,Nr)
        Rspecr = fltarr(Nspec,Nth,Nr)
        Cspecr = fltarr(Nspec,Nr)
        spec = fltarr(Nspec)
        rdata = fltarr(3,Nr)
        spec_s = dblarr(6,Nspec,Nth)
        spec2 = dblarr(Nspec,Nth)
        Qspec2 = dblarr(Nspec,Nth)
        Uspec2 = dblarr(Nspec,Nth)

        Mdotr = fltarr(N1)
        for ir=0,Nr2-1 do begin
            if (rdata2(0,ir) gt r(0)) then begin
                inrr = where(r le rdata2(0,ir))
                rdex = max(inrr)
                Mdotr(rdex) = rdata2(1,ir)
            endif
        endfor
        if ((aa eq 0.9)or(aa eq 0.99)) then Mdotr(0:Nr/2-1)=rdata2(1,*)
        Mdot0 = mean(Mdotr(where((r gt 2)and(r lt 6))))
        Mdot0 = mean(Mdotr(where((r gt Rhor)and(r lt 2*Rhor))))

        if ((aa eq 0.9)or(aa eq 0.99)) then mdotr = smooth(mdotr,21)
        Mdotr(where(r gt R_attach)) = Mdot0

        rdatafile = 'data/scat_spec.0000.dat'
        dumpstr = string(spec_runid,format='(I4.4)')
        strput,rdatafile,dumpstr,15
        openr,1,rdatafile
        readf,1,Ispec,Qspec,Uspec
        for isort=0,5 do begin
            readf,1,spec2,Qspec2,Uspec2
            spec_s(isort,*,*)=spec2
        endfor
        close,1
        rdatafile = 'data/scat_spcr.0000.dat'
        dumpstr = string(spec_runid,format='(I4.4)')
        strput,rdatafile,dumpstr,15
        openr,1,rdatafile
        readf,1,Ispecr,rdata,Qspecr,Uspecr,Rspecr,Cspecr
        close,1
        
        nu = fltarr(Nspec)
        dnu = nu
        e_min = 0.001
        e_max = 1000.0
        nu_min = e_min*2.
        nu_max = e_max/2.
        nu(0) = e_min
        for j=1.,Nspec-1 do begin
            nu(j) = e_min*10.0^(j/(Nspec-1.)*alog10(e_max/e_min))
            dnu(j-1) = nu(j)-nu(j-1)
        endfor
        dnu = deriv(nu)
        dnurr = dnu#(fltarr(Nr)+1)

        Ispec(*,*)=0
        Qspec(*,*)=0
        Uspec(*,*)=0
        for ir=0,Nr-1 do begin
            for i_dex = 0,Nth-1 do begin
                spec(*)=Ispecr(*,i_dex,ir)
                spec0 = spec
                scale_spec,spec,nu,dnu,nu_min,nu_max,1.,1.,mdot0/mdotr(ir)
                Ispecr(*,i_dex,ir)=spec
                Qspecr(*,i_dex,ir)=Qspecr(*,i_dex,ir)*spec/(spec0+1.)
                Uspecr(*,i_dex,ir)=Uspecr(*,i_dex,ir)*spec/(spec0+1.)
            endfor
            Ispec = Ispec+Ispecr(*,*,ir)
            Qspec = Qspec+Qspecr(*,*,ir)
            Uspec = Uspec+Uspecr(*,*,ir)
        endfor

        idex75 = 14
        deg_spec = sqrt(Qspec^2+Uspec^2)/Ispec
        ang_spec = atan(Uspec,Qspec)/2.*!radeg
        neg_ang = where(ang_spec lt -5)
        if (neg_ang(0) ge 0) then ang_spec(neg_ang)=ang_spec(neg_ang)+180.

        if (ifig eq 2.1) then begin
            if (iaa eq 0) then begin
                plot_oi,nu,deg_spec(*,idex75)*100.,$
                  xrange=[0.1,20.],yrange=[0,10],xstyle=1,ystyle=1,$
                  xtitle='E!Lobs!N (keV)',ytitle='polarization degree (%)',thick=thk
                plots,[0.15,0.3],[1,1]*7.5,thick=thk,color= 60
                plots,[0.15,0.3],[1,1]*7.0,thick=thk,color= 90
                plots,[0.15,0.3],[1,1]*6.5,thick=thk,color= 120
                plots,[0.15,0.3],[1,1]*6.0,thick=thk,color= 150
                xyouts,0.32,7.4,'a/M=0'
                xyouts,0.32,6.9,'     0.5'
                xyouts,0.32,6.4,'     0.9'
                xyouts,0.32,5.9,'   0.99'
            endif
            oplot,nu,deg_spec(*,idex75)*100,thick=thk,color=60+iaa*30
        endif
        min_deg = 20
        if (ifig eq 2.2) then begin
            if (iaa eq 0) then begin
                plot_oi,nu,smooth(ang_spec(*,idex75),1),$
                  xrange=[0.1,20.],yrange=[-20,150],xstyle=1,ystyle=1,$
                  xtitle='E!Lobs!N (keV)',ytitle='polarization angle (!Uo!N)'
                plots,[0.15,0.3],[1,1]*60+min_deg,thick=thk,color=60
                plots,[0.15,0.3],[1,1]*50+min_deg,thick=thk,color=90
                plots,[0.15,0.3],[1,1]*40+min_deg,thick=thk,color=120
                plots,[0.15,0.3],[1,1]*30+min_deg,thick=thk,color=150
                xyouts,0.32,58+min_deg,'a/M=0'
                xyouts,0.32,48+min_deg,'     0.5'
                xyouts,0.32,38+min_deg,'     0.9'
                xyouts,0.32,28+min_deg,'   0.99'
            endif
            oplot,nu,smooth(ang_spec(*,idex75),1),thick=thk,color=60+iaa*30
        endif

;        stop
    endfor
endif

if ((ifig ge 3)and(ifig lt 5)) then begin

    sigma_SB = 5.6705d-5
    N1 = 0
    N2 = 0
    N3 = 0
    Nth = 41
    Nspec = 101
    R_attach = 15.
    nteff0 = 0.06
    if (ifig eq 3.1) then begin
        runid = 8151
        aa = 0.
        Nr2 = 912
        mdotfile = 'mdot_tavg_00.dat'
        phrunid = 6050
    endif
    if (ifig eq 3.3) then begin
        runid = 8251
        aa = 0.5
        Nr2 = 960
        nteff = 0.09
        mdotfile = 'mdot_tavg_05.dat'
        phrunid = 6060
    endif
    if (ifig eq 3.2) then begin
        runid = 8352
        aa = 0.9
;        Nr2 = 1020
        Nr2 = 204
        nteff = 0.163
;        mdotfile = 'mdot_tavg_09.dat'
        mdotfile = 'new_mdot_tavg_09.dat'
        phrunid = 6068
    endif
    if (ifig eq 4.1) then begin
        runid = 8151
        aa = 0.
        Nr2 = 912
        mdotfile = 'mdot_tavg_00.dat'
        nteff = 0.06
        ntrunid = 6000
        phrunid = 6050
    endif
    if (ifig eq 4.2) then begin
        runid = 8251
        aa = 0.5
        Nr2 = 960
        nteff = 0.086
        mdotfile = 'mdot_tavg_05.dat'
        ntrunid = 6010
        phrunid = 6060
    endif
    if (ifig eq 4.3) then begin
        runid = 8351
        aa = 0.9
        Nr2 = 1020
        nteff = 0.163
        mdotfile = 'mdot_tavg_09.dat'
        ntrunid = 6018
        phrunid = 6068
    endif
    if (ifig eq 4.4) then begin
        runid = 8451
        aa = 0.99
        Nr2 = 1056
        nteff = 0.277
        mdotfile = 'mdot_tavg_99.dat'
        ntrunid = 6027
        phrunid = 6077
    endif
    a2 = aa*aa
    Rhor = 1.+sqrt(1.-a2)
    rdata2 = fltarr(4,Nr2)
    if (ifig eq 3.2) then rdata2 = fltarr(2,Nr2)
    openr,1,mdotfile
    readf,1,rdata2
    close,1

    spec_runid = runid
    rdatafile = 'data/gr_0000.dat'
    dumpstr = string(runid,format='(I4.4)')
    strput,rdatafile,dumpstr,8
    openr,1,rdatafile
    readf,1,N1,N2,N3
    r = fltarr(N1)
    t = fltarr(N2)
    p = fltarr(N3)
    readf,1,r,t,p
    close,1
    Nr = N1
    drr = deriv(r)

    Rspec = fltarr(Nspec,Nr)
    Ispec = fltarr(Nspec,Nth)
    Ispecr = fltarr(Nspec,Nth,Nr)
    Qspec = fltarr(Nspec,Nth)
    Uspec = fltarr(Nspec,Nth)
    spec = fltarr(Nspec)
    rdata = fltarr(3,Nr)
    spec_s = dblarr(6,Nspec,Nth)
    spec2 = dblarr(Nspec,Nth)
    Qspec2 = dblarr(Nspec,Nth)
    Uspec2 = dblarr(Nspec,Nth)
    
    Mdotr = fltarr(N1)
    for ir=0,Nr2-1 do begin
        if (rdata2(0,ir) gt r(0)) then begin
            inrr = where(r le rdata2(0,ir))
            rdex = max(inrr)
            Mdotr(rdex) = rdata2(1,ir)
        endif
    endfor
    if (ifig eq 3.2) then Mdotr(0:Nr2-1)=rdata2(1,*)
    Mdot0 = mean(Mdotr(where((r gt 2)and(r lt 6))))
    Mdot0 = mean(Mdotr(where((r gt Rhor)and(r lt 2*Rhor))))
    if (ifig eq 3.2) then mdotr = smooth(mdotr,21)
    Mdotr(where(r gt R_attach)) = Mdot0
    
    rdatafile = 'data/scat_spec.0000.dat'
    dumpstr = string(spec_runid,format='(I4.4)')
    strput,rdatafile,dumpstr,15
    openr,1,rdatafile
    readf,1,Ispec,Qspec,Uspec
    for isort=0,5 do begin
        readf,1,spec2,Qspec2,Uspec2
        spec_s(isort,*,*)=spec2
    endfor
    close,1
    rdatafile = 'data/scat_spcr.0000.dat'
    dumpstr = string(spec_runid,format='(I4.4)')
    strput,rdatafile,dumpstr,15
    openr,1,rdatafile
    readf,1,Ispecr
    close,1
    
    nu = fltarr(Nspec)
    dnu = nu
    e_min = 0.001
    e_max = 1000.0
    nu_min = e_min*2.
    nu_max = e_max/2.
    nu(0) = e_min
    for j=1.,Nspec-1 do begin
        nu(j) = e_min*10.0^(j/(Nspec-1.)*alog10(e_max/e_min))
        dnu(j-1) = nu(j)-nu(j-1)
    endfor
    dnu = deriv(nu)
    dnurr = dnu#(fltarr(Nr)+1)

    Ispec(*,*)=0
    Qspec(*,*)=0
    Uspec(*,*)=0
    for ir=0,Nr-1 do begin
        for i_dex = 0,Nth-1 do begin
            spec(*)=Ispecr(*,i_dex,ir)
            spec0 = spec
            scale_spec,spec,nu,dnu,nu_min,nu_max,1.,1.,mdot0/mdotr(ir)
            Ispecr(*,i_dex,ir)=spec
        endfor
        Ispec = Ispec+Ispecr(*,*,ir)
    endfor

    Ispec = Ispec*2.42d17
    idex75 = 15
    idex75a = 9
    idex75b = 10
    idex60 = 10
    idex60a = 19
    idex60b = 20

    if (ifig eq 3.1) then begin
        ntids = [6000,6004,6008]
        plot_oo,nu,Ispec(*,idex60)*nu,$
          xrange=[1,20.],yrange=[1d35,1d37],xstyle=1,ystyle=1,$
          xtitle='E!Lobs!N (keV)',ytitle='E F!LE!N',thick=thk
        oplot,nu,Ispec(*,idex60)*nu,color=120,thick=thk 
        plots,[1.15,1.5],[1,1]*10.^(36.15),thick=thk,color= 120
        plots,[1.15,1.5],[1,1]*10.^(36.0),thick=thk,linestyle=1
        plots,[1.15,1.5],[1,1]*10.^(35.85),thick=2,linestyle=2
        plots,[1.15,1.5],[1,1]*10.^(35.7),thick=2,linestyle=2
        plots,[1.15,1.5],[1,1]*10.^(35.55),thick=2,linestyle=2
        xyouts,1.52,10.^(36.1),'Harm3d, a/M=0.0'
        xyouts,1.52,10.^(35.95),'SBPL, a/M=0.0'
        xyouts,1.52,10.^(35.8),'NT, a/M=0.4'
        xyouts,1.52,10.^(35.65),'NT, a/M=0.2'
        xyouts,1.52,10.^(35.5),'NT, a/M=0.0'
        for idex = 0,2 do begin
            runid = ntids(idex)
            rdatafile = 'data/scat_spec.0000.dat'
            dumpstr = string(runid,format='(I4.4)')
            strput,rdatafile,dumpstr,15
            openr,1,rdatafile
            readf,1,Ispec
            close,1
;scale up to 0.3 L_edd
            spec(*)=0.5*(Ispec(*,idex60b)+Ispec(*,idex60a))*2.42d17
            scale_spec,spec,nu,dnu,nu_min,nu_max,1.,1.,3.
            oplot,nu,spec*nu,linestyle=2,thick=2
        endfor
        rdatafile = 'data/scat_spec.0000.dat'
        dumpstr = string(phrunid,format='(I4.4)')
        strput,rdatafile,dumpstr,15
        openr,1,rdatafile
        readf,1,Ispec
        close,1
;scale up to 0.3 L_edd
        spec(*)=0.5*(Ispec(*,idex60b)+Ispec(*,idex60a))*2.42d17
        scale_spec,spec,nu,dnu,nu_min,nu_max,1.,1.,2.55
        oplot,nu,spec*nu,linestyle=1,thick=thk
    endif

    if (ifig eq 3.2) then begin
        ntids = [6018,6023,6026]
        plot_oo,nu,Ispec(*,idex60)*nu,$
          xrange=[0.1,30.],yrange=[1d34,1d37],xstyle=1,ystyle=1,$
          xtitle='E!Lobs!N (keV)',ytitle='E F!LE!N',thick=thk
        oplot,nu,Ispec(*,idex60)*nu,color=120,thick=thk
        plots,[0.25,0.5],[1,1]*10.^(35.55),thick=thk,color= 120
        plots,[0.25,0.5],[1,1]*10.^(35.3),thick=thk,linestyle=1
        plots,[0.25,0.5],[1,1]*10.^(35.05),thick=thk,linestyle=2
        plots,[0.25,0.5],[1,1]*10.^(34.8),thick=thk,linestyle=2
        plots,[0.25,0.5],[1,1]*10.^(34.55),thick=thk,linestyle=2
        xyouts,0.52,10.^(35.5),'Harm3d, a/M=0.9'
        xyouts,0.52,10.^(35.25),'SBPL, a/M=0.9'
        xyouts,0.52,10.^(35.00),'NT, a/M=0.98'
        xyouts,0.52,10.^(34.75),'NT, a/M=0.95'
        xyouts,0.52,10.^(34.5),'NT, a/M=0.9'
        for idex = 0,2 do begin
            runid = ntids(idex)
            rdatafile = 'data/scat_spec.0000.dat'
            dumpstr = string(runid,format='(I4.4)')
            strput,rdatafile,dumpstr,15
            openr,1,rdatafile
            readf,1,Ispec
            close,1
;scale up to 0.3 L_edd
            spec(*)=0.5*(Ispec(*,idex60b)+Ispec(*,idex60a))*2.42d17
            scale_spec,spec,nu,dnu,nu_min,nu_max,1.,1.,3./(nteff/nteff0)
            oplot,nu,spec*nu,linestyle=2,thick=thk
        endfor
        rdatafile = 'data/scat_spec.0000.dat'
        dumpstr = string(phrunid,format='(I4.4)')
        strput,rdatafile,dumpstr,15
        openr,1,rdatafile
        readf,1,Ispec
        close,1
;scale up to 0.3 L_edd
        spec(*)=0.5*(Ispec(*,idex60b)+Ispec(*,idex60a))*2.42d17
        scale_spec,spec,nu,dnu,nu_min,nu_max,1.,1.,3.
        oplot,nu,spec*nu,linestyle=1,thick=thk
    endif
    if (ifig eq 3.3) then begin
        ntids = [6010,6012,6014]
        plot_oo,nu,Ispec(*,idex60)*nu,$
          xrange=[1,20.],yrange=[1d35,1d37],xstyle=1,ystyle=1,$
          xtitle='E!Lobs!N (keV)',ytitle='E F!LE!N',thick=thk
        oplot,nu,Ispec(*,idex60)*nu,color=120,thick=thk
        plots,[1.15,1.5],[1,1]*10.^(36.15),thick=thk,color= 120
        plots,[1.15,1.5],[1,1]*10.^(36.0),thick=thk,linestyle=1
        plots,[1.15,1.5],[1,1]*10.^(35.85),thick=2,linestyle=2
        plots,[1.15,1.5],[1,1]*10.^(35.7),thick=2,linestyle=2
        plots,[1.15,1.5],[1,1]*10.^(35.55),thick=2,linestyle=2
        xyouts,1.52,10.^(36.1),'Harm3d, a/M=0.5'
        xyouts,1.52,10.^(35.95),'SBPL, a/M=0.5'
        xyouts,1.52,10.^(35.8),'NT, a/M=0.7'
        xyouts,1.52,10.^(35.65),'NT, a/M=0.6'
        xyouts,1.52,10.^(35.5),'NT, a/M=0.5'
        for idex = 0,2 do begin
            runid = ntids(idex)
            rdatafile = 'data/scat_spec.0000.dat'
            dumpstr = string(runid,format='(I4.4)')
            strput,rdatafile,dumpstr,15
            openr,1,rdatafile
            readf,1,Ispec
            close,1
;scale up to 0.3 L_edd
            spec(*)=0.5*(Ispec(*,idex60b)+Ispec(*,idex60a))*2.42d17
            scale_spec,spec,nu,dnu,nu_min,nu_max,1.,1.,3./(nteff/nteff0)
            oplot,nu,spec*nu,linestyle=2,thick=2
        endfor
        rdatafile = 'data/scat_spec.0000.dat'
        dumpstr = string(phrunid,format='(I4.4)')
        strput,rdatafile,dumpstr,15
        openr,1,rdatafile
        readf,1,Ispec
        close,1
;scale up to 0.3 L_edd
        spec(*)=0.5*(Ispec(*,idex60b)+Ispec(*,idex60a))*2.42d17
        scale_spec,spec,nu,dnu,nu_min,nu_max,1.,1.,2.5
        oplot,nu,spec*nu,linestyle=1,thick=thk
    endif
stop
    if (fix(ifig) eq 4) then begin
        plot_oo,nu,Ispec(*,idex60)*nu,$
          xrange=[0.1,20.],yrange=[1d34,1d38],xstyle=1,ystyle=1,$
          xtitle='E!Lobs!N (keV)',ytitle='E F!LE!N',thick=thk
        oplot,nu,Ispec(*,idex60)*nu,color=120,thick=thk
        plots,[0.15,0.3],[1,1]*10.^(37.3),thick=thk,color= 120
        plots,[0.15,0.3],[1,1]*10.^(37.05),thick=thk,linestyle=2
        plots,[0.15,0.3],[1,1]*10.^(36.8),thick=thk,linestyle=1
        xyouts,0.32,10.^(37.25),'Harm3d'
        xyouts,0.32,10.^(37.00),'NT'
        xyouts,0.32,10.^(36.75),'SBPL'

        rdatafile = 'data/scat_spec.0000.dat'
        dumpstr = string(ntrunid,format='(I4.4)')
        strput,rdatafile,dumpstr,15
        openr,1,rdatafile
        readf,1,Ispec
        close,1
;scale up to 0.3 L_edd
        spec(*)=0.5*(Ispec(*,idex60b)+Ispec(*,idex60a))*2.42d17
        scale_spec,spec,nu,dnu,nu_min,nu_max,1.,1.,3./(nteff/nteff0)
        oplot,nu,spec*nu,linestyle=2,thick=thk

        rdatafile = 'data/scat_spec.0000.dat'
        dumpstr = string(phrunid,format='(I4.4)')
        strput,rdatafile,dumpstr,15
        openr,1,rdatafile
        readf,1,Ispec
        close,1
;scale up to 0.3 L_edd
        spec(*)=0.5*(Ispec(*,idex60b)+Ispec(*,idex60a))*2.42d17
        scale_spec,spec,nu,dnu,nu_min,nu_max,1.,1.,3.
        oplot,nu,spec*nu,linestyle=1,thick=thk

    endif

endif

if (ifig eq 5) then begin
    runids = [8151,8251,8352,8452]
    aa_ = [0.0,0.5,0.9,0.99]
    sigma_SB = 5.6705d-5
    cc = 3d10
    N1 = 0
    N2 = 0
    N3 = 0
    for iaa = 0,3 do begin
        aa = aa_(iaa)
        a2 = aa*aa
        Z1 = 1.+(1-aa*aa)^(1./3.)*((1.+aa)^(1./3.)+(1-aa)^(1./3.))
        Z2 = sqrt(3.*aa*aa+Z1^2.)
        Risco = (3.+Z2-((3.-Z1)*(3+Z1+2.*Z2))^0.5)
        Eisco = (Risco*Risco-2*Risco+aa*sqrt(Risco))/ $
          (Risco*sqrt(Risco*Risco-3*Risco+2*aa*sqrt(Risco))) 
        eta = 1.-Eisco
        Rhor = 1.+sqrt(1.-a2)

        R_attach = 15.
        if (aa eq 0.0) then begin
            mdotfile = 'mdot_tavg_00.dat'
            Nr2 = 912
        endif
        if (aa eq 0.5) then begin
            mdotfile = 'mdot_tavg_05.dat'
            Nr2 = 960
        endif
        if (aa eq 0.9) then begin
            mdotfile = 'new_mdot_tavg_09.dat'
;            Nr2 = 1020
            Nr2 = 204
        endif
        if (aa eq 0.99) then begin
            mdotfile = 'new_mdot_tavg_99.dat'
;            Nr2 = 1056
            Nr2 = 210
            R_attach = 10.
        endif
        if ((aa eq 0.0)or(aa eq 0.5)) then rdata2 = fltarr(4,Nr2)
        if ((aa eq 0.9)or(aa eq 0.99)) then rdata2 = fltarr(2,Nr2)
        openr,1,mdotfile
        readf,1,rdata2
        close,1

        runid = runids(iaa)
        rdatafile = 'data/gr_0000.dat'
        dumpstr = string(runid,format='(I4.4)')
        strput,rdatafile,dumpstr,8
        openr,1,rdatafile
        readf,1,N1,N2,N3
        r = fltarr(N1)
        t = fltarr(N2)
        p = fltarr(N3)
        readf,1,r,t,p
        close,1
        Tdisk = fltarr(N1,N3)
        blnk = fltarr(N1,N3)

        Mdotr = fltarr(N1)
        for ir=0,Nr2-1 do begin
            if (rdata2(0,ir) gt r(0)) then begin
                inrr = where(r le rdata2(0,ir))
                rdex = max(inrr)
                Mdotr(rdex) = rdata2(1,ir)
            endif
        endfor
        if ((aa eq 0.9)or(aa eq 0.99)) then Mdotr(0:Nr2-1)=rdata2(1,0:Nr2-1)
        Mdot0 = mean(Mdotr(where((r gt 2)and(r lt 6))))
        Mdot0 = mean(Mdotr(where((r gt Rhor)and(r lt 2*Rhor))))

        Mdotr(where(r gt R_attach)) = Mdot0
        if ((aa eq 0.9)or(aa eq 0.99)) then mdotr = smooth(mdotr,21)

        drr = deriv(r)
        rrp = r
        Rmatch = Risco*1.5-0.5*Rhor
        Rmatch = Risco*1.5
        if (aa eq 0) then Rmatch0 = Rmatch
;        innerdisk = where(r lt Risco*1.5)
        innerdisk = where(r lt Rmatch)
        innerir = max(innerdisk)
        for ir=innerir,0,-1 do $
          rrp(ir)=rrp(ir+1)-0.5*(drr(ir)+drr(ir+1))/ $
          sqrt(1.-4./(r(ir)+r(ir+1))+4.*aa*aa/(r(ir)+r(ir+1))^2.)
        for ir=innerir,N1-1 do $
          rrp(ir)=rrp(ir-1)+0.5*(drr(ir)+drr(ir-1))/ $
          sqrt(1.-4./(r(ir)+r(ir-1))+4.*aa*aa/(r(ir)+r(ir-1))^2.)
        irout = min(where(r gt Rhor))

        rdatafile = 'data/ph_0000.dat'
        dumpstr = string(runid,format='(I4.4)')
        strput,rdatafile,dumpstr,8
        openr,1,rdatafile
        readf,1,blnk,blnk,Tdisk,blnk,blnk,blnk,blnk
        close,1

        llr = fltarr(N1)
        llr(*) = 2.*sigma_SB*Tdisk(*,0)^4.
        llr = llr*(mdot0/mdotr)
        Tdr = (llr/(2.*sigma_SB))^0.25

        x0 = sqrt(Risco)                
        x1 = 2.*cos(acos(aa)/3.-!PI/3.)
        x2 = 2.*cos(acos(aa)/3.+!PI/3.)
        x3 = -2.*cos(acos(aa)/3.)
        x = sqrt(r)
        ntB = 1.+aa/r^1.5
        ntC = 1.-3./r+2.*aa/r^1.5
        M_BH = 2d34
        Mdot_cgs = 0.1*1.3d38*(M_BH/2d33)/(eta*cc^2.)
        Mstar = M_BH/(3.*2d33)
        Mdstar = Mdot_cgs/(1.0d17)

        ptQ = ((1.+aa/(x*x*x))/(x*sqrt(1.-3./(x*x)+2.*aa/(x*x*x))))* $
              (x-x0-3./2.*aa*alog(x/x0) $
               -3.*(x1-aa)^2./(x1*(x1-x2)*(x1-x3))*alog((x-x1)/(x0-x1)) $
               -3.*(x2-aa)^2./(x2*(x2-x1)*(x2-x3))*alog((x-x2)/(x0-x2)) $
               -3.*(x3-aa)^2./(x3*(x3-x1)*(x3-x2))*alog((x-x3)/(x0-x3)))
        ntQ = ptQ
        nt_flux = 5.4885d25*(Mdstar/Mstar^2.)*ntQ/(r^3.*ntB*sqrt(ntC))
        nt_flux(where(r lt Risco))=0.1
        Tnt = (nt_flux/(2.*sigma_SB))^0.25
        Tnt = Tnt*(Tdr(200)/Tnt(200))
        llr_nt = fltarr(N1)
        llr_nt(*) = 2.*sigma_SB*Tnt^4.
;stop
        Rpisco = rrp(max(where(r lt Risco*1.5)))
        Rpisco = rrp(max(where(r lt Rmatch)))
        Rfact = (Rmatch/Rmatch0)
        drp = deriv(rrp)
        rstar = rrp/Rpisco

        alpha = 1.7
        beta = -2.0
        Rn = 1.0
        R0 = 1.6
        delta = 0.35
        C = 1.e23

        alpha = 1.73
        beta = -2.0
        Rn = 1.0
        R0 = 1.68
        delta = 0.4*alog(10.)
        C = 0.8e23

        xi = (beta-alpha)/2.
        phi = (beta+alpha)/2.
;        Lstar2 = C*(rstar/Rn)^(phi-1.)* $
;                 (cosh(alog10(rstar/R0)/delta)/ $
;                  cosh(alog10(Rn/R0)/delta))^(xi*delta*alog(10.))
        Lstar2 = C*(rstar/Rn)^(phi-1.)*(cosh(alog(rstar/R0)/delta)/ $
                                        cosh(alog(Rn/R0)/delta))^(xi*delta)
        LLstar = Lstar2*(drp/drr)*(6./Risco)^2.
        Tstar = (LLstar/(2.*sigma_SB))^0.25
        rstar_01 = min(rstar(where(rstar gt 0.1)))
        Tstar_01 = Tstar(where(rstar eq rstar_01))
        Tstar(where(rstar lt rstar_01)) = Tstar_01
;        for i=0,100 do $
;           print,i,drr(i)/drp(i),rstar(i),Tstar(i)

        llr = llr*2.4d14
        llr_nt = llr_nt*2.4d14
        if (iaa eq 0) then $
           plot_oo,rrp/Rpisco,llr*(Risco/6.)^2.*(drr/drp)*(rrp/Rpisco),$
                   xrange=[0.1,200],yrange=[1d34,1d38],xstyle=1,xtitle='r*', ytitle = 'dL/dr*'
        oplot,rrp/Rpisco,llr*(Risco/6.)^2.*(drr/drp)*(rrp/Rpisco),color=60+iaa*30, $
              thick=thk
        if (iaa eq 0) then $
           oplot,rrp/Rpisco,llr_nt*(Risco/6.)^2.*(drr/drp)*(rrp/Rpisco),linestyle=1,thick=2
        if (iaa eq 2) then $
           oplot,rrp/Rpisco,llr_nt*(Risco/6.)^2.*(drr/drp)*(rrp/Rpisco),linestyle=3,thick=2
;        if (iaa eq 0) then $
;          plot_oo,rrp/Rpisco,llr*Rfact^2.*(drr/drp)*(rrp/Rmatch),$
;          xrange=[0.1,100],xtitle='R*', ytitle = 'dL/dR*'
;        oplot,rrp/Rpisco,llr*Rfact^2.*(drr/drp)*(rrp/Rmatch),color=60+iaa*30, $
;          thick=thk
;        plot_oi,r,Tdisk(*,0)*(mdot0/mdotr)^0.25,xrange=[2,100]
;        oplot,r,Tstar,linestyle=2
;        stop
    endfor
    rstar = rrp/Rpisco
    alpha = 1.7
    beta = -2.0
    Rn = 1.0
    R0 = 1.6
    C = 1.e23
    delta = 0.35

    alpha = 1.73
    beta = -2.0
    Rn = 1.0
    R0 = 1.68
    delta = 0.4*alog(10.)
    C = 0.8e23

    xi = (beta-alpha)/2.
    phi = (beta+alpha)/2.
;    Lstar2 = C*(rstar/Rn)^(phi-1.)*(cosh(alog10(rstar/R0)/delta)/ $
;                               cosh(alog10(Rn/R0)/delta))^(xi*delta*alog(10.))
    Lstar2 = C*(rstar/Rn)^(phi-1.)*(cosh(alog(rstar/R0)/delta)/ $
                               cosh(alog(Rn/R0)/delta))^(xi*delta)
;    oplot,rstar,Lstar*rstar,thick=thk
    Lstar2 = Lstar2*2.4d14
    oplot,rstar,Lstar2*rstar,thick=thk,linestyle=2
    !P.charsize=1.25
    plots,[10,25],[1,1]*10.^37.7,thick=thk,color=60
    plots,[10,25],[1,1]*10.^37.5,thick=thk,color=90
    plots,[10,25],[1,1]*10.^37.3,thick=thk,color=120
    plots,[10,25],[1,1]*10.^37.1,thick=thk,color=150
    plots,[10,25],[1,1]*10.^36.9,thick=thk,linestyle=2
    plots,[10,25],[1,1]*10.^36.7,thick=2,linestyle=1
    plots,[10,25],[1,1]*10.^36.5,thick=2,linestyle=3
    xyouts,27,10.^37.65,'a/M=0'
    xyouts,27,10.^37.45,'a/M=0.5'
    xyouts,27,10.^37.25,'a/M=0.9'
    xyouts,27,10.^37.05,'a/M=0.99'
    xyouts,27,10.^36.85,'SBPL fit'
    xyouts,27,10.^36.65,'NT (a/M=0)'
    xyouts,27,10.^36.45,'NT (a/M=0.9)'
    stop
endif

if (fix(ifig) eq 6) then begin
    runids = [8151,8251,8351,8451]
    aa_ = [0.0,0.5,0.9,0.99]
    sigma_SB = 5.6705d-5
    N1 = 0
    N2 = 0
    N3 = 0
    for iaa = 1,1 do begin

        aa = aa_(iaa)
        a2 = aa*aa
        Z1 = 1.+(1-aa*aa)^(1./3.)*((1.+aa)^(1./3.)+(1-aa)^(1./3.))
        Z2 = sqrt(3.*aa*aa+Z1^2.)
        Risco = (3.+Z2-((3.-Z1)*(3+Z1+2.*Z2))^0.5)
        Eisco = (Risco*Risco-2*Risco+aa*sqrt(Risco))/ $
          (Risco*sqrt(Risco*Risco-3*Risco+2*aa*sqrt(Risco))) 
        eta = 1.-Eisco
        Rhor = 1.+sqrt(1.-a2)

        runid = runids(iaa)
        rdatafile = 'data/gr_0000.dat'
        dumpstr = string(runid,format='(I4.4)')
        strput,rdatafile,dumpstr,8
        openr,1,rdatafile
        readf,1,N1,N2,N3
        r = fltarr(N1)
        t = fltarr(N2)
        p = fltarr(N3)
        readf,1,r,t,p
        close,1
        u0 = fltarr(N1,N2,N3)
        u1 = fltarr(N1,N2,N3)
        u2 = fltarr(N1,N2,N3)
        u3 = fltarr(N1,N2,N3)
        ll = fltarr(N1,N2,N3)  
        rr = fltarr(N1,N2,N3)
        tt = fltarr(N1,N2,N3)
        pp = fltarr(N1,N2,N3)
        dV = fltarr(N1,N2,N3)
        g_ = fltarr(N1,N2,4,4)
        g_up = fltarr(N1,N2,4,4)
        Tdisk = fltarr(N1,N3)
        blnk = fltarr(N1,N3)

        drr = deriv(r)
        rrp = r
        innerdisk = where(r lt Risco*1.5)
        innerir = max(innerdisk)
        for ir=innerir,0,-1 do $
          rrp(ir)=rrp(ir+1)-0.5*(drr(ir)+drr(ir+1))/ $
          sqrt(1.-4./(r(ir)+r(ir+1))+4.*aa*aa/(r(ir)+r(ir+1))^2.)
        for ir=innerir,N1-1 do $
          rrp(ir)=rrp(ir-1)+0.5*(drr(ir)+drr(ir-1))/ $
          sqrt(1.-4./(r(ir)+r(ir-1))+4.*aa*aa/(r(ir)+r(ir-1))^2.)
        irout = min(where(r gt Rhor))

        rdatafile = 'data/u0_0000.dat'
        dumpstr = string(runid,format='(I4.4)')
        strput,rdatafile,dumpstr,8
        openr,1,rdatafile
        readf,1,u0
        close,1
        rdatafile = 'data/u1_0000.dat'
        dumpstr = string(runid,format='(I4.4)')
        strput,rdatafile,dumpstr,8
        openr,1,rdatafile
        readf,1,u1
        close,1
        rdatafile = 'data/u2_0000.dat'
        dumpstr = string(runid,format='(I4.4)')
        strput,rdatafile,dumpstr,8
        openr,1,rdatafile
        readf,1,u2
        close,1
        rdatafile = 'data/u3_0000.dat'
        dumpstr = string(runid,format='(I4.4)')
        strput,rdatafile,dumpstr,8
        openr,1,rdatafile
        readf,1,u3
        close,1
        rdatafile = 'data/ll_0000.dat'
        dumpstr = string(runid,format='(I4.4)')
        strput,rdatafile,dumpstr,8
        openr,1,rdatafile
        readf,1,ll
        close,1

        u0r = fltarr(N1)
        u1r = fltarr(N1)
        u2r = fltarr(N1)
        u3r = fltarr(N1)
        llr = fltarr(N1)

        u0r(*) = 0.5*(u0(*,N2/2-1,0)+u0(*,N2/2,0))
        u1r(*) = 0.5*(u1(*,N2/2-1,0)+u1(*,N2/2,0))
        u2r(*) = 0.5*(u2(*,N2/2-1,0)+u2(*,N2/2,0))
        u3r(*) = 0.5*(u3(*,N2/2-1,0)+u3(*,N2/2,0))
        llr(*) = 2.*sigma_SB*Tdisk(*,0)^4.

        Rplunge = Risco
        Eplunge = (Rplunge*Rplunge-2*Rplunge+aa*sqrt(Rplunge))/ $
          (Rplunge*sqrt(Rplunge*Rplunge-3*Rplunge+2*aa*sqrt(Rplunge))) 
        Lplunge = sqrt(Rplunge)*(Rplunge*Rplunge-2.*aa*sqrt(Rplunge)+aa^2.)/ $
          (Rplunge*sqrt(Rplunge*Rplunge-3.*Rplunge+2.*aa*sqrt(Rplunge)))
        
        u0a = fltarr(N1)
        u1a = fltarr(N1)
        u2a = fltarr(N1)
        u3a = fltarr(N1)
        p_ = fltarr(4)
        v_ = fltarr(4)
        for ir = irout,N1-1 do begin
            x_ = [0.,r(ir),!PI/2.,0]
            calc_g,g_dn,g_up,x_,aa
            if (r(ir) ge Rplunge) then begin
                Ecirc = (r(ir)^2.-2.*r(ir)+aa*sqrt(r(ir)))/ $
                  (r(ir)*sqrt(r(ir)^2.-3.*r(ir)+2.*aa*sqrt(r(ir)))) 
                Lcirc = sqrt(r(ir))*(r(ir)^2.-2.*aa*sqrt(r(ir))+a2)/ $
		  (r(ir)*sqrt(r(ir)^2.-3.*r(ir)+2.*aa*sqrt(r(ir))))
		p_(0) = -Ecirc
		p_(1) = 0
		p_(2) = 0
                p_(3) = Lcirc
            endif
            if (r(ir) lt Rplunge) then begin
                Ecirc = Eplunge 
                Lcirc = Lplunge
		p_(0) = -Ecirc
		p_(1) = -sqrt(abs((1.+g_up(0,0)*Ecirc*Ecirc $
                                   -2.*g_up(0,3)*Ecirc*Lcirc+ $
                                   g_up(3,3)*Lcirc*Lcirc)/g_up(1,1)))
		p_(2) = 0
		p_(3) = Lcirc
            endif
;            print,ir,r(ir),Ecirc,Lcirc,p_
            v_(0) = g_up(0,0)*p_(0)+g_up(0,3)*p_(3)
            v_(1) = g_up(1,1)*p_(1)
            v_(2) = g_up(2,2)*p_(2)
            v_(3) = g_up(3,0)*p_(0)+g_up(3,3)*p_(3)
            u0a(ir) = v_(0)
            u1a(ir) = v_(1)
            u2a(ir) = v_(2)
            u3a(ir) = v_(3)
        endfor
        Rpisco = rrp(max(where(u1a lt 0)))

        if (ifig eq 6.1) then begin
            plot_oo,r-Rhor,u0a,xrange=[0.01,100],yrange=[1,100],$
              xtitle='R-R!Lhor!N (M)',ytitle = 'u!Ut!N',psym=4,thick=thk
            oplot,r-Rhor,u0r,color=120,thick=thk
            plots,[10.^0.,10.^0.6],[1,1]*10.^(1.8),color=120,thick=thk
            plots,[10.^0.,10.^0.2,10.^0.4,10.^0.6],[1,1,1,1]*10.^(1.6),psym=4,thick=thk
            xyouts,10.^0.7,10.^1.75,'Harm3d data'
            xyouts,10.^0.7,10.^1.55,'geodesics'
        endif
        if (ifig eq 6.2) then begin
            plot_oi,r-Rhor,u1a,xrange=[0.01,100],yrange=[-0.6,0.2],$
              xtitle='R-R!Lhor!N (M)',ytitle = 'u!Ur!N',psym=4,thick=thk
            oplot,r-Rhor,u1r,color=120,thick=thk
            plots,[10.^0.,10.^0.6],[1,1]*0.15,color=120,thick=thk
            plots,[10.^0.,10.^0.2,10^0.4,10.^0.6],[1,1,1,1]*0.08,psym=4,thick=thk
            xyouts,10.^0.7,0.13,'Harm3d data'
            xyouts,10.^0.7,0.06,'geodesics'
        endif
        if (ifig eq 6.3) then begin
            plot_oi,r-Rhor,u2a,xrange=[0.01,100],yrange=[-0.01,0.01],$
              xtitle='R-R!Lhor!N (M)',ytitle = 'u!U!9q!3!N',psym=4,thick=thk
            oplot,r-Rhor,u2r,color=120,thick=thk
            plots,[10.^0.,10.^0.6],[1,1]*0.008,color=120,thick=thk
            plots,[10.^0.,10.^0.2,10.^0.4,10.^0.6],[1,1,1,1]*0.006,psym=4,thick=thk
            xyouts,10.^0.7,0.0075,'Harm3d data'
            xyouts,10.^0.7,0.0055,'geodesics'
        endif
        if (ifig eq 6.4) then begin
            plot_oo,r-Rhor,u3a,xrange=[0.01,100],yrange=[0.001,10],$
              xtitle='R-R!Lhor!N (M)',ytitle = 'u!U!9f!3!N',psym=4,thick=thk
            oplot,r-Rhor,u3r,color=120,thick=thk
            plots,[10.^0.,10.^0.6],[1,1]*10.^(0.8),color=120,thick=thk
            plots,[10.^0.,10.^0.2,10.^0.4,10.^0.6],[1,1,1,1]*10.^(0.5),psym=4,thick=thk
            xyouts,10.^0.7,10.^0.75,'Harm3d data'
            xyouts,10.^0.7,10.^0.45,'geodesics'
        endif
    endfor
endif

if (ifig eq 7) then begin
    runids = [8151,8251,8351,8451]
    spec_runids = [8199,8299,8399,8499]
    aa_ = [0.0,0.5,0.9,0.99]
    sigma_SB = 5.6705d-5
    N1 = 0
    N2 = 0
    N3 = 0
    Nth = 41
    Nspec = 101
    for iaa = 0,3 do begin
        aa = aa_(iaa)
        a2 = aa*aa
        Z1 = 1.+(1-aa*aa)^(1./3.)*((1.+aa)^(1./3.)+(1-aa)^(1./3.))
        Z2 = sqrt(3.*aa*aa+Z1^2.)
        Risco = (3.+Z2-((3.-Z1)*(3+Z1+2.*Z2))^0.5)
        Eisco = (Risco*Risco-2*Risco+aa*sqrt(Risco))/ $
          (Risco*sqrt(Risco*Risco-3*Risco+2*aa*sqrt(Risco))) 
        eta = 1.-Eisco
        Rhor = 1.+sqrt(1.-a2)

        R_attach = 15.
        if (aa eq 0.0) then begin
            subid = [7100,7101,7102,7103,7104,7105,7106,7107,7108,7109,$
                     7110,7111,7112,7113,7114,7115,7116,7117,7118,7119,$
                     7120,7121,7122,7123,7124,7125,7126,7127,7128,7129,$
                     7130,7131,7132,7133,7134,7135,7136,7137,7138,7139,$
                     7140,7141,7142,7143,7144,7145,7146,7147,7148,7149,7150]
            Nframes = 51.
        endif
        if (aa eq 0.5) then begin
            subid = [7200,7201,7202,7203,7204,7205,7206,7207,7208,7209,$
                     7210,7211,7212,7213,7214,7215,7216,7217,7218,7219,$
                     7220,7221,7222,7223,7224,7225,7226,7227,7228,7229,$
                     7230,7231,7232,7233,7234,7235,7236,7237,7238,7239,$
                     7240,7241,7242,7243,7244,7245,7246,7247,7248,7249,7250]
            Nframes = 51.
        endif
        if (aa eq 0.9) then begin
            subid = [7300,7301,7302,7303,7304,7305,7306,7307,7308,7309,$
                     7310,7311,7312,7313,7314,7315,7316,7317,7318,7319,$
                     7320,7321,7322,7323,7324,7325,7326,7327,7328,7329,$
                     7330,7331,7332,7333,7334,7335,7336,7337,7338,7339,$
                     7340,7341,7342,7343,7344,7345,7346,7347,7348,7349,7350]
            Nframes = 51.
        endif
        if (aa eq 0.99) then begin
            subid = [7420,7421,7422,7423,7424,7425,7426,7427,7428,7429,$
                     7430,7431,7432,7433,7434,7435,7436,7437,7438,7439,$
                     7440,7441,7442,7443,7444,7445,7446,7447,7448,7449,7450]
            Nframes = 31.
        endif

        nu = fltarr(Nspec)
        dnu = nu
        e_min = 0.001
        e_max = 1000.0
        nu_min = e_min*2.
        nu_max = e_max/2.
        nu(0) = e_min
        for j=1.,Nspec-1 do begin
            nu(j) = e_min*10.0^(j/(Nspec-1.)*alog10(e_max/e_min))
            dnu(j-1) = nu(j)-nu(j-1)
        endfor
        dnu = deriv(nu)
        I_t = dblarr(Nframes)
        I_tt = dblarr(Nframes,Nth)
        sigma_th = dblarr(Nth)
        for it=0,Nframes-1 do begin
            runid = subid(it)
            spec_runid = spec_runids(iaa)
            Ispec = dblarr(Nspec,Nth)
            rdatafile = 'data/scat_spec.0000.dat'
            dumpstr = string(runid,format='(I4.4)')
            strput,rdatafile,dumpstr,15
            openr,1,rdatafile
            readf,1,Ispec
            close,1
            for ith=0,Nth-1 do $
              I_tt(it,ith)=total(Ispec(*,ith)*dnu)
            I_t(it) = total(I_tt(it,*))
        endfor
        xx = dindgen(Nframes)
        coeffs = poly_fit(xx,I_t,1)
        y_fit = coeffs(0)+coeffs(1)*xx;+coeffs(2)*xx^2.*0.
        for ith=0,Nth-1 do $
          I_tt(*,ith)=I_tt(*,ith)/y_fit(*)
        I_norm = total(I_tt,2)
        for ith=0,Nth-1 do $
          sigma_th(ith) = variance(I_tt(*,ith))
        stop
    endfor
endif

if (fix(ifig) eq 8) then begin
    Nspec = 101
    id0 = 7120
    Nt = 41
    Nph = 40
    Nframes = 1
    spec = dblarr(Nspec,Nph,Nt)
    spec2 = dblarr(Nspec,Nph,Nt)

    Nbands = 4
    b1lo = 0.1
    b1hi = 0.3
    b2lo = 0.3
    b2hi = 1.0
    b3lo = 1.0
    b3hi = 3.0
    b4lo = 3.
    b4hi = 10.

    amp = dblarr(Nbands,Nt)
    amp2 = dblarr(Nbands,Nt)

    lcurve = dblarr(Nbands,Nph)
    emin = 0.001              
    emax = 1000.
    nu = emin*10.^(findgen(Nspec)/(Nspec-1.)*alog10(emax/emin))
    dnu = deriv(nu)
    b1dex = where((nu gt b1lo)and(nu le b1hi))
    b2dex = where((nu gt b2lo)and(nu le b2hi))
    b3dex = where((nu gt b3lo)and(nu le b3hi))
    b4dex = where((nu gt b4lo)and(nu le b4hi))
    pdex1 = findgen(Nph)
    pdex2 = (pdex1+Nph/4.)mod Nph
    pdex3 = (pdex1+Nph*2/4.)mod Nph
    pdex4 = (pdex1+Nph*3/4.)mod Nph
    tframe = findgen(Nframes)*0.5+10.
    inc = acos(abs(1.-(findgen(Nt)+0.5)/Nt*2.))*!radeg
    Nframes = 1
    datapts = findgen(4,Nt,2*Nframes)
    rdatafile = 'data/scat_spcp.0000.dat'
    iframes = 0
    for run_id = id0,id0+Nframes-1,1 do begin
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,15
        openr,1,rdatafile
        readf,1,spec2
        close,1
        spec=spec2
        iframes=iframes+1
        for it = 1,Nt-2 do begin
            idex = [0,0]+it
            idex = [-1,0,1]+it
            ;idex = [-2,-1,0,1,2]+it
            for ip = 0,Nph-1 do begin
                lcurve(0,ip)=$
                  total(total(spec(b1dex,ip,idex),3))
                lcurve(1,ip)=$
                  total(total(spec(b2dex,ip,idex),3))
                lcurve(2,ip)=$
                  total(total(spec(b3dex,ip,idex),3))
                lcurve(3,ip)=$
                  total(total(spec(b4dex,ip,idex),3))
            endfor
            lcurve(*,*)=0.25*(lcurve(*,pdex1)+lcurve(*,pdex2)+$
                              lcurve(*,pdex3)+lcurve(*,pdex4))
            plot,lcurve(2,*)
            wait,1.0
;            for iband=0,Nbands-1 do lcurve(iband,*)=smooth(lcurve(iband,*),1)
            for iband=0,Nbands-1 do begin
                if (it lt Nt/2) then $
                  datapts(iband,it,iframes*2-2)= $
                  sqrt((variance(lcurve(iband,*)))/mean(lcurve(iband,*))^2.)
                if (it gt Nt/2) then $
                  datapts(iband,Nt-it-1,iframes*2-1)= $
                  sqrt((variance(lcurve(iband,*)))/mean(lcurve(iband,*))^2.)
                amp2(iband,it)=amp2(iband,it)+$
                  (variance(lcurve(iband,*)))/mean(lcurve(iband,*))^2.
            endfor
            print,it,amp2(2,it)
        endfor
    endfor
    datamean = fltarr(4,Nt/2)
    datasig = fltarr(4,Nt/2)
    for iband = 0,Nbands-1 do begin
        for it=0,Nt/2-1 do begin
            datamean(iband,it)=mean(datapts(iband,it,*))
            datasig(iband,it)=sqrt(variance(datapts(iband,it,*)))
            if (datasig(iband,it) gt datamean(iband,it)) then $
              datasig(iband,it)=0.7*datamean(iband,it)
        endfor
    endfor

    amp = sqrt(amp2/Nframes)
    tframe = findgen(Nph)/Nph
    inc = inc(0:Nt/2-1)
    odex = findgen(Nt/4)*2
    edex = odex+1
    plot_io,inc(edex),datamean(1,edex)*100,yrange=[0.1,100],$
      psym=2,xrange=[10,80],xstyle=1,xtitle='inclination (deg)',ytitle='RMS (%)'
    oplot,inc(edex),datamean(0,edex)*100,psym=2,color=60,thick=thk
    oplot,inc(edex),datamean(1,edex)*100,psym=2,color=90,thick=thk
    oplot,inc(edex),datamean(2,edex)*100,psym=2,color=120,thick=thk
    oplot,inc(edex),datamean(3,edex)*100,psym=2,color=150,thick=thk
    for it = 1,Nt/2-5,2 do begin
        oplot,[1,1]*inc(it),datamean(0,it)*100+[-datasig(0,it),datasig(0,it)]*100,$
          thick=thk,color=60
        oplot,[1,1]*inc(it),datamean(1,it)*100+[-datasig(1,it),datasig(1,it)]*100,$
          thick=thk,color=90
        oplot,[1,1]*inc(it),datamean(2,it)*100+[-datasig(2,it),datasig(2,it)]*100,$
          thick=thk,color=120
        oplot,[1,1]*inc(it),datamean(3,it)*100+[-datasig(3,it),datasig(3,it)]*100,$
          thick=thk,color=150
        oplot,[inc(it-1),inc(it+1)],datamean(0,it)*100*[1,1],thick=thk,color=60
        oplot,[inc(it-1),inc(it+1)],datamean(1,it)*100*[1,1],thick=thk,color=90
        oplot,[inc(it-1),inc(it+1)],datamean(2,it)*100*[1,1],thick=thk,color=120
        oplot,[inc(it-1),inc(it+1)],datamean(3,it)*100*[1,1],thick=thk,color=150
    endfor
stop
    plot_io,inc,amp(1,*)*100,yrange=[0.1,100],psym=2,xrange=[20,80],$
      xstyle=1,xtitle='inclination',ytitle='rms (%)'
    oplot,inc,amp(0,*)*100,psym=2,color=60,thick=4
    oplot,inc,amp(1,*)*100,psym=2,color=90,thick=4
    oplot,inc,amp(2,*)*100,psym=2,color=120,thick=4
    oplot,inc,amp(3,*)*100,psym=2,color=150,thick=4
    
stop
endif


END

PRO calc_g, g_dn, g_up, x_,aa

  g_dn = fltarr(4,4)
  g_up = fltarr(4,4)
  r = x_(1)
  r2 = r*r
  a2 = aa*aa
  cth = cos(x_(2))
  cth2 = cth*cth
  sth2 = 1.-cth2
  Sig = r2+a2*cth2
  Del = r2-2*r+a2

  g_dn(0,0) = -(1.-2.*r/Sig)
  g_dn(1,1) = Sig/Del
  g_dn(2,2) = Sig
  g_dn(3,3) = (r2+a2+2.*a2*r*sth2/Sig)*sth2
  g_dn(0,3) = -2.*aa*r*sth2/Sig
  g_dn(3,0) = g_dn(0,3)

  g_up(0,0) = -((r2+a2)*(r2+a2)-a2*Del*sth2)/(Sig*Del)
  g_up(1,1) = Del/Sig
  g_up(2,2) = 1./Sig
  if (sth2 ne 0) then $
    g_up(3,3) = (Del-a2*sth2)/(Sig*Del*sth2) else $
    g_up(3,3) = 0.
  g_up(0,3) = -2.*aa*r/(Sig*Del)
  g_up(3,0) = g_up(0,3)

END
