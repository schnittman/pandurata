;Collection of plotting scripts for generating figures from ApJ 769,
;156 (Schnittman, Krolik, & Noble 2013)

PRO figures_pan,ifig

!P.font = 0
!P.charsize = 1.5
thk = 5
if ((ifig eq 0.1)or(ifig eq 0.2)) then begin
    N1 = 0
    N2 = 0
    N3 = 0
    aa = 0.0
    t_frame = 0.

    run_id = 1250
    rdatafile = 'data/gr_0000.dat'
    dumpstr = string(run_id,format='(I4.4)')
    strput,rdatafile,dumpstr,8
    openr,1,rdatafile

    readf,1,N1,N2,N3
    r = fltarr(N1)
    t = fltarr(N2)
    t_ = fltarr(N2+1)
    p = fltarr(N3)
    readf,1,r,t,p
    close,1
    
    dr = deriv(r)
    dt = deriv(t)
    dp = deriv(p)
    t_(0:N2-1)=t-0.5*dt(0)
    t_(N2)=t(N2-1)+0.5*dt(0)

    rr = fltarr(N1,N2,N3)
    tt = fltarr(N1,N2,N3)
    pp = fltarr(N1,N2,N3)
    dV = fltarr(N1,N2,N3)
    g_ = fltarr(N1,N2,4,4)
    ll_cgs = fltarr(N1,N2,N3)
    ll_tmp = fltarr(N1,N2,N3)
    T_cor = fltarr(N1,N2,N3)
    diskbody_ijk = intarr(N1,N2,N3)

    for ir=0,N1-1 do rr(ir,*,*)=r(ir)
    for it=0,N2-1 do tt(*,it,*)=t(it)
    for ip=0,N3-1 do pp(*,*,ip)=p(ip)

    R_hor = 1.+sqrt(1.-aa*aa)
    Sigma = rr^2.+aa^2.*cos(tt)^2.
    Delta = rr^2.-2.*rr+aa^2.
    alpha = sqrt(Sigma*Delta/(Sigma*Delta+2.*rr*(aa^2.+rr^2.)))
    omega = 2.*rr*aa/(Sigma*Delta+2.*rr*(aa^2.+rr^2.))
    varph = sqrt((Sigma*Delta+2.*rr*(aa^2.+rr^2.))/Sigma*sin(tt)^2.)

    g_(*,*,0,0) = -alpha(*,*,0)^2.+omega(*,*,0)^2.*varph(*,*,0)^2.
    g_(*,*,0,3) = -omega(*,*,0)*varph(*,*,0)^2.
    g_(*,*,1,1) = Sigma(*,*,0)/Delta(*,*,0)
    g_(*,*,2,2) = Sigma(*,*,0)
    g_(*,*,3,3) = varph(*,*,0)^2.
    g_(*,*,3,0) = g_(*,*,0,3)

    detg = fltarr(N1,N2)
    detg(*,*) = g_(*,*,0,0)*g_(*,*,1,1)*g_(*,*,2,2)*g_(*,*,3,3) - $
      g_(*,*,0,3)*g_(*,*,1,1)*g_(*,*,2,2)*g_(*,*,0,3)
    detg2 = (r^4.)#(sin(t)^2.)
;    for ip=0,N3-1 do dV(*,*,ip) = sqrt(-detg)*(dr#dt)*dp(0)
    for ip=0,N3-1 do dV(*,*,ip) = sqrt(detg2)*(dr#dt)*dp(0)
;stop
    G_N = 6.6726d-8
    M_BH = 10.*(2.d33)
    cc = 3d10
    t_cgs = G_N*M_BH/cc^3.
    r_cgs = G_N*M_BH/cc^2.
    dV_cgs = dV*r_cgs^3.

    if (ifig eq 0.2) then begin
        rdatafile = 'data/ll_0000.dat'
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,8
        openr,1,rdatafile
        readf,1,ll_cgs
        close,1
        rdatafile = 'data/te_0000.dat'
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,8
        openr,1,rdatafile
        readf,1,T_cor
        close,1
        rdatafile = 'data/db_0000.dat'
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,8
        openr,1,rdatafile
        readf,1,diskbody_ijk
        close,1
        incor = where((diskbody_ijk eq 0)and(rr gt 2))
        L_max = max (ll_cgs(incor))
        L_min = 1d-12*L_max
        T_max = max (T_cor(incor))
        T_min = 1d-6*T_max
        Nlbin = 101.
        llbin = L_min*10.^(findgen(Nlbin)/(Nlbin-1.)*alog10(L_max/L_min))
        dllbin = deriv(llbin)
        Tbin = T_min*10.^(findgen(Nlbin)/(Nlbin-1.)*alog10(T_max/T_min))
        dTbin = deriv(Tbin)
        Vlbin = fltarr(Nlbin)
        VTbin = fltarr(Nlbin)
        for ir=0,N1-1 do begin
            for ith=0,N2-1 do begin
                for iph=0,N3-1 do begin
                    if ((rr(ir,ith,iph) gt 2)and(diskbody_ijk(ir,ith,iph) eq 0)) then begin
                        jL = fix(Nlbin*alog10(ll_cgs(ir,ith,iph)/L_min)/alog10(L_max/L_min))
                        if (jL lt 0) then jL = 0
                        if (jL gt Nlbin-1) then jL=Nlbin-1
                        Vlbin(jL) += dV_cgs(ir,ith,iph)/dllbin(jL)
                        jL = fix(Nlbin*alog10(T_cor(ir,ith,iph)/T_min)/alog10(T_max/T_min))
                        if (jL lt 0) then jL = 0
                        if (jL gt Nlbin-1) then jL=Nlbin-1
                        VTbin(jL) += dV_cgs(ir,ith,iph)/dTbin(jL)
                    endif
                endfor
            endfor
        endfor
        V_cum = fltarr(Nlbin)
        L_cum = fltarr(Nlbin)
        for ib=1,Nlbin-1 do begin
            V_cum(ib)=V_cum(ib-1)+Vlbin(ib)*dllbin(ib)
            L_cum(ib)=L_cum(ib-1)+Vlbin(ib)*llbin(ib)*dllbin(ib)
        endfor
        plot_oi,L_cum/max(L_cum),V_cum/max(V_cum),xrange=[0.001,1],yrange=[0.0,1],$
          xtitle='fraction of L!Lcor!N',ytitle='fraction of V!Lcor!N',thick=thk

        stop
    endif
        
    ll_cgs(*,*,*)=0.
    run_id_start = 1250
    run_id_stop = 1250
    Nframes = 0.
    for run_id = run_id_start,run_id_stop,50 do begin
        rdatafile = 'data/ll_0000.dat'
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,8
        openr,1,rdatafile
        readf,1,ll_tmp
        close,1
        ll_cgs = ll_cgs+ll_tmp
        Nframes = Nframes+1.
    endfor
    ll_cgs = ll_cgs/Nframes

    insim = where(rr gt R_hor)
    harm_tot = ll_cgs*0.
    harm_tot(insim)=ll_cgs(insim)*dv_cgs(insim)*4.
    dldrharm = total(total(harm_tot,3),2)

    plot_oo,r,dldrharm/dr*r,xrange=[2,60],yrange=[1d36,1d38],xstyle=1,$
      xtitle='r/M',ytitle='dL/d(log r)',thick=thk

    aa = 0.
    Nr = 1001.
    Mstar = 3.33
    Mdstar = 25.56
    Z1 = 1.+(1-aa*aa)^(1./3.)*((1.+aa)^(1./3.)+(1-aa)^(1./3.))
    Z2 = sqrt(3.*aa*aa+Z1^2.)
    Risco = (3.+Z2-((3.-Z1)*(3+Z1+2.*Z2))^0.5)
    Revnt = 1.001*(1.+sqrt(1.-aa*aa))
    Rout = 200.
    dr = (Rout-Risco)/(Nr-1.)
    rr = findgen(Nr)*dr+Risco
    r_cm = rr*Mstar*4.4d5

    ntA = 1.+aa^2./rr^2.+2.*aa^2./rr^3.
    ntB = 1.+aa/rr^1.5
    ntC = 1.-3./rr+2.*aa/rr^1.5
    ntD = 1.-2./rr+aa^2./rr^2.
    ntE = 1.+4.*aa^2./rr^2.-4.*aa^2./rr^3.+3.*aa^4./rr^4.
    ntE2 = 2.+5.*aa^2./rr^2.-2.*aa^2./rr^3.+3.*aa^4./rr^4.
    ntF = 1.-2.*aa/rr^1.5+aa^2./rr^2.
    ntJ = dblarr(Nr)
    Nu = 101.
    for i=0,Nr-1 do begin
        u_star = 1./rr(i)
        du = u_star/(Nu-1.)
        subu = findgen(Nu)*du
        ntJ(i) = exp(1.5*du*total((1.-2.*aa*subu^1.5+aa^2.*subu^2.)/ $
                                  ((1.+aa*subu^1.5)*(1.-3.*subu+2.*aa*subu^1.5))))
    endfor
    
    L_ms = 2./sqrt(3.)*(3.*sqrt(Risco)-2.*aa)/sqrt(Risco)
    ntL = ntF/sqrt(ntC)-L_ms/sqrt(rr)
    ntQ = dblarr(Nr)
    for i=0,Nr-1 do begin
        Q_int = ntL(0:i)*ntF(0:i)/(ntB(0:i)*ntC(0:i)*ntJ(0:i)*rr(0:i)^1.5)
        ntQ(i) = ntL(i)-3.*ntJ(i)/(2.*sqrt(rr(i)))* $
          dr*(total(Q_int(0:i))-0.5*(Q_int(0)+Q_int(i)))
    endfor
    nt_flux = 6.d25*(Mdstar/Mstar^2.)*ntQ/(rr^3.*ntB*sqrt(ntC))
    oplot,rr,nt_flux*4.*!PI*r_cm^2.,thick=thk,linestyle=2
    oplot,[2.5,4],[1,1]*10^37.8,thick=thk
    oplot,[2.5,4],[1,1]*10^37.6,thick=thk,linestyle=2
    xyouts,4.1,10^37.75,'Harm3d'
    xyouts,4.1,10^37.55,'N-T'
    stop
endif

if (fix(ifig) eq 1) then begin
    N1 = 0
    N2 = 0
    N3 = 0
    aa = 0.0
    t_frame = 0.
    phi_avg = 0.

    run_id = 1250
    rdatafile = 'data/gr_0000.dat'
    dumpstr = string(run_id,format='(I4.4)')
    strput,rdatafile,dumpstr,8
    openr,1,rdatafile

    readf,1,N1,N2,N3
    r = fltarr(N1)
    t = fltarr(N2)
    t_ = fltarr(N2+1)
    p = fltarr(N3)
    readf,1,r,t,p
    close,1
    
    dr = deriv(r)
    dt = deriv(t)
    dp = deriv(p)
    t_(0:N2-1)=t-0.5*dt(0)
    t_(N2)=t(N2-1)+0.5*dt(0)

    T_cor = fltarr(N1,N2,N3)
    rho = fltarr(N1,N2,N3)
    LL = fltarr(N1,N2,N3)
    rr = fltarr(N1,N2,N3)
    tt = fltarr(N1,N2,N3)
    tt_ = fltarr(N1,N2+1,N3)
    pp = fltarr(N1,N2,N3)
    dV = fltarr(N1,N2,N3)
    g_ = fltarr(N1,N2,4,4)

    for ir=0,N1-1 do rr(ir,*,*)=r(ir)
    for it=0,N2-1 do tt(*,it,*)=t(it)
    for it=0,N2 do tt_(*,it,*)=t_(it)
    for ip=0,N3-1 do pp(*,*,ip)=p(ip)

    rho_cgs = fltarr(N1,N2,N3)
    ll_cgs = fltarr(N1,N2,N3)
    bb_cgs = fltarr(N1,N2,N3)
    diskbody_ijk = intarr(N1,N2,N3)
    tau_es = fltarr(N1,N2+1,N3)

    sigtau_es = fltarr(N1,N3)
    netflux = fltarr(N1,N3)
    diskflux = fltarr(N1,N3)
    em_top = fltarr(N1,N3)
    em_bot = fltarr(N1,N3)
    ref_top = fltarr(N1,N3)
    ref_bot = fltarr(N1,N3)
    diskbody = fltarr(N1,N3)
    Tdisk = fltarr(N1,N3)

    rdatafile = 'data/ph_0000.dat'
    dumpstr = string(run_id,format='(I4.4)')
    strput,rdatafile,dumpstr,8
    openr,1,rdatafile
    readf,1,diskbody,sigtau_es,Tdisk,em_top,em_bot,ref_top,ref_bot
    close,1
    rdatafile = 'data/rh_0000.dat'
    dumpstr = string(run_id,format='(I4.4)')
    strput,rdatafile,dumpstr,8
    openr,1,rdatafile
    readf,1,rho_cgs
    close,1
    rdatafile = 'data/te_0000.dat'
    dumpstr = string(run_id,format='(I4.4)')
    strput,rdatafile,dumpstr,8
    openr,1,rdatafile
    readf,1,T_cor
    close,1
    rdatafile = 'data/ta_0000.dat'
    dumpstr = string(run_id,format='(I4.4)')
    strput,rdatafile,dumpstr,8
    openr,1,rdatafile
    readf,1,tau_es
    close,1
    rdatafile = 'data/db_0000.dat'
    dumpstr = string(run_id,format='(I4.4)')
    strput,rdatafile,dumpstr,8
    openr,1,rdatafile
    readf,1,diskbody_ijk
    close,1
    rdatafile = 'data/ll_0000.dat'
    dumpstr = string(run_id,format='(I4.4)')
    strput,rdatafile,dumpstr,8
    openr,1,rdatafile
    readf,1,ll_cgs
    close,1
    rdatafile = 'data/bb_0000.dat'
    dumpstr = string(run_id,format='(I4.4)')
    strput,rdatafile,dumpstr,8
    openr,1,rdatafile
    readf,1,bb_cgs
    close,1
    bb_cgs = bb_cgs/(8.*!PI) ;gives energy density U_b

    ;calculate radiation pressure in corona and disk
    Prad = ll_cgs*(9.11d-28*3d10^2)/(12.*1.38d-16*(T_cor+100))/(3d10*6.65d-25*rho_cgs/1.67d-24)
    Praddisk = 1./3.*4.*5.67d-5/3d10*T_cor^4.
    indisk = where(diskbody_ijk ge 1)
    incor = where(diskbody_ijk eq 0)
    Prad(indisk) = Praddisk(indisk)
    subi = 5
    Prad2 = Prad
    inhotcor = where((diskbody_ijk eq 0)and(Prad gt 0))
    incoldcor = where((diskbody_ijk eq 0)and(Prad eq 0))
    for ir=subi,N1-1-subi do begin
        for jth=subi,N2-1-subi do begin
            if ((diskbody_ijk(ir,jth,0) eq 0)and(Prad(ir,jth,0) eq 0)) then begin
                subP = Prad(ir-subi:ir+subi,jth-subi:jth+subi,0:subi)
                hotgas = where(subP gt 0)
                if (hotgas(0) ge 0) then Prad2(ir,jth,0)=mean(subP(hotgas))
                if (hotgas(0) lt 0) then Prad2(ir,jth,0)=mean(Prad(inhotcor))
            endif
        endfor
    endfor
;    Prad(incoldcor)=Prad2(incoldcor)
;stop

    if (phi_avg eq 1) then begin
       rho_cgs(*,*,0) = total(rho_cgs,3)/N3
       ll_cgs(*,*,0) = total(ll_cgs,3)/N3
       tau_es(*,*,0) = total(tau_es,3)/N3
       T_cor(*,*,0) = total(T_cor,3)/N3
       diskbody_ijk(*,*,0) = total(diskbody_ijk,3)/N3
       Tdisk(*,0) = total(Tdisk,2)/N3
       for ip = 1,N3-1 do begin
          rho_cgs(*,*,ip)=rho_cgs(*,*,0)
          ll_cgs(*,*,ip)=ll_cgs(*,*,0)
          tau_es(*,*,ip)=tau_es(*,*,0)
          T_cor(*,*,ip)=T_cor(*,*,0)
          diskbody_ijk(*,*,ip)=diskbody_ijk(*,*,0)
          Tdisk(*,ip)=Tdisk(*,0)
       endfor
    endif

    for ir=0,N1-1 do begin
        if (diskbody(ir,0) ge 1) then begin
            indisk = where(diskbody_ijk(ir,*,0) ge 0.5)
            T_cor(ir,indisk)=Tdisk(ir,0)
        endif
    endfor
    
    dA = (r*dr)#(fltarr(N3)+!PI/2./N3)
    for ir=0,N1-1 do rr(ir,*,*)=r(ir)
    for it=0,N2-1 do tt(*,it,*)=t(it)
    for ip=0,N3-1 do pp(*,*,ip)=p(ip)
    
;    rho_cgs(where(rho_cgs lt 1d-9))=1d-9
    rho_cgs(where(rho_cgs lt 1d-10))=1d-10
    if (ifig eq 1.1) then begin
        contour,rho_cgs(*,*,0),r#sin(t),r#cos(t),xrange=[0,60],yrange=[-30,30],$
          levels = 10.^(findgen(100)/99.*6.-11.),/isotropic,/fill,ystyle=1
;better for mdot=0.01:
;        contour,rho_cgs(*,*,0),r#sin(t),r#cos(t),xrange=[0,60],yrange=[-30,30],$
;          levels = 10.^(findgen(100)/99.*5.-10.),/isotropic,/fill,ystyle=1
        contour,tau_es(*,*,0),r#sin(t_),r#cos(t_),xrange=[0,60],yrange=[-30,30],$
          levels = [0.01,0.1,1],/isotropic,/noerase,thick=4,c_labels=[1,1,1],$
          xtitle = 'x/M', ytitle = 'z/M',ystyle=1,c_colors=[0,0,0]
        phi = (findgen(101)/100.-0.5)*!PI
        oplot,2.*cos(phi),2.*sin(phi),thick=5
    endif
    if (ifig eq 1.11) then begin
        Npx = 12600
        Nx = 100
        Ny = 5
        xx = findgen(Nx)
        yy = findgen(Ny)
        xscale = findgen(Nx)/(Nx-1.)*6.-11.
;        xscale = findgen(Nx)/(Nx-1.)*5.-10.
        xscale = 10.^xscale
        yscale = fltarr(Nx)+1
        movscl = fltarr(Ny,Nx)
        for i=0,Nx-1 do movscl(*,i)=xscale(i)
;    data = byte(255*(alog10(movscl/max(movscl)+1e-5)+5.01)/5.)
        enlarge = rebin(movscl,50,1000)
;      tvscl,enlarge,0,1000,/Device
        contour,enlarge,Position=[2600,6000,3100,12000],/Device,$
          yticks=1,yminor=1,ytickname=[' ',' '],$
          xticks=1,xminor=1,xtickname=[' ',' '],$
          levels = 10.^(findgen(100)/99.*6.-11.),/fill,/Noerase
;          levels = 10.^(findgen(100)/99.*5.-10.),/fill,/Noerase
        plot_oi,yscale,xscale,Position=[2600,6000,3100,12000],xticks=1,$
          xtickname=[' ',' '],yticks=3,yminor=1,$
          ytickname=['10!U-11!N','10!U-9!N','10!U-7!N','10!U-5!N'],$
          ytitle='!9r!3 (g/cm!U3!N)',/Noerase,/Device
;        plot_oi,yscale,xscale,Position=[2600,6000,3100,12000],xticks=1,$
;          xtickname=[' ',' '],yticks=5,yminor=1,$
;          ytickname=['10!U-10!N','10!U-9!N','10!U-8!N','10!U-7!N','10!U-6!N','10!U-5!N'],$
;          ytitle='!9r!3 (g/cm!U3!N)',/Noerase,/Device
        
    endif
    if (ifig eq 1.2) then begin
        ll_cgs(where(ll_cgs lt 1d10))=1d10
        contour,ll_cgs(*,*,0),r#sin(t),r#cos(t),xrange=[0,60],yrange=[-30,30],$
          levels = 10.^(findgen(100)/99.*6.+10.),/isotropic,/fill,ystyle=1
        contour,tau_es(*,*,0),r#sin(t_),r#cos(t_),xrange=[0,60],yrange=[-30,30],$
          levels = [0.01,0.1,1],/isotropic,/noerase,thick=4,c_labels=[1,1,1],$
          xtitle = 'x/M', ytitle = 'z/M',ystyle=1,c_colors=[255,255,255]
    endif
    if (ifig eq 1.21) then begin
        Npx = 12600
        Nx = 100
        Ny = 5
        xx = findgen(Nx)
        yy = findgen(Ny)
        xscale = findgen(Nx)/(Nx-1.)*6.+10.
        xscale = 10.^xscale
        yscale = fltarr(Nx)+1
        movscl = fltarr(Ny,Nx)
        for i=0,Nx-1 do movscl(*,i)=xscale(i)
        enlarge = rebin(movscl,50,1000)
        contour,enlarge,Position=[2600,6000,3100,12000],/Device,$
          yticks=1,yminor=1,ytickname=[' ',' '],$
          xticks=1,xminor=1,xtickname=[' ',' '],$
          levels = 10.^(findgen(100)/99.*6.+10.),/fill,/Noerase
        plot_oi,yscale,xscale,Position=[2600,6000,3100,12000],xticks=1,$
          xtickname=[' ',' '],yticks=3,yminor=1,$
          ytickname=['10!U10!N','10!U12!N','10!U14!N','10!U16!N'],$
          ytitle='!8L!3 (erg/s/cm!U3!N)',/Noerase,/Device
        
    endif
    if (ifig eq 1.3) then begin
        T_cor(where(T_cor lt 1d5)) = 1.1d5
        contour,T_cor(*,*,0),r#sin(t),r#cos(t),xrange=[0,40],yrange=[-20,20],$
          levels = 10.^(findgen(100)/99.*6.+5.),/isotropic,/fill,ystyle=1
        contour,tau_es(*,*,0),r#sin(t_),r#cos(t_),xrange=[0,40],yrange=[-20,20],$
          levels = [0.01,0.1,1],/isotropic,/noerase,thick=4,c_labels=[1,1,1],$
          xtitle = 'x/M', ytitle = 'z/M',ystyle=1
        phi = (findgen(101)/100.-0.5)*!PI
;        phi(100)=phi(0)
;        polyfill,2.*cos(phi),2.*sin(phi),color=0,/device
        oplot,2.*cos(phi),2.*sin(phi),thick=5
    endif
    if (ifig eq 1.31) then begin
        
        Npx = 12600
        Nx = 100
        Ny = 5
        xx = findgen(Nx)
        yy = findgen(Ny)
        xscale = findgen(Nx)/(Nx-1.)*6.+5.
        xscale = 10.^xscale/11000000. ;temp scale in keV
        yscale = fltarr(Nx)+1
        movscl = fltarr(Ny,Nx)
        for i=0,Nx-1 do movscl(*,i)=xscale(i)
        enlarge = rebin(movscl,50,1000)
        contour,enlarge,Position=[2600,6000,3100,12000],/Device,$
          yticks=1,yminor=1,ytickname=[' ',' '],$
          xticks=1,xminor=1,xtickname=[' ',' '],$
          levels = 10.^(findgen(100)/99.*6.+5.)/11000000.,/fill,/Noerase
        plot_oi,yscale,xscale,Position=[2600,6000,3100,12000],xticks=1,$
          xtickname=[' ',' '],yticks=3,yminor=1,$
          ytickname=['10!U-2!N','10!U0!N','10!U2!N','10!U4!N'],$
          ytitle='T (keV)',/Noerase,/Device
        
    endif
    if (ifig eq 1.4) then begin
        T_cor(where(T_cor lt 1d3)) = 1d3
        contour,bb_cgs(*,*,0),r#sin(t),r#cos(t),xrange=[0,60],yrange=[-30,30],$
          levels = 10.^(findgen(100)/99.*6.+8.),/isotropic,/fill,ystyle=1
        contour,tau_es(*,*,0),r#sin(t_),r#cos(t_),xrange=[0,60],yrange=[-30,30],$
          levels = [0.01,0.1,1],/isotropic,/noerase,thick=4,c_labels=[1,1,1],$
          xtitle = 'x/M', ytitle = 'z/M',ystyle=1
    endif
    if (ifig eq 1.41) then begin
        
        Npx = 12600
        Nx = 100
        Ny = 5
        xx = findgen(Nx)
        yy = findgen(Ny)
        xscale = findgen(Nx)/(Nx-1.)*6.+8.
        xscale = 10.^xscale/11000000. ;temp scale in keV
        yscale = fltarr(Nx)+1
        movscl = fltarr(Ny,Nx)
        for i=0,Nx-1 do movscl(*,i)=xscale(i)
        enlarge = rebin(movscl,50,1000)
        contour,enlarge,Position=[2600,6000,3100,12000],/Device,$
          yticks=1,yminor=1,ytickname=[' ',' '],$
          xticks=1,xminor=1,xtickname=[' ',' '],$
          levels = 10.^(findgen(100)/99.*6.+8.)/11000000.,/fill,/Noerase
        plot_oi,yscale,xscale,Position=[2600,6000,3100,12000],xticks=1,$
          xtickname=[' ',' '],yticks=3,yminor=1,$
          ytickname=['10!U8!N','10!U10!N','10!U12!N','10!U14!N'],$
          ytitle='U!LB!N (erg/cm!U3!N)',/Noerase,/Device
        
    endif
    stop
endif

if (fix(ifig) eq 2) then begin
    N1 = 0
    N2 = 0
    N3 = 0
    aa = 0.0
    t_frame = 0.
    R_hor = 2.0

    run_id = 1250
    rdatafile = 'data/gr_0000.dat'
    dumpstr = string(run_id,format='(I4.4)')
    strput,rdatafile,dumpstr,8
    openr,1,rdatafile

    readf,1,N1,N2,N3
    r = fltarr(N1)
    t = fltarr(N2)
    t_ = fltarr(N2+1)
    p = fltarr(N3)
    readf,1,r,t,p
    close,1
    
    dr = deriv(r)
    dt = deriv(t)
    dp = deriv(p)

    rho = fltarr(N1,N2,N3)
    LL = fltarr(N1,N2,N3)
    rr = fltarr(N1,N2,N3)
    tt = fltarr(N1,N2,N3)
    pp = fltarr(N1,N2,N3)
    dV = fltarr(N1,N2,N3)
    g_ = fltarr(N1,N2,4,4)

    for ir=0,N1-1 do rr(ir,*,*)=r(ir)
    for it=0,N2-1 do tt(*,it,*)=t(it)
    for ip=0,N3-1 do pp(*,*,ip)=p(ip)

    rho_cgs = fltarr(N1,N2,N3)
    ll_cgs = fltarr(N1,N2,N3)
    diskbody_ijk = intarr(N1,N2,N3)

    R_hor = 1.+sqrt(1.-aa*aa)
    Sigma = rr^2.+aa^2.*cos(tt)^2.
    Delta = rr^2.-2.*rr+aa^2.
    alpha = sqrt(Sigma*Delta/(Sigma*Delta+2.*rr*(aa^2.+rr^2.)))
    omega = 2.*rr*aa/(Sigma*Delta+2.*rr*(aa^2.+rr^2.))
    varph = sqrt((Sigma*Delta+2.*rr*(aa^2.+rr^2.))/Sigma*sin(tt)^2.)

    g_(*,*,0,0) = -alpha(*,*,0)^2.+omega(*,*,0)^2.*varph(*,*,0)^2.
    g_(*,*,0,3) = -omega(*,*,0)*varph(*,*,0)^2.
    g_(*,*,1,1) = Sigma(*,*,0)/Delta(*,*,0)
    g_(*,*,2,2) = Sigma(*,*,0)
    g_(*,*,3,3) = varph(*,*,0)^2.
    g_(*,*,3,0) = g_(*,*,0,3)

    detg = fltarr(N1,N2)
    detg(*,*) = g_(*,*,0,0)*g_(*,*,1,1)*g_(*,*,2,2)*g_(*,*,3,3) - $
      g_(*,*,0,3)*g_(*,*,1,1)*g_(*,*,2,2)*g_(*,*,0,3)
    for ip=0,N3-1 do dV(*,*,ip) = sqrt(-detg)*(dr#dt)*dp(0)

    G_N = 6.6726d-8
    M_BH = 10.*(2.d33)
    cc = 3d10
    t_cgs = G_N*M_BH/cc^3.
    r_cgs = G_N*M_BH/cc^2.
    dV_cgs = dV*r_cgs^3.

    sigtau_es = fltarr(N1,N3)
    netflux = fltarr(N1,N3)
    diskflux = fltarr(N1,N3)
    em_top = fltarr(N1,N3)
    em_bot = fltarr(N1,N3)
    ref_top = fltarr(N1,N3)
    ref_bot = fltarr(N1,N3)
    diskbody = fltarr(N1,N3)
    Tdisk = fltarr(N1,N3)

    rdatafile = 'data/ph_0000.dat'
    dumpstr = string(run_id,format='(I4.4)')
    strput,rdatafile,dumpstr,8
    openr,1,rdatafile
    readf,1,diskbody,sigtau_es,Tdisk,em_top,em_bot,ref_top,ref_bot
    close,1
    rdatafile = 'data/rh_0000.dat'
    dumpstr = string(run_id,format='(I4.4)')
    strput,rdatafile,dumpstr,8
    openr,1,rdatafile
    readf,1,rho_cgs
    close,1
    rdatafile = 'data/db_0000.dat'
    dumpstr = string(run_id,format='(I4.4)')
    strput,rdatafile,dumpstr,8
    openr,1,rdatafile
    readf,1,diskbody_ijk
    close,1
    rdatafile = 'data/ll_0000.dat'
    dumpstr = string(run_id,format='(I4.4)')
    strput,rdatafile,dumpstr,8
    openr,1,rdatafile
    readf,1,ll_cgs
    close,1

    corpow_ijk = fltarr(N1,N2,N3)
    readcp = fltarr(N3,N2,N1)
    rdatafile = 'data/scat_cpow.0000.dat'
    dumpstr = string(run_id,format='(I4.4)')
    strput,rdatafile,dumpstr,15
    openr,1,rdatafile
    readf,1,readcp
    close,1
    for ir=0,N1-1 do corpow_ijk(ir,*,*) = transpose(readcp(*,*,ir))

  incor = where((diskbody_ijk eq 0)and(rr gt R_hor))
  inatm = where((diskbody_ijk eq 1)and(rr gt R_hor))
  indisk = where((diskbody_ijk eq 2)and(rr gt R_hor))
  insim = where((diskbody_ijk le 2)and(rr gt R_hor))
  botcor = where((diskbody_ijk eq 0)and(rr gt R_hor)and(tt ge !PI/2))
  topcor = where((diskbody_ijk eq 0)and(rr gt R_hor)and(tt lt !PI/2))

  harm_cor = ll_cgs
  harm_cor(*,*,*)=0.
  harm_cor(incor)=ll_cgs(incor)*dv_cgs(incor)*4.
  rt_cor = ll_cgs
  rt_cor(*,*,*)=0.
  rt_cor(incor)=corpow_ijk(incor)*2.42d17
  nan_problem = where(finite(rt_cor,/NAN))
  
  if (nan_problem(0) ge 0) then rt_cor(nan_problem)=harm_cor(nan_problem)
  dldrharm = total(total(harm_cor,3),2)
  dldrrt = total(total(rt_cor,3),2)
  if (ifig eq 2.1) then begin
      plot_oo,r,dldrharm/dr*r,xrange=[2,60],yrange=[1d35,1d38],xstyle=1,$
        xtitle='r/M',ytitle='dL/d(log r)',thick=thk
      oplot,r,dldrrt/dr*r,color=120,thick=thk
      plots,[2.5,4],[1,1]*10^(37.7),thick=thk
      plots,[2.5,4],[1,1]*10^(37.4),color=120,thick=thk
      xyouts,4.1,10^(37.65),'Harm3d'
      xyouts,4.1,10^(37.35),'Pandurata'
  endif

  dldtharm = total(total(harm_cor,3),1)
  dldtrt = total(total(rt_cor,3),1)
  meantop = mean(cos(em_top))*0.5
  meanbot = mean(cos(em_bot))*0.5
  meandisk = where((cos(t) lt meantop)and(cos(t) gt meanbot))
  dldtharm(meandisk)=0
  dldtrt(meandisk)=0
  if (ifig eq 2.2) then begin
      plot_io,cos(t),dldtharm/abs(deriv(cos(t))),$
        xrange=[-1,1],yrange=[1d35,1d38],$
        xtitle='cos(!9q!3)',ytitle='dL/d(cos !9q!3)',thick=thk
      oplot,cos(t),dldtrt/abs(deriv(cos(t))),color=120,thick=thk
      plots,[-0.9,-0.7],[1,1]*10^(37.8),thick=thk
      plots,[-0.9,-0.7],[1,1]*10^(37.5),color=120,thick=thk
      xyouts,-0.65,10^(37.75),'Harm3d'
      xyouts,-0.65,10^(37.45),'Pandurata'
  endif
  stop
endif

if (fix(ifig) eq 3) then begin
    N1 = 0
    N2 = 0
    N3 = 0
    aa = 0.0
    t_frame = 0.
    R_hor = 2.0

    run_id = 1251
    rdatafile = 'data/gr_0000.dat'
    dumpstr = string(run_id,format='(I4.4)')
    strput,rdatafile,dumpstr,8
    openr,1,rdatafile

    readf,1,N1,N2,N3
    r = fltarr(N1)
    t = fltarr(N2)
    p = fltarr(N3)
    readf,1,r,t,p
    close,1
    
    dr = deriv(r)
    dt = deriv(t)
    dp = deriv(p)

    ll_cgs = fltarr(N1,N2,N3)
    wgt = fltarr(N1,N2,N3)
    rho_cgs = fltarr(N1,N2,N3)
    diskbody_ijk = intarr(N1,N2,N3)
    Te = fltarr(N1,N2,N3)
    tau_es = fltarr(N1,N2+1,N3)
    tauc = fltarr(N1,N2,N3)
    rr = fltarr(N1,N2,N3)
    tt = fltarr(N1,N2,N3)
    pp = fltarr(N1,N2,N3)
    dV = fltarr(N1,N2,N3)
    g_ = fltarr(N1,N2,4,4)

    for ir=0,N1-1 do rr(ir,*,*)=r(ir)
    for it=0,N2-1 do tt(*,it,*)=t(it)
    for ip=0,N3-1 do pp(*,*,ip)=p(ip)

    R_hor = 1.+sqrt(1.-aa*aa)
    Sigma = rr^2.+aa^2.*cos(tt)^2.
    Delta = rr^2.-2.*rr+aa^2.
    alpha = sqrt(Sigma*Delta/(Sigma*Delta+2.*rr*(aa^2.+rr^2.)))
    omega = 2.*rr*aa/(Sigma*Delta+2.*rr*(aa^2.+rr^2.))
    varph = sqrt((Sigma*Delta+2.*rr*(aa^2.+rr^2.))/Sigma*sin(tt)^2.)

    g_(*,*,0,0) = -alpha(*,*,0)^2.+omega(*,*,0)^2.*varph(*,*,0)^2.
    g_(*,*,0,3) = -omega(*,*,0)*varph(*,*,0)^2.
    g_(*,*,1,1) = Sigma(*,*,0)/Delta(*,*,0)
    g_(*,*,2,2) = Sigma(*,*,0)
    g_(*,*,3,3) = varph(*,*,0)^2.
    g_(*,*,3,0) = g_(*,*,0,3)

    detg = fltarr(N1,N2)
    detg(*,*) = g_(*,*,0,0)*g_(*,*,1,1)*g_(*,*,2,2)*g_(*,*,3,3) - $
      g_(*,*,0,3)*g_(*,*,1,1)*g_(*,*,2,2)*g_(*,*,0,3)
    for ip=0,N3-1 do dV(*,*,ip) = sqrt(-detg)*(dr#dt)*dp(0)

    G_N = 6.6726d-8
    M_BH = 10.*(2.d33)
    cc = 3d10
    t_cgs = G_N*M_BH/cc^3.
    r_cgs = G_N*M_BH/cc^2.
    dV_cgs = dV*r_cgs^3.
    
    sigtau_es = fltarr(N1,N3)
    netflux = fltarr(N1,N3)
    diskflux = fltarr(N1,N3)
    em_top = fltarr(N1,N3)
    em_bot = fltarr(N1,N3)
    ref_top = fltarr(N1,N3)
    ref_bot = fltarr(N1,N3)
    diskbody = fltarr(N1,N3)
    Tdisk = fltarr(N1,N3)
    ir_start = min(where(r gt R_hor))

    if ((ifig eq 3.1)or(ifig eq 3.2)or(ifig eq 3.21)) then begin
        for ilum = 0,4 do begin
            Tcor_r = fltarr(N1)
            Tcor_rtot = fltarr(N1)
            fcor_r = fltarr(N1)
            fcor_rt = fltarr(N1)
            Lcor = 0.
            Latm = 0.
            Ltot = 0.
            Nt = 0.
            Rin = 0.
            Rtrans = 0.
            HR = 0.
            for tslice=1001,1501,50 do begin
                Nt = Nt+1.
                run_id = tslice+ilum
                rdatafile = 'data/ph_0000.dat'
                dumpstr = string(run_id,format='(I4.4)')
                strput,rdatafile,dumpstr,8
                openr,1,rdatafile
                readf,1,diskbody,sigtau_es,Tdisk,em_top,em_bot,ref_top,ref_bot
                close,1
                dumpstr = string(run_id,format='(I4.4)')
                rdatafile = 'data/rh_0000.dat'
                strput,rdatafile,dumpstr,8
                openr,1,rdatafile
                readf,1,rho_cgs
                close,1
                rdatafile = 'data/ll_0000.dat'
                strput,rdatafile,dumpstr,8
                openr,1,rdatafile
                readf,1,ll_cgs
                close,1
                rdatafile = 'data/db_0000.dat'
                strput,rdatafile,dumpstr,8
                openr,1,rdatafile
                readf,1,diskbody_ijk
                close,1
                rdatafile = 'data/te_0000.dat'
                strput,rdatafile,dumpstr,8
                openr,1,rdatafile
                readf,1,Te
                close,1
                rdatafile = 'data/ta_0000.dat'
                strput,rdatafile,dumpstr,8
                openr,1,rdatafile
                readf,1,tau_es
                close,1
                tauc=0.5*(tau_es(*,0:N2-1,*)+tau_es(*,1:N2,*))
                incor = where((diskbody_ijk eq 0)and(rr gt R_hor))
                indisk = where((diskbody_ijk gt 0)and(rr gt R_hor))
                inatm = where((tauc gt 0.1)and(tauc lt 1.0)and(rr gt R_hor))
                insim = where((diskbody_ijk ge 0)and(rr gt R_hor))
                Latm = Latm+total(ll_cgs(inatm)*dV_cgs(inatm))
                Lcor = Lcor+total(ll_cgs(incor)*dV_cgs(incor))
                Ltot = Ltot+total(ll_cgs(insim)*dV_cgs(insim))
                for ir=0,N1-1 do begin
                    incor_r = where((diskbody_ijk eq 0)and(rr eq r(ir)))
                    insim_r = where(rr eq r(ir))
                    fcor_r(ir)  =total(ll_cgs(incor_r)*dV_cgs(incor_r))/ $
                      total(ll_cgs(insim_r)*dV_cgs(insim_r))
                endfor
                fcor_rt = fcor_rt+fcor_r
                harm_cor = ll_cgs*0
                harm_cor(incor) = ll_cgs(incor)*dV_cgs(incor)
                dldrcor = total(total(harm_cor,3),2)
                harm_disk = ll_cgs*0
                harm_disk(indisk) = ll_cgs(indisk)*dV_cgs(indisk)
                dldrdisk = total(total(harm_disk,3),2)
                rdex = min(where(dldrdisk gt dldrcor))
                Rtrans = Rtrans+r(rdex)
                middex = 100
                for iph=0,N3-1 do begin
                    rdex = min(where(diskbody(*,iph) gt 0))
                    Rin = Rin+r(rdex)
                    HR = HR+0.5*(abs(cos(em_top(middex,iph)))+ $
                                 abs(cos(em_bot(middex,iph))))
                endfor
                wgt(*,*,*)=0.
                if (ifig eq 3.1) then wgt(incor)=ll_cgs(incor)
                if (ifig eq 3.2) then wgt(incor)=rho_cgs(incor)
                for ir=ir_start,N1-1 do $
                  Tcor_r(ir) = total(Te(ir,*,*)*dV(ir,*,*)*wgt(ir,*,*))/$
                  total(dV(ir,*,*)*wgt(ir,*,*))
                Tcor_r = Tcor_r/(1.1d7) ;convert from K to keV
                Tcor_rtot = Tcor_rtot+Tcor_r
            endfor
            fcor_rt = fcor_rt/Nt
            Rtrans = Rtrans/Nt
            Rin = Rin/(Nt*N3)
            HR = HR/(Nt*N3)
            Latm = Latm/Nt
            Lcor = Lcor/Nt
            Ltot = Ltot/Nt
            print,Latm/Lcor,Lcor/Ltot,Rin,Rtrans,HR
            Tcor_r = Tcor_rtot/Nt
            if ((ifig eq 3.1)or(ifig eq 3.2)) then begin
                if (ilum eq 0) then begin
                    plot_oo,r,Tcor_r,xrange=[2,60],xstyle=1,yrange=[1,10000],ystyle=1,$
                      xtitle='r/M',ytitle = 'coronal T (keV)',thick=thk
                    plots,[13,20],[1,1]*10^3.6,color=60,thick=thk
                    plots,[13,20],[1,1]*10^3.4,color=90,thick=thk
                    plots,[13,20],[1,1]*10^3.2,color=120,thick=thk
                    plots,[13,20],[1,1]*10^3.0,color=150,thick=thk
                    plots,[13,20],[1,1]*10^2.8,color=180,thick=thk
                    !P.charsize=1.25
                    xyouts,21,10^3.55,'L/L!LEdd!N=0.01'
                    xyouts,32,10^3.35,'0.03'
                    xyouts,32,10^3.15,'0.1'
                    xyouts,32,10^2.95,'0.3'
                    xyouts,32,10^2.75,'1.0'
                    !P.charsize=1.5
                endif
                oplot,r,Tcor_r,color=60+ilum*30,thick=thk
            endif
            if (ifig eq 3.21) then begin
                if (ilum eq 0) then begin
                    plot_oi,r,fcor_rt,xrange=[2,60],xstyle=1,yrange=[0,1],ystyle=1,$
                      xtitle='r/M',ytitle = 'coronal fraction',thick=thk
                    plots,[13,20],[1,1]*0.9,color=60,thick=thk
                    plots,[13,20],[1,1]*0.85,color=90,thick=thk
                    plots,[13,20],[1,1]*0.8,color=120,thick=thk
                    plots,[13,20],[1,1]*0.75,color=150,thick=thk
                    plots,[13,20],[1,1]*0.7,color=180,thick=thk
                    xyouts,21,0.88,'L/L!LEdd!N=0.01'
                    xyouts,34,0.83,'0.03'
                    xyouts,34,0.78,'0.1'
                    xyouts,34,0.73,'0.3'
                    xyouts,34,0.68,'1.0'
                endif
                oplot,r,fcor_rt,color=60+ilum*30,thick=thk
            endif
        endfor
    endif

    if (ifig eq 3.3) then begin
        Ntau = 41.
        tau = 10.^(findgen(Ntau)/(Ntau-1.)*4.-4.)
        np_tau = fltarr(Ntau)
        Nt = 0.
        Te_r_tau = fltarr(N1,Ntau)
        np_r_tau = fltarr(N1,Ntau)
        ilum = 2
        for tslice=1001,1501,50 do begin
            Nt = Nt+1.
            run_id = tslice+ilum
            dumpstr = string(run_id,format='(I4.4)')
            rdatafile = 'data/te_0000.dat'
            strput,rdatafile,dumpstr,8
            openr,1,rdatafile
            readf,1,Te
            close,1
            rdatafile = 'data/ta_0000.dat'
            strput,rdatafile,dumpstr,8
            openr,1,rdatafile
            readf,1,tau_es
            close,1
            for it=0,N2-1 do tauc(*,it,*)=0.5*(tau_es(*,it,*)+tau_es(*,it+1,*))

            taudex = fix((alog10(tauc)+4.)/4.*Ntau)
            for ir=0,N1-1,20 do begin
                for ith=0,N2-1 do begin
                    for iph=0,N3-1 do begin
                        td = taudex(ir,ith,iph)
                        if ((td ge 0)and(td lt Ntau)) then begin
                            Te_r_tau(ir,td)=Te_r_tau(ir,td)+Te(ir,ith,iph)
                            np_r_tau(ir,td)=np_r_tau(ir,td)+1.
                        endif
                    endfor
                endfor
            endfor
        endfor
        for ir=0,N1-1,20 do Te_r_tau(ir,*)=Te_r_tau(ir,*)/np_r_tau(ir,*)
        Te_r_tau = Te_r_tau/1.1d7
        plot_oo,tau,Te_r_tau(40,*),yrange=[1,1000],xtitle='optical depth',$
          ytitle='T!Lcor!N (keV)',thick=thk
        for ir=40,140,20 do oplot,tau,Te_r_tau(ir,*),thick=thk,color=ir+30
        for ir=40.,140.,20 do $
          plots,[2d-4,5d-4],[1,1]*10^(ir/100.),thick=thk,color=(180-ir)+30
        xyouts,6d-4,10^(1.35),'R/M=4'
        xyouts,18d-4,10^(1.15),'6'
        xyouts,18d-4,10^(0.95),'9'
        xyouts,18d-4,10^(0.75),'13'
        xyouts,18d-4,10^(0.55),'20'
        xyouts,18d-4,10^(0.35),'30'
    endif
    
    if (ifig eq 3.4) then begin
        hrho = fltarr(N1,N3)
        hrad = fltarr(N1,N3)
        sig_rho = fltarr(N1,N3)
        sig_rad = fltarr(N1,N3)
        dthV = fltarr(N1,N2,N3)
        rho_i = fltarr(N1,N2,N3)
        ll_i = fltarr(N1,N2,N3)
        for ip=0,N3-1 do dthV(*,*,ip)=r#dt
        Nt = 0.
        for tslice=1001,1501,100 do begin
            hrho2 = fltarr(N1,N3)
            hrad2 = fltarr(N1,N3)
            Nt = Nt+1.
            run_id = tslice
            dumpstr = string(run_id,format='(I4.4)')
            rdatafile = 'data/rh_0000.dat'
            strput,rdatafile,dumpstr,8
            openr,1,rdatafile
            readf,1,rho_cgs
            close,1
            rdatafile = 'data/ll_0000.dat'
            strput,rdatafile,dumpstr,8
            openr,1,rdatafile
            readf,1,ll_cgs
            close,1
            for it=1,N2-1 do begin
                rho_i(*,it,*)=rho_i(*,it-1,*)+rho_cgs(*,it,*)*dthV(*,it,*)
                ll_i(*,it,*)=ll_i(*,it-1,*)+ll_cgs(*,it,*)*dthV(*,it,*)
            endfor
            sig_rho(*,*) = rho_i(*,N2-1,*)
            sig_rad(*,*) = ll_i(*,N2-1,*)
            for ir=0,N1-1 do begin
                for iph=0,N3-1 do begin
                    indisk = where(rho_i(ir,*,iph) gt 0.25*sig_rho(ir,iph))
                    hrho2(ir,iph)=t(indisk(0))
                    indisk = where(ll_i(ir,*,iph) gt 0.25*sig_rad(ir,iph))
                    hrad2(ir,iph)=t(indisk(0))
                endfor
            endfor
            hrho = hrho+hrho2
            hrad = hrad+hrad2
        endfor
        hrho = hrho/Nt
        hrad = hrad/Nt
        stop
    endif
endif    

if (ifig eq 4) then begin
    !P.charsize=1.25
    Nspec = 1001
    id0 = 2021
    Nframes = 0
    Nt = 41
    spec = dblarr(Nspec,Nt)
    spec2 = dblarr(Nspec,Nt)

    emin = 0.001              
    emax = 10000.
    nu = emin*10.^(findgen(Nspec)/(Nspec-1.)*alog10(emax/emin))
    dnu = deriv(nu)

    rdatafile = 'data/scat_spec.0000.dat'
    for run_id = id0,id0+500,50 do begin
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,15
        openr,1,rdatafile
        readf,1,spec2
        close,1
        spec=spec+spec2
        Nframes=Nframes+1
    endfor
    spec = spec/Nframes
    plot_oo,nu,nu*total(spec,2)*2.4d17,xrange=[0.1,1000],yrange=[1d35,1d39], $
      xstyle=1,xtitle='E!Lobs!N (keV)',ytitle='Luminosity !9n!3L!L!9n!N!3', $
      thick=thk,/isotropic
    oplot,nu,nu*total(spec,2)*2.4d17,thick=thk,color=60 
    print,total(dnu*total(spec,2))*2.4d17
   ;inclination dependence
;    Ispec = 0.5*(spec(*,0)+spec(*,Nt-1))
;    Ispec = Ispec/total(Ispec*dnu)*1.3d37
;    plot_oo,nu,nu*Ispec,xrange=[0.1,1000],yrange=[1d35,1d37], $
;      xstyle=1,xtitle='E!Lobs!N (keV)',ytitle='Luminosity !9n!3L!L!9n!N!3', $
;      thick=thk
;    for inc=0,Nt/2,5 do begin
;        Ispec = 0.5*(spec(*,inc)+spec(*,Nt-1-inc))
;        Ispec = Ispec/total(Ispec*dnu)*1.3d37
;        oplot,nu,nu*Ispec,thick=thk,color=60+inc*4
;    endfor
;stop

    for id0=2022,2025 do begin
        spec(*,*)=0.
        Nframes = 0
        for run_id = id0,id0+500,50 do begin
            dumpstr = string(run_id,format='(I4.4)')
            strput,rdatafile,dumpstr,15
            openr,1,rdatafile
            readf,1,spec2
            close,1
            spec=spec+spec2
            Nframes=Nframes+1
        endfor
        spec = spec/Nframes
        oplot,nu,nu*total(spec,2)*2.4d17,thick=thk,color=60+30*(id0-2021)
        print,total(dnu*total(spec,2))*2.4d17
;stop        
    endfor
    plots,[20,60],[1.,1.]*10.^38.5,thick=thk,color=60
    plots,[20,60],[1.,1.]*10.^38.25,thick=thk,color=90
    plots,[20,60],[1.,1.]*10.^38.0,thick=thk,color=120
    plots,[20,60],[1.,1.]*10.^37.75,thick=thk,color=150
    plots,[20,60],[1.,1.]*10.^37.5,thick=thk,color=180
    xyouts,70,10.^38.4,'L/L!LEdd!N=0.01'
    xyouts,340,10.^38.15,'0.03'
    xyouts,340,10.^37.9,'0.1'
    xyouts,340,10.^37.65,'0.3'
    xyouts,340,10.^37.4,'1.0'
    stop
endif

if (ifig eq 4.01) then begin
    Nspec = 101
    id0 = 2003
    Nframes = 0
    Nt = 41
    spec = dblarr(Nspec,Nt)
    spec2 = dblarr(Nspec,Nt)

    emin = 0.001              
    emax = 10000.
    nu = emin*10.^(findgen(Nspec)/(Nspec-1.)*alog10(emax/emin))
    dnu = deriv(nu)

    rdatafile = 'data/scat_spec.0000.dat'
    for run_id = id0,id0+500,100 do begin
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,15
        openr,1,rdatafile
        readf,1,spec2
        close,1
        spec=spec+spec2
        Nframes=Nframes+1
    endfor
    spec = spec/Nframes
    plot_oo,nu,nu*total(spec,2)*2.4d17,xrange=[0.1,1000],yrange=[1d34,1d39], $
      xstyle=1,xtitle='E!Lobs!N (keV)',ytitle='Luminosity !9n!3L!L!9n!N!3', $
      thick=thk
    oplot,nu,nu*total(spec,2)*2.4d17,thick=thk,color=60

    for id0=2003,2005,2 do begin
        spec(*,*)=0.
        Nframes = 0
        for run_id = id0,id0+500,100 do begin
            dumpstr = string(run_id,format='(I4.4)')
            strput,rdatafile,dumpstr,15
            openr,1,rdatafile
            readf,1,spec2
            close,1
            spec=spec+spec2
            Nframes=Nframes+1
        endfor
        spec = spec/Nframes
        oplot,nu,nu*total(spec,2)*2.4d17,thick=thk,color=60+30*(id0-2001)
    endfor
    plots,[40,100],[1.,1.]*10.^38.5,thick=thk,color=60
    plots,[40,100],[1.,1.]*10.^38.2,thick=thk,color=120
    plots,[40,100],[1.,1.]*10.^37.9,thick=thk,color=180
    xyouts,110,10.^38.4,'L/L!LEdd!N=0.01'
    xyouts,350,10.^38.1,'0.1'
    xyouts,350,10.^37.8,'1.0'
    stop
endif

if (ifig eq 4.1) then begin
    Nspec = 1001
    id0 = 1023
    Nframes = 0
    Nt = 41
    idex15 = 0
    idex30 = 3
    idex45 = 6
    idex60 = 10
    idex75 = 15
    spec = dblarr(Nspec,Nt)
    spec2 = dblarr(Nspec,Nt)

    emin = 0.001              
    emax = 10000.
    nu = emin*10.^(findgen(Nspec)/(Nspec-1.)*alog10(emax/emin))
    dnu = deriv(nu)

    rdatafile = 'data/scat_spec.0000.dat'
    for run_id = id0,id0+500,50 do begin
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,15
        openr,1,rdatafile
        readf,1,spec2
        close,1
        spec=spec+spec2
        Nframes=Nframes+1
    endfor
    for it = 0,(Nt-1)/2. do $
      spec(*,it)=spec(*,it)+spec(*,Nt-it-1)
    spec = spec/Nframes*Nt/2.
    plot_oo,nu,nu*spec(*,idex15)*2.4d17,xrange=[0.1,1000],yrange=[1d34,1d39], $
      xstyle=1,xtitle='E!Lobs!N (keV)',ytitle='Luminosity !9n!3L!L!9n!N!3', $
      thick=thk
    oplot,nu,nu*spec(*,idex15)*2.4d17,thick=thk,color=60
    oplot,nu,nu*spec(*,idex30)*2.4d17,thick=thk,color=90
    oplot,nu,nu*spec(*,idex45)*2.4d17,thick=thk,color=120
    oplot,nu,nu*spec(*,idex60)*2.4d17,thick=thk,color=150
    oplot,nu,nu*spec(*,idex75)*2.4d17,thick=thk,color=180
    stop
endif
    
if ((ifig eq 4.21)or(ifig eq 4.22)or(ifig eq 4.23)) then begin
    Nspec = 1001
    if (ifig eq 4.21) then id0 = 1021
    if (ifig eq 4.22) then id0 = 1023
    if (ifig eq 4.23) then id0 = 1025
    Nframes = 0
    Nt = 41
    Nr = 180
    spec = dblarr(Nspec,Nt)
    spec2 = dblarr(Nspec,Nt)
    specr = dblarr(Nspec,Nt,Nr)
    specr2 = dblarr(Nspec,Nt,Nr)
    cspecr = dblarr(Nspec,Nr)
    cspecr2 = dblarr(Nspec,Nr)
    N1 = 0
    N2 = 0
    N3 = 0

    r0 = 2.
    r1 = 6.
    r2 = 15.
    r3 = 30.
    r4 = 60.

    emin = 0.001              
    emax = 10000.
    nu = emin*10.^(findgen(Nspec)/(Nspec-1.)*alog10(emax/emin))
    dnu = deriv(nu)

    rdatafile = 'data/gr_0000.dat'
    dumpstr = string(id0,format='(I4.4)')
    strput,rdatafile,dumpstr,8
    openr,1,rdatafile
    readf,1,N1,N2,N3
    r = fltarr(N1)
    t = fltarr(N2)
    p = fltarr(N3)
    readf,1,r,t,p
    close,1
    rdex1 = where((r ge r0)and(r lt r1))
    rdex2 = where((r ge r1)and(r lt r2))
    rdex3 = where((r ge r2)and(r lt r3))
    rdex4 = where((r ge r3)and(r lt r4))

    if (ifig eq 4.21) then ymax = 1d38
    if (ifig eq 4.22) then ymax = 1d39
    if (ifig eq 4.23) then ymax = 1d40
    ymin = ymax/1d6
    tmp = fltarr(N1)
    qspecr = specr
    uspecr = specr
    rspecr = specr
    rdatafile = 'data/scat_spcr.0000.dat'
    for run_id = id0,id0+500,100 do begin
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,15
        openr,1,rdatafile
        readf,1,specr2,tmp,tmp,tmp,qspecr,uspecr,rspecr,cspecr2
        close,1
        specr=specr+specr2
        cspecr=cspecr+cspecr2
        Nframes=Nframes+1
    endfor
    specr = specr/Nframes
    cspecr = cspecr/Nframes
    dnur = dnu#(fltarr(N1)+1)
    C_r = total(cspecr*dnur,1)
    E_r = total(total(specr,2)*dnur,1)

    spec = total(specr,3)
    plot_oo,nu,nu*total(spec,2)*2.4d17,xrange=[0.1,1000],yrange=[ymin,ymax], $
      xstyle=1,xtitle='E!Lobs!N (keV)',ytitle='Luminosity !9n!3L!L!9n!N!3', $
      thick=thk
    oplot,nu,nu*total(total(specr(*,*,rdex1),3),2)*2.4d17,color=150,thick=thk
    oplot,nu,nu*total(total(specr(*,*,rdex2),3),2)*2.4d17,color=120,thick=thk
    oplot,nu,nu*total(total(specr(*,*,rdex3),3),2)*2.4d17,color=90,thick=thk
    oplot,nu,nu*total(total(specr(*,*,rdex4),3),2)*2.4d17,color=60,thick=thk
    plots,[30,100],[1,1]*ymax*10.^(-0.4),thick=thk
    plots,[30,100],[1,1]*ymax*10.^(-0.7),thick=thk,color=150
    plots,[30,100],[1,1]*ymax*10.^(-1.0),thick=thk,color=120
    plots,[30,100],[1,1]*ymax*10.^(-1.3),thick=thk,color=90
    plots,[30,100],[1,1]*ymax*10.^(-1.6),thick=thk,color=60
    !P.charsize=1.25
    if (ifig eq 4.21) then xyouts, 0.2,ymax*10.^(-0.5),'L/L!LEdd!N=0.01'
    if (ifig eq 4.22) then xyouts, 0.2,ymax*10.^(-0.5),'L/L!LEdd!N=0.1'
    if (ifig eq 4.23) then xyouts, 0.2,ymax*10.^(-0.5),'L/L!LEdd!N=1.0'
    xyouts,110,ymax*10.^(-0.44),'total'
    xyouts,110,ymax*10.^(-0.74),'r/M=2-6'
    xyouts,110,ymax*10.^(-1.04),'     6-15'
    xyouts,110,ymax*10.^(-1.34),'     15-30'
    xyouts,110,ymax*10.^(-1.64),'     30-60'
    !P.charsize=1.5
    stop
endif
    
if (ifig eq 4.3) then begin
    Nspec = 1001
    id0 = 1023
    Nframes = 0
    Nt = 41
    spec = dblarr(Nspec,Nt)
    spec2 = dblarr(Nspec,Nt)
    N1 = 0
    N2 = 0
    N3 = 0

    emin = 0.001              
    emax = 10000.
    nu = emin*10.^(findgen(Nspec)/(Nspec-1.)*alog10(emax/emin))
    dnu = deriv(nu)

    ymax = 1d39
    ymin = ymax/1d5
    rdatafile = 'data/scat_spec.0000.dat'
    for run_id = id0,id0+500,50 do begin
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,15
        openr,1,rdatafile
        readf,1,spec2
        close,1
        spec=spec+spec2
        Nframes=Nframes+1
    endfor
    spec = spec/Nframes
    plot_oo,nu,nu*total(spec,2)*2.4d17,xrange=[0.1,1000],yrange=[ymin,ymax], $
      xstyle=1,xtitle='E!Lobs!N (keV)',ytitle='Luminosity !9n!3L!L!9n!N!3', $
      thick=thk
    print,total(dnu*total(spec,2)*2.4d17)

    Nframes = 0
    Nspec = 1001.
    nu = emin*10.^(findgen(Nspec)/(Nspec-1.)*alog10(emax/emin))
    dnu = deriv(nu)
    spec = dblarr(Nspec,Nt)
    spec2 = dblarr(Nspec,Nt)
    id0 = 1007
    for run_id = id0,id0+500,100 do begin
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,15
        openr,1,rdatafile
        readf,1,spec2
        close,1
        spec=spec+spec2
        Nframes=Nframes+1
    endfor
    spec = spec/Nframes
    oplot,nu,nu*total(spec,2)*2.4d17,thick=thk,color=120
    print,total(dnu*total(spec,2)*2.4d17)

    Nframes = 0
    Nspec = 1001.
    nu = emin*10.^(findgen(Nspec)/(Nspec-1.)*alog10(emax/emin))
    dnu = deriv(nu)
    spec = dblarr(Nspec,Nt)
    spec2 = dblarr(Nspec,Nt)
    id0 = 1017
    for run_id = id0,id0+500,100 do begin
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,15
        openr,1,rdatafile
        readf,1,spec2
        close,1
        spec=spec+spec2
        Nframes=Nframes+1
    endfor
    spec = spec/Nframes
    oplot,nu,nu*total(spec,2)*2.4d17,thick=thk,color=60
    print,total(dnu*total(spec,2)*2.4d17)

    plots,[30,100],[1,1]*ymax*10.^(-0.4),thick=thk,color=60
    plots,[30,100],[1,1]*ymax*10.^(-0.8),thick=thk
    plots,[30,100],[1,1]*ymax*10.^(-1.2),thick=thk,color=120
    xyouts,110,ymax*10.^(-0.44),'!9t!3!Lphot!N=2.0'
    xyouts,110,ymax*10.^(-0.84),'!9t!3!Lphot!N=1.0'
    xyouts,110,ymax*10.^(-1.24),'!9t!3!Lphot!N=0.5'


    stop
endif
    

if (fix(ifig) eq 5) then begin
    Nspec = 1001
    id0 = 2012
    Nt = 41
    idex15 = 0
    idex30 = 2
    idex45 = 5
    idex60 = 10
    idex75 = 13
    spec = dblarr(Nspec,Nt)
    specl = dblarr(Nspec,Nt)
    spect = dblarr(Nspec,Nt)
    spec2 = dblarr(Nspec,Nt)
    spec3 = dblarr(Nspec,Nt)

    emin = 0.001              
    emax = 10000.
    nu = emin*10.^(findgen(Nspec)/(Nspec-1.)*alog10(emax/emin))
    dnu = deriv(nu)

    idfin = 500
    rdatafile = 'data/scat_spec.0000.dat'
    for run_id = id0,id0+idfin,50 do begin
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,15
        openr,1,rdatafile
        readf,1,spec2
        close,1
        spec=spec+spec2
;        dumpstr = string(run_id+30,format='(I4.4)')
;        strput,rdatafile,dumpstr,15
;        openr,1,rdatafile
;        readf,1,spec3
;        close,1
;        plot_oi,nu,spec3(*,idex45)/spec2(*,idex45),xrange=[1,100],yrange=[0.5,1.5],$
;          xstyle=1,xtitle='E!Lobs!N (keV)',ytitle='line ratio'
;        stop
    endfor
    for run_id = id0+30,id0+idfin+30,50 do begin
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,15
        openr,1,rdatafile
        readf,1,spec2
        close,1
        specl=specl+spec2
    endfor
    for run_id = id0+10,id0+idfin+10,50 do begin
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,15
        openr,1,rdatafile
        readf,1,spec2
        close,1
        spect=spect+spec2
    endfor
    for it = 0,(Nt-1)/2. do begin
        spec(*,it)=spec(*,it)+spec(*,Nt-it-1)
        specl(*,it)=specl(*,it)+specl(*,Nt-it-1)
        spect(*,it)=spect(*,it)+spect(*,Nt-it-1)
    endfor
    speca = spect-specl+spec
    if (ifig eq 5.1) then begin
        plot,nu,specl(*,idex15)/spec(*,idex15),xrange=[2,10],yrange=[1.0,1.4],ystyle=1,$
          xstyle=1,xtitle='E!Lobs!N (keV)',ytitle='line ratio',thick=thk
        oplot,nu,specl(*,idex15)/spec(*,idex15),thick=thk,color=60
        oplot,nu,specl(*,idex30)/spec(*,idex30),thick=thk,color=90
        oplot,nu,specl(*,idex45)/spec(*,idex45),thick=thk,color=120
        oplot,nu,specl(*,idex60)/spec(*,idex60),thick=thk,color=150
        oplot,nu,specl(*,idex75)/spec(*,idex75),thick=thk,color=180
        plots,[2.5,3.5],[1,1]*1.35,thick=thk,color=60
        plots,[2.5,3.5],[1,1]*1.33,thick=thk,color=90
        plots,[2.5,3.5],[1,1]*1.31,thick=thk,color=120
        plots,[2.5,3.5],[1,1]*1.29,thick=thk,color=150
        plots,[2.5,3.5],[1,1]*1.27,thick=thk,color=180
        xyouts,3.6,1.34,'i=15!Uo!N'
        xyouts,3.9,1.32,'30!Uo!N'
        xyouts,3.9,1.30,'45!Uo!N'
        xyouts,3.9,1.28,'60!Uo!N'
        xyouts,3.9,1.26,'75!Uo!N'
stop
        Elow = 4.0
        Ecut = 8.0
        lineratio = specl(*,idex60)/spec(*,idex60)
        inline = where((nu gt Elow)and(nu le Ecut))
        aslope = (lineratio(max(inline))-1.)/(Ecut-Elow)
        bint = -Elow*aslope
        lineratio(inline)=lineratio(inline)-(aslope*nu(inline)+bint)
        lineratio(where(nu gt Ecut)) = 1.
        print,total((lineratio-1.)*dnu)
        oplot,nu,lineratio
        lineratio = specl(*,idex30)/spec(*,idex30)
        lineratio(where(nu gt 7)) = 1.
        print,total((lineratio-1.)*dnu)
        lineratio = specl(*,idex45)/spec(*,idex45)
        lineratio(where(nu gt 7.5)) = 1.
        print,total((lineratio-1.)*dnu)
        lineratio = specl(*,idex60)/spec(*,idex60)
        lineratio(where(nu gt 8)) = 1.
        print,total((lineratio-1.)*dnu)
        lineratio = specl(*,idex75)/spec(*,idex75)
        lineratio(where(nu gt 8)) = 1.
        print,total((lineratio-1.)*dnu)
    endif
    if (ifig eq 5.2) then begin
        plot_oi,nu,specl(*,idex30)/spec(*,idex30),xrange=[1,100],yrange=[0.5,1.4],$
          xstyle=1,ystyle=1,xtitle='E!Lobs!N (keV)',ytitle='line ratio'
        oplot,nu,specl(*,idex30)/spec(*,idex30),color=60,thick=thk
        oplot,nu,spect(*,idex30)/spec(*,idex30),thick=thk
        oplot,nu,speca(*,idex30)/spec(*,idex30),color=120,thick=thk
    endif
    if (ifig eq 5.3) then begin
        P_ = fltarr(5)
        spect = spect/22.
        Elow = 4.0
        Ecut = 8.0
        inline = where((nu gt Elow)and(nu le Ecut))
        openr,1,'bfp_001_i15.dat'
        readf,1,P_
        close,1
        fitspec=calc_diskpl_spec(nu,P_)
        plot,nu,spect(*,idex15)/fitspec,xrange=[2,10],yrange=[0.9,1.5],ystyle=1,$
          xstyle=1,xtitle='E!Lobs!N (keV)',ytitle='line ratio',thick=thk
        oplot,nu,spect(*,idex15)/fitspec,color=60,thick=thk
        lineratio = spect(*,idex15)/fitspec-1.
        print,total(lineratio(inline)*dnu(inline))
        openr,1,'bfp_001_i30.dat'
        readf,1,P_
        close,1
        fitspec=calc_diskpl_spec(nu,P_)
        oplot,nu,spect(*,idex30)/fitspec,color=90,thick=thk
        lineratio = spect(*,idex30)/fitspec-1.
        print,total(lineratio(inline)*dnu(inline))
        openr,1,'bfp_001_i45.dat'
        readf,1,P_
        close,1
        fitspec=calc_diskpl_spec(nu,P_)
        oplot,nu,spect(*,idex45)/fitspec,color=120,thick=thk
        lineratio = spect(*,idex45)/fitspec-1.
        print,total(lineratio(inline)*dnu(inline))
        openr,1,'bfp_001_i60.dat'
        readf,1,P_
        close,1
        fitspec=calc_diskpl_spec(nu,P_)
        oplot,nu,spect(*,idex60)/fitspec,color=150,thick=thk
        lineratio = spect(*,idex60)/fitspec-1.
        print,total(lineratio(inline)*dnu(inline))
        openr,1,'bfp_001_i75.dat'
        readf,1,P_
        close,1
        fitspec=calc_diskpl_spec(nu,P_)
        oplot,nu,spect(*,idex75)/fitspec,color=180,thick=thk
        lineratio = spect(*,idex75)/fitspec-1.
        print,total(lineratio(inline)*dnu(inline))
        plots,[2.5,3.5],[1,1]*1.45,thick=thk,color=60
        plots,[2.5,3.5],[1,1]*1.42,thick=thk,color=90
        plots,[2.5,3.5],[1,1]*1.39,thick=thk,color=120
        plots,[2.5,3.5],[1,1]*1.36,thick=thk,color=150
        plots,[2.5,3.5],[1,1]*1.33,thick=thk,color=180
        xyouts,3.6,1.44,'i=15!Uo!N'
        xyouts,3.9,1.41,'30!Uo!N'
        xyouts,3.9,1.38,'45!Uo!N'
        xyouts,3.9,1.35,'60!Uo!N'
        xyouts,3.9,1.32,'75!Uo!N'
    endif
    stop
endif

if (ifig eq 6) then begin
    Nspec = 1001
    Nt = 41
    id0 = 2011
    idfin = 500
    idex15 = 0
    idex30 = 1
    idex45 = 5
    idex60 = 10
    idex75 = 15
    spec = dblarr(Nspec,Nt)
    specl = dblarr(Nspec,Nt)
    spect = dblarr(Nspec,Nt)
    spec2 = dblarr(Nspec,Nt)
    
    spec_lum = dblarr(5,Nspec,Nt)
    specl_lum = dblarr(5,Nspec,Nt)
    spect_lum = dblarr(5,Nspec,Nt)

    emin = 0.001              
    emax = 10000.
    nu = emin*10.^(findgen(Nspec)/(Nspec-1.)*alog10(emax/emin))
    dnu = deriv(nu)

    rdatafile = 'data/scat_spec.0000.dat'
    for ilum = 0,2 do begin
        for run_id = id0+ilum,id0+ilum+idfin,50 do begin
            dumpstr = string(run_id,format='(I4.4)')
            strput,rdatafile,dumpstr,15
            openr,1,rdatafile
            readf,1,spec2
            close,1
            spec=spec+spec2
        endfor
        for run_id = id0+30+ilum,id0+idfin+30+ilum,50 do begin
            dumpstr = string(run_id,format='(I4.4)')
            strput,rdatafile,dumpstr,15
            openr,1,rdatafile
            readf,1,spec2
            close,1
            specl=specl+spec2
        endfor
;        for run_id = 1021+ilum,1521+ilum,50 do begin
;            dumpstr = string(run_id,format='(I4.4)')
;            strput,rdatafile,dumpstr,15
;            openr,1,rdatafile
;            readf,1,spec2
;            close,1
;            spect=spect+spec2
;        endfor
        spec_lum(ilum,*,*)=spec
        specl_lum(ilum,*,*)=specl
;        spect_lum(ilum,*,*)=spect
    endfor
    for it = 0,(Nt-1)/2. do begin
        spec_lum(*,*,it)=spec_lum(*,*,it)+spec_lum(*,*,Nt-it-1)
        specl_lum(*,*,it)=specl_lum(*,*,it)+specl_lum(*,*,Nt-it-1)
;        spect_lum(*,*,it)=spect_lum(*,*,it)+spect_lum(*,*,Nt-it-1)
    endfor
    plot,nu,specl_lum(0,*,idex30)/spec_lum(0,*,idex30),xrange=[2,10],yrange=[1,1.4],$
      xstyle=1,xtitle='E!Lobs!N (keV)',ytitle='line ratio'
    oplot,nu,specl_lum(0,*,idex30)/spec_lum(0,*,idex30),thick=thk,color=60
    oplot,nu,specl_lum(1,*,idex30)/spec_lum(1,*,idex30),thick=thk,color=90
    oplot,nu,specl_lum(2,*,idex30)/spec_lum(2,*,idex30),thick=thk,color=120
;    oplot,nu,specl_lum(3,*,idex30)/spec_lum(3,*,idex30),thick=thk,color=150
;    oplot,nu,specl_lum(4,*,idex30)/spec_lum(4,*,idex30),thick=thk,color=180
    plots,[2.5,3.5],[1,1]*1.35,thick=thk,color=60
    plots,[2.5,3.5],[1,1]*1.33,thick=thk,color=90
    plots,[2.5,3.5],[1,1]*1.31,thick=thk,color=120
;    plots,[2.5,3.5],[1,1]*1.36,thick=thk,color=150
;    plots,[2.5,3.5],[1,1]*1.33,thick=thk,color=180
    xyouts,3.6,1.34,'L/L!LEdd!N=0.01'
    xyouts,4.7,1.32,'0.03'
    xyouts,4.7,1.30,'0.1'
;    xyouts,4.7,1.35,'0.3'
;    xyouts,4.7,1.32,'1.0'
endif

if (ifig eq 7) then begin
    N1 = 0
    N2 = 0
    N3 = 0

    run_id = 1225
    rdatafile = 'data/gr_0000.dat'
    dumpstr = string(run_id,format='(I4.4)')
    strput,rdatafile,dumpstr,8
    openr,1,rdatafile

    readf,1,N1,N2,N3
    r = fltarr(N1)
    t = fltarr(N2)
    p = fltarr(N3)
    readf,1,r,t,p
    close,1
    dr = deriv(r)
    dt = deriv(t)
    Nspec = 1001
    emin = 0.001              
    emax = 10000.
    nu = emin*10.^(findgen(Nspec)/(Nspec-1.)*alog10(emax/emin))
    dnu = deriv(nu)
    I_r = fltarr(N1)
    nu_lo = 6.4
    nu_hi = 30.0
    erange = where((nu gt nu_lo)and(nu lt nu_hi))
    erange0 = where((nu gt 2)and(nu lt nu_hi))
    ll = fltarr(N1,N2,N3)
    ll2 = ll
    dthV = fltarr(N1,N2,N3)
    diskbody_ijk = intarr(N1,N2,N3)
    for ip=0,N3-1 do dthV(*,*,ip)=r#dt
    Nframes = 0.
    id0=1024
    for run_id = id0+20,id0+520,50 do begin
        rdatafile = 'data/ll_0000.dat'
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,8
        openr,1,rdatafile
        readf,1,ll2
        close,1
        rdatafile = 'data/db_0000.dat'
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,8
        openr,1,rdatafile
        readf,1,diskbody_ijk
        close,1
        ll2(where(diskbody_ijk eq 0))=0.
        ll = ll+ll2
        Nframes=Nframes+1.
    endfor
    ll = ll/Nframes
    Ir_em = total(total(ll*dthV,2),2)/N3*1.45d6
    T_r = 1.8*(Ir_em/5.67d-5)^0.25
    Nfe_r = fltarr(N1)
    for ir=0,N1-1 do begin
        x = 1.d3*nu/(T_r(ir)*8.6173d-5)
        Inu = nu^3./(exp(x)-1.)
        nupeak = nu(where(Inu eq max(Inu)))
;        erange0 = where((nu gt 2*nupeak(0))and(nu lt nu_hi))
        Nfe_r(ir) = total(Inu(erange0)*dnu(erange0)*2.4d17/(nu(erange0)*1.6d-12))
;        Nfe_r(ir) = total(Inu*dnu*2.4d17/(nu*1.6d-12))
    endfor
;stop
    Nfe_r = Nfe_r/Nfe_r(100)*1d36

    Rspec = fltarr(Nspec,N1)
    Rspec2 = fltarr(Nspec,N1)
    rdatafile = 'data/scat_disk.0000.dat'
    for id0=1021,1025,1 do begin
        Nframes = 0.
        for run_id = id0+20,id0+520,50 do begin
            dumpstr = string(run_id,format='(I4.4)')
            strput,rdatafile,dumpstr,15
            openr,1,rdatafile
            readf,1,Rspec2
            close,1
            Rspec=Rspec+Rspec2
            Nframes=Nframes+1.
        endfor
        Rspec = Rspec/Nframes
        for ir = 0,N1-1 do $
          I_r(ir)=total(Rspec(erange,ir)*dnu(erange)*2.4d17/(nu(erange)*1.6d-12))
        if (id0 eq 1021) then begin
            plot_oo,r,I_r,xrange=[2,60],xstyle=1,thick=thk, $
              xtitle='r/M', ytitle='Fe K absorption [cm!U-2!N s!U-1!N]',yrange=[1d31,1d39],$
              ystyle=1
            plots,[12,20],[1.,1.]*10.^38.5,thick=thk,color=60
            plots,[12,20],[1.,1.]*10.^38.1,thick=thk,color=90
            plots,[12,20],[1.,1.]*10.^37.7,thick=thk,color=120
            plots,[12,20],[1.,1.]*10.^37.3,thick=thk,color=150
            plots,[12,20],[1.,1.]*10.^36.9,thick=thk,color=180
            xyouts,22,10.^38.4,'L/L!LEdd!N=0.01'
            xyouts,37,10.^38.0,'0.03'
            xyouts,37,10.^37.6,'0.1'
            xyouts,37,10.^37.2,'0.3'
            xyouts,37,10.^36.8,'1.0'
            plots,[2.5,4],[1,1]*10.^38.5,thick=thk,linestyle=2
            xyouts,4.1,10.^38.4,'F!Lem!N'

        endif
        oplot,r,I_r,thick=thk,color=60+30*(id0-1021)
    endfor
    oplot,r(where(r gt 2)),Nfe_r(where(r gt 2)),thick=thk,linestyle=2
stop
endif

if (fix(ifig) eq 8) then begin
    N1 = 0
    N2 = 0
    N3 = 0
    aa = 0.0
    t_frame = 0.
    R_hor = 2.0

    id0 = 2003
    run_id = id0
    rdatafile = 'data/gr_0000.dat'
    dumpstr = string(run_id,format='(I4.4)')
    strput,rdatafile,dumpstr,8
    openr,1,rdatafile

    readf,1,N1,N2,N3
    r = fltarr(N1)
    t = fltarr(N2)
    t_ = fltarr(N2+1)
    p = fltarr(N3)
    readf,1,r,t,p
    close,1
    
    dr = deriv(r)
    dt = deriv(t)
    dp = deriv(p)

    rho = fltarr(N1,N2,N3)
    LL = fltarr(N1,N2,N3)
    rr = fltarr(N1,N2,N3)
    tt = fltarr(N1,N2,N3)
    pp = fltarr(N1,N2,N3)
    dV = fltarr(N1,N2,N3)
    g_ = fltarr(N1,N2,4,4)

    for ir=0,N1-1 do rr(ir,*,*)=r(ir)
    for it=0,N2-1 do tt(*,it,*)=t(it)
    for ip=0,N3-1 do pp(*,*,ip)=p(ip)

    ll_cgs = fltarr(N1,N2,N3)
    diskbody_ijk = intarr(N1,N2,N3)

    
    G_N = 6.6726d-8
    M_BH = 10.*(2.d33)
    cc = 3d10
    t_cgs = G_N*M_BH/cc^3.
    r_cgs = G_N*M_BH/cc^2.
    dV_cgs = dV*r_cgs^3.

    sigtau_es = fltarr(N1,N3)
    netflux = fltarr(N1,N3)
    diskflux = fltarr(N1,N3)
    em_top = fltarr(N1,N3)
    em_bot = fltarr(N1,N3)
    ref_top = fltarr(N1,N3)
    ref_bot = fltarr(N1,N3)
    diskbody = fltarr(N1,N3)
    Tdisk = fltarr(N1,N3)
    dldr = fltarr(N1)
    dldr0 = fltarr(N1)
    dldr001 = fltarr(N1)
    dldr01 = fltarr(N1)
    dldr1 = fltarr(N1)
    L001 = 7.9d36               ;seed photon luminosity
    L01 = 1.07d38
    L1 = 1.2d39

    Nframes = 0
    id0 = 2001
    for run_id=id0,id0+500,100 do begin
        rdatafile = 'data/ph_0000.dat'
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,8
        openr,1,rdatafile
        readf,1,diskbody,sigtau_es,Tdisk,em_top,em_bot,ref_top,ref_bot
        close,1
        rdatafile = 'data/db_0000.dat'
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,8
        openr,1,rdatafile
        readf,1,diskbody_ijk
        close,1
        corpow_ijk = fltarr(N1,N2,N3)
        readcp = fltarr(N3,N2,N1)
        rdatafile = 'data/scat_cpow.0000.dat'
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,15
        openr,1,rdatafile
        readf,1,readcp
        close,1
        for ir=0,N1-1 do corpow_ijk(ir,*,*) = transpose(readcp(*,*,ir))
        
        incor = where((diskbody_ijk eq 0)and(rr gt R_hor))
        rt_cor = ll_cgs
        rt_cor(*,*,*)=0.
        rt_cor(incor)=corpow_ijk(incor)*2.42d17
        nan_problem = where(finite(rt_cor,/NAN))
        if (nan_problem(0) ge 0) then rt_cor(nan_problem)=0.
        dldr0 = total(total(rt_cor,3),2)
        dldr = dldr+dldr0
        Nframes = Nframes+1
    endfor
    dldr = dldr/Nframes
    dldr001 = dldr/L001

    dldr(*)=0.
    Nframes = 0
    id0 = 2003
    for run_id=id0,id0+500,100 do begin
        rdatafile = 'data/ph_0000.dat'
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,8
        openr,1,rdatafile
        readf,1,diskbody,sigtau_es,Tdisk,em_top,em_bot,ref_top,ref_bot
        close,1
        rdatafile = 'data/db_0000.dat'
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,8
        openr,1,rdatafile
        readf,1,diskbody_ijk
        close,1
        corpow_ijk = fltarr(N1,N2,N3)
        readcp = fltarr(N3,N2,N1)
        rdatafile = 'data/scat_cpow.0000.dat'
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,15
        openr,1,rdatafile
        readf,1,readcp
        close,1
        for ir=0,N1-1 do corpow_ijk(ir,*,*) = transpose(readcp(*,*,ir))
        
        incor = where((diskbody_ijk eq 0)and(rr gt R_hor))
        rt_cor = ll_cgs
        rt_cor(*,*,*)=0.
        rt_cor(incor)=corpow_ijk(incor)*2.42d17
        nan_problem = where(finite(rt_cor,/NAN))
        if (nan_problem(0) ge 0) then rt_cor(nan_problem)=0.
        dldr0 = total(total(rt_cor,3),2)
        dldr = dldr+dldr0
        Nframes = Nframes+1
    endfor
    dldr = dldr/Nframes
    dldr01 = dldr/L01

    dldr(*)=0.
    Nframes = 0
    id0 = 2005
    for run_id=id0,id0+500,100 do begin
        rdatafile = 'data/ph_0000.dat'
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,8
        openr,1,rdatafile
        readf,1,diskbody,sigtau_es,Tdisk,em_top,em_bot,ref_top,ref_bot
        close,1
        rdatafile = 'data/db_0000.dat'
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,8
        openr,1,rdatafile
        readf,1,diskbody_ijk
        close,1
        corpow_ijk = fltarr(N1,N2,N3)
        readcp = fltarr(N3,N2,N1)
        rdatafile = 'data/scat_cpow.0000.dat'
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,15
        openr,1,rdatafile
        readf,1,readcp
        close,1
        for ir=0,N1-1 do corpow_ijk(ir,*,*) = transpose(readcp(*,*,ir))
        
        incor = where((diskbody_ijk eq 0)and(rr gt R_hor))
        rt_cor = ll_cgs
        rt_cor(*,*,*)=0.
        rt_cor(incor)=corpow_ijk(incor)*2.42d17
        nan_problem = where(finite(rt_cor,/NAN))
        if (nan_problem(0) ge 0) then rt_cor(nan_problem)=0.
        dldr0 = total(total(rt_cor,3),2)
        dldr = dldr+dldr0
        Nframes = Nframes+1
    endfor
    dldr = dldr/Nframes
    dldr1 = dldr/L1

    dldr1 = smooth(dldr1,7)*r(1)/dr(1)
    dldr1(where(dldr1 lt 0))=-1d-8/dldr1(where(dldr1 lt 0))
    dldr01 = smooth(dldr01,7)*r(1)/dr(1)
    dldr01(where(dldr01 lt 0))=-1d-8/dldr01(where(dldr01 lt 0))
    dldr001 = smooth(dldr001,7)*r(1)/dr(1)
    dldr001(where(dldr001 lt 0))=-1d-8/dldr001(where(dldr001 lt 0))

    plot_oi,r,alog10(dldr1),xrange=[2,60],xstyle=1,yrange=[-8,0],$
      ytickname=['-1.0','-10!U-2!N','0','10!U-2','1.0'],$
      xtitle='r/M',ytitle='dL/d(log r)/L!Lseed!N',thick=thk
    oplot,r,alog10(dldr001),thick=thk,color=60
    oplot,r,alog10(dldr01),thick=thk,color=120
    oplot,r,alog10(dldr1),thick=thk,color=180
    oplot,r,alog10(r*0+1d-4),thick=thk,linestyle=2
    plots,[15,25],[1,1]*(-0.9),thick=thk,color=60
    plots,[15,25],[1,1]*(-1.5),thick=thk,color=120
    plots,[15,25],[1,1]*(-2.1),thick=thk,color=180
    xyouts,27,-1.0,'L/L!LEdd!N=0.01'
    xyouts,45,-1.6,'0.1'
    xyouts,45,-2.2,'1.0'

    stop
endif


if (ifig eq 9) then begin
    N1 = 0
    N2 = 0
    N3 = 0
    aa = 0.0
    t_frame = 0.
    R_hor = 2.0

    run_id = 2301
    rdatafile = 'data/gr_0000.dat'
    dumpstr = string(run_id,format='(I4.4)')
    strput,rdatafile,dumpstr,8
    openr,1,rdatafile

    readf,1,N1,N2,N3
    r = fltarr(N1)
    t = fltarr(N2)
    p = fltarr(N3)
    readf,1,r,t,p
    close,1
    
    dr = deriv(r)
    dt = deriv(t)
    dp = deriv(p)

    ll_cgs = fltarr(N1,N2,N3)
    wgt = fltarr(N1,N2,N3)
    rho_cgs = fltarr(N1,N2,N3)
    diskbody_ijk = intarr(N1,N2,N3)
    Te = fltarr(N1,N2,N3)
    tau_es = fltarr(N1,N2+1,N3)
    tauc = fltarr(N1,N2,N3)
    rr = fltarr(N1,N2,N3)
    drr = fltarr(N1,N2,N3)
    tt = fltarr(N1,N2,N3)
    pp = fltarr(N1,N2,N3)
    u1 = fltarr(N1,N2,N3)
    u2 = fltarr(N1,N2,N3)
    u3 = fltarr(N1,N2,N3)
    dV = fltarr(N1,N2,N3)
    g_ = fltarr(N1,N2,4,4)

    for ir=0,N1-1 do rr(ir,*,*)=r(ir)
    for ir=0,N1-1 do drr(ir,*,*)=dr(ir)
    for it=0,N2-1 do tt(*,it,*)=t(it)
    for ip=0,N3-1 do pp(*,*,ip)=p(ip)
    dV = rr^2*sin(tt)*drr*dt(0)*dp(0)

    sigtau_es = fltarr(N1,N3)
    netflux = fltarr(N1,N3)
    diskflux = fltarr(N1,N3)
    em_top = fltarr(N1,N3)
    em_bot = fltarr(N1,N3)
    ref_top = fltarr(N1,N3)
    ref_bot = fltarr(N1,N3)
    diskbody = fltarr(N1,N3)
    Tdisk = fltarr(N1,N3)
    ir_start = min(where(r gt R_hor))

    vr_r = fltarr(N1)
    vt_r = fltarr(N1)
    vp_r = fltarr(N1)
    vr_t = fltarr(N1)
    vt_t = fltarr(N1)
    vp_t = fltarr(N1)
    vr2_r = fltarr(N1)
    vt2_r = fltarr(N1)
    vp2_r = fltarr(N1)
    Tdisk_r = fltarr(N1)
    delr = 10

    for id0 = 2001,2005,2 do begin
        Nframes=0.
        Tdisk_r(*)=0.
        for run_id=id0,id0+500,100 do begin
            rdatafile = 'data/ph_0000.dat'
            dumpstr = string(run_id,format='(I4.4)')
            strput,rdatafile,dumpstr,8
            openr,1,rdatafile
            readf,1,diskbody,sigtau_es,Tdisk,em_top,em_bot,ref_top,ref_bot
            close,1
            rdatafile = 'data/u1_0000.dat'
            strput,rdatafile,dumpstr,8
            openr,1,rdatafile
            readf,1,u1
            close,1
            rdatafile = 'data/u2_0000.dat'
            strput,rdatafile,dumpstr,8
            openr,1,rdatafile
            readf,1,u2
            close,1
            rdatafile = 'data/u3_0000.dat'
            strput,rdatafile,dumpstr,8
            openr,1,rdatafile
            readf,1,u3
            close,1
            rdatafile = 'data/rh_0000.dat'
            strput,rdatafile,dumpstr,8
            openr,1,rdatafile
            readf,1,rho_cgs
            close,1
            rdatafile = 'data/ll_0000.dat'
            strput,rdatafile,dumpstr,8
            openr,1,rdatafile
            readf,1,ll_cgs
            close,1
            rdatafile = 'data/ta_0000.dat'
            strput,rdatafile,dumpstr,8
            openr,1,rdatafile
            readf,1,tau_es
            close,1
            for it=0,N2-1 do tauc(*,it,*)=0.5*(tau_es(*,it,*)+tau_es(*,it+1,*))
            zz = abs(rr*cos(tt))
            r_plane = rr(*,0,*)
            diskavgr = where((r_plane gt 6)and(r_plane lt 60))
            H_dens = (total(rho_cgs*zz*dV,2)/total(rho_cgs*dV,2))/rr(*,0,*)
            H_dens = mean(H_dens(diskavgr))
            H_diss = (total(ll_cgs*zz*dV,2)/total(ll_cgs*dV,2))/rr(*,0,*)
            H_diss = mean(H_diss(diskavgr))
;            H_diss = total(ll_cgs*zz*dV)/total(ll_cgs*dV)
            print,H_dens,H_diss,H_diss/H_dens
;            for ir=delr,N1-1-delr do begin
;                atmr = where((rr gt r(ir-delr))and(rr lt r(ir+delr))and $
;                             (tauc gt 0.2)and(tauc lt 1.0))
;                if (atmr(0) ge 0) then begin
;                    vv = u1(atmr)
;                    vavg = mean(vv)
;                    vr2_r(ir) = (mean((vv-vavg)^2.))
;                    vv = u2(atmr)*rr(atmr)
;                    vavg = mean(vv)
;                    vt2_r(ir) = (mean((vv-vavg)^2.))
;                    vv = u3(atmr)*rr(atmr)*sin(tt(atmr))
;                    vavg = mean(vv)
;                    vp2_r(ir) = (mean((vv-vavg)^2.))
;                endif
;            endfor
            vr_t = vr_t+vr2_r
            vt_t = vt_t+vt2_r
            vp_t = vp_t+vp2_r
            Nframes = Nframes+1.
            print,run_id
            Tdisk_r = Tdisk_r+total(Tdisk,2)/N3
        endfor
        vr_t = vr_t/Nframes
        vt_t = vt_t/Nframes
        vp_t = vp_t/Nframes
        Tdisk_r = Tdisk_r/Nframes
        vtot_r = sqrt(vr_t+vt_t+vp_t)
        E_r = 0.5*vtot_r^2.
        if (id0 eq 2001) then begin
            plot_oo,r,E_r,xrange=[2,60],xstyle=1,xtitle='r/M',ytitle='E/m!Le!Nc!U2!N',$
              thick=thk,yrange=[1d-4,0.1]
            plots,[15,25],[1,1]*10.^(-1.2),thick=thk,color=60
            plots,[15,25],[1,1]*10.^(-1.5),thick=thk,color=120
            plots,[15,25],[1,1]*10.^(-1.8),thick=thk,color=180
            xyouts,27,10.^(-1.3),'L/L!LEdd!N=0.01'
            xyouts,45,10.^(-1.6),'0.1'
            xyouts,45,10.^(-1.9),'1.0'
        endif
        oplot,r,E_r,thick=thk,color=60+30*(id0-2001)
        Eph_r = Tdisk_r/1.1d7*2.8/511.
        oplot,r,Eph_r,thick=thk,linestyle=2,color=60+30*(id0-2001)
    endfor
  stop
endif    

if (fix(ifig) eq 10) then begin
    Nspec = 101
    id0 = 3003
    Nt = 41
    Nframes = 11
    idex15 = 0
    idex30 = 3
    idex45 = 6
    idex60 = [8,9,10,11,12]
    idex75 = 15
    spec = dblarr(Nspec,Nt,Nframes)
    spec2 = dblarr(Nspec,Nt)

    Nbands = 4
    b1lo = 0.1
    b1hi = 1.
    b2lo = 1.
    b2hi = 10.
    b3lo = 10.
    b3hi = 100.
    b4lo = 100.
    b4hi = 10000.
    b1lo = 0.1
    b1hi = 3.
    b2lo = 3.
    b2hi = 10.
    b3lo = 10.
    b3hi = 30.
    b4lo = 30.
    b4hi = 10000.

    lcurve = dblarr(Nbands,Nframes)
    emin = 0.001              
    emax = 10000.
    nu = emin*10.^(findgen(Nspec)/(Nspec-1.)*alog10(emax/emin))
    dnu = deriv(nu)
    b1dex = where((nu gt b1lo)and(nu le b1hi))
    b2dex = where((nu gt b2lo)and(nu le b2hi))
    b3dex = where((nu gt b3lo)and(nu le b3hi))
    b4dex = where((nu gt b4lo)and(nu le b4hi))
    tframe = findgen(Nframes)*0.5+10.
    rdatafile = 'data/scat_spec.0000.dat'
    iframe=0
    for run_id = id0,id0+500,50 do begin
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,15
        openr,1,rdatafile
        readf,1,spec2
        close,1
        spec(*,*,iframe)=spec2
        iframe=iframe+1
    endfor
    for iframe = 0,Nframes-1 do begin
        lcurve(0,iframe)=total(total(spec(b1dex,idex60,iframe),2)*dnu(b1dex))
        lcurve(1,iframe)=total(total(spec(b2dex,idex60,iframe),2)*dnu(b2dex))
        lcurve(2,iframe)=total(total(spec(b3dex,idex60,iframe),2)*dnu(b3dex))
        lcurve(3,iframe)=total(total(spec(b4dex,idex60,iframe),2)*dnu(b4dex))
    endfor
    lcurve = lcurve(*,0:Nframes-1)
    ltot = total(lcurve(0:1,*),1)
    print,'total var',sqrt(variance(ltot))/mean(ltot)
    ltot = total(lcurve,1)
    tframe = tframe(0:Nframes-1)
    coeff = poly_fit(tframe,ltot,1)
    lsec = coeff(0)+tframe*coeff(1)
    for ibin = 0,3 do lcurve(ibin,*)=lcurve(ibin,*)/lsec
    plot,tframe,lcurve(0,*)/mean(lcurve(0,*)),$
      xtitle='t (1000 M)',ytitle='Intensity',xrange=[10,15],xstyle=1,$
      yrange=[0,2.2],ystyle=1
    for iband=0,Nbands-1 do begin
        oplot,tframe,lcurve(iband,*)/mean(lcurve(iband,*)),thick=thk,color=60+30*iband
        print,sqrt(variance(lcurve(iband,*)))/mean(lcurve(iband,*))
    endfor
    plots,[11.2,11.7],[1,1]*2.05,thick=thk,color=60
    plots,[11.2,11.7],[1,1]*1.9,thick=thk,color=90
    plots,[11.2,11.7],[1,1]*1.75,thick=thk,color=120
    plots,[11.2,11.7],[1,1]*1.6,thick=thk,color=150
    xyouts,11.8,2.0,'0.1-3 keV'
    xyouts,11.8,1.85,'3-10 keV'
    xyouts,11.8,1.7,'10-30 keV'
    xyouts,11.8,1.55,'>30 keV'
stop
endif

if (fix(ifig) eq 11) then begin
    Nspec = 101
    id0 = 3003
    Nt = 41
    Nph = 40
    Nframes = 1
    spec = dblarr(Nspec,Nph,Nt)
    spec2 = dblarr(Nspec,Nph,Nt)

    Nbands = 4
    b1lo = 0.1
    b1hi = 3.
    b2lo = 3.
    b2hi = 6.
    b3lo = 6.
    b3hi = 25.
    b4lo = 30.
    b4hi = 10000.

    amp = dblarr(Nbands,Nt)
    amp2 = dblarr(Nbands,Nt)

    lcurve = dblarr(Nbands,Nph)
    emin = 0.001              
    emax = 10000.
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
    Nframes = 11
    datapts = findgen(4,Nt,2*Nframes)
    rdatafile = 'data/scat_spcp.0000.dat'
    Nframes = 0
    for run_id = id0,id0+500,50 do begin
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,15
        openr,1,rdatafile
        readf,1,spec2
        close,1
        spec=spec2
        Nframes=Nframes+1
        for it = 1,Nt-2 do begin
            idex = [0,0]+it
            idex = [-1,0,1]+it
            ;idex = [-2,-1,0,1,2]+it
            for iframe = 0,Nph-1 do begin
                lcurve(0,iframe)=$
                  total(total(spec(b1dex,iframe,idex),3))
                lcurve(1,iframe)=$
                  total(total(spec(b2dex,iframe,idex),3))
                lcurve(2,iframe)=$
                  total(total(spec(b3dex,iframe,idex),3))
                lcurve(3,iframe)=$
                  total(total(spec(b4dex,iframe,idex),3))
            endfor
            lcurve(*,*)=0.25*(lcurve(*,pdex1)+lcurve(*,pdex2)+$
                              lcurve(*,pdex3)+lcurve(*,pdex4))
            for iband=0,Nbands-1 do lcurve(iband,*)=smooth(lcurve(iband,*),1)
            for iband=0,Nbands-1 do begin
                if (it lt Nt/2) then $
                  datapts(iband,it,Nframes*2-2)= $
                  sqrt((variance(lcurve(iband,*)))/mean(lcurve(iband,*))^2.)
                if (it gt Nt/2) then $
                  datapts(iband,Nt-it-1,Nframes*2-1)= $
                  sqrt((variance(lcurve(iband,*)))/mean(lcurve(iband,*))^2.)
                amp2(iband,it)=amp2(iband,it)+$
                  (variance(lcurve(iband,*)))/mean(lcurve(iband,*))^2.
            endfor
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

  
