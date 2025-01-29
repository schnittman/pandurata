PRO figures5a, ifig
!P.font = 0
!P.charsize = 1.5
thk = 4

if (ifig eq 0) then begin
  NN = [26,35,51,76,106,165,262,420,680,5000]
  err_N = [1d-5,1d-6,1d-7,1d-8,1d-9,1d-10,1d-11,1d-12,0.9d-13]
  plot_oo,NN,err_N,psym=1,thick=thk,xrange=[10,1000],xstyle=1,$
    xtitle = 'N steps',ytitle='fractional error !9d!3Q/Q'
  oplot,NN/5,(NN/5)^(-5.)*40.,thick=thk
  xyouts, 100,1.5d-8,'!9e!8 ~ N!U-5!N!3'
stop
endif

if (ifig eq 1) then begin
  N1 = 180
  N2 = 160
  N3 = 64
  run_id = 1204
  rdatafile = 'harm3d/data/gr_1203.dat'
  dumpstr = string(run_id,format='(I4.4)')
  strput,rdatafile,dumpstr,15
  openr,1,rdatafile
  readf,1,N1,N2,N3
  r = fltarr(N1)
  t = fltarr(N2)
  p = fltarr(N3)
  readf,1,r,t,p
  close,1

  em_top = fltarr(N1,N3)
  em_bot = fltarr(N1,N3)
  vertices_top = fltarr(3,N1,N3)
  vertices_bot = fltarr(3,N1,N3)
  sigtau_es = fltarr(N1,N3)
  diskbody = fltarr(N1,N3)
  Tdisk = fltarr(N1,N3)
  rdatafile = 'harm3d/data/ph_0000.dat'
  strput,rdatafile,dumpstr,15
  openr,1,rdatafile
  readf,1,diskbody,sigtau_es,Tdisk,em_top,em_bot
  close,1
  for ir=0,N1-1 do begin
      for ip=0,N3-1 do begin
          vertices_top(*,ir,ip)=r(ir)*[sin(em_top(ir,ip))*cos(p(ip)), $
                                       sin(em_top(ir,ip))*sin(p(ip)), $
                                       cos(em_top(ir,ip))]
          vertices_bot(*,ir,ip)=r(ir)*[sin(em_bot(ir,ip))*cos(p(ip)), $
                                       sin(em_bot(ir,ip))*sin(p(ip)), $
                                       cos(em_bot(ir,ip))]
      endfor
  endfor

;  win_size = 600
;  window, 0, xsize = win_size, ysize = win_size
  t3d, /reset, rot = [-70, 10, 0]
  scale3d
  bounds= 25.
  warp = bounds * 1.0
;  warp = bounds * 0.71
  zbounds = bounds * 0.71
  !X.S = [-(-bounds), 1.0]/(bounds-(-bounds))
  !Y.S = [-(-warp), 1.0]/(warp-(-warp))
  !Z.S = [-(-zbounds), 1.0]/(zbounds-(-zbounds))
  
  dx = -20.
  dz = -10.
  vertices_top(0,*,*)=vertices_top(0,*,*)+dx
  vertices_top(2,*,*)=vertices_top(2,*,*)+dz
;  vertices_bot(0,*,*)=vertices_bot(0,*,*)+dx
;  vertices_bot(*,*,*)=vertices_top(*,*,*)
  Ngrid = 40
  dgrid = 2.
  pphi = findgen(100)/99.*!PI/2
  for rg=2,Ngrid,dgrid do plots,rg*cos(pphi)+dx,rg*sin(pphi),0*sin(pphi)+dz,/t3d,/data
  for pg=0.,Ngrid,dgrid do plots,[0,Ngrid]*cos(pg/Ngrid*!PI/2.)+dx,$
    [0,Ngrid]*sin(pg/Ngrid*!PI/2.),[0,Ngrid]*0+dz,/t3d,/data
  vt = vertices_top
  maxT = max(Tdisk)
  for ir =1,N1-25 do begin
      for ip=0,N3-2 do begin
          polyfill,[vt(0,ir,ip),vt(0,ir-1,ip),vt(0,ir-1,ip+1),vt(0,ir,ip+1)], $
            [vt(1,ir,ip),vt(1,ir-1,ip),vt(1,ir-1,ip+1),vt(1,ir,ip+1)], $
            [vt(2,ir,ip),vt(2,ir-1,ip),vt(2,ir-1,ip+1),vt(2,ir,ip+1)], $
            /t3d,/data,color=Tdisk(ir,ip)/maxT*255
      endfor
  endfor
  dr = 4
  for ir =dr,N1-70,dr do begin
      for ip=0,N3-2 do begin
          plots,[vertices_top(0,ir,ip),vertices_top(0,ir,ip+1)],$
            [vertices_top(1,ir,ip),vertices_top(1,ir,ip+1)],$
            [vertices_top(2,ir,ip),vertices_top(2,ir,ip+1)],/t3d,/data,color=0
          plots,[vertices_top(0,ir,ip),vertices_top(0,ir-dr,ip)],$
            [vertices_top(1,ir,ip),vertices_top(1,ir-dr,ip)],$
            [vertices_top(2,ir,ip),vertices_top(2,ir-dr,ip)],/t3d,/data,color=0
      endfor
      plots,[vertices_top(0,ir,ip),vertices_top(0,ir-dr,ip)],$
        [vertices_top(1,ir,ip),vertices_top(1,ir-dr,ip)],$
        [vertices_top(2,ir,ip),vertices_top(2,ir-dr,ip)],/t3d,/data,color=0
  endfor
  for ir =N1-70,N1-25 do begin
      for ip=0,N3-2 do begin
          plots,[vertices_top(0,ir,ip),vertices_top(0,ir,ip+1)],$
            [vertices_top(1,ir,ip),vertices_top(1,ir,ip+1)],$
            [vertices_top(2,ir,ip),vertices_top(2,ir,ip+1)],/t3d,/data,color=0
          plots,[vertices_top(0,ir,ip),vertices_top(0,ir-1,ip)],$
            [vertices_top(1,ir,ip),vertices_top(1,ir-1,ip)],$
            [vertices_top(2,ir,ip),vertices_top(2,ir-1,ip)],/t3d,/data,color=0
      endfor
      plots,[vertices_top(0,ir,ip),vertices_top(0,ir-1,ip)],$
        [vertices_top(1,ir,ip),vertices_top(1,ir-1,ip)],$
        [vertices_top(2,ir,ip),vertices_top(2,ir-1,ip)],/t3d,/data,color=0
  endfor
  plots,[0,40]+dx,[0,0],[0,0]+dz,/t3d,/data,thick=thk
  plots,[0,0]+dx,[0,40],[0,0]+dz,/t3d,/data,thick=thk
;  plots,[0,0]+dx,[0,0],[0,40]+dz,/t3d,/data,thick=thk
  stop
END


if ((ifig eq 8.1)or(ifig eq 8.2)) then begin
  Nspec = 501.
  Nth = 81.
  Ispec = fltarr(Nspec,Nth)

  if (ifig eq 8.1) then run_id = 1002
  if (ifig eq 8.2) then run_id = 1003
  rdatafile = 'scat_spec.0000.dat'
  dumpstr = string(run_id,format='(I4.4)')
  strput,rdatafile,dumpstr,10
  openr,1,rdatafile
  readf,1,Ispec
  close,1

  rdata = fltarr(2,Nspec)
  rdatafile = 'kerrspec_1001_i15.dat'
  strput,rdatafile,dumpstr,9
  openr,1,rdatafile
  readf,1,rdata
  close,1
  rtspec15 = rdata(1,*)
  rdatafile = 'kerrspec_1001_i45.dat'
  strput,rdatafile,dumpstr,9
  openr,1,rdatafile
  readf,1,rdata
  close,1
  rtspec45 = rdata(1,*)
  rdatafile = 'kerrspec_1001_i75.dat'
  strput,rdatafile,dumpstr,9
  openr,1,rdatafile
  readf,1,rdata
  close,1
  rtspec75 = rdata(1,*)

  nu = fltarr(Nspec)
  dnu = nu
  ;linear scale
  e_min = 0.0
  e_max = 2.0
  nu(0) = e_min
  for j=1.,Nspec-1 do begin
    nu(j) = e_min+j/(Nspec-1)*(e_max-e_min) ;
    dnu(j-1) = nu(j)-nu(j-1)
  endfor

  idex0 = 80
  idex15 = 78
  idex30 = 71
  idex45 = 57
  idex60 = 40
  idex75 = 21
  Ispec15 = Ispec(*,idex15)/total(Ispec(*,idex15)*dnu)
  Ispec45 = Ispec(*,idex45)/total(Ispec(*,idex45)*dnu)
  Ispec75 = Ispec(*,idex75)/total(Ispec(*,idex75)*dnu)
  peak15 = where(Ispec15 eq max(Ispec15))
  peak45 = where(Ispec45 eq max(Ispec45))
  peak75 = where(Ispec75 eq max(Ispec75))
  plot,nu,Ispec(*,idex15)/$
    total(Ispec(*,idex15)*dnu),$
    xtitle='E!Lobs!N/E!Lem!N',ytitle='Normalized Intensity',$
    xrange=[0.,1.5],yrange=[0,6]
  oplot,nu,Ispec15,thick=thk,color = 120
  oplot,nu,Ispec45,thick=thk,color = 120
  oplot,nu,Ispec75,thick=thk,color = 120
  oplot,rdata(0,*),rtspec15/total(rtspec15*deriv(rdata(0,*))),thick=thk
  oplot,rdata(0,*),rtspec45/total(rtspec45*deriv(rdata(0,*))),thick=thk
  oplot,rdata(0,*),rtspec75/total(rtspec75*deriv(rdata(0,*))),thick=thk
  
  xyouts,nu(peak15),Ispec15(peak15)+0.3,'i=15!Uo!N'
  xyouts,nu(peak45),Ispec45(peak45)+0.3,'45!Uo!N'
  xyouts,nu(peak75),Ispec75(peak75)+0.3,'75!Uo!N'
  if (ifig eq 8.1) then xyouts,0.1,5.4,'a/M=0'
  if (ifig eq 8.2) then xyouts,0.1,5.4,'a/M=0.99'
  plots,[0.1,0.3],[1,1]*5,thick=thk
  plots,[0.1,0.3],[1,1]*4.6,thick=thk,color=120
  xyouts,0.32,4.9,'obs-to-em'
  xyouts,0.32,4.5,'em-to-obs'
  stop
endif

if (ifig eq 9) then begin
  Nspec = 501.
  Nth = 41.
  spec = fltarr(Nspec,Nth)
  I_lo = fltarr(Nspec,Nth)
  I_hi = fltarr(Nspec,Nth)
  nu = fltarr(Nspec)
  dnu = fltarr(Nspec)
  Nruns = 6
  runids = [1009,1004,1005,1006,1007,1008]
  NN = [3.,10.,30.,100.,300.,1000.]
  err_N = fltarr(Nruns)
  NN = (NN+1.)^2.*501.
  rdatafile = 'scat_spec.0000.dat'
  ;linear scale
  e_min = 0.0
  e_max = 2.0
  nu(0) = e_min
  for j=1.,Nspec-1 do begin
    nu(j) = e_min+j/(Nspec-1)*(e_max-e_min) ;
    dnu(j-1) = nu(j)-nu(j-1)
  endfor

  dnuth = dnu#(fltarr(Nth)+1.)
  dumpstr = string(runids(Nruns-1),format='(I4.4)')
  strput,rdatafile,dumpstr,10
  openr,1,rdatafile
  readf,1,spec
  close,1
  I_hi = spec
;  plot,nu,spec(*,Nth/10.)
  for rundex = 0,Nruns-2 do begin
      dumpstr = string(runids(rundex),format='(I4.4)')
      strput,rdatafile,dumpstr,10
      openr,1,rdatafile
      readf,1,spec
      close,1
      I_lo = spec
      err = sqrt(total((I_lo-I_hi)^2.*dnuth))/sqrt(total((I_hi)^2.*dnuth))
      print,err
      err_N(rundex)=err
;      oplot,nu,spec(*,Nth/10),psym=rundex+1
  endfor
;stop
  plot_oo,NN(0:Nruns-2),err_N(0:Nruns-2),psym=4,thick=thk,$
    xtitle = 'N(photons)',ytitle='fractional error'
  oplot,NN/5,1./sqrt(NN/5)*220,thick=thk
  xyouts, 1e4,3,'!9e!8 ~ N!U-1/2!N!3'
  oplot,[1,1]*3d6,[1,1]*10^(0.5),psym=2,thick=thk
  oplot,[1,1]*3d6,[1,1]*10^(0.7),psym=4,thick=thk
  xyouts,4d6,10^(0.45),'obs-to-emit'
  xyouts,4d6,10^(0.65),'emit-to-obs'
  stop
endif

if (ifig eq 9.1) then begin
  Nspec = 501.
  Ninc = 3
  spec = fltarr(Nspec,Ninc)
  I_lo = fltarr(Nspec,Ninc)
  I_hi = fltarr(Nspec,Ninc)
  nu = fltarr(Nspec)
  dnu = fltarr(Nspec)
  rdata = fltarr(2,Nspec)
  Nruns = 6
  runids = [0010,0030,0100,0300,1000,3000]
  runinc = [15,45,75]
  NN = [10.,30.,100.,300.,1000.,3000.]
  err_N = fltarr(Nruns)
  NN = (NN+1.)^2.*40
  rdatafile = 'kerrspec_N0000_i15.dat'
  ;linear scale
  e_min = 0.0
  e_max = 2.0
  nu(0) = e_min
  for j=1.,Nspec-1 do begin
    nu(j) = e_min+j/(Nspec-1)*(e_max-e_min) ;
    dnu(j-1) = nu(j)-nu(j-1)
  endfor
  dnuth = dnu#(fltarr(Ninc)+1.)

  for idex = 0,Ninc-1 do begin
      dumpstr = string(runids(Nruns-1),format='(I4.4)')
      strput,rdatafile,dumpstr,10
      dumpstr = string(runinc(idex),format='(I2.2)')
      strput,rdatafile,dumpstr,16
      openr,1,rdatafile
      readf,1,rdata
      close,1
      spec(*,idex)=rdata(1,*)/total(rdata(1,*))
  endfor
  I_hi = spec
;  plot,nu,I_hi(*,0)
  for rundex = 0,Nruns-2 do begin
      for idex = 0,Ninc-1 do begin
          dumpstr = string(runids(rundex),format='(I4.4)')
          strput,rdatafile,dumpstr,10
          dumpstr = string(runinc(idex),format='(I2.2)')
          strput,rdatafile,dumpstr,16
          openr,1,rdatafile
          readf,1,rdata
          close,1
          spec(*,idex)=rdata(1,*)/total(rdata(1,*))
      endfor
      I_lo = spec
      err = sqrt(total((I_lo-I_hi)^2.*dnuth))/sqrt(total((I_hi)^2.*dnuth))
      print,err
      err_N(rundex)=err
      oplot,nu,I_lo(*,0),psym=rundex+1
  endfor
;stop
;  plot_oo,NN(0:Nruns-2),err_N(0:Nruns-2),psym=1,thick=thk,$
;    xtitle = 'N(photons)',ytitle='fractional error'
  oplot,NN(0:Nruns-2),err_N(0:Nruns-2),psym=2,thick=thk
;  oplot,NN/5,1./sqrt(NN/5)*100,thick=thk
;  xyouts, 2e4,1.5,'!9e!8 ~ N!U-1/2!N!3'
  
  stop
endif


if (fix(ifig) eq 10) then begin
    N = 401.
    aa = 0.99
    Ixy = dblarr(N,N)
    wght = dblarr(N,N)
    rdata = dblarr(10,N,N)

    if (aa eq 0.0) then openr,1,'pol400_00_75.dat'
    if (aa eq 0.99) then begin
        if ((ifig eq 10.4)or(ifig eq 10.5)or(ifig eq 10.1)or(ifig eq 10.2)) then $
          openr,1,'pol400_99_75.dat'
        if ((ifig eq 10.41)or(ifig eq 10.51)) then openr,1,'pol400_99_60.dat'
        if ((ifig eq 10.42)or(ifig eq 10.52)) then openr,1,'pol400_99_45.dat'
    endif
    readf,1,rdata
    close,1

    if (aa eq 0.0) then run_id = 1010
    if (aa eq 0.99) then run_id = 1011
    if ((ifig eq 10.4)or(ifig eq 10.5)or(ifig eq 10.2)) then incdex = 10
    if ((ifig eq 10.41)or(ifig eq 10.51)) then incdex = 20
    if ((ifig eq 10.42)or(ifig eq 10.52)) then incdex = 29
    Nr = 500.
    Rmax = max(rdata(1,*,*))*1.01
    Rmax = 20.
;    aa = 0.9
    deleff = 0.0
    M = 10.
;    Mdot = 9.26
    Mdot = 5.47
    calc_nt_tr,M,Mdot,aa,deleff,Nr,Rmin,Rmax,T_r
    dr = (Rmax-Rmin)/(Nr-1)
    rr = findgen(Nr)*dr+Rmin

    ;reverse images left-to-right
    idex = indgen(N)
    Ixy = dblarr(N,N)
    for i=0,9 do rdata(i,*,*)=rdata(i,idex,N-idex-1)

                                ;convert from line intensity to bb flux 
    I_NT = (rdata(3,*,*))^(4./3.)
    D_fact = fltarr(N,N)
    D_fact(*,*)=rdata(4,*,*)
    for ix=0,N-1 do begin
        for iy=0,N-1 do begin
            ir = (rdata(1,ix,iy)-Rmin)/dr
            if ((ir ge 0)and(ir lt Nr)) then I_NT(0,ix,iy)=I_NT(0,ix,iy)*T_r(ir)^4.*D_fact(ix,iy)
        endfor
    endfor
    D_fact=transpose(D_fact)
    Ixy(*,*)=transpose(I_NT(0,*,*))
    movmax = max(Ixy)
    data = byte(255*(alog10(Ixy/movmax+1e-4)+4.01)/4.)
    N2 = (fix(600./N))*N
    Eph = rdata(3,*,*)^(1./3.)
    psi = atan(rdata(9,*,*),rdata(8,*,*))
    deg = sqrt(rdata(8,*,*)^2.+rdata(9,*,*)^2.)

    if (ifig eq 10.1) then begin
        enlarge = rebin(data,N2,N2)
        erase
        tvscl,enlarge,0,0,/Device
        dd = 150.
;  maxIxy = max(Ixy)
        step =25.
        for i=N/step,N-1,N/step do begin
            for j=N/step,N-1,N/step do begin
                x0 = j/(N-1.)*12600.
                y0 = i/(N-1.)*12600.
                dff = 100.*dd*([rdata(8,i,j),rdata(9,i,j),0])
                if (rdata(3,i,j) gt 0) then $
                  plots,[x0-dff(0)/2.,x0+dff(0)/2.],[y0-dff(1)/2.,y0+dff(1)/2.],$
                  color=0,/Device,thick=thk
            endfor
        endfor
        plots,[10000,10000+5.*dd],[10500,10500],color=255,thick=thk,/Device
        xyouts,10000,9800,'deg=5%',color=255,/Device
        xyouts,1300,10000,'obs-to-em',color=255,/Device
        plots,[1300,11400],[4500,4500],color=255,thick=thk,/Device
        xyouts,6000,3900,'30M',color=255,/Device
    endif

    Nnu = 101.
    h = 6.63d-27
    kB = 1.38d-16
    mp = 1.67d-27
    f_hard = 1.8
    emax = 1000.
    emin = 0.001
    nu_ev = emin*10.^(dindgen(Nnu)/(Nnu-1.)*alog10(emax/emin))*1000.
    nu_hz = nu_ev*1.6d-12/h
    
    spec = dblarr(Nnu)
    specx = dblarr(Nnu)
    specy = dblarr(Nnu)
    spec(*)=0.
    for ix=0,N-1 do begin
        for iy=0,N-1 do begin
            if ((rdata(3,ix,iy) gt 0)and(rdata(1,ix,iy) lt Rmax)and(rdata(1,ix,iy) gt Rmin)) then begin
                rdex = round((rdata(1,ix,iy)-Rmin)/dr)
                Tr = T_r(rdex)*Eph(0,ix,iy)*f_hard
                xnu = nu_ev/(8.62d-5*Tr)
                f_nu = Tr^3.*xnu^3/(exp(xnu)-1.)/f_hard^4.*D_fact(iy,ix)
                spec = spec+f_nu
                specx = specx+f_nu*deg(0,ix,iy)*cos(2.*psi(0,ix,iy))
                specy = specy+f_nu*deg(0,ix,iy)*sin(2.*psi(0,ix,iy))
            endif
        endfor
    endfor
    
    if (ifig eq 10.4) then begin
        plot_oi,nu_ev/1d3,sqrt(specx^2+specy^2)/spec*100.,thick=thk,xrange=[0.1,50],xstyle=1,$
          xtitle='photon energy (keV)',ytitle='Polarization degree (%)',yrange=[0,4]
        xyouts,0.2,2.5,'i=75!Uo!N'
        xyouts,0.2,1.0,'60!Uo!N'
        xyouts,0.2,0.3,'45!Uo!N'
        plots,[10^0.,10^0.3],[1,1]*3.55,thick=thk,color=120
        plots,[10^0.05,10^0.15,10^0.25],[1,1,1]*3.25,thick=thk,psym=4
        xyouts,10^0.35,3.5,'obs-to-em'
        xyouts,10^0.35,3.2,'em-to-obs'
    endif
    if ((ifig eq 10.4)or(ifig eq 10.41)or(ifig eq 10.42)) then $
      oplot,nu_ev/1d3,sqrt(specx^2+specy^2)/spec*100.,thick=thk,color=120
    
                                ;compare with emitter5
    Nspec = 101.
    Nt = 41.
;    incdex = 20
    spec = dblarr(Nspec,Nt)
    Qspec = dblarr(Nspec,Nt)
    Uspec = dblarr(Nspec,Nt)
    deg_spec = dblarr(Nspec,Nt)
    ang_spec = dblarr(Nspec,Nt)
    emin = 0.001                ;AGN scale
    emax = 1000.
    nu = emin*10.^(findgen(Nspec)/(Nspec-1.)*alog10(emax/emin))
    
    rdatafile = 'scat_spec.0000.dat'
    dumpstr = string(run_id,format='(I4.4)')
    strput,rdatafile,dumpstr,10
    openr,1,rdatafile
    readf,1,spec,Qspec,Uspec
    close,1
    deg_spec = sqrt(Qspec^2+Uspec^2)/spec
    ang_spec = atan(Uspec,Qspec)/2.*!radeg
    if ((ifig eq 10.4)or(ifig eq 10.41)or(ifig eq 10.42)) then begin
        oplot,nu,deg_spec(*,incdex)*100.,thick=thk,psym=4
    endif
    
    if (ifig eq 10.5) then begin
        plot_oi,nu_ev/1d3,(atan(specy,specx)/2.*!radeg),thick=thk,xrange=[0.1,50],$
          xstyle=1,xtitle='photon energy (keV)',ytitle='Polarization angle (deg)',$
          yrange=[-90,90]
        xyouts,20,-20,'i=75!Uo!N'
        xyouts,20,-60,'60!Uo!N'
        xyouts,20,-90,'45!Uo!N'
        plots,[10^(-0.6),10^(-0.3)],[1,1]*80,thick=thk,color=120
        plots,[10^(-0.55),10^(-0.45),10^(-0.34)],[1,1,1]*70,thick=thk,psym=4
;        plots,[0.2,0.5],[1,1]*80,thick=thk
;        plots,[0.2,0.5],[1,1]*70,thick=thk,psym=3
        xyouts,0.6,78,'obs-to-em'
        xyouts,0.6,68,'em-to-obs'
        oplot,nu,ang_spec(*,incdex),psym=4,thick=thk
        oplot,nu_ev/1d3,atan(specy,specx)/2.*!radeg,thick=thk,color=120
    endif
    if ((ifig eq 10.51)or(ifig eq 10.52)) then begin
        oplot,nu_ev/1d3,atan(specy,specx)/2.*!radeg,thick=thk,color=120
        oplot,nu,ang_spec(*,incdex),psym=4,thick=thk
    endif
    if (ifig eq 10.2) then begin
        N = 201.
        Nt = 41.
        rdata = dblarr(2,N,N)
        Ixy = dblarr(N,N)
        wght = dblarr(N,N)
        mov = dblarr(N,N,Nt)
        movx = dblarr(N,N,Nt)
        movy = dblarr(N,N,Nt)
        mov2 = dblarr(N,N,Nt)
        movx2 = dblarr(N,N,Nt)
        movy2 = dblarr(N,N,Nt)

        rdatafile = 'scat_imag.0000.dat'
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,10
        openr,1,rdatafile
        readf,1,mov
        close,1

        rdatafile = 'scat_ipol.0000.dat'
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,10
        openr,1,rdatafile
        readf,1,movx,movy
        close,1

        for it=0,Nt-1 do begin
            mov(*,*,it)=transpose(mov(*,*,it))
            movx(*,*,it)=transpose(movx(*,*,it))
            movy(*,*,it)=transpose(movy(*,*,it))
        endfor
        sortmov = mov(sort(mov))
        movmax = max(mov)
        movmax=sortmov(1.*N*N*Nt-1000.)
        outliers = where(mov gt movmax)
        if (outliers(0) lt 0) then outliers = 0
        movx(outliers)=movx(outliers)/mov(outliers)*movmax
        movy(outliers)=movy(outliers)/mov(outliers)*movmax
        mov(outliers)=movmax

        it = incdex
        Ixy = dblarr(N,N)
        Xpol = Ixy
        Ypol = Ixy
        Ixy(*,*) = mov(*,*,it)
        Xpol(*,*) = movx(*,*,it)
        Ypol(*,*) = movy(*,*,it)
        Ixy = Ixy+1d-10
        psi = atan(Ypol/Ixy,Xpol/Ixy)/2.
        deg = sqrt((Ypol/Ixy)^2.+(Xpol/Ixy)^2.)
        Xpol = deg*cos(psi)
        Ypol = deg*sin(psi)
        data = byte(255*(alog10(Ixy/movmax+10^(-3.5))+3.51)/3.5)
        N2 = (fix(600./N))*N
        enlarge = rebin(data,N2,N2)
        erase
        tvscl,enlarge,0,0,/Device
        dd = 150.
        step=25.
        maxIxy = max(Ixy)
        for i=N/step,N-1,N/step do begin
            for j=N/step,N-1,N/step do begin
                x0 = j/(N-1.)*12600.
                y0 = i/(N-1.)*12600.
                dff = 100.*dd*([Xpol(j,i),Ypol(j,i),0])
                if (Ixy(j,i) gt 1d-4*maxIxy) then $
                  plots,[x0-dff(0)/2.,x0+dff(0)/2.],[y0-dff(1)/2.,y0+dff(1)/2.],$
                  color=0,/Device,thick=thk
            endfor
        endfor
;        plots,[10000,10000+5.*dd],[11500,11500],color=255,thick=thk,/Device
;        xyouts,10000,10900,'deg=5%',color=255,/Device
        plots,[10000,10000+5.*dd],[10500,10500],color=255,thick=thk,/Device
        xyouts,10000,9800,'deg=5%',color=255,/Device
        xyouts,1300,10000,'em-to-obs',color=255,/Device
        plots,[1300,11400],[4500,4500],color=255,thick=thk,/Device
        xyouts,6000,3900,'30M',color=255,/Device
    endif

endif

if ((ifig eq 11.1)or(ifig eq 11.2)) then begin
    Nspec = 101.
    Nt = 21.
    spec = dblarr(Nspec,Nt)
    Qspec = dblarr(Nspec,Nt)
    Uspec = dblarr(Nspec,Nt)
    deg_spec = dblarr(Nspec,Nt)
    ang_spec = dblarr(Nspec,Nt)
    shortx = acos((findgen(Nt)+0.5)/Nt)*!radeg
    chandra_I = [0.41441,0.47490,0.52397,0.57001,0.61439,0.65770,$
                 0.70029,0.74234,0.78398,0.82530,0.86637,$
                 0.90722,0.94789,0.98842,1.02882,1.06911,$
                 1.10931,1.14943,1.18947,1.22945,1.26938]
    chandra_d = [0.11713,0.08979,0.07448,0.06311,0.05410,0.04667,$
                 0.04041,0.03502,0.03033,0.02619,0.02252,$
                 0.01923,0.01627,0.01358,0.011123,0.008880,$
                 0.006818,0.004919,0.003155,0.001522,0.0]
    ccth = findgen(21)/20.
    cth = findgen(21)/21.+0.5/21.
    dcth = deriv(cth)
;LOG ENERGY SPACING
    emin = 0.001
    emax = 10000.
    nu = emin*10.^(findgen(Nspec)/(Nspec-1.)*alog10(emax/emin))
    dnu = deriv(nu)
    dnuth = dnu#(fltarr(Nt)+1.)

    runids = [20,21,22,23]
    psm = [1,2,4,5]
    if (ifig eq 11.1) then begin
        plot,cth,chandra_I,xtitle='cos !8i!3',ytitle='Intensity',thick=thk
        plots,[1,1.5,2]*0.1,[1,1,1]*1.4,psym=psm(3),thick=thk
        plots,[1,1.5,2]*0.1,[1,1,1]*1.33,psym=psm(2),thick=thk
        plots,[1,1.5,2]*0.1,[1,1,1]*1.26,psym=psm(1),thick=thk
        plots,[1,1.5,2]*0.1,[1,1,1]*1.19,psym=psm(0),thick=thk
        plots,[1,1.5,2]*0.1,[1,1,1]*1.12,thick=thk
        xyouts,0.22,1.38,'!9t!3=5'
        xyouts,0.22,1.31,'!9t!3=3'
        xyouts,0.22,1.24,'!9t!3=2'
        xyouts,0.22,1.17,'!9t!3=1'
        xyouts,0.22,1.10,'Chandra limit'
    endif
    if (ifig eq 11.2) then begin
        plot,cth,chandra_d*100,xtitle='cos !8i!3',ytitle='Polarization (%)',$
          yrange=[0,15],thick=thk
        plots,[6,6.5,7]*0.1,[1,1,1]*14,psym=psm(3),thick=thk
        plots,[6,6.5,7]*0.1,[1,1,1]*13.3,psym=psm(2),thick=thk
        plots,[6,6.5,7]*0.1,[1,1,1]*12.6,psym=psm(1),thick=thk
        plots,[6,6.5,7]*0.1,[1,1,1]*11.9,psym=psm(0),thick=thk
        plots,[6,6.5,7]*0.1,[1,1,1]*11.2,thick=thk
        xyouts,0.72,13.8,'!9t!3=5'
        xyouts,0.72,13.1,'!9t!3=3'
        xyouts,0.72,12.4,'!9t!3=2'
        xyouts,0.72,11.7,'!9t!3=1'
        xyouts,0.72,11.0,'Chandra limit'
    endif
    for rdex = 0,3 do begin
        run_id = runids(rdex)
        rdatafile = 'scat_spec.0000.dat'
        dumpstr = string(run_id,format='(I4.4)')
        strput,rdatafile,dumpstr,10
        openr,1,rdatafile
        readf,1,spec,Qspec,Uspec
        close,1
        deg_spec = sqrt(Qspec^2+Uspec^2)/spec
        ang_spec = atan(Uspec,Qspec)/2.
        
        spec_th = total(spec*dnuth,1)/cth
        spec_th = spec_th/spec_th(10)*chandra_I(10)
        if (ifig eq 11.1) then $
          oplot,cth,spec_th,psym=psm(rdex),thick=thk
        if (ifig eq 11.2) then $
;          oplot,cth,deg_spec(Nspec/2,*)*100,psym=psm(rdex),thick=thk
        oplot,cth,sqrt(total(Qspec*dnuth,1)^2+total(Uspec*dnuth,1)^2)/total(spec*dnuth,1)*100,psym=psm(rdex),thick=thk
    endfor
    stop
endif

if ((ifig eq 12.1)or(ifig eq 12.2)or(ifig eq 12.3)or(ifig eq 12.4)) then begin
    Nt = 21.
    Nspec = 101.
    run_id = 99
    run_sort = -1

    spec = dblarr(Nspec,Nt)
    Qspec = dblarr(Nspec,Nt)
    Uspec = dblarr(Nspec,Nt)
    spec_s = dblarr(6,Nspec,Nt)
    Qspec_s = dblarr(6,Nspec,Nt)
    Uspec_s = dblarr(6,Nspec,Nt)
    spec2 = dblarr(Nspec,Nt)
    Qspec2 = dblarr(Nspec,Nt)
    Uspec2 = dblarr(Nspec,Nt)
    deg_spec = dblarr(Nspec,Nt)
    ang_spec = dblarr(Nspec,Nt)

    rdatafile = 'scat_spec.0000.dat'
    dumpstr = string(run_id,format='(I4.4)')
    strput,rdatafile,dumpstr,10
    openr,1,rdatafile
    readf,1,spec,Qspec,Uspec
    for isort=0,5 do begin
        readf,1,spec2,Qspec2,Uspec2
        spec_s(isort,*,*)=spec2
        Qspec_s(isort,*,*)=Qspec2
        Uspec_s(isort,*,*)=Uspec2
    endfor
    close,1

    emin = 0.001
    emax = 10000.
    nu = emin*10.^(findgen(Nspec)/(Nspec-1.)*alog10(emax/emin))
    x = alog10(nu/511.)
    dnu = deriv(nu)
    if (ifig eq 12.1) then begin
        plot_io,x,nu*spec(*,1),xrange=[-5,1],xstyle=1,yrange=[1d-14,1d-8],$
          xtitle='log E/(m!Le!Nc!U2!N)',ytitle = 'E F!LE!N',thick=thk
        oplot,x,nu*spec_s(0,*,1)*10,thick=thk,linestyle=1
        oplot,x,nu*spec_s(1,*,1),thick=thk,linestyle=2
        oplot,x,nu*spec_s(2,*,1),thick=thk,linestyle=3
        oplot,x,nu*spec_s(3,*,1),thick=thk,linestyle=4
        oplot,x,nu*spec_s(5,*,1),thick=thk,linestyle=5
        xyouts,-3,10^(-8.5),'cos !8i!3 = 0.11'
        xyouts,-3,10^(-9.),'T!Lbb!N=10 eV'
        xyouts,-1,10^(-8.5),'!9t!3=0.50'
        xyouts,-1,10^(-9.),'T!Le!N=56 keV'
        
    endif
    if (ifig eq 12.2) then begin
        plot_io,x,nu*spec(*,10),xrange=[-5,1],xstyle=1,yrange=[1d-14,1d-8],$
          xtitle='log E/(m!Le!Nc!U2!N)',ytitle = 'E F!LE!N',thick=thk
        oplot,x,nu*spec_s(0,*,10),thick=thk,linestyle=1
        oplot,x,nu*spec_s(1,*,10),thick=thk,linestyle=2
        oplot,x,nu*spec_s(2,*,10),thick=thk,linestyle=3
        oplot,x,nu*spec_s(3,*,10),thick=thk,linestyle=4
        oplot,x,nu*spec_s(5,*,10),thick=thk,linestyle=5
        xyouts,-3,10^(-8.5),'cos !8i!3 = 0.50'
    endif
    if (ifig eq 12.3) then begin
        plot,x,-Qspec(*,1)/spec(*,1)*100.,xrange=[-5,1],yrange=[-15,30],$
          xstyle=1,ystyle=1,$
          xtitle='log E/(m!Le!Nc!U2!N)',ytitle = '!9d!3 (%)',thick=thk
;        oplot,x,-Qspec_s(0,*,1)/spec_s(0,*,1)*100,thick=thk,linestyle=1
;        oplot,x,-Qspec_s(1,*,1)/spec_s(1,*,1)*100,thick=thk,linestyle=2
;        oplot,x,-Qspec_s(2,*,1)/spec_s(2,*,1)*100,thick=thk,linestyle=3
;        oplot,x,-Qspec_s(3,*,1)/spec_s(3,*,1)*100,thick=thk,linestyle=4
;        oplot,x,-Qspec_s(5,*,1)/spec_s(5,*,1)*100,thick=thk,linestyle=5
    endif
    if (ifig eq 12.4) then begin
        plot,x,-Qspec(*,10)/spec(*,10)*100.,xrange=[-5,1],yrange=[-15,30],$
          xstyle=1,ystyle=1,$
          xtitle='log E/(m!Le!Nc!U2!N)',ytitle = '!9d!3 (%)',thick=thk
;        oplot,x,-Qspec_s(0,*,10)/spec_s(0,*,10)*100,thick=thk,linestyle=1
;        oplot,x,-Qspec_s(1,*,10)/spec_s(1,*,10)*100,thick=thk,linestyle=2
;        oplot,x,-Qspec_s(2,*,10)/spec_s(2,*,10)*100,thick=thk,linestyle=3
;        oplot,x,-Qspec_s(3,*,10)/spec_s(3,*,10)*100,thick=thk,linestyle=4
;        oplot,x,-Qspec_s(5,*,10)/spec_s(5,*,10)*100,thick=thk,linestyle=5
    endif
endif

if ((ifig eq 13.1)or(ifig eq 13.2)or(ifig eq 13.3)or(ifig eq 13.4)) then begin
    Nt = 21.
    Nspec = 101.
    run_id = 199
    run_sort = -1

    spec = dblarr(Nspec,Nt)
    Qspec = dblarr(Nspec,Nt)
    Uspec = dblarr(Nspec,Nt)
    spec_s = dblarr(6,Nspec,Nt)
    Qspec_s = dblarr(6,Nspec,Nt)
    Uspec_s = dblarr(6,Nspec,Nt)
    spec2 = dblarr(Nspec,Nt)
    Qspec2 = dblarr(Nspec,Nt)
    Uspec2 = dblarr(Nspec,Nt)
    deg_spec = dblarr(Nspec,Nt)
    ang_spec = dblarr(Nspec,Nt)

    rdatafile = 'scat_spec.0000.dat'
    dumpstr = string(run_id,format='(I4.4)')
    strput,rdatafile,dumpstr,10
    openr,1,rdatafile
    readf,1,spec,Qspec,Uspec
    for isort=0,5 do begin
        readf,1,spec2,Qspec2,Uspec2
        spec_s(isort,*,*)=spec2
        Qspec_s(isort,*,*)=Qspec2
        Uspec_s(isort,*,*)=Uspec2
    endfor
    close,1

    emin = 0.001
    emax = 10000.
    nu = emin*10.^(findgen(Nspec)/(Nspec-1.)*alog10(emax/emin))
    x = alog10(nu/511.)
    dnu = deriv(nu)
    if (ifig eq 13.1) then begin
        plot_io,x,nu*spec(*,1),xrange=[-5,1],xstyle=1,yrange=[1d-11,1d-8],$
          xtitle='log E/(m!Le!Nc!U2!N)',ytitle = 'E F!LE!N',thick=thk
        oplot,x,nu*spec_s(0,*,1),thick=thk,linestyle=1
        oplot,x,nu*spec_s(1,*,1),thick=thk,linestyle=2
        oplot,x,nu*spec_s(2,*,1),thick=thk,linestyle=3
        oplot,x,nu*spec_s(3,*,1),thick=thk,linestyle=4
        oplot,x,nu*spec_s(5,*,1),thick=thk,linestyle=5
        xyouts,-3,10^(-8.5),'cos !8i!3 = 0.11'
        xyouts,-3,10^(-8.75),'T!Lbb!N=10 eV'
        xyouts,-1,10^(-8.5),'!9t!3=0.05'
        xyouts,-1,10^(-8.75),'T!Le!N=352 keV'
        
    endif
    if (ifig eq 13.2) then begin
        plot_io,x,nu*spec(*,10),xrange=[-5,1],xstyle=1,yrange=[1d-11,1d-8],$
          xtitle='log E/(m!Le!Nc!U2!N)',ytitle = 'E F!LE!N',thick=thk
        oplot,x,nu*spec_s(0,*,10),thick=thk,linestyle=1
        oplot,x,nu*spec_s(1,*,10),thick=thk,linestyle=2
        oplot,x,nu*spec_s(2,*,10),thick=thk,linestyle=3
        oplot,x,nu*spec_s(3,*,10),thick=thk,linestyle=4
        oplot,x,nu*spec_s(5,*,10),thick=thk,linestyle=5
        xyouts,-3,10^(-8.5),'cos !8i!3 = 0.50'
    endif
    if (ifig eq 13.3) then begin
        plot,x,-Qspec(*,1)/spec(*,1)*100.,xrange=[-5,1],yrange=[-5,10],$
          xstyle=1,ystyle=1,$
          xtitle='log E/(m!Le!Nc!U2!N)',ytitle = '!9d!3 (%)',thick=thk
;        oplot,x,-Qspec_s(0,*,1)/spec_s(0,*,1)*100,thick=thk,linestyle=1
;        oplot,x,-Qspec_s(1,*,1)/spec_s(1,*,1)*100,thick=thk,linestyle=2
;        oplot,x,-Qspec_s(2,*,1)/spec_s(2,*,1)*100,thick=thk,linestyle=3
;        oplot,x,-Qspec_s(3,*,1)/spec_s(3,*,1)*100,thick=thk,linestyle=4
;        oplot,x,-Qspec_s(5,*,1)/spec_s(5,*,1)*100,thick=thk,linestyle=5
    endif
    if (ifig eq 13.4) then begin
        plot,x,-Qspec(*,10)/spec(*,10)*100.,xrange=[-5,1],yrange=[-5,10],$
          xstyle=1,ystyle=1,$
          xtitle='log E/(m!Le!Nc!U2!N)',ytitle = '!9d!3 (%)',thick=thk
;        oplot,x,-Qspec_s(0,*,10)/spec_s(0,*,10)*100,thick=thk,linestyle=1
;        oplot,x,-Qspec_s(1,*,10)/spec_s(1,*,10)*100,thick=thk,linestyle=2
;        oplot,x,-Qspec_s(2,*,10)/spec_s(2,*,10)*100,thick=thk,linestyle=3
;        oplot,x,-Qspec_s(3,*,10)/spec_s(3,*,10)*100,thick=thk,linestyle=4
;        oplot,x,-Qspec_s(5,*,10)/spec_s(5,*,10)*100,thick=thk,linestyle=5
    endif
endif

if (ifig eq 14) then begin
    Nt = 21.
    Nspec = 101.
    run_id = 30
    run_sort = -1

    spec = dblarr(Nspec,Nt)
    Qspec = dblarr(Nspec,Nt)
    Uspec = dblarr(Nspec,Nt)
    spec_s = dblarr(6,Nspec,Nt)
    Qspec_s = dblarr(6,Nspec,Nt)
    Uspec_s = dblarr(6,Nspec,Nt)
    spec2 = dblarr(Nspec,Nt)
    Qspec2 = dblarr(Nspec,Nt)
    Uspec2 = dblarr(Nspec,Nt)
    deg_spec = dblarr(Nspec,Nt)
    ang_spec = dblarr(Nspec,Nt)

    rdatafile = 'scat_spec.0000.dat'
    dumpstr = string(run_id,format='(I4.4)')
    strput,rdatafile,dumpstr,10
    openr,1,rdatafile
    readf,1,spec,Qspec,Uspec
    for isort=0,5 do begin
        readf,1,spec2,Qspec2,Uspec2
        spec_s(isort,*,*)=spec2
        Qspec_s(isort,*,*)=Qspec2
        Uspec_s(isort,*,*)=Uspec2
    endfor
    close,1

    emin = 0.001
    emax = 10000.
    nu = emin*10.^(findgen(Nspec)/(Nspec-1.)*alog10(emax/emin))
    x = alog10(nu/511.)
    dnu = deriv(nu)
    if (ifig eq 14) then begin
        plot_io,x,nu*total(spec,2),xrange=[-5,1],xstyle=1,yrange=[1d-14,1d-8],$
          xtitle='log E/(m!Le!Nc!U2!N)',ytitle = 'E F!LE!N',thick=thk
    endif
    stop
endif


END

PRO calc_nt_tr,M,Mdot,aa,deleff,Nr,Rin,Rout,T_r
  kapp = 0.4
  ak = 7.566d-15
  kB = 1.38d-16
  c = 3.d10
  G = 6.67d-8
  mp = 1.67d-24
  sigma = 5.6705d-5
;  deleff=0.00
;  Mstar = 3.33e0
  Mstar = M/3.
  Mbh = Mstar*3.*(1.99d33)
  Mdstar = 10.e0
  Mdstar = Mdot
  Mdot = Mdstar*(1.0d17) 
;Mdot of 1e17 g/s is ~0.1 L_edd for a/M=1, Mbh = 3Msun (0.01 L_edd for a/M=0)
  alph = 1.e-1

  Z1 = 1.+(1-aa*aa)^(1./3.)*((1.+aa)^(1./3.)+(1-aa)^(1./3.))
  Z2 = sqrt(3.*aa*aa+Z1^2.)
  Risco = (3.+Z2-((3.-Z1)*(3+Z1+2.*Z2))^0.5)
  Eisco = (Risco*Risco-2*Risco+aa*sqrt(Risco))/ $
    (Risco*sqrt(Risco*Risco-3*Risco+2*aa*sqrt(Risco)))
  Lisco = sqrt(Risco)*(Risco*Risco-2.*aa*sqrt(Risco)+aa^2)/ $
    (Risco*sqrt(Risco*Risco-3.*Risco+2.*aa*sqrt(Risco)))
  Revnt = 1.001*(1.+sqrt(1.-aa*aa))
  Rin = Risco
;  Rout = 50.
  dr = (Rout-Rin)/(Nr-1.)
  rr = findgen(Nr)*dr+Rin
  r_cm = rr*Mstar*4.4d5
  Nrp = 101.
  drp = (Rin-Revnt)/(Nrp-1.)
  rrp = findgen(Nrp)*drp+Revnt
  
  ntA = 1.+aa^2./rr^2.+2.*aa^2./rr^3.
  ntB = 1.+aa/rr^1.5
  ntC = 1.-3./rr+2.*aa/rr^1.5
  ntD = 1.-2./rr+aa^2./rr^2.
  ntE = 1.+4.*aa^2./rr^2.-4.*aa^2./rr^3.+3.*aa^4./rr^4.
  ntE2 = 2.+5.*aa^2./rr^2.-2.*aa^2./rr^3.+3.*aa^4./rr^4.
  ntF = 1.-2.*aa/rr^1.5+aa^2./rr^2.
  ntJ = dblarr(Nr)
  ntSig = 0.75*sqrt(Mbh*G/r_cm^3.)*ntD/ntC
  Nu = 101.
  for i=0,Nr-1 do begin
      u_star = 1./rr(i)
      du = u_star/(Nu-1.)
      subu = findgen(Nu)*du
      ntJ(i) = exp(1.5*du*total((1.-2.*aa*subu^1.5+aa^2.*subu^2.)/ $
                                ((1.+aa*subu^1.5)*(1.-3.*subu+2.*aa*subu^1.5))))
  endfor
  
  rhC = 1.-4.*aa/rr^1.5+3.*aa^2/rr^2
  rhB = 1.-3./rr+2.*aa/rr^1.5
  
  L_ms = 2./sqrt(3.)*(3.*sqrt(Risco)-2.*aa)/sqrt(Risco)
  ntL = ntF/sqrt(ntC)-L_ms/sqrt(rr)
  ntQ = dblarr(Nr)
  for i=0,Nr-1 do begin
      Q_int = ntL(0:i)*ntF(0:i)/(ntB(0:i)*ntC(0:i)*ntJ(0:i)*rr(0:i)^1.5)
      ntQ(i) = ntL(i)-3.*ntJ(i)/(2.*sqrt(rr(i)))* $
        dr*(total(Q_int(0:i))-0.5*(Q_int(0)+Q_int(i)))
  endfor
  ntQ0 = ntQ
  
  Rak1 = ntQ0/(ntB*sqrt(ntC))
  Rak2 = Risco^1.5*sqrt(ntC(0))*deleff/ntC/sqrt(rr)
  ntQ = ntB*sqrt(ntC)*(Rak1+Rak2)
  nt_flux = 6.d25*(Mdstar/Mstar^2.)*ntQ/(rr^3.*ntB*sqrt(ntC))
  nt_sdens = 20.*(Mstar/Mdstar/alph)*rr^1.5*ntB^3*sqrt(ntC)*ntE/(ntA^2.*ntQ)
  vr = -Mdot/(2.*!PI*r_cm*nt_sdens*sqrt(ntD))/c
  nt_Ts = 6.8d8*((alph^2.*Mstar^(-10.)*Mdstar^8.)* $
                 (ntA^8.*ntQ^8.)/(rr^15.*ntB^16.*ntC^2.*ntD^2.*ntE^4.))^(1./9.)
  nt_Ts(0)=1d-10
  Lum = 4.*!PI*total(nt_flux*rr*sqrt(ntD/ntA-0.*aa^2/rr^4/ntA))*dr*(4.4d5*Mstar)^2
  Ledd = 3.6d38*Mstar
  vr0=vr(0)
  T_r = (nt_flux/sigma)^0.25+1.
;stop
END

PRO calc_plunge,aa,rr,hz_rp,A_rp,vr0

  sz = size(rr)
  Nr = sz(1)
  M = 1.0d
  Re = 1.+sqrt(1.-aa*aa)
  g_dn00 = dblarr(Nr)
  g_dn03 = dblarr(Nr)
  g_dn33 = dblarr(Nr)
  g_up00 = dblarr(Nr)
  g_up03 = dblarr(Nr)
  g_up33 = dblarr(Nr)
  v_coor = dblarr(4,Nr)
  v_zamo = dblarr(4,Nr)

  Z1 = 1.+(1-aa*aa)^(1./3.)*((1.+aa)^(1./3.)+(1-aa)^(1./3.))
  Z2 = sqrt(3.*aa*aa+Z1^2.)
  Risco = (3.+Z2-((3.-Z1)*(3+Z1+2.*Z2))^0.5)
  print,aa,Risco
  Rout = Risco
  dr = (Rout-Re)/(Nr-1.)
  rr = dindgen(Nr)*dr+Re
  ntA = 1.+aa^2./rr^2.+2.*aa^2./rr^3.
  ntD = 1.-2./rr+aa^2./rr^2.
  vr0 = -4.296329d-4
  vr0 = -0.0007d
;  vr0=0.
  Eisco = (Risco*Risco-2*M*Risco+aa*sqrt(M*Risco))/ $
    (Risco*sqrt(Risco*Risco-3*M*Risco+2*aa*sqrt(M*Risco)))
  Lisco = sqrt(M*Risco)*(Risco*Risco-2*aa*sqrt(M*Risco)+aa*aa)/ $
    (Risco*sqrt(Risco*Risco-3*M*Risco+2*aa*sqrt(M*Risco)))
;  Lisco=Lisco*0.99
  g_dn00 = -ntD/ntA+(2.*M*aa/rr^2.)^2/ntA
  g_dn03 = -2.*M*aa/rr
  g_dn33 = rr^2*ntA
  g_up00 = -ntA/ntD
  g_up03 = -2.*M*aa/(rr^3*ntD)
  g_up33 = 1./(rr^2*ntA)-(2.*M*aa/rr^3)^2/(ntA*ntD)

  D0 = ntD(Nr-1)
  A0 = ntA(Nr-1)
  r0 = rr(nr-1)
  E0 = D0/A0*(2.*M*aa/(r0^3.*D0)*Lisco+ $
    sqrt(4.*M*M*aa*aa/(r0^6.*D0^2)*Lisco^2+A0/D0*(vr0^2/D0+ $
         (1./(r0^2*A0)-4.*aa^2/(r0^6.*A0*D0))*Lisco^2 + 1.)))

  vr = -sqrt(-ntD*(1.+g_up00*E0^2-2.*g_up03*E0*Lisco+g_up33*Lisco^2))
  dvdr = 1./(2.*vr)*(-2./rr^2+2.*aa^2/rr^3+E0^2*(-2.*aa^2/rr^3-6.*aa^2/rr^4) $
                    +12.*M*aa*E0*Lisco/rr^4 $
         -Lisco^2*((-2./rr^3+6./rr^4-4.*aa^2/rr^5+24.*M^2*aa^2/rr^7)/ntA $
           +(ntD/rr^2-(2.*M*aa/rr^3)^2)*(2.*aa^2/rr^3+6.*aa^2/rr^4)/ntA^2))
  thet_r = dvdr+2.*vr/rr
  area_r = thet_r*0.
  int_r = thet_r/vr
  for i=0,Nr-1 do begin
    area_r(i) = exp(-(total(int_r(i:Nr-1))-0.5*int_r(Nr-1)-0.5*int_r(i))*dr)
  endfor
  v_coor(0,*) = -g_up00*E0+g_up03*Lisco
  v_coor(1,*) = vr
  v_coor(3,*) = -g_up03*E0+g_up33*Lisco
  v_zamo(0,*) = sqrt(ntD/ntA)*v_coor(0,*)
  v_zamo(1,*) = v_coor(1,*)/sqrt(ntD)
  v_zamo(3,*) = -2.*M*aa/(rr^2*sqrt(ntA))*v_coor(0,*)+rr*sqrt(ntA)*v_coor(3,*)
  norm_coor = v_coor(0,*)^2*g_dn00+v_coor(1,*)^2/ntD+v_coor(3,*)^2*g_dn33+$
              2.*v_coor(0,*)*v_coor(3,*)*g_dn03
  norm_zamo = -v_zamo(0,*)^2+v_zamo(1,*)^2+v_zamo(3,*)^2
  zz = rr^2*aa^2*ntD/(rr^2+aa^2)^2
  S_r = 3.*sqrt(zz)/(1.-zz)
  reim_r = M/rr^3*(v_zamo(0,*)^2*(1.+2.*zz)/(1.-zz) $
    -2.*v_zamo(0,*)*v_zamo(3,*)*S_r-v_zamo(1,*)^2 $
    +v_zamo(3,*)^2*(2.+zz)/(1.-zz))
  reim_r = reim_r/reim_r(Nr-1)
  hz_r = (1./(area_r*reim_r^3.))^(1./7.)
  sig_r = 1./area_r
  hz_rp = hz_r
  a_rp = area_r
;  stop
END

