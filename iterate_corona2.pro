PRO iterate_corona2,run_id,nxt_id,ilevel,rmsw,aa
;  aa = 0.
  N1 = 0
  N2 = 0
  N3 = 0

;  rdatafile = 'datt/gr_0000.dat'
  rdatafile = 'data/gr_0000.dat'
  dumpstr = string(run_id,format='(I4.4)')
  strput,rdatafile,dumpstr,8
  openr,1,rdatafile
  readf,1,N1,N2,N3
  r = dblarr(N1)
  t = dblarr(N2)
  t_ = dblarr(N2+1)
  p = dblarr(N3)
  readf,1,r,t,p
  close,1

  T_cor = dblarr(N1,N2,N3)
  rr = dblarr(N1,N2,N3)
  tt = dblarr(N1,N2,N3)
  pp = dblarr(N1,N2,N3)
  dV = dblarr(N1,N2,N3)
  rdsh_wght = dblarr(N1,N2,N3)
  g_ = dblarr(N1,N2,4,4)
  g_up = dblarr(N1,N2,4,4)
  b2_cgs = dblarr(N1,N2,N3)
  ll_cgs = dblarr(N1,N2,N3)
  diskbody_ijk = intarr(N1,N2,N3)

  sigtau_es = dblarr(N1,N3)
  netflux = dblarr(N1,N3)
  diskflux = dblarr(N1,N3)
  em_top = dblarr(N1,N3)
  em_bot = dblarr(N1,N3)
  ref_top = dblarr(N1,N3)
  ref_bot = dblarr(N1,N3)
  diskbody = dblarr(N1,N3)
  Tdisk = dblarr(N1,N3)

  rdatafile = 'data/ph_0000.dat'
  dumpstr = string(run_id,format='(I4.4)')
  strput,rdatafile,dumpstr,8
  openr,1,rdatafile
  readf,1,diskbody,sigtau_es,Tdisk,em_top,em_bot,ref_top,ref_bot
  close,1

  Tavg = total(Tdisk,2)/N3
  sigavg = total(sigtau_es,2)/N3
  wdata = dblarr(3,N1)
  wdata(0,*)=r
  wdata(1,*)=Tavg
  wdata(2,*)=sigavg
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
  rdatafile = 'data/ll_0000.dat'
  dumpstr = string(run_id,format='(I4.4)')
  strput,rdatafile,dumpstr,8
  openr,1,rdatafile
  readf,1,ll_cgs
  close,1

  T_cor0 = T_cor

  dr = deriv(r)
  dt = deriv(t)
  dp = deriv(p)
  dA = (r*dr)#(dblarr(N3)+!PI/2./N3)
  for ir=0,N1-1 do rr(ir,*,*)=r(ir)
  for it=0,N2-1 do tt(*,it,*)=t(it)
  for ip=0,N3-1 do pp(*,*,ip)=p(ip)

  R_hor = 1.+sqrt(1.-aa*aa)
  Sigma = rr^2.+aa^2.*cos(tt)^2.
  Delta = rr^2.-2.*rr+aa^2.
  alpha = sqrt(Sigma*Delta/(Sigma*Delta+2.*rr*(aa^2.+rr^2.)))
  omega = 2.*rr*aa/(Sigma*Delta+2.*rr*(aa^2.+rr^2.))
  varph = sqrt((Sigma*Delta+2.*rr*(aa^2.+rr^2.))/Sigma*sin(tt)^2.)
  rdsh_wght = alpha
  rdsh_wght(where(rr le 2))=0.

  g_(*,*,0,0) = -alpha(*,*,0)^2.+omega(*,*,0)^2.*varph(*,*,0)^2.
  g_(*,*,0,3) = -omega(*,*,0)*varph(*,*,0)^2.
  g_(*,*,1,1) = Sigma(*,*,0)/Delta(*,*,0)
  g_(*,*,2,2) = Sigma(*,*,0)
  g_(*,*,3,3) = varph(*,*,0)^2.
  g_(*,*,3,0) = g_(*,*,0,3)

  g_up(*,*,0,0) = -((rr(*,*,0)^2.+aa^2.)^2.-aa^2.*Delta(*,*,0)*sin(tt(*,*,0))^2.)/ $
    (Sigma(*,*,0)*Delta(*,*,0))                   
  g_up(*,*,1,1) = Delta(*,*,0)/Sigma(*,*,0)
  g_up(*,*,2,2) = 1./Sigma(*,*,0)
  g_up(*,*,3,3) = (Delta(*,*,0)-aa^2.*sin(tt)^2.)/ $
    (Sigma(*,*,0)*Delta(*,*,0)*sin(tt(*,*,0))^2.)
  g_up(*,*,0,3) = -2*aa*rr(*,*,0)/(Sigma(*,*,0)*Delta(*,*,0))
  g_up(*,*,3,0) = g_up(*,*,0,3)

  t_(0:N2-1)=t-0.5*dt(0)
  t_(N2)=t(N2-1)+0.5*dt(0)
  for ip=0,N3-1 do dV(*,*,ip)=(dr#dt)*rr(*,*,ip)^2.*sin(tt(*,*,ip))*dp(ip)
  detg = dblarr(N1,N2)
  detg(*,*) = g_(*,*,0,0)*g_(*,*,1,1)*g_(*,*,2,2)*g_(*,*,3,3) - $
    g_(*,*,0,3)*g_(*,*,1,1)*g_(*,*,2,2)*g_(*,*,0,3)
  for ip=0,N3-1 do dV(*,*,ip) = sqrt(-detg)*(dr#dt)*dp(0)

;  L_Ledd = 0.1
  kappa = 0.4
  G_N = 6.6726d-8
  M_BH = 10.*(2.d33)
;  M_BH = 1.d7*(2.d33)
  cc = 3d10
  sigma_SB = 5.6705d-5
  sigma_T = 6.6525d-25
  mp = 1.6726d-24
  me = 9.1094d-28
  za = 7.5657d-15
  kB = 1.3807d-16
  kB_ev = 8.6173d-5
  qe = 4.8032d-10

  gam20 = (3./2.*kB*T_cor0/(me*cc*cc)+1.)^2.

  t_cgs = G_N*M_BH/cc^3.
  r_cgs = G_N*M_BH/cc^2.
  dV_cgs = dV*r_cgs^3.
  cor_scaler = dblarr(N1)
;  rdata = dblarr(2,N1)
;  openr,1,'scale_corona_power.dat'
;  readf,1,rdata
;  close,1
;  cor_scaler(*)=rdata(1,*)

  corpow_ijk = dblarr(N1,N2,N3)
  readcp = dblarr(N3,N2,N1)
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
;stop
  f_step = 0.5
  if (ilevel eq -3) then begin
      R_s = 0.2
      th_s = 0.5
      ph_s = 0.5
  endif
  if (ilevel eq -2) then begin
      R_s = 0.3
      th_s = 1.0
      ph_s = 1.0
  endif
  if (ilevel eq -1) then begin
      R_s = 0.5
      th_s = 1.5
      ph_s = 1.5
  endif
  if (ilevel eq 0) then begin
      R_s = 0.5
      th_s = 1.5
      ph_s = 1.5
  endif
  if (ilevel eq 1) then begin
      R_s = 0.3
      th_s = 1.0
      ph_s = 1.0
  endif
  if (ilevel eq 2) then begin
      R_s = 0.2
      th_s = 0.5
      ph_s = 0.5
  endif
  if (ilevel eq 5) then begin
      R_s = 0.1
      th_s = 0.1
      ph_s = 0.1
  endif

  rms_threshold = 0.00
;stop
  cor_scale = dblarr(N1,N2,N3)+1.
  nan_problem = where(finite(rt_cor,/NAN))
  
  if (nan_problem(0) ge 0) then rt_cor(nan_problem)=harm_cor(nan_problem)
  harm_norml = harm_cor
  if (ilevel le 5) then begin
      smooth_corona,harm_cor,harm_norml,r,t,p,0.2,0.5,0.5
      smooth_corona,harm_cor,harm_smooth,r,t,p,R_s,th_s,ph_s
      smooth_corona,rt_cor,rt_smooth,r,t,p,R_s,th_s,ph_s
      harm_smooth = double(harm_smooth)
      harm_smooth = harm_smooth+1d-10*mean(harm_smooth)
      rt_smooth = double(rt_smooth)
      rt_smooth = rt_smooth+1d-10*mean(rt_smooth)
      rms = total((harm_smooth-rt_smooth)^2./harm_norml)/total(harm_smooth)
      rmsw = total((harm_smooth-rt_smooth)^2./harm_norml*rdsh_wght)/ $
        total(harm_smooth*rdsh_wght)
      dldrharm = total(total(harm_smooth,3),2)
      dldrrt = total(total(rt_smooth,3),2)
      cor_scaler = dldrharm/dldrrt
      cor_scaler(where(r lt 3)) = mean(cor_scaler(where((r gt 3)and(r lt 4))))
      cor_scaler(where(r gt 50)) = mean(cor_scaler(where((r gt 45)and(r lt 50))))
      if (ilevel ge 0) then cor_scale = harm_smooth/rt_smooth
      if (ilevel lt 0) then $
        cor_scale = (harm_smooth*f_step+rt_smooth*(1.-f_step))/rt_smooth
  endif
stop
  off_scale_hi = where(cor_scale gt 2)
  off_scale_lo = where(cor_scale lt 0)
  if (off_scale_hi(0) ge 0) then cor_scale(off_scale_hi)=2.
  if (off_scale_lo(0) ge 0) then cor_scale(off_scale_lo)=0.1
  gam21 = gam20
  T_cor1 = T_cor0
  ;for ir=0,N1-1 do gam21(ir,*,*)=(gam20(ir,*,*)-1.)*cor_scaler(ir)+1.

  gam21=(gam20-1.)*cor_scale+1.
  T_cor1 = (sqrt(gam21)-1.)*(2./3.*me*cc^2/kB)
;  plot_oo,rr(incor),cor_scale(incor),psym=3
;  if (rms lt 0.2) then stop
;stop
  v_bulk = dblarr(N1)
  eph_r = dblarr(N1)

  harm_cor = harm_smooth
  rt_cor = rt_smooth
  dldrharm = total(total(harm_cor,3),2)
  dldrrt = total(total(rt_cor,3),2)
  ymax = max(dldrharm/dr*r)
  plot_oo,r,dldrharm/dr*r,xrange=[2,60],yrange=[ymax*0.1,ymax*10.],xstyle=1,$
    xtitle='R/M',ytitle='dL/d(log R)'
  oplot,r,dldrrt/dr*r,color=120
  print,run_id,rms,rmsw

  if (run_id ne nxt_id) then begin
      dumpstr_in = string(run_id,format='(I4.4)')
      dumpstr_out = string(nxt_id,format='(I4.4)')
      cpstr = 'cp data/ph_0000.dat data/ph_0000.dat'
      strput,cpstr,dumpstr_in,11
      strput,cpstr,dumpstr_out,28
      spawn,cpstr
      cpstr = 'cp data/gr_0000.dat data/gr_0000.dat'
      strput,cpstr,dumpstr_in,11
      strput,cpstr,dumpstr_out,28
      spawn,cpstr
      cpstr = 'cp data/rh_0000.dat data/rh_0000.dat'
      strput,cpstr,dumpstr_in,11
      strput,cpstr,dumpstr_out,28
      spawn,cpstr
      cpstr = 'cp data/u0_0000.dat data/u0_0000.dat'
      strput,cpstr,dumpstr_in,11
      strput,cpstr,dumpstr_out,28
      spawn,cpstr
      cpstr = 'cp data/u1_0000.dat data/u1_0000.dat'
      strput,cpstr,dumpstr_in,11
      strput,cpstr,dumpstr_out,28
      spawn,cpstr
      cpstr = 'cp data/u2_0000.dat data/u2_0000.dat'
      strput,cpstr,dumpstr_in,11
      strput,cpstr,dumpstr_out,28
      spawn,cpstr
      cpstr = 'cp data/u3_0000.dat data/u3_0000.dat'
      strput,cpstr,dumpstr_in,11
      strput,cpstr,dumpstr_out,28
      spawn,cpstr
      cpstr = 'cp data/ta_0000.dat data/ta_0000.dat'
      strput,cpstr,dumpstr_in,11
      strput,cpstr,dumpstr_out,28
      spawn,cpstr
      cpstr = 'cp data/db_0000.dat data/db_0000.dat'
      strput,cpstr,dumpstr_in,11
      strput,cpstr,dumpstr_out,28
      spawn,cpstr
      cpstr = 'cp data/bb_0000.dat data/bb_0000.dat'
      strput,cpstr,dumpstr_in,11
      strput,cpstr,dumpstr_out,28
      spawn,cpstr
      cpstr = 'cp data/ll_0000.dat data/ll_0000.dat'
      strput,cpstr,dumpstr_in,11
      strput,cpstr,dumpstr_out,28
      spawn,cpstr
  endif

  if (run_id eq nxt_id) then begin
      wdatafile = 'data/te_0000.dat'
      dumpstr = string(nxt_id,format='(I4.4)')
      strput,wdatafile,dumpstr,8
      openw,1,wdatafile
      printf,1,float(T_cor1)
      close,1
  endif

stop
END

PRO smooth_corona,cor,corsm,r,t,p,R_s,th_s,ph_s

  sz = size(r)
  N1 = sz(1)
  sz = size(t)
  N2 = sz(1)
  sz = size(p)
  N3 = sz(1)

  k_r = dblarr(N1)
  k_th = dblarr(N2)
  k_ph = dblarr(N3)
  corsm_r = cor*0.
  corsm_t = cor*0.
  corsm_p = cor*0.
  eta = alog(r/min(r))
  for ir=0,N1-1 do begin
      eta0 = eta(ir)
      k_r = exp(-(eta-eta0)^2./(2.*R_s^2.))
      k_r = k_r/total(k_r)
;      if (ir eq 0) then plot,eta,k_r
;      if ((ir gt 0)and((ir mod 10) eq 0)) then oplot,eta,k_r
      for irp = 0,N1-1 do $
        if (k_r(irp) gt 1d-6) then corsm_r(ir,*,*)=corsm_r(ir,*,*)+k_r(irp)*cor(irp,*,*)
  endfor

  for it=0,N2-1 do begin
      k_th(*) = 0.
      for itp=0,N2-1 do begin
          if ((abs(t(it)-t(itp)) lt th_s)and(cos(t(it))*cos(t(itp)) gt 0)) then $
            k_th(itp)=cos(!PI*(t(it)-t(itp))/th_s)+1.
      endfor
      k_th = k_th/total(k_th)
      for itp=0,N2-1 do $
        if (k_th(itp) gt 1d-6) then $
        corsm_t(*,it,*)=corsm_t(*,it,*)+k_th(itp)*corsm_r(*,itp,*)
;      if (it eq 0) then plot,t/(!PI),k_th
;      if ((it gt 0)and((it mod 8) eq 0)) then oplot,t/(!PI),k_th
;      if (it eq N2-1) then oplot,t/(!PI),k_th
  endfor
;stop

  for ip=0,N3-1 do begin
      k_ph(*) = 0.
      for ipp=0,N3-1 do begin
          if (abs(p(ip)-p(ipp)) lt ph_s) then $
            k_ph(ipp)=cos(!PI*(p(ip)-p(ipp))/ph_s)+1.
          if ((!PI/2-p(ipp)+p(ip)) lt ph_s) then $
            k_ph(ipp)=k_ph(ipp)+cos(!PI*(p(ip)-(p(ipp)-!PI/2))/ph_s)+1.
          if ((!PI/2-p(ip)+p(ipp)) lt ph_s) then $
            k_ph(ipp)=k_ph(ipp)+cos(!PI*(p(ipp)-(p(ip)-!PI/2))/ph_s)+1.
      endfor
      k_ph = k_ph/total(k_ph)
      for ipp=0,N3-1 do $
        if (k_ph(ipp) gt 1d-6) then $
        corsm_p(*,*,ip)=corsm_p(*,*,ip)+k_ph(ipp)*corsm_t(*,*,ipp)
;      if (ip eq 0) then plot,p/(!PI/2),k_ph
;      if ((ip gt 0)and((ip mod 10) eq 0)) then oplot,p/(!PI/2),k_ph
  endfor
  corsm=corsm_p
;  stop
END
