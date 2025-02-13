PRO read_harm3d_data2,run_id_in,run_id_out,aa,L_Ledd,corona_tau
;like read_harm3d_data, but with new raw data format from scott
;Reads in raw data in ASCII format from Harm3d, interpolates it onto a
;new grid with uniform dtheta, and calculates the location of the
;photosphere, disk surface temperature, and estimates the coronal
;temperature based on inverse-Compton power balance. 

  N1 = 0
  N2 = 0
  N3 = 0
  Nsym = 4 ;Nsym-fold symmetry, copy sim data into 2PI domain
  M_sun = 14.8 ;BH mass in solar mass units
  t_frame = 0. 
  attach_NT = 0.;used for adding a Novikov-Thorne disk at large radius
  aa_NT = aa
  if (aa lt 0) then begin
;      aa_NT = -aa
      print,aa
      aa = 0.
  endif
  blnk = 'ABC'
  bdatafile = 'rawdata/RADFLUX_ASCII2_012500.dat'
  if (aa eq 0.0) then rdatafile = 'rawdata/RADFLUX_ASCII_012000.dat'
  if (aa eq 0.5) then rdatafile = 'a05data/RADFLUX_ASCII_012000.dat'
  if (aa eq 0.9) then rdatafile = 'a09data/RADFLUX_ASCII_012000.dat'
  if (aa eq 0.99) then rdatafile ='a099data/RADFLUX_ASCII_012000.dat'
  run_id = run_id_in
  dumpstr = string(run_id*10,format='(I5.5)')
  if (aa eq 0.0) then strput,rdatafile,dumpstr,23
  if (aa eq 0.5) then strput,rdatafile,dumpstr,23
  if (aa eq 0.9) then strput,rdatafile,dumpstr,23
  if (aa eq 0.99) then strput,rdatafile,dumpstr,24
  print,rdatafile
  openr,1,rdatafile
  readf,1,blnk
  readf,1,N1,N2,N3,t_frame
  
  ;set up grid points of coordinates. r,t,p are cell-centered,
  ;t_ is at cell boundary

  r = fltarr(N1)
  t = fltarr(N2)
  t_ = fltarr(N2+1)
  p = fltarr(N3)
  readf,1,blnk
  readf,1,r
  readf,1,blnk
  readf,1,t
  readf,1,blnk
  readf,1,p
  readf,1,blnk

  ;For thin H/R, re-sample r so that N1 ~ 200

  N1old = N1
  i_samp = round(N1old/200.)
;  i_samp = 1                    ;no resample
  N1 = fix(N1old/i_samp)
  N1 = fix(N1/6.)*6
  rho = fltarr(N1,N2,N3)
  uu = fltarr(N1,N2,N3)
  u0 = fltarr(N1,N2,N3)
  u1 = fltarr(N1,N2,N3)
  u2 = fltarr(N1,N2,N3)
  u3 = fltarr(N1,N2,N3)
  LL = fltarr(N1,N2,N3)  
  rr = fltarr(N1,N2,N3)
  tt = fltarr(N1,N2,N3)
  pp = fltarr(N1,N2,N3)
  dV = fltarr(N1,N2,N3)
  g_ = fltarr(N1,N2,4,4)
  g_up = fltarr(N1,N2,4,4)
  rdata = fltarr(6,N3,N2,N1old)
  readf,1,rdata
  close,1

;magnetic field data. generally omit

  bdata = fltarr(5,N3,N2,N1old)
;  openr,1,bdatafile
;  readf,1,blnk
;  readf,1,N1,N2,N3,t_frame
;  readf,1,blnk
;  readf,1,r
;  readf,1,blnk
;  readf,1,t
;  readf,1,blnk
;  readf,1,p
;  readf,1,blnk
  N1 = fix(N1old/i_samp)
  N1 = fix(N1/6.)*6
  bb = fltarr(N1,N2,N3)
;  LLu0 = fltarr(N1,N2,N3)
;  LLu1 = fltarr(N1,N2,N3)
;  LLu2 = fltarr(N1,N2,N3)
;  LLu3 = fltarr(N1,N2,N3)
;  readf,1,bdata
;  close,1

;down-sample data from N1old to N1 points in r

  rnew = fltarr(N1)
  for ir=0,N1-1 do begin
      rnew(ir)=r(ir*i_samp)
      for it=0,N2-1 do begin
          LL(ir,it,*)=rdata(0,*,it,ir*i_samp)
          rho(ir,it,*)=rdata(1,*,it,ir*i_samp)
          u0(ir,it,*)=rdata(2,*,it,ir*i_samp)
          u1(ir,it,*)=rdata(3,*,it,ir*i_samp)
          u2(ir,it,*)=rdata(4,*,it,ir*i_samp)
          u3(ir,it,*)=rdata(5,*,it,ir*i_samp)
;          LLu0(ir,it,*)=bdata(0,*,it,ir*i_samp)
;          LLu1(ir,it,*)=bdata(1,*,it,ir*i_samp)
;          LLu2(ir,it,*)=bdata(2,*,it,ir*i_samp)
;          LLu3(ir,it,*)=bdata(3,*,it,ir*i_samp)
          bb(ir,it,*)=bdata(4,*,it,ir*i_samp)
      endfor
  endfor
  r = rnew

;set up 3d coordinate arrays

  dr = deriv(r)
  dt = deriv(t)
  dp = deriv(p)
  dA = (r*dr)#(dblarr(N3)+!PI/2./N3)
  for ir=0,N1-1 do rr(ir,*,*)=r(ir)
  for it=0,N2-1 do tt(*,it,*)=t(it)
  for ip=0,N3-1 do pp(*,*,ip)=p(ip)

;define metric in Boyer-Lindquist

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

;check to see if velocity is normalized g_\mu\nu u^\mu u^\nu=-1

  nrm = g_(*,*,0,0)*u0(*,*,0)*u0(*,*,0)+ $
    2.*g_(*,*,0,3)*u0(*,*,0)*u3(*,*,0)+ $
    g_(*,*,1,1)*u1(*,*,0)*u1(*,*,0)+ $
    g_(*,*,2,2)*u2(*,*,0)*u2(*,*,0)+ $
    g_(*,*,3,3)*u3(*,*,0)*u3(*,*,0)
;stop
;interpolate onto uniform grid in theta
;inside the funnel region, set the density to floor and velocity to
;ZAMO

  newdt = !PI/N2
  newt = newdt/2.+findgen(N2)*newdt
  newtt = tt
  for it=0,N2-1 do newtt(*,it,*)=newt(it)
  Sigma = rr^2.+aa^2.*cos(newtt)^2.
  Delta = rr^2.-2.*rr+aa^2.
  alpha = sqrt(Sigma*Delta/(Sigma*Delta+2.*rr*(aa^2.+rr^2.)))
  omega = 2.*rr*aa/(Sigma*Delta+2.*rr*(aa^2.+rr^2.))
  varph = sqrt((Sigma*Delta+2.*rr*(aa^2.+rr^2.))/Sigma*sin(newtt)^2.)
  
  g_(*,*,0,0) = -alpha(*,*,0)^2.+omega(*,*,0)^2.*varph(*,*,0)^2.
  g_(*,*,0,3) = -omega(*,*,0)*varph(*,*,0)^2.
  g_(*,*,1,1) = Sigma(*,*,0)/Delta(*,*,0)
  g_(*,*,2,2) = Sigma(*,*,0)
  g_(*,*,3,3) = varph(*,*,0)^2.
  g_(*,*,3,0) = g_(*,*,0,3)

  g_up(*,*,0,0) = -((rr(*,*,0)^2.+aa^2.)^2.-aa^2.*Delta(*,*,0)*sin(newtt(*,*,0))^2.)/ $
    (Sigma(*,*,0)*Delta(*,*,0))                   
  g_up(*,*,1,1) = Delta(*,*,0)/Sigma(*,*,0)
  g_up(*,*,2,2) = 1./Sigma(*,*,0)
  g_up(*,*,3,3) = (Delta(*,*,0)-aa^2.*sin(newtt)^2.)/ $
    (Sigma(*,*,0)*Delta(*,*,0)*sin(newtt(*,*,0))^2.)
  g_up(*,*,0,3) = -2*aa*rr(*,*,0)/(Sigma(*,*,0)*Delta(*,*,0))
  g_up(*,*,3,0) = g_up(*,*,0,3)

  v_ZAMO = fltarr(N1,N2,4)
  v_ZAMO(*,*,0) = g_up(*,*,0,0)*(-alpha(*,*,0))
  v_ZAMO(*,*,3) = g_up(*,*,3,0)*(-alpha(*,*,0))
  
;interpolate onto new theta grid

  rho_new = rho
  uu_new = uu
  u0_new = u0
  u1_new = u1
  u2_new = u2
  u3_new = u3
  ll_new = ll
  bb_new = bb
  rho_new(*,0,*)=rho(*,0,*)
  uu_new(*,0,*)=uu(*,0,*)
  u0_new(*,0,*)=u0(*,0,*)
  u1_new(*,0,*)=u1(*,0,*)
  u2_new(*,0,*)=u2(*,0,*)
  u3_new(*,0,*)=u3(*,0,*)
  ll_new(*,0,*)=ll(*,0,*)
  bb_new(*,0,*)=bb(*,0,*)
  rho_new(*,N2-1,*)=rho(*,N2-1,*)
  uu_new(*,N2-1,*)=uu(*,N2-1,*)
  u0_new(*,N2-1,*)=u0(*,N2-1,*)
  u1_new(*,N2-1,*)=u1(*,N2-1,*)
  u2_new(*,N2-1,*)=u2(*,N2-1,*)
  u3_new(*,N2-1,*)=u3(*,N2-1,*)
  ll_new(*,N2-1,*)=ll(*,N2-1,*)
  bb_new(*,N2-1,*)=bb(*,N2-1,*)
  wlo_it = fltarr(N2)
  whi_it = fltarr(N2)
  rho_min = min(rho)
  uu_min = uu(where(rho eq rho_min))
  uu_min = uu_min(0)
  for it=0,N2-1 do begin
      if (newt(it) lt t(0)) then begin
          rho_new(*,it,*) = rho_min
          uu_new(*,it,*) = uu_min
          ll_new(*,it,*) = 0.
          bb_new(*,it,*) = 0.
          u0_new(*,it,*) = v_ZAMO(*,it,0)#(fltarr(N3)+1.)
          u1_new(*,it,*) = v_ZAMO(*,it,1)#(fltarr(N3)+1.)
          u2_new(*,it,*) = v_ZAMO(*,it,2)#(fltarr(N3)+1.)
          u3_new(*,it,*) = v_ZAMO(*,it,3)#(fltarr(N3)+1.)
      endif
      if (newt(it) gt t(N2-1)) then begin
          rho_new(*,it,*) = rho_min
          uu_new(*,it,*) = uu_min
          ll_new(*,it,*) = 0.
          bb_new(*,it,*) = 0.
          u0_new(*,it,*) = v_ZAMO(*,it,0)#(fltarr(N3)+1.)
          u1_new(*,it,*) = v_ZAMO(*,it,1)#(fltarr(N3)+1.)
          u2_new(*,it,*) = v_ZAMO(*,it,2)#(fltarr(N3)+1.)
          u3_new(*,it,*) = v_ZAMO(*,it,3)#(fltarr(N3)+1.)
      endif
      if ((newt(it) ge t(0))and(newt(it) le t(N2-1))) then begin
          hidex = where(t ge newt(it))
          ihi = hidex(0)
          ilo = ihi-1
          wlo = (t(ihi)-newt(it))/(t(ihi)-t(ilo))
          whi = (newt(it)-t(ilo))/(t(ihi)-t(ilo))
          wlo_it(it)=wlo
          whi_it(it)=whi
          rho_new(*,it,*)=wlo*rho(*,ilo,*)+whi*rho(*,ihi,*)
          uu_new(*,it,*)=wlo*uu(*,ilo,*)+whi*uu(*,ihi,*)
          u0_new(*,it,*)=wlo*u0(*,ilo,*)+whi*u0(*,ihi,*)
          u1_new(*,it,*)=wlo*u1(*,ilo,*)+whi*u1(*,ihi,*)
          u2_new(*,it,*)=wlo*u2(*,ilo,*)+whi*u2(*,ihi,*)
          u3_new(*,it,*)=wlo*u3(*,ilo,*)+whi*u3(*,ihi,*)
          ll_new(*,it,*)=wlo*ll(*,ilo,*)+whi*ll(*,ihi,*)
          bb_new(*,it,*)=wlo*bb(*,ilo,*)+whi*bb(*,ihi,*)
      endif
  endfor

  y = fltarr(N2)
  newy = fltarr(N2)

  rho = rho_new
  uu = uu_new
  u0 = u0_new
  u1 = u1_new
  u2 = u2_new
  u3 = u3_new
  ll = ll_new
  bb = bb_new

  ;free up memory
  rho_new = 0
  uu_new = 0
  u1_new = 0
  u2_new = 0
  u3_new = 0
  ll_new = 0
  bb_new = 0

;calculate metric on new coord grid

  t(*) = newt(*)
  dt = deriv(t)
  t_(0:N2-1)=t-0.5*dt(0)
  t_(N2)=t(N2-1)+0.5*dt(0)
  for it=0,N2-1 do tt(*,it,*)=t(it)
  for ip=0,N3-1 do dV(*,*,ip)=(dr#dt)*rr(*,*,ip)^2.*sin(tt(*,*,ip))*dp(ip)
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
  detg = dblarr(N1,N2)
  detg(*,*) = g_(*,*,0,0)*g_(*,*,1,1)*g_(*,*,2,2)*g_(*,*,3,3) - $
    g_(*,*,0,3)*g_(*,*,1,1)*g_(*,*,2,2)*g_(*,*,0,3)
  for ip=0,N3-1 do dV(*,*,ip) = sqrt(-detg)*(dr#dt)*dp(0)

;  stop
;  for ip=0,N3-1 do dV(*,*,ip)=(dr#dt)*rr(*,*,ip)^2.*sin(tt(*,*,ip))*dp(ip)

  ;fix u^t so that g_munu u^mu u^nu = -1
  za = g_(*,*,0,0)
  for ip = 0,N3-1 do begin
      zb = 2.*g_(*,*,0,3)*u3(*,*,ip)
      zc = g_(*,*,1,1)*u1(*,*,ip)^2.+g_(*,*,2,2)*u2(*,*,ip)^2.+$
        g_(*,*,3,3)*u3(*,*,ip)^2.+1.
      u0_new(*,*,ip)=(-zb-sqrt(zb*zb-4.*za*zc))/(2.*za)
  endfor
  ;set u0 inside the horizon to the value just outside
  ir_out = where(r gt R_hor)
  ir_out = ir_out(0)
  for ir = 0,N1-1 do $
    if (r(ir) lt R_hor) then $
    u0_new(ir,*,*)=u0_new(ir_out,*,*)
          
  u0 = u0_new
  u0_new = 0
  ph0 = 12
  newnrm = g_(*,*,0,0)*u0(*,*,ph0)*u0(*,*,ph0)+ $
    2.*g_(*,*,0,3)*u0(*,*,ph0)*u3(*,*,ph0)+ $
    g_(*,*,1,1)*u1(*,*,ph0)*u1(*,*,ph0)+ $
    g_(*,*,2,2)*u2(*,*,ph0)*u2(*,*,ph0)+ $
    g_(*,*,3,3)*u3(*,*,ph0)*u3(*,*,ph0)
  E_0 = -(g_(*,*,0,0)*u0(*,*,0)+g_(*,*,0,3)*u3(*,*,0))
;stop

;Mdotc is the accretion rate in harm3d units. M_BH is black hole mass
;(in gm)

  if (aa eq 0) then begin
      Mdotc = 3e-4              
      eta = 0.06                
  endif
  if (aa eq 0.5) then begin
      Mdotc = 1.03e-4              
      eta = 0.08 
  endif 
  if (aa eq 0.9) then begin
      Mdotc = 2.38e-4              
      eta = 0.15
  endif 
  if (aa eq 0.99) then begin
      Mdotc = 6.8e-5
      eta = 0.25                   
  endif 
  M_BH = M_sun*(2.d33)
  Z1 = 1.+(1-aa*aa)^(1./3.)*((1.+aa)^(1./3.)+(1-aa)^(1./3.))
  Z2 = sqrt(3.*aa*aa+Z1^2.)
  Risco = (3.+Z2-((3.-Z1)*(3+Z1+2.*Z2))^0.5)
  Eisco = (Risco*Risco-2*Risco+aa*sqrt(Risco))/ $
    (Risco*sqrt(Risco*Risco-3*Risco+2*aa*sqrt(Risco))) 
  eta = 1.-Eisco

;define bunch of constants
  kappa = 0.4
  G_N = 6.6726d-8
  cc = 3d10
  sigma_SB = 5.6705d-5
  sigma_T = 6.6525d-25
  mp = 1.6726d-24
  me = 9.1094d-28
  za = 7.5657d-15 ;=4\sigma/c
  kB = 1.3807d-16
  kB_ev = 8.6173d-5
  qe = 4.8032d-10

;convert code units to cgs
  Mdot_cgs = 1.38d38*(M_BH/(2.d33))*L_Ledd/(eta*cc^2.)
;  Mdot_cgs = 4.772d17
  t_cgs = G_N*M_BH/cc^3.
  r_cgs = G_N*M_BH/cc^2.
  rho_cgs = rho*(4.*!PI*cc^2.)/(kappa*G_N*M_BH*Mdotc)*(L_Ledd/eta)
  b2_cgs = bb*(4.*!PI*cc^2.)/(kappa*G_N*M_BH*Mdotc)*(L_Ledd/eta)*cc^2.
;  b2_cgs = rho_cgs*0.
;  openr,1,'data/bb_1001.dat'
;  readf,1,b2_cgs
;  close,1
  Ub_cgs = b2_cgs/(8.*!PI)
;  Ub_cgs = Ub_cgs*1d-8 ;estimate of reduced synchroton power, after accounting 
;                       ;for self-absorption
  ll_cgs = ll*(4.*!PI*cc^2.)/(kappa*G_N*M_BH*Mdotc)*(L_Ledd/eta)*cc^2/t_cgs
  dV_cgs = dV*r_cgs^3.
;stop
;calculate optical depth tau and flux arrays

  tau_es = fltarr(N1,N2+1,N3)
  tau_es_top = fltarr(N1,N2+1,N3)
  tau_es_bot = fltarr(N1,N2+1,N3)
  flux_cgs = fltarr(N1,N2+1,N3)
  flux_cgs_bot = fltarr(N1,N2+1,N3)
  T_cor = ll_cgs*0
  for it=1,N2 do begin
      tau_es_top(*,it,*)=tau_es_top(*,it-1,*)+ $
        dt(it-1)*rr(*,it-1,*)*kappa*rho_cgs(*,it-1,*)*r_cgs
      flux_cgs(*,it,*)=flux_cgs(*,it-1,*)+ $
        dt(it-1)*sin(t(it-1))*rr(*,it-1,*)*ll_cgs(*,it-1,*)*r_cgs
  endfor
  for it=N2-1,0,-1 do begin
      tau_es_bot(*,it,*)=tau_es_bot(*,it+1,*)+ $
        dt(it)*rr(*,it,*)*kappa*rho_cgs(*,it,*)*r_cgs
      flux_cgs_bot(*,it,*)=flux_cgs_bot(*,it+1,*)+ $
        dt(it)*sin(t(it))*rr(*,it,*)*ll_cgs(*,it,*)*r_cgs
  endfor
  sigtau_es = fltarr(N1,N3)
  netflux = fltarr(N1,N3)
  diskflux = fltarr(N1,N3)
  em_top = fltarr(N1,N3)
  em_bot = fltarr(N1,N3)
  ref_top = fltarr(N1,N3)
  ref_bot = fltarr(N1,N3)
  diskbody = fltarr(N1,N3)
  diskbody_ijk = intarr(N1,N2,N3)
  vertices_top = fltarr(3,N1,N3)
  vertices_bot = fltarr(3,N1,N3)
  ccx = fltarr(N1,N3)
  ccy = fltarr(N1,N3)
  ccx(*,*)=rr(*,0,*)*cos(pp(*,0,*))
  ccy(*,*)=rr(*,0,*)*sin(pp(*,0,*))
  sigtau_es(*,*)=tau_es_top(*,N2,*)
  netflux(*,*) = flux_cgs(*,N2,*)
  Tdisk = (0.5*netflux/sigma_SB)^0.25 ;surface temperature of disk
  Uph = 0.5*za*Tdisk^4.         ;radiation density
  Uph_ijk = fltarr(N1,N2,N3)
;  cor_scaler = fltarr(N1)
;  rdata = fltarr(2,N1)
;  openr,1,'scale_corona_power.dat'
;  readf,1,rdata
;  close,1
;  cor_scaler(*)=rdata(1,*)
;  for ir=0,N1-1 do Uph(ir,*)=Uph(ir,*)/cor_scaler(ir)
  for i=0,N1-1 do begin
      Uph(i,*)=mean(Uph(i,*)) ;smooth out rad field in azimuth
      Uph_ijk(i,*,*)=(fltarr(N2)+1)#Uph(i,*)
  endfor
  Utot_ijk = Uph_ijk+Ub_cgs
  llsyn_cgs = Ub_cgs/Utot_ijk*ll_cgs
  llcomp_cgs = ll_cgs-llsyn_cgs
  gam2 = T_cor
  T_cor2 = T_cor
  Tdisk(*,*)=0.
  n_e = rho_cgs/mp
;  corona_tau = 2.0
  photo_tau = corona_tau+1.
;stop

;calculate location of reflection surfaces (photosphere with optical
;depth =2)
  for ir=0,N1-1 do begin
      for ip=0,N3-1 do begin
          subtau = tau_es_top(ir,*,ip)
          if (sigtau_es(ir,ip) gt 2.*photo_tau) then begin
              diskbody(ir,ip) = 2.
              irtopdex = where(subtau gt photo_tau)
              irtopdex = irtopdex(0)
              ihit = irtopdex
              ilot = ihit-1
              wlot = (subtau(ihit)-photo_tau)/(subtau(ihit)-subtau(ilot))
              whit = (photo_tau-subtau(ilot))/(subtau(ihit)-subtau(ilot))
              ref_top(ir,ip)=t_(ilot)*wlot+t_(ihit)*whit
              subtau = tau_es_bot(ir,*,ip)
              irbotdex = where(subtau gt photo_tau)
              irbotdex = max(irbotdex)
              ilob = irbotdex
              ihib = ilob+1
              wlob = (subtau(ihib)-photo_tau)/(subtau(ihib)-subtau(ilob))
              whib = (photo_tau-subtau(ilob))/(subtau(ihib)-subtau(ilob))
              ref_bot(ir,ip)=t_(ilob)*wlob+t_(ihib)*whib
              if (ilot le ilob) then diskbody_ijk(ir,ilot:ilob,ip)=2.
          endif else begin
              ref_top(ir,ip)=!PI/2.
              ref_bot(ir,ip)=!PI/2.
          endelse
      endfor
  endfor

;calculate location of emission surfaces (photosphere with optical
;depth = corona_tau) This is what pandurata uses as photosphere
;Use linear interpolation to get the exact value of \theta for the
;photosphere and the flux between the corona_tau surfaces.

  for ir=0,N1-1 do begin
      for ip=0,N3-1 do begin
          subtau = tau_es_top(ir,*,ip)
          neg_z = where(subtau gt 0.5*sigtau_es(ir,ip))
          if (sigtau_es(ir,ip) gt 2.*corona_tau) then begin
              if (sigtau_es(ir,ip) le 2.*photo_tau) then diskbody(ir,ip)=1.
              irtopdex = where(subtau gt corona_tau)
              irtopdex = irtopdex(0)
              ihit = irtopdex
              ilot = ihit-1
              wlot = (subtau(ihit)-corona_tau)/(subtau(ihit)-subtau(ilot))
              whit = (corona_tau-subtau(ilot))/(subtau(ihit)-subtau(ilot))
              em_top(ir,ip)=t_(ilot)*wlot+t_(ihit)*whit
              subtau = tau_es_bot(ir,*,ip)
              irbotdex = where(subtau gt corona_tau)
              irbotdex = max(irbotdex)
              ilob = irbotdex
              ihib = ilob+1
              wlob = (subtau(ihib)-corona_tau)/(subtau(ihib)-subtau(ilob))
              whib = (corona_tau-subtau(ilob))/(subtau(ihib)-subtau(ilob))
              em_bot(ir,ip)=t_(ilob)*wlob+t_(ihib)*whib
              diskflux(ir,ip)=(flux_cgs(ir,ilob,ip)*wlob+flux_cgs(ir,ihib,ip)*whib) - $
                (flux_cgs(ir,ilot,ip)*wlot+flux_cgs(ir,ihit,ip)*whit)

              ;dump corona energy into thermal disk
              ;diskflux(ir,ip)=netflux(ir,ip)

              if (diskflux(ir,ip) lt 0) then diskflux(ir,ip)=0.
              Tdisk(ir,ip)=(0.5*diskflux(ir,ip)/sigma_SB)^0.25
              ;print,ir,ip,ilot,ilob
              if (ilot le ilob) then begin
                  ;print,ir,ip
                  ;diskbody_ijk(ir,ilot:ilob,ip)=1
                  for it=ilot,ilob do $
                    if (diskbody_ijk(ir,it,ip) ne 2) then $
                    diskbody_ijk(ir,it,ip) = 1
              endif
          endif else begin
              em_top(ir,ip)=!PI/2.
              em_bot(ir,ip)=!PI/2.
          endelse
          vertices_top(*,ir,ip)=r(ir)*[sin(em_top(ir,ip))*cos(p(ip)), $
                                       sin(em_top(ir,ip))*sin(p(ip)), $
                                       cos(em_top(ir,ip))]
          vertices_bot(*,ir,ip)=r(ir)*[sin(em_bot(ir,ip))*cos(p(ip)), $
                                       sin(em_bot(ir,ip))*sin(p(ip)), $
                                       cos(em_bot(ir,ip))]
          subtau = tau_es_top(ir,*,ip)
          neg_z = where(subtau gt 0.5*sigtau_es(ir,ip))
          if (neg_z(0) ge 0) then subtau(neg_z)=tau_es_bot(ir,neg_z,ip)
          tau_es(ir,*,ip)=subtau
          subflx = flux_cgs(ir,*,ip)
          neg_z = where(subflx gt 0.5*netflux(ir,ip))
          if (neg_z(0) ge 0) then subflx(neg_z)=flux_cgs_bot(ir,neg_z,ip)
          flux_cgs(ir,*,ip)=subflx
          T_cor(ir,*,ip)=0

;following Schnittman, Krolik, & Noble (2013), estimate coronal
;temperature by balancing IC cooling with harm3d emissivity
          for ith=0,N2-1 do begin
              if (Utot_ijk(ir,ith,ip) gt 0) then begin
                  gam2 = 3./4.*ll_cgs(ir,ith,ip)/ $
                    (cc*sigma_T*n_e(ir,ith,ip)*Utot_ijk(ir,ith,ip))+1.
                  ;for low temp
                  T_cor(ir,ith,ip)=(sqrt(gam2)-1.)*(2./3.*me*cc^2/kB)
              endif
          endfor
          if (sigtau_es(ir,ip) gt 2.*corona_tau) then $
            if (ilob ge ilot) then T_cor(ir,ilot:ilob,ip)=Tdisk(ir,ip)
      endfor
  endfor
  ;plot_surfaces,N1,N3,sigtau_es,vertices_top,vertices_bot

  incor = where((diskbody_ijk eq 0)and(rr gt R_hor))
  inatm = where((diskbody_ijk eq 1)and(rr gt R_hor))
  indisk = where((diskbody_ijk eq 2)and(rr gt R_hor))
  insim = where((diskbody_ijk le 2)and(rr gt R_hor))
  if (incor(0) eq -1) then incor = 0
  brem_pow = 1.4d-27*sqrt(T_cor)*n_e^2.
  syn_pow = 4./3.*((kB_ev*T_cor/5.11d5+1.)^2.-1.)*cc*sigma_T*Ub_cgs*n_e
  syn_pow2 = 4./3.*((kB_ev*T_cor2/5.11d5+1.)^2.-1.)*cc*sigma_T*Ub_cgs*n_e
  ;T_cor = T_cor2
  Uph_ij = total(Uph_ijk,3)/N3
  Ubb_ij = total(Ub_cgs,3)/N3
  print,total(ll_cgs(incor)*dv_cgs(incor)),total(ll_cgs(indisk)*dv_cgs(indisk)),$
    total(ll_cgs(incor)*dv_cgs(incor))/total(ll_cgs(insim)*dv_cgs(insim))

;stop
  ;THIN-DISK LIMIT
  scale_factor = 0.01
  scale_factor = 1.0
;  scale_factor = 0.5
  for ir=0,N1-1 do begin
      for ip=0,N3-1 do begin
          if (diskbody(ir,ip) ge 1) then begin
              em_top(ir,ip) = !PI/2.-scale_factor*(!PI/2.-em_top(ir,ip))
              em_bot(ir,ip) = !PI/2.+scale_factor*(em_bot(ir,ip)-!PI/2.)
;              em_top(ir,ip) = !PI/2.-scale_factor
;              em_bot(ir,ip) = !PI/2.+scale_factor
          endif
          if (diskbody(ir,ip) ge 2) then begin
              ref_top(ir,ip) = !PI/2.-scale_factor*(!PI/2.-ref_top(ir,ip))
              ref_bot(ir,ip) = !PI/2.+scale_factor*(ref_bot(ir,ip)-!PI/2.)
;              ref_top(ir,ip) = !PI/2.-scale_factor
;              ref_bot(ir,ip) = !PI/2.+scale_factor
          endif

          ;force optically thin regions to emit thermal flux
          ;if ((diskbody(ir,ip) eq 0)and(netflux(ir,ip) gt 0)) then $
          ;  diskbody(ir,ip) = 1
      endfor
  endfor
;stop
;;NOVIKOV-THORNE LIMIT
  scale_height = 0.01
  R_attach = max(r)
  M_BH = M_BH
  Mdot_cgs = Mdot_cgs
  if (attach_NT eq 1) then begin
      scale_height = 0.01
      R_attach = 0.
;      R_attach = 15.
;      calc_JS_data,rr,tt,pp,diskbody,sigtau_es,Tdisk,$
;        em_top,em_bot,ref_top,ref_bot, $
;        rho,T_cor,u0,u1,u2,u3,tau_es,diskbody_ijk,ll_cgs,dV_cgs,$
;        aa_NT,L_Ledd,M_BH,Mdot_cgs,scale_height,N1,N2,N3,R_attach
      calc_NT_data2,rr,tt,pp,diskbody,sigtau_es,Tdisk,$
        em_top,em_bot,ref_top,ref_bot, $
        rho,T_cor,u0,u1,u2,u3,tau_es,diskbody_ijk,ll_cgs,dV_cgs,$
        aa_NT,L_Ledd,M_BH,Mdot_cgs,scale_height,N1,N2,N3,R_attach
      rho_cgs = rho
      Nr2 = 2*N1
      r = dblarr(Nr2)
      r(*) = rr(*,N2/2,0)
      N1=Nr2
  endif

  corpow_ijk = fltarr(N1,N2,N3)
  readcp = fltarr(N3,N2,N1)
;  rdatafile = 'data/scat_cpow.0000.dat'
;  dumpstr = string(run_id,format='(I4.4)')
;  strput,rdatafile,dumpstr,15
;  openr,1,rdatafile
;  readf,1,readcp
;  close,1
;  for ir=0,N1-1 do corpow_ijk(ir,*,*) = transpose(readcp(*,*,ir))
;stop
  ;SMOOTH DISK SURFACE
;  em_top_avg = mean(em_top(50:150,*))
;  em_bot_avg = mean(em_bot(50:150,*))
;  ref_top_avg = mean(ref_top(50:150,*))
;  ref_bot_avg = mean(ref_bot(50:150,*))
;  Tdisk_avg = Tdisk
;  for ir=0,N1-1 do begin
;      for ip=0,N3-1 do begin
;          if (diskbody(ir,ip) ge 1) then begin
;              em_top(ir,ip) = em_top_avg
;              em_bot(ir,ip) = em_bot_avg
;          endif
;          if (diskbody(ir,ip) ge 2) then begin
;              ref_top(ir,ip) = ref_top_avg
;              ref_bot(ir,ip) = ref_bot_avg
;          endif
;      endfor
;      Tdisk_avg(ir,*) = (mean(Tdisk(ir,*)^4.))^0.25
;  endfor
;  Tdisk = Tdisk_avg

  incor = where((diskbody_ijk eq 0)and(rr gt R_hor)and(rr lt R_attach))
  inatm = where((diskbody_ijk eq 1)and(rr gt R_hor))
  indisk = where((diskbody_ijk eq 2)and(rr gt R_hor))
  insim = where((diskbody_ijk le 2)and(rr gt R_hor))
  if (incor(0) eq -1) then incor = 0
  contour,T_cor(*,*,0),r#sin(t),r#cos(t),xrange=[0,60],yrange=[-30,30],$
    levels = 10.^(findgen(100)/99.*6.+5.),/isotropic,/fill
  contour,tau_es(*,*,0),r#sin(t_),r#cos(t_),xrange=[0,60],yrange=[-30,30],$
    levels = [0.01,0.1,1],/isotropic,/noerase,thick=4,c_labels=[1,1,1],$
    xtitle = 'x/M', ytitle = 'z/M',title='gas temperature'
  harm_cor = ll_cgs
  harm_cor(*,*,*)=0.
  harm_cor(incor)=ll_cgs(incor)*dv_cgs(incor)*4.
  harm_cor(incor)=llcomp_cgs(incor)*dv_cgs(incor)*4.
  harm_syn = harm_cor*0
  harm_syn(incor)=llsyn_cgs(incor)*dv_cgs(incor)*4.
  harm_tot = harm_cor*0
  harm_tot(incor)=ll_cgs(incor)*dv_cgs(incor)*4.
  rt_cor = ll_cgs
  rt_cor(*,*,*)=0.
  rt_cor(incor)=corpow_ijk(incor)*2.42d17
  rt_tot = rt_cor+harm_syn
  dldrcor = total(total(harm_cor,3),2)
  dldrsyn = total(total(harm_syn,3),2)
  dldrtot = total(total(harm_tot,3),2)
  dldrrt = total(total(rt_cor,3),2)
  dldrrttot = total(total(rt_tot,3),2)
;stop
  print,total(ll_cgs(insim)*dV_cgs(insim))*4./(Mdot_cgs*cc*cc)
  contour,Tdisk,r#cos(p),r#sin(p),xrange=[0,40],yrange=[0,40],$
    levels = 5d6*findgen(100)/99.,/isotropic,/fill
;  stop

  pdex = indgen(N3)
  if (Nsym gt 1) then begin
     diskbody2 = fltarr(N1,Nsym*N3)
     sigtau_es2 = fltarr(N1,Nsym*N3)
     Tdisk2 = fltarr(N1,Nsym*N3)
     em_top2 = fltarr(N1,Nsym*N3)
     em_bot2 = fltarr(N1,Nsym*N3)
     ref_top2 = fltarr(N1,Nsym*N3)
     ref_bot2 = fltarr(N1,Nsym*N3)
     p2 = fltarr(Nsym*N3)
     rho_cgs2 = fltarr(N1,N2,Nsym*N3)
     T_cor2 = fltarr(N1,N2,Nsym*N3)
     u0_2 = fltarr(N1,N2,Nsym*N3)
     u1_2 = fltarr(N1,N2,Nsym*N3)
     u2_2 = fltarr(N1,N2,Nsym*N3)
     u3_2 = fltarr(N1,N2,Nsym*N3)
     diskbody_ijk2 = intarr(N1,N2,Nsym*N3)
     tau_es2 = fltarr(N1,N2+1,Nsym*N3)
     b2_cgs2 = fltarr(N1,N2,Nsym*N3)
     ll_cgs2 = fltarr(N1,N2,Nsym*N3)
     for jquad =0,Nsym-1 do begin
        pdex2 = N3*jquad+pdex
        diskbody2(*,pdex2)=diskbody(*,pdex)
        sigtau_es2(*,pdex2)=sigtau_es(*,pdex)
        Tdisk2(*,pdex2)=Tdisk(*,pdex)
        em_top2(*,pdex2)=em_top(*,pdex)
        em_bot2(*,pdex2)=em_bot(*,pdex)
        ref_top2(*,pdex2)=ref_top(*,pdex)
        ref_bot2(*,pdex2)=ref_bot(*,pdex)
        p2(pdex2)=p(pdex)+jquad*2.*!PI/Nsym
        rho_cgs2(*,*,pdex2)=rho_cgs(*,*,pdex)
        T_cor2(*,*,pdex2)=T_cor(*,*,pdex)
        u0_2(*,*,pdex2)=u0(*,*,pdex)
        u1_2(*,*,pdex2)=u1(*,*,pdex)
        u2_2(*,*,pdex2)=u2(*,*,pdex)
        u3_2(*,*,pdex2)=u3(*,*,pdex)
        diskbody_ijk2(*,*,pdex2)=diskbody_ijk(*,*,pdex)
        tau_es2(*,*,pdex2)=tau_es(*,*,pdex)
        b2_cgs2(*,*,pdex2)=b2_cgs(*,*,pdex)
        ll_cgs2(*,*,pdex2)=ll_cgs(*,*,pdex)
     endfor
     diskbody = diskbody2
     sigtau_es = sigtau_es2
     Tdisk = Tdisk2
     em_top = em_top2
     em_bot = em_bot2
     ref_top = ref_top2
     ref_bot = ref_bot2
     p = p2
     rho_cgs = rho_cgs2
     b2_cgs = b2_cgs2
     ll_cgs = ll_cgs2
     tau_es = tau_es2
     T_cor = T_cor2
     u0 = u0_2
     u1 = u1_2
     u2 = u2_2
     u3 = u3_2
     diskbody_ijk = diskbody_ijk2
  endif
     
     
  run_id = run_id_out
  wdatafile = 'data/ph_0000.dat'
  dumpstr = string(run_id,format='(I4.4)')
  strput,wdatafile,dumpstr,8
  openw,1,wdatafile
  printf,1,float(diskbody),float(sigtau_es),float(Tdisk),float(em_top),$
    float(em_bot),float(ref_top),float(ref_bot)
  close,1
stop
  wdatafile = 'data/gr_0000.dat'
  dumpstr = string(run_id,format='(I4.4)')
  strput,wdatafile,dumpstr,8
  print,wdatafile
  openw,1,wdatafile
  printf,1,N1,N2,N3*Nsym
  printf,1,r,t,p
  close,1
  wdatafile = 'data/rh_0000.dat'
  dumpstr = string(run_id,format='(I4.4)')
  strput,wdatafile,dumpstr,8
  openw,1,wdatafile
  printf,1,float(rho_cgs)
  close,1
  wdatafile = 'data/te_0000.dat'
  dumpstr = string(run_id,format='(I4.4)')
  strput,wdatafile,dumpstr,8
  openw,1,wdatafile
  printf,1,float(T_cor)
  close,1
  wdatafile = 'data/u0_0000.dat'
  dumpstr = string(run_id,format='(I4.4)')
  strput,wdatafile,dumpstr,8
  openw,1,wdatafile
  printf,1,float(u0)
  close,1
  wdatafile = 'data/u1_0000.dat'
  dumpstr = string(run_id,format='(I4.4)')
  strput,wdatafile,dumpstr,8
  openw,1,wdatafile
  printf,1,float(u1)
  close,1
  wdatafile = 'data/u2_0000.dat'
  dumpstr = string(run_id,format='(I4.4)')
  strput,wdatafile,dumpstr,8
  openw,1,wdatafile
  printf,1,float(u2)
  close,1
  wdatafile = 'data/u3_0000.dat'
  dumpstr = string(run_id,format='(I4.4)')
  strput,wdatafile,dumpstr,8
  openw,1,wdatafile
  printf,1,float(u3)
  close,1
  wdatafile = 'data/ta_0000.dat'
  dumpstr = string(run_id,format='(I4.4)')
  strput,wdatafile,dumpstr,8
  openw,1,wdatafile
  printf,1,float(tau_es)
  close,1
  wdatafile = 'data/db_0000.dat'
  dumpstr = string(run_id,format='(I4.4)')
  strput,wdatafile,dumpstr,8
  openw,1,wdatafile
  printf,1,diskbody_ijk
  close,1
  wdatafile = 'data/bb_0000.dat'
  dumpstr = string(run_id,format='(I4.4)')
  strput,wdatafile,dumpstr,8
  openw,1,wdatafile
  printf,1,float(b2_cgs)
  close,1
  wdatafile = 'data/ll_0000.dat'
  dumpstr = string(run_id,format='(I4.4)')
  strput,wdatafile,dumpstr,8
  openw,1,wdatafile
  printf,1,float(ll_cgs)
  close,1
  ;oplot,rr(*,0,0),Tdisk(*,0)  
stop  
END

PRO plot_surfaces,N1,N3,sigtau_es,vertices_top,vertices_bot
  win_size = 600
  window, 0, xsize = win_size, ysize = win_size
  t3d, /reset, rot = [-80, 0, 0]
  scale3d
  bounds= 30.
  warp = bounds * 1.0
;   warp = bounds * 0.71
  
  !X.S = [-(-bounds), 1.0]/(bounds-(-bounds))
  !Y.S = [-(-warp), 1.0]/(warp-(-warp))
  !Z.S = [-(-bounds), 1.0]/(bounds-(-bounds))
  
  vertices_top(0,*,*)=vertices_top(0,*,*)-30.
  vertices_bot(0,*,*)=vertices_bot(0,*,*)-30.
  for ir =0,N1-2 do begin
      for ip=0,N3-2 do begin
          if (sigtau_es(ir,ip) gt 2.) then begin          
              if (sigtau_es(ir,ip+1) gt 2.) then begin
                  plots,[vertices_top(0,ir,ip),vertices_top(0,ir,ip+1)],$
                    [vertices_top(1,ir,ip),vertices_top(1,ir,ip+1)],$
                    [vertices_top(2,ir,ip),vertices_top(2,ir,ip+1)],/t3d,/data
                  plots,[vertices_bot(0,ir,ip),vertices_bot(0,ir,ip+1)],$
                    [vertices_bot(1,ir,ip),vertices_bot(1,ir,ip+1)],$
                    [vertices_bot(2,ir,ip),vertices_bot(2,ir,ip+1)],/t3d,/data
              endif
              if (sigtau_es(ir+1,ip) gt 2.) then begin
                  plots,[vertices_top(0,ir,ip),vertices_top(0,ir+1,ip)],$
                    [vertices_top(1,ir,ip),vertices_top(1,ir+1,ip)],$
                    [vertices_top(2,ir,ip),vertices_top(2,ir+1,ip)],/t3d,/data
                  plots,[vertices_bot(0,ir,ip),vertices_bot(0,ir+1,ip)],$
                    [vertices_bot(1,ir,ip),vertices_bot(1,ir+1,ip)],$
                    [vertices_bot(2,ir,ip),vertices_bot(2,ir+1,ip)],/t3d,/data
              endif
;              if (sigtau_es(ir+1,ip+1) gt 2.) then begin
;                  plots,[vertices_top(0,ir,ip),vertices_top(0,ir+1,ip+1)],$
;                    [vertices_top(1,ir,ip),vertices_top(1,ir+1,ip+1)],$
;                    [vertices_top(2,ir,ip),vertices_top(2,ir+1,ip+1)],/t3d,/data
;                  plots,[vertices_bot(0,ir,ip),vertices_bot(0,ir+1,ip+1)],$
;                    [vertices_bot(1,ir,ip),vertices_bot(1,ir+1,ip+1)],$
;                    [vertices_bot(2,ir,ip),vertices_bot(2,ir+1,ip+1)],/t3d,/data
;              endif
          endif
      endfor
  endfor
;  stop
;  rr_c = fltarr(N1,N2)
;  zz_c = fltarr(N1,N2)
;  rr_c(*,*)=rr(*,*,0)
;  zz_c(*,*)=rr_c*cos(tt(*,*,0))
;  stop
END

PRO color_bar,loglo,loghi
    Npx = 12600
    Nx = 100
    Ny = 5
    xx = findgen(Nx)
    yy = findgen(Ny)
    xscale = findgen(Nx)/(Nx-1.)*(loghi-loglo)+loglo
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
      levels = xscale,/fill,/Noerase
    plot_io,yscale,xscale,Position=[2600,6000,3100,12000],xticks=1,$
      xtickname=[' ',' '],/Noerase,/Device,yrange=[xscale(0),xscale(Nx-1)],$
      ystyle=1
stop
;    ,yticks=5,yminor=1,$
;      ytickname=['10!U-2!N','10!U-1!N','10!U0!N','10!U1!N',$
;                 '10!U2!N','10!U3!N'],$
;      ytitle='T (keV)',/Noerase,/Device
END
