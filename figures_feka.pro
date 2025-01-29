PRO figures_feka,ifig,run_id
!P.font = 0
!P.charsize=1.5
thk = 5.

if (ifig eq 1) then begin
    Nspec = 1001
;    run_id = 3256
    Nt = 41
    Nph = 40
    Nframes = 1
    spec = dblarr(Nspec,Nph,Nt)
    spec2 = dblarr(Nspec,Nph,Nt)
    mu = (findgen(Nt)+0.5)/Nt*2.-1.

    rdatafile = 'data/gr_0000.dat'
    dumpstr = string(1256,format='(I4.4)')
    strput,rdatafile,dumpstr,8
    openr,1,rdatafile
    readf,1,N1,N2,N3
    r = fltarr(N1)
    t = fltarr(N2)
    p = fltarr(N3)
    readf,1,r,t,p
    close,1
    dr = deriv(r)
    Rspec = dblarr(101,N1)
    Rspecp = dblarr(101,N3,N1)
    Rspecp_top = dblarr(101,N3,N1)
    Rspecp_bot = dblarr(101,N3,N1)
    Ispec = dblarr(101,Nt)
    Espec_top = dblarr(Nspec,N3,N1)
    Espec_bot = dblarr(Nspec,N3,N1)
    n64 = 0.
    n67 = 0.
    n69 = 0.
    e69 = 64
    e_64 = 634
    e_67 = 638
    e_69 = 640
    openr,1,'brooks/scat_disk.1256.dat'
    readf,1,Rspec
    readf,1,Rspecp
    readf,1,Rspecp_top
    readf,1,Rspecp_bot
    close,1
    openr,1,'data/scat_spec.1256.dat'
    readf,1,Ispec
    close,1
stop
    emin = 0.001              
    emax = 1000.
    nu = emin*10.^(findgen(Nspec)/(Nspec-1.)*alog10(emax/emin))
    dnu = deriv(nu)
    shorte = emin*10.^(findgen(101.)/(101.-1.)*alog10(emax/emin))
    de = deriv(shorte)
    blnkr = 1.0
    blnkp = 1.0
;    openw,1,'brooks/obs_specta.dat'
;    for it=0,Nt-1 do begin
;       for inu=0,100 do begin
;          printf,1,mu(it),shorte(inu),Ispec(inu,it)*2.42d17/(mu(1)-mu(0))
;       endfor
;    endfor
;    close,1

    openr,1,'feka_top-3.1256.dat'
    for ir=0,N1-1 do begin
       for ip=0,N3-1 do begin
          readf,1,blnkr,blnkp,n64,n67,n69
;          Espec_top(e69,ip,ir)=(n64*6.4+n67*6.7+n69*6.9)*1.6d-9/de(e69)
          Espec_top(e_64,ip,ir)=(n64*6.4)*1.6d-9/dnu(e_64)
          Espec_top(e_67,ip,ir)=(n67*6.7)*1.6d-9/dnu(e_67)
          Espec_top(e_69,ip,ir)=(n69*6.9)*1.6d-9/dnu(e_69)
       endfor
    endfor
    close,1
    openr,1,'feka_bot-3.1256.dat'
    for ir=0,N1-1 do begin
       for ip=0,N3-1 do begin
          readf,1,blnkr,blnkp,n64,n67,n69
;          print,r(ir),p(ip),n64,n67,n69
;          Espec_bot(e69,ip,ir)=(n64*6.4+n67*6.7+n69*6.9)*1.6d-9/de(e69)
          Espec_bot(e_64,ip,ir)=(n64*6.4)*1.6d-9/dnu(e_64)
          Espec_bot(e_67,ip,ir)=(n67*6.7)*1.6d-9/dnu(e_67)
          Espec_bot(e_69,ip,ir)=(n69*6.9)*1.6d-9/dnu(e_69)
       endfor
    endfor
    close,1
    Espec = (Espec_top+Espec_bot)/2.
    Especr = dblarr(Nspec,N1)
    Especr(*,*) = total(Espec,2)/N3
    EspecrN = dblarr(N1)
    for ir=0,N1-1 do EspecrN(ir)=total(Especr(*,ir)/(nu*1.6d-9)*dnu)
    EspecF = dblarr(Nspec)
    for ir=0,N1-1 do EspecF(*) = EspecF(*)+Especr(*,ir)*4.*!PI*r(ir)*dr(ir)*(1.5d6)^2.
    Rspecr = dblarr(101,N1)
    Rspecr(*,*)=total(Rspecp,2)/N3/2.*2.42d17
    RspecrN = dblarr(N1)
    for ir=0,N1-1 do RspecrN(ir)= $
;       total(Rspecr(e69+1:100,ir)/(shorte(e69+1:100)*1.6d-9)*de(e69+1:100))
       total(Rspecr(e69:100,ir)/(shorte(e69:100)*1.6d-9)*de(e69:100))
    RspecF = dblarr(101)
    for ir=0,N1-1 do RspecF(*) = RspecF(*)+Rspecr(*,ir)*4.*!PI*r(ir)*dr(ir)*(1.5d6)^2.
    
stop
    idex0=0
    idex30=2
    idex45=6
    idex60=10
    idex75=15

    pdex1 = findgen(Nph)
    pdex2 = (pdex1+Nph/4.)mod Nph
    pdex3 = (pdex1+Nph*2/4.)mod Nph
    pdex4 = (pdex1+Nph*3/4.)mod Nph
    inc = acos(abs(1.-(findgen(Nt)+0.5)/Nt*2.))*!radeg
    rdatafile = 'data/scat_spcp.0000.dat'
    dumpstr = string(run_id,format='(I4.4)')
    strput,rdatafile,dumpstr,15
    openr,1,rdatafile
    readf,1,spec
    close,1
    spec2(*,*,*)=0.25*(spec(*,pdex1,*)+spec(*,pdex2,*)+$
                       spec(*,pdex3,*)+spec(*,pdex4,*))
    spec_nuth = dblarr(Nspec,Nt)
    spec_nuth_sig = spec_nuth
    spec_nuth(*,*)=total(spec2,2)
    negthdex = Nt-1-findgen(Nt)
    spec_mean = spec_nuth(*,*)/Nph
;    spec_mean = (spec_nuth(*,*)+spec_nuth(*,negthdex))/(2.*Nph)
    for inu=0,Nspec-1 do begin
       for ith=0,(Nt-1)/2 do begin
          for iph=0,Nph/4-1 do begin
             spec_nuth_sig(inu,ith)=spec_nuth_sig(inu,ith)+ $
                                    (spec2(inu,iph,ith)-spec_mean(inu,ith))^2.
             spec_nuth_sig(inu,ith)=spec_nuth_sig(inu,ith)+ $
                                    (spec2(inu,iph,Nt-1-ith)-spec_mean(inu,ith))^2.
          endfor
       endfor
    endfor
;stop
    spec_nuth_sig = sqrt(spec_nuth_sig/(Nph/2.))
    ymax = max(spec_mean(*,idex0))
    plot,nu,spec_mean(*,idex0),xrange=[0,10],yrange=[0,ymax*1.2],ystyle=1,$
         xtitle='!3E!Lobs!N (keV)',ytitle='!3flux (a.u.)',thick=thk
    oplot,nu,spec_mean(*,idex0)+spec_nuth_sig(*,idex0),linestyle=2,thick=thk
    oplot,nu,spec_mean(*,idex0)-spec_nuth_sig(*,idex0),linestyle=2,thick=thk
    xyouts,1,ymax,'!3i=0!Uo!N'
;stop
    ymax = max(spec_mean(*,idex30))
    plot,nu,spec_mean(*,idex30),xrange=[0,10],yrange=[0,ymax*1.2],ystyle=1,$
         xtitle='!3E!Lobs!N (keV)',ytitle='!3flux (a.u.)',thick=thk
    oplot,nu,spec_mean(*,idex30)+spec_nuth_sig(*,idex30),linestyle=2,thick=thk
    oplot,nu,spec_mean(*,idex30)-spec_nuth_sig(*,idex30),linestyle=2,thick=thk
    xyouts,1,ymax,'!3i=30!Uo!N'
;stop
    ymax = max(spec_mean(*,idex60))
    plot,nu,spec_mean(*,idex60),xrange=[0,10],yrange=[0,ymax*1.2],ystyle=1,$
         xtitle='!3E!Lobs!N (keV)',ytitle='!3flux (a.u.)',thick=thk
    oplot,nu,spec_mean(*,idex60)+spec_nuth_sig(*,idex60),linestyle=2,thick=thk
    oplot,nu,spec_mean(*,idex60)-spec_nuth_sig(*,idex60),linestyle=2,thick=thk
    xyouts,1,ymax,'!3i=60!Uo!N'
;stop
    ymax = max(spec_mean(*,idex75))
    plot,nu,spec_mean(*,idex75),xrange=[0,10],yrange=[0,ymax*1.2],ystyle=1,$
         xtitle='!3E!Lobs!N (keV)',ytitle='!3flux (a.u.)',thick=thk
    oplot,nu,spec_mean(*,idex75)+spec_nuth_sig(*,idex75),linestyle=2,thick=thk
    oplot,nu,spec_mean(*,idex75)-spec_nuth_sig(*,idex75),linestyle=2,thick=thk
    xyouts,1,ymax,'!3i=75!Uo!N'
stop
    ymax = max(spec_mean(*,idex0))
    plot,nu,spec_mean(*,idex0),xrange=[2,10],yrange=[0,ymax],ystyle=1,$
         xtitle='!3E!Lobs!N (keV)',ytitle='!3flux (a.u.)',thick=thk
    oplot,nu,spec_mean(*,idex30),thick=thk
    oplot,nu,spec_mean(*,idex45),thick=thk
    oplot,nu,spec_mean(*,idex60),thick=thk
    oplot,nu,spec_mean(*,idex75),thick=thk
    xyouts,2.5,0.9*ymax,'PL profile'
    xyouts,2.5,0.8*ymax,'6.4keV line'

    wdata = fltarr(6,Nspec)
    wdata(0,*) = nu
    wdata(1,*) = spec_mean(*,idex0)
    wdata(2,*) = spec_mean(*,idex30)
    wdata(3,*) = spec_mean(*,idex45)
    wdata(4,*) = spec_mean(*,idex60)
    wdata(5,*) = spec_mean(*,idex75)
stop
    EW_th = dblarr(Nt)
    cont_th = dblarr(Nt)
    line_th = dblarr(Nt)
    e_lo = 2.
    e_hi = 20.
    inband_lo = where((shorte gt e_lo)and(shorte lt e_hi))
    inband_hi = where((nu gt e_lo)and(nu lt e_hi))
    costh = findgen(41)/40.*2.-1
    costh = -costh
    for ith=0,Nt-1 do begin
       pl = (nu/shorte(67))^(-1.5)*Ispec(67,ith)*2.42d17
       line = spec_mean(*,ith)*40.*2.42d17
       EW_th(ith) = total(((line+pl)/pl-1.)*dnu)*1000.
       cont_th(ith) = total(Ispec(inband_lo,ith)*de(inband_lo))*2.42d17
       line_th(ith) = total(spec_mean(inband_hi,ith)*dnu(inband_hi))*2.42d17*40.
    endfor
    plot,costh,line_th/cont_th*(e_hi-e_lo)*1000.,thick=thk,$
         xtitle='cos(i)',ytitle='EW (eV)'
    oplot,costh,EW_th,thick=thk,linestyle=2
stop
    wdata = fltarr(6,Nspec)
    wdata(0,*) = nu
    wdata(1,*) = spec_mean(*,idex0)
    wdata(2,*) = spec_mean(*,idex30)
    wdata(3,*) = spec_mean(*,idex45)
    wdata(4,*) = spec_mean(*,idex60)
    wdata(5,*) = spec_mean(*,idex75)
    openw,1,'brooks/line_ICprofiles_rawdata.dat'
    printf,1,wdata
    close,1
stop

endif


END
