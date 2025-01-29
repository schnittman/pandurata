PRO figures17b, ifig
!P.font = 0
!P.charsize=1.5
  thk = 4.

  if ((ifig eq 1)or(fix(ifig) eq 2)or(ifig eq 3)) then begin
      i_plot = 1
      N1 = 0
      N2 = 0
      N3 = 0
      Nth = 41
      Nspec = 101
      run_id = 1000
      run_sort = -1
      first_read = 0
      target_inc=60.
      rdatafile = 'data/gr_0000.dat'
      dumpstr = string(run_id,format='(I4.4)')
      strput,rdatafile,dumpstr,8
      openr,1,rdatafile
      readf,1,N1,N2,N3
      rr = fltarr(N1)
      Nr = N1
      readf,1,rr
      close,1
      drr = deriv(rr)
      Rspec = fltarr(Nspec,Nr)
      Ispec = fltarr(Nspec,Nth)
      Ispecr = fltarr(Nspec,Nth,Nr)
      Ispecr2 = fltarr(Nspec,Nth,Nr)
      Rspecr = fltarr(Nspec,Nth,Nr)
      Rspecr2 = fltarr(Nspec,Nth,Nr)
      Cspecr = fltarr(Nspec,Nr)
      Cspecr2 = fltarr(Nspec,Nr)
      spec = fltarr(Nspec)
      spec_r = fltarr(Nspec)
      spec_d = fltarr(Nspec)
      rdata = fltarr(3,Nr)
      Lfactr = fltarr(Nr)
      Gfactr = fltarr(Nr)
      I_r = fltarr(Nr)
      mean_nu_t = fltarr(Nr)
      mean_nu_d = fltarr(Nr)
      mean_nu_r = fltarr(Nr)
      Mdotr = fltarr(Nr)
      spec_s = dblarr(6,Nspec,Nth)
      spec_s2 = dblarr(6,Nspec,Nth)
      spec2 = dblarr(Nspec,Nth)
      spec0 = fltarr(Nspec,Nth)
      Qspec2 = dblarr(Nspec,Nth)
      Uspec2 = dblarr(Nspec,Nth)

      Nr2 = 912
      rdata2 = fltarr(4,Nr2)
      openr,1,'mdot_tavg.dat'
      readf,1,rdata2
      close,1
      for ir=0,Nr2-1 do begin
          if (rdata2(ir) gt rr(0)) then begin
              inrr = where(rr lt rdata2(0,ir))
              rdex = max(inrr)
              Mdotr(rdex) = rdata2(1,ir)
          endif
      endfor
      Mdot0 = mean(Mdotr(where((rr gt 5)and(rr lt 7))))
      Mdotr(where(rr gt 20)) = Mdot0
      nu = fltarr(Nspec)
      dnu = nu
      e_min = 0.001
      e_max = 1000.0
      nu_min = e_min*2.
      nu_max = e_max/2.
      nu(0) = e_min
      for j=1.,Nspec-1 do $
          nu(j) = e_min*10.0^(j/(Nspec-1.)*alog10(e_max/e_min))
      dnu = deriv(nu)
      dnurr = dnu#(fltarr(Nr)+1)
      dnuinc = dnu#(fltarr(Nth)+1)
;stop
      aa=0.0
      Nframes = 51.
      if ((first_read eq 1)and((ifig eq 1)or(ifig eq 2.1)or(ifig eq 3))) then $
        subid = [1000,1010,1020,1030,1040,1050,1060,1070,1080,1090, $
                 1100,1110,1120,1130,1140,1150,1160,1170,1180,1190, $
                 1200,1210,1220,1230,1240,1250,1260,1270,1280,1290, $
                 1300,1310,1320,1330,1340,1350,1360,1370,1380,1390, $
                 1400,1410,1420,1430,1440,1450,1460,1470,1480,1490,1500]
      if ((first_read eq 1)and(ifig eq 2.2)) then begin
          subid = [1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500]
;          subid = [900,901,902,903,904,905,906,907,908,909,910]
          subid=subid+1
          Nframes=11.
      endif
      if (first_read eq 1) then begin
          for id_dex = 0,Nframes-1 do begin
              rdatafile = 'data/scat_spec.0000.dat'
              dumpstr = string(subid(id_dex),format='(I4.4)')
              strput,rdatafile,dumpstr,15
              openr,1,rdatafile
              readf,1,spec0,Qspec2,Uspec2
              for isort=0,5 do begin
                  readf,1,spec2,Qspec2,Uspec2
                  spec_s(isort,*,*)=spec2
              endfor
              close,1
              spec_s2=spec_s2+spec_s
              rdatafile = 'data/scat_spcr.0000.dat'
              dumpstr = string(subid(id_dex),format='(I4.4)')
              strput,rdatafile,dumpstr,15
              openr,1,rdatafile
              readf,1,Ispecr,rdata,Rspecr,Rspecr,Rspecr,Cspecr
              close,1
              Ispecr2=Ispecr2+Ispecr
              Rspecr2=Rspecr2+Rspecr
              Cspecr2=Cspecr2+Cspecr
              Ispec = total(Ispecr,3)
              print,total(total(Ispecr,3)*dnuinc),total(total(Rspecr,3)*dnuinc),$
                total(total(Cspecr,2)*dnu)
              
;          stop
          endfor
          spec_s = spec_s2/Nframes
          Ispecr = Ispecr2/Nframes
          Rspecr = Rspecr2/Nframes
          Cspecr = Cspecr2/Nframes
          if ((ifig eq 1)or(ifig eq 2.1)or(ifig eq 3)) then $
            openw,1,'tint_spcr.dat'
          if (ifig eq 2.2) then $
            openw,1,'tint_spcr_thin.dat'
          printf,1,rdata,Ispecr,Rspecr,Cspecr,spec_s
          close,1
      endif
      if ((first_read eq 0)and((ifig eq 1)or(ifig eq 2.1)or(ifig eq 3))) then begin
          openr,1,'tint_spcr.dat'
          readf,1,rdata,Ispecr,Rspecr,Cspecr,spec_s
          close,1
      endif          
      if ((first_read eq 0)and(ifig eq 2.2)) then begin
          openr,1,'tint_spcr_thin.dat'
          readf,1,rdata,Ispecr,Rspecr,Cspecr,spec_s
          close,1
      endif          
;      scn_data = fltarr(2,5000)
;      openr,1,'plots/scn/Inu_i65.dat'
;      readf,1,scn_data
;      close,1
      rr =rdata(0,*)
      drr = deriv(rr)

;scale all the spectra by Mdot(r)
;      mdotr(*)=mdot0
      for ir=0,Nr-1 do begin
          for i_dex = 0,Nth-1 do begin
              spec(*)=Ispecr(*,i_dex,ir)
              scale_spec,spec,nu,dnu,nu_min,nu_max,mdot0/mdotr(ir)
              Ispecr(*,i_dex,ir)=spec
              spec(*)=Rspecr(*,i_dex,ir)
              scale_spec,spec,nu,dnu,nu_min,nu_max,mdot0/mdotr(ir)
              Rspecr(*,i_dex,ir)=spec
          endfor
          spec(*)=Cspecr(*,ir)
          scale_spec,spec,nu,dnu,nu_min,nu_max,mdot0/mdotr(ir)
          Cspecr(*,ir)=spec
      endfor

      Nstop = Nr
;      Nstop = 130
      Ispec = Ispec*0
      for ir = 0,Nstop-1 do begin
          Ispec = Ispec+Ispecr(*,*,ir)
      endfor

      for ir=0,Nr-1 do begin
          spec(*) = total(Ispecr(*,*,ir),2)
          mean_nu_t(ir) = total(spec*nu*dnu)/total(spec*dnu)
          spec_r(*) = total(Rspecr(*,*,ir),2)
          spec_d = spec-spec_r
          mean_nu_r(ir) = total(spec_r*nu*dnu)/total(spec_r*dnu)
          mean_nu_d(ir) = total(spec_d*nu*dnu)/total(spec_d*dnu)
      endfor
      Mdot_cgs = 2.556e18       ;a/m=0 HARM3D data for eta=0.06
      harm_eff = total(total(Ispec,2)*dnu)/(Mdot_cgs*3d10^2.)*2.42d17
      spec2(*,*)=spec_s(0,*,*)
      Ispec_d = spec2
      harm_eff_d =  total(total(spec2,2)*dnu)/(Mdot_cgs*3d10^2.)*2.42d17
      spec2(*,*)=spec_s(1,*,*)
      Ispec_r = spec2
      harm_eff_r =  total(total(spec2,2)*dnu)/(Mdot_cgs*3d10^2.)*2.42d17
      harm_eff_c =  total(total(Cspecr,2)*dnu)/(Mdot_cgs*3d10^2.)*2.42d17
      
      inc_dex = round(cos(target_inc/!radeg)*(Nth-1.))
      dLdR = dblarr(Nr)
      dLdR_th = dblarr(Nr,Nth)
      dLCdR = dblarr(Nr)
      for ir=0,Nr-1 do begin
          for i_dex = 0,Nth-1 do begin
              dLdR(ir)=dLdr(ir)+total(Ispecr(*,i_dex,ir)*dnu)/drr(ir)
              dLdR_th(ir,i_dex)=dLdr_th(ir,i_dex)+$
                total(Ispecr(*,i_dex,ir)*dnu)/drr(ir)
          endfor
          dLCdR(ir)=total(Cspecr(*,ir)*dnu)/drr(ir)
      endfor

      if (ifig eq 1) then begin
          plot_oo,rr,(dLdR)*2.42d17,xrange=[2,40],xstyle=1,$
;            yrange=[dLdR(Nr/2),max(dLdR)]*2.42d17,$
            yrange=[1d35,4d37],ystyle=1,$
            xtitle = 'r/M', ytitle='dL/dr  [erg/s/(GM/c!U2!N)]',thick=thk
          oplot,rr,(dLCdR+dLdR)*2.42d17,linestyle=1,thick=thk,color=120
          oplot,rr,dLdR*2.42d17,thick=thk,color=120

          dLdR_nocorona = dLdR
          openr,1,'dldr_nocorona.dat'
          readf,1,dLdR_nocorona
          close,1
          oplot,rr,dLdR_nocorona*2.42d17,thick=thk,color=60
          plots,[2.5,4],[1,1]*10^37.4,thick=thk,color=120
          plots,[2.5,4],[1,1]*10^37.2,thick=thk,linestyle=1,color=120
          plots,[2.5,4],[1,1]*10^37.0,thick=thk,linestyle=2
          !P.charsize=1.25
          xyouts,4.2,10^37.35,'ThinHR observed'
          xyouts,4.2,10^37.15,'ThinHR emitted'
          xyouts,4.2,10^36.95,'NT observed'
          xyouts,15,10^37.35,'optically thick'
          !P.charsize=1.5
      endif
stop
      if (ifig eq 2.1) then begin
          inc_dex0 = round(cos(0./!radeg)*(Nth-1.))
          inc_dex = round(cos(30./!radeg)*(Nth-1.))
          plot_oo,nu,nu*Ispec(*,inc_dex)*2.42d17*Nth,xrange=[0.2,10],xstyle=1,$
            yrange=[0.3,10]*2d37,ystyle=1,$
            xtitle='E (keV)', ytitle='!9n !3L!9!Ln!N !3',thick=thk
          inc_dex = round(cos(30./!radeg)*(Nth-1.))
          oplot,nu,nu*Ispec(*,inc_dex)*2.42d17*Nth,thick=thk,color=60
          inc_dex = round(cos(45./!radeg)*(Nth-1.))
          oplot,nu,nu*Ispec(*,inc_dex)*2.42d17*Nth,thick=thk,color=90
          inc_dex = round(cos(60./!radeg)*(Nth-1.))
          oplot,nu,nu*Ispec(*,inc_dex)*2.42d17*Nth,thick=thk,color=120
          inc_dex = round(cos(75./!radeg)*(Nth-1.))
          oplot,nu,nu*Ispec(*,inc_dex)*2.42d17*Nth,thick=thk,color=150
          !P.charsize=1.25
          xyouts,3,10^38.1,'optically thick'
          !P.charsize=1.5
      endif
      if (ifig eq 2.2) then begin
          inc_dex0 = round(cos(0./!radeg)*(Nth-1.))
          inc_dex = round(cos(30./!radeg)*(Nth-1.))
          plot_oo,nu,nu*Ispec(*,inc_dex)*2.42d17*Nth,xrange=[0.2,10],xstyle=1,$
            yrange=[0.3,10]*2d37,ystyle=1,$
            xtitle='E (keV)', ytitle='!9n !3F!9!Ln!N !3',thick=thk
          inc_dex = round(cos(30./!radeg)*(Nth-1.))
          oplot,nu,nu*Ispec(*,inc_dex)*2.42d17*Nth,thick=thk,color=60
          inc_dex = round(cos(45./!radeg)*(Nth-1.))
          oplot,nu,nu*Ispec(*,inc_dex)*2.42d17*Nth,thick=thk,color=90
          inc_dex = round(cos(60./!radeg)*(Nth-1.))
          oplot,nu,nu*Ispec(*,inc_dex)*2.42d17*Nth,thick=thk,color=120
          inc_dex = round(cos(75./!radeg)*(Nth-1.))
          oplot,nu,nu*Ispec(*,inc_dex)*2.42d17*Nth,thick=thk,color=150
          plots,[0.25,0.5],[1,1]*10.^38.15,thick=thk,color=60
          plots,[0.25,0.5],[1,1]*10.^38.05,thick=thk,color=90
          plots,[0.25,0.5],[1,1]*10.^37.95,thick=thk,color=120
          plots,[0.25,0.5],[1,1]*10.^37.85,thick=thk,color=150
          !P.charsize=1.25
          xyouts,0.52,10^38.13,'i = 30!Uo!N'
          xyouts,0.52,10^38.03,'i = 45!Uo!N'
          xyouts,0.52,10^37.93,'i = 60!Uo!N'
          xyouts,0.52,10^37.83,'i = 75!Uo!N'
          xyouts,3,10^38.1,'optically thin'
          !P.charsize=1.5
      endif
      if (ifig eq 3) then begin
          plot_oo,nu,nu*Ispec(*,inc_dex)*2.42d17*Nth,xrange=[0.2,10],xstyle=1,$
            yrange=[0.3,10]*2d37,ystyle=1,$
            xtitle='E (keV)', ytitle='!9n !3L!9!Ln!N !3',thick=thk
          oplot,nu,nu*Ispec(*,inc_dex)*2.42d17*Nth,color=120,thick=thk
          plots,[0.25,0.5],[1,1]*10.^38.15,thick=thk,color=120
          plots,[0.25,0.5],[1,1]*10.^38.0,thick=thk,linestyle=2 
          !P.charsize=1.25
          xyouts,0.52,10^38.13,'ThinHR'
          xyouts,0.52,10^37.98,'NT'
          xyouts,3,10^38.1,'optically thick'
          !P.charsize=1.5
;          oplot,nu,nu*Ispec_d(*,inc_dex)*2.42d17,linestyle=2,thick=thk
;          oplot,nu,nu*Ispec_r(*,inc_dex)*2.42d17,linestyle=1,thick=thk
      endif
      Naa = 5
      aa = fltarr(Naa)
      aa_all = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99,0.995,0.999]
      aaid = [900,901,902,903,904,905,906,907,908,909,910]
      aa(*) = aa_all(0:Naa-1)
      Ispecr_fit = Ispecr
      rdatafile = 'data/scat_spcr.0000.dat'
      for iaa = 0,Naa-1,2 do begin
          dumpstr = string(aaid(iaa),format='(I4.4)')
          strput,rdatafile,dumpstr,15
          openr,1,rdatafile
          readf,1,Ispecr_fit,rdata,Rspecr,Rspecr,Rspecr,Cspecr
          close,1
          rr = rdata(0,*)
          Ispecr_fit(*,*,Nr-1)=0.
          Ispec_fit = total(Ispecr_fit,3)
          NT_eff = total(total(Ispec_fit,2)*dnu)*2.42d17/(2.55d18*3d10^2.)
          spec2(*,*)=spec_s(0,*,*)
          Ispec_fit_d = spec2
          NT_eff_d =  total(total(spec2,2)*dnu)*2.42d17/(2.55d18*3d10^2.)
          spec2(*,*)=spec_s(1,*,*)
          Ispec_fit_r = spec2
          NT_eff_r =  total(total(spec2,2)*dnu)*2.42d17/(2.55d18*3d10^2.)
          dLdR = dblarr(Nr)
          dLdR_th = dblarr(Nr,Nth)
          dLCdR = dblarr(Nr)
          if ((ifig eq 1) and (iaa ge 0)) then begin
              for ir=0,Nr-1 do begin
                  for i_dex = 0,Nth-1 do begin
                      dLdR(ir)=dLdr(ir)+total(Ispecr_fit(*,i_dex,ir)*dnu)/drr(ir)
                      dLdR_th(ir,i_dex)=dLdr_th(ir,i_dex)+$
                        total(Ispecr_fit(*,i_dex,ir)*dnu)/drr(ir)
                  endfor
                  dLCdR(ir)=total(Cspecr(*,ir)*dnu)/drr(ir)
              endfor
;              oplot,rr,dLCdR*2.42d17,linestyle=1,thick=thk,color=120
              oplot,rr,dLdR*2.42d17,linestyle=2,thick=thk
;              oplot,rr,dLdR_th(*,5)*2.42d17*20,linestyle=1,thick=thk
;              oplot,rr,dLdR_th(*,20)*2.42d17*20,linestyle=3,thick=thk
          endif
          if (ifig eq 3) then $
            oplot,nu,Ispec_fit(*,inc_dex)*nu*2.418d17*Nth,linestyle=2,thick=thk
      endfor
  endif

  if (ifig eq 4) then begin
      i_plot = 1
      N1 = 0
      N2 = 0
      N3 = 0
      Nr = 360
      Nth = 41
      Nspec = 101
      run_id = 1000
      L_eta_fact = Nth/(3d10^2.*2.55d18)
      rdatafile = 'data/gr_0000.dat'
      dumpstr = string(run_id,format='(I4.4)')
      strput,rdatafile,dumpstr,8
      openr,1,rdatafile
      readf,1,N1,N2,N3
      rr = fltarr(N1)
      Nr = N1
      readf,1,rr
      close,1
      drr = deriv(rr)
      Ispecr = fltarr(Nspec,Nth,Nr)
      Ispecr_nt = fltarr(Nspec,Nth,Nr)
      Rspecr = fltarr(Nspec,Nth,Nr)
      Cspecr = fltarr(Nspec,Nr)
      Ispec = fltarr(Nspec,Nth)
      sub_spec = fltarr(Nspec)
      spec_s = dblarr(6,Nspec,Nth)
      Mdotr = fltarr(Nr)

      Nr2 = 912
      rdata2 = fltarr(4,Nr2)
      openr,1,'mdot_tavg.dat'
      readf,1,rdata2
      close,1
      for ir=0,Nr2-1 do begin
          if (rdata2(ir) gt rr(0)) then begin
              inrr = where(rr lt rdata2(0,ir))
              rdex = max(inrr)
              Mdotr(rdex) = rdata2(1,ir)
          endif
      endfor
      Mdot0 = mean(Mdotr(where((rr gt 5)and(rr lt 7))))
      Mdotr(where(rr gt 20)) = Mdot0
;stop
      rdata = fltarr(3,Nr)
      openr,1,'tint_spcr_thin.dat'
      readf,1,rdata,Ispecr,Rspecr,Cspecr,spec_s
      close,1
      rr =rdata(0,*)
      drr = deriv(rr)
      nu = fltarr(Nspec)
      dnu = nu
      e_min = 0.001
      e_max = 1000.0
      nu_min = e_min*2.
      nu_max = e_max/2.
      nu(0) = e_min
      for j=1.,Nspec-1 do $
          nu(j) = e_min*10.0^(j/(Nspec-1.)*alog10(e_max/e_min))
      dnu = deriv(nu)
      dnurr = dnu#(fltarr(Nr)+1)
      dnuinc = dnu#(fltarr(Nth)+1)

;scale all the spectra by Mdot(r)
      for ir=0,Nr-1 do begin
          for i_dex = 0,Nth-1 do begin
              sub_spec(*)=Ispecr(*,i_dex,ir)
              scale_spec,sub_spec,nu,dnu,nu_min,nu_max,mdot0/mdotr(ir)
              Ispecr(*,i_dex,ir)=sub_spec
          endfor
      endfor

      Nstop = Nr
;      Nstop = 130
      Ispec = Ispec*0
      for ir = 0,Nstop-1 do begin
          Ispec = Ispec+Ispecr(*,*,ir)
      endfor
      spec = total(Ispec*dnuinc,1)*2.42d17
      th = acos((findgen(Nth)+0.5)/Nth)*!radeg
      th_new = [th(3:Nth-1),0]
      spec_new = [spec(3:Nth-1),spec(Nth-1)+0.5*(spec(Nth-1)-spec(Nth-2))]
      plot,th_new,spec_new*L_eta_fact,xrange=[0,90],xstyle=1,xtitle='observer inc (deg)',$
        ytitle='d!9h!3/d(cos!9q!3)',thick=thk,yrange=[0,0.2]
      oplot,th_new,spec_new*L_eta_fact,color=120,thick=thk
      plots,[5,15],[1,1]*0.18,thick=thk,color=120
      plots,[5,15],[1,1]*0.16,thick=thk
      !P.charsize=1.25
      xyouts,16,0.178,'ThinHR'
      xyouts,16,0.158,'NT'
      plots,[50,60],0.18,thick=thk
      plots,[50,60],0.16,thick=thk,linestyle=2
      xyouts,61,0.178,'Optically thin'
      xyouts,61,0.158,'Optically thick'
      !P.charsize=1.5
      
      rdatafile = 'data/scat_spcr.0899.dat'
      openr,1,rdatafile
      readf,1,Ispecr_nt,rdata
      close,1
      rr = rdata(0,*)
      Ispecr_nt(*,*,Nr-1)=0.
      Ispec_nt = total(Ispecr_nt,3)
      spec_nt = total(Ispec_nt*dnuinc,1)*2.42d17
      spec_nt_new = [spec_nt(3:Nth-1),spec_nt(Nth-1)+ $
                     0.5*(spec_nt(Nth-1)-spec_nt(Nth-2))]
      oplot,th_new,spec_nt_new*L_eta_fact,thick=thk


      openr,1,'tint_spcr.dat'
      readf,1,rdata,Ispecr,Rspecr,Cspecr,spec_s
      close,1
;scale all the spectra by Mdot(r)
      for ir=0,Nr-1 do begin
          for i_dex = 0,Nth-1 do begin
              sub_spec(*)=Ispecr(*,i_dex,ir)
              scale_spec,sub_spec,nu,dnu,nu_min,nu_max,mdot0/mdotr(ir)
              Ispecr(*,i_dex,ir)=sub_spec
          endfor
      endfor

      Ispec = Ispec*0
      for ir = 0,Nstop-1 do begin
          Ispec = Ispec+Ispecr(*,*,ir)
      endfor
      spec = total(Ispec*dnuinc,1)*2.42d17
      spec_new = [spec(3:Nth-1),spec(Nth-1)+0.5*(spec(Nth-1)-spec(Nth-2))]
      oplot,th_new,spec_new*L_eta_fact,thick=thk,linestyle=2,color=120

      rdatafile = 'data/scat_spcr.0900.dat'
      openr,1,rdatafile
      readf,1,Ispecr_nt,rdata
      close,1
      rr = rdata(0,*)
      Ispecr_nt(*,*,Nr-1)=0.
      Ispec_nt = total(Ispecr_nt,3)
      spec_nt = total(Ispec_nt*dnuinc,1)*2.428d17
      spec_nt_new = [spec_nt(3:Nth-1),spec_nt(Nth-1)+ $
                     0.5*(spec_nt(Nth-1)-spec_nt(Nth-2))]
      oplot,th_new,spec_nt_new*L_eta_fact,thick=thk,linestyle=2

      stop
  endif

  if (fix(ifig) eq 5) then begin
      Naa = 20
      Nth = 41
      aa = fltarr(Naa)
;      aa_all = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
      aa_all = [0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,$
                0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95]
      aa(*) = aa_all(0:Naa-1)
      th = acos((findgen(Nth)+0.5)/(Nth))*!radeg
      newth = [th,0]
      chi2grid = dblarr(Naa,Nth)+1.
      aagrid = aa#(dblarr(Nth+1)+1.)
      thgrid = (dblarr(Naa)+1.)#newth
      newchi2 = dblarr(Naa,Nth+1)
      if ifig eq 5.1 then openr,1,'chi2grid_i15.dat'
      if ifig eq 5.2 then openr,1,'chi2grid_i30.dat'
;      if ifig eq 5.3 then openr,1,'chi2grid_i45.dat'
      if ifig eq 5.3 then openr,1,'chi2grid_i45_new.dat'
      if ifig eq 5.4 then openr,1,'chi2grid_i60.dat'
;      if ifig eq 5.5 then openr,1,'chi2grid_i75.dat'
      if ifig eq 5.5 then openr,1,'chi2grid_i75_new.dat'
      readf,1,chi2grid
      close,1
      newchi2(*,0:Nth-1)=chi2grid
      newchi2(*,Nth)=newchi2(*,Nth-1)+0.5*(newchi2(*,Nth-1)-newchi2(*,Nth-2))
      newchi2 = abs(newchi2)
      newchi2 = newchi2/min(newchi2)
      contour,newchi2,thgrid,aagrid, $
        levels=10.^(findgen(101)/100.*4),$
        xrange=[10,85],xstyle=1,yrange=[0,0.9],ystyle=1,$
        xtitle='inclination (deg)',ytitle='spin (a/M)',/fill
  endif
  stop
END


PRO scale_spec,Ispec_fit,nu,dnu,nu_min,nu_max,best_mdot
  inrange = where((nu ge nu_min)and(nu le nu_max))
  inrange = [min(inrange)-1,inrange,max(inrange)+1]
  Ispec_fit0 = Ispec_fit
  flux0 = total(Ispec_fit0(inrange)*dnu(inrange))
  best_dist = 1.
  best_mass = 1.
  dist_scale = best_dist
  mass_scale = best_mass
  mdot_scale = best_mdot
  dE = (mass_scale^(-0.5))*(mdot_scale^0.25)
  logdE = alog(dE)/alog(nu(1)/nu(0))
  ilo = floor(logdE)
  ihi = ilo+1
  wlo = ihi-logdE
  whi = 1.-wlo
  Ispec_fit0 = Ispec_fit*dist_scale^(-2.)*mdot_scale/dE
  Ispec_fit(inrange)=(Ispec_fit0(inrange-ilo)*wlo + $
    Ispec_fit0(inrange-ihi)*whi)
  return
END
