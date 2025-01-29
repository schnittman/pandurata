PRO figures5b,ifig,line_yes
  set_plot,'ps'
  !P.Charsize=1.5
  !P.font= 0
  thk = 4.0
  
  if ((ifig eq 1)or(ifig eq 4)) then begin
      !P.Charsize=1
      N = 81.
      Nt = 21.
      Nspec = 201

      if (ifig eq 1) then run_id = 153
      if (ifig eq 4) then run_id = 145
      run_sort = 0
      rdata = dblarr(2,N,N)
      Ixy = dblarr(N,N)
      mov = dblarr(N,N,Nt)
      movx = dblarr(N,N,Nt)
      movy = dblarr(N,N,Nt)
      I_tot = dblarr(Nt)

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
      mov(N-1,*,*)=0
      movx(N-1,*,*)=0
      movy(N-1,*,*)=0
      mov(0,*,*)=0
      movx(0,*,*)=0
      movy(0,*,*)=0
      mov(*,N-1,*)=0
      movx(*,N-1,*)=0
      movy(*,N-1,*)=0
      sortmov = mov(sort(mov))
      movmax = max(mov)
      movmax=sortmov(1.*N*N*Nt-100.)
      outliers = where(mov gt movmax)
      movx(outliers)=movx(outliers)/mov(outliers)*movmax
      movy(outliers)=movy(outliers)/mov(outliers)*movmax
      mov(outliers)=movmax

      if (Nt eq 11) then $
        inc_lbl = ['i=87!Uo!N','i=82!Uo!N','i=77!Uo!N','i=71!Uo!N','i=66!Uo!N',$
                   'i=60!Uo!N','i=54!Uo!N','i=47!Uo!N','i=39!Uo!N','i=30!Uo!N','i=17!Uo!N']
      if (Nt eq 21) then $
        inc_lbl = ['i=89!Uo!N','i=86!Uo!N','i=83!Uo!N','i=80!Uo!N','i=78!Uo!N',$
                   'i=75!Uo!N','i=72!Uo!N','i=69!Uo!N','i=66!Uo!N','i=63!Uo!N',$
                   'i=60!Uo!N','i=57!Uo!N','i=53!Uo!N','i=50!Uo!N','i=46!Uo!N',$
                   'i=42!Uo!N','i=38!Uo!N','i=33!Uo!N','i=28!Uo!N','i=22!Uo!N',$
                   'i=12!Uo!N']
      it = 5
      if (N eq 81) then Ixy = dblarr(121,121)
      ishift = N/2
      Xpol = Ixy
      Ypol = Ixy
      Ixy(ishift/2:ishift/2+N-1,ishift:N+ishift-1) = mov(*,*,it)
      Xpol(ishift/2:ishift/2+N-1,ishift:N+ishift-1) = movx(*,*,it)
      Ypol(ishift/2:ishift/2+N-1,ishift:N+ishift-1) = movy(*,*,it)
      Ixy = Ixy+1d-10
      psi = atan(Ypol/Ixy,Xpol/Ixy)/2.
      tot_ang = atan(total(Ypol),total(Xpol))/2.*!radeg
      deg = sqrt((Ypol/Ixy)^2.+(Xpol/Ixy)^2.)
      tot_deg = sqrt(total(Xpol)^2+total(Ypol)^2)/total(Ixy)
      Xpol = deg*cos(psi)
      Ypol = deg*sin(psi)
      I_tot(it)=total(Ixy)
      data = byte(255*(alog10(Ixy/movmax+1e-4)+4.01)/4.)
      N15 = (N-1)*1.5+1
      N2 = (fix(600./N15))*N15
      enlarge = rebin(data,N2,N2)
      enlarge(*,0:N2/4)=255
      erase
      tvscl,enlarge,0,0,/Device
      dd = 100.
      pstep=4
      for i=N15/6.,N15*(5./6.),pstep do begin
          for j=N15/3.,N15,pstep do begin
              x0 = i/(N15-1.)*12600.
              y0 = j/(N15-1.)*12600.
              dff = 5d3*dd*([Xpol(i,j),Ypol(i,j),0])*Ixy(i,j)/movmax
              dff = 100.*dd*([Xpol(i,j),Ypol(i,j),0]);*data(i,j)
              if (Ixy(i,j) gt 3d-5) then $
                plots,[x0-dff(0)/2.,x0+dff(0)/2.],[y0-dff(1)/2.,y0+dff(1)/2.],$
                color=0,/Device,thick=thk
          endfor
      endfor
      plots,[10900,10900+5.*dd],[12000,12000],color=255,thick=thk,/Device
      xyouts,10900,11500,'deg=5%',color=255,/Device
      plots,[7333,10500],[3800,3800],color=255,thick=thk,/Device
      plots,[2200,5367],[3800,3800],color=255,thick=thk,/Device
      plots,[2200,2200],[3700,3900],color=255,thick=thk,/Device
      plots,[10500,10500],[3700,3900],color=255,thick=thk,/Device
      xyouts,5600,3700,'D = 40M',color=255,/Device
;stop
      Npx = 12600
      Nx = 10
      Ny = 5
      xx = findgen(Nx)
      yy = findgen(Ny)
      xscale = findgen(Nx)/(Nx-1.)*5.-5.
      xscale = 10.^xscale
      yscale = fltarr(Nx)+1
      movscl = fltarr(Ny,Nx)
      for i=0,Nx-1 do movscl(*,i)=xscale(i)
      data = byte(255*(alog10(movscl/max(movscl)+1e-5)+5.01)/5.)
      enlarge = rebin(data,50,1000)
;      tvscl,enlarge,0,1000,/Device
      contour,enlarge,Position=[1600,6000,2100,12000],/Device,$
        yticks=1,yminor=1,ytickname=[' ',' '],$
        xticks=1,xminor=1,xtickname=[' ',' '],$
        levels=indgen(255),c_color=indgen(255),/fill,/Noerase
      plot_oi,yscale,xscale,Position=[1600,6000,2100,12000],xticks=1,$
        xtickname=[' ',' '],yticks=5,yminor=1,$
        ytickname=['10!U-5!N','10!U-4!N','10!U-3!N','10!U-2!N','10!U-1!N','1'],$
        ytitle='I/I!Lmax!N',/Noerase,/Device,color=255
      xyouts,11000,10500,inc_lbl(it),/Device,color=255
  endif

  if ((ifig eq 2) or (ifig eq 3)) then begin
      Nt = 21
      Nspec = 201
      
      if (Nt eq 21) then begin
          idex75 = 5
          idex60 = 10
          idex45 = 14
          
          inc_lbl = ['i=89!Uo!N','i=86!Uo!N','i=83!Uo!N','i=80!Uo!N','i=78!Uo!N',$
                     'i=75!Uo!N','i=72!Uo!N','i=69!Uo!N','i=66!Uo!N','i=63!Uo!N',$
                     'i=60!Uo!N','i=57!Uo!N','i=53!Uo!N','i=50!Uo!N','i=46!Uo!N',$
                     'i=42!Uo!N','i=38!Uo!N','i=33!Uo!N','i=28!Uo!N','i=22!Uo!N',$
                     'i=12!Uo!N']
      endif

      spec = dblarr(Nspec,Nt)
      Qspec = dblarr(Nspec,Nt)
      Uspec = dblarr(Nspec,Nt)
      spec2 = dblarr(Nspec,Nt)
      Qspec2 = dblarr(Nspec,Nt)
      Uspec2 = dblarr(Nspec,Nt)
      spec_s = dblarr(6,Nspec,Nt)
      Qspec_s = dblarr(6,Nspec,Nt)
      Uspec_s = dblarr(6,Nspec,Nt)
;LOG ENERGY SPACING
      emin = 0.1
      emax = 1000
      shorte = emin*10.^(findgen(Nspec)/(Nspec-1.)*alog10(emax/emin))
      de = deriv(shorte)
      
      run_id = 143
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
      spec(*,*)=spec_s(0,*,*)
      Qspec(*,*)=Qspec_s(0,*,*)
      Uspec(*,*)=Uspec_s(0,*,*)
      deg_spec = sqrt(Qspec^2+Uspec^2)/spec
      ang_spec = atan(Uspec,Qspec)/2.

      if (ifig eq 2) then begin
          plot_oo,shorte,deg_spec(*,idex45)*100.,$
            xrange=[emin,emax/10.],yrange=[1d-1,10],xstyle=1,ystyle=1,$
            xtitle='E!Lobs!N (keV)',ytitle='polarization degree (%)',thick=thk
          oplot,shorte,deg_spec(*,idex45)*100,thick=thk,color=60
          oplot,shorte,deg_spec(*,idex60)*100,thick=thk,color=120
          oplot,shorte,deg_spec(*,idex75)*100,thick=thk
          xyouts,0.2,deg_spec(0,idex45)*100,'i=45!Uo!N',color=60
          xyouts,0.2,deg_spec(0,idex60)*100,'i=60!Uo!N',color=120
          xyouts,0.2,deg_spec(0,idex75)*100,'i=75!Uo!N'
          plots,[10,25],[1,1]*(7),thick=thk
          plots,[10,25],[1,1]*(5),thick=thk,linestyle=2
          xyouts,27,7,'a/m = 0'
          xyouts,27,5,'a/m = 0.9'
      endif

      if (ifig eq 3) then begin
          ang_spec = ang_spec*!radeg
          plot_oi,shorte,smooth(ang_spec(*,idex45),1),$
            xrange=[emin,emax/10.],yrange=[-90,10],xstyle=1,ystyle=1,$
            xtitle='E!Lobs!N (keV)',ytitle='polarization angle (deg)'
          oplot,shorte,smooth(ang_spec(*,idex45),1),thick=thk,color=60
          oplot,shorte,smooth(ang_spec(*,idex60),1),thick=thk,color=120
          oplot,shorte,smooth(ang_spec(*,idex75),1),thick=thk
          plots,[0.15,0.3],[1,1]*(-20.),thick=thk
          plots,[0.15,0.3],[1,1]*(-25.),thick=thk,linestyle=2
          xyouts,0.32,-20,'a/m = 0'
          xyouts,0.32,-25,'a/m = 0.9'
      endif

      run_id = 145
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
      spec(*,*)=spec_s(0,*,*)
      Qspec(*,*)=Qspec_s(0,*,*)
      Uspec(*,*)=Uspec_s(0,*,*)
      deg_spec = sqrt(Qspec^2+Uspec^2)/spec
      ang_spec = atan(Uspec,Qspec)/2.

      if (ifig eq 2) then begin
          oplot,shorte,deg_spec(*,idex45)*100,thick=thk,color=60,linestyle=2
          oplot,shorte,deg_spec(*,idex60)*100,thick=thk,color=120,linestyle=2
          oplot,shorte,deg_spec(*,idex75)*100,thick=thk,linestyle=2
      endif
      if (ifig eq 3) then begin
          ang_spec = ang_spec*!radeg
          oplot,shorte,ang_spec(*,idex45),thick=thk,color=60,linestyle=2
          oplot,shorte,ang_spec(*,idex60),thick=thk,color=120,linestyle=2
          oplot,shorte,ang_spec(*,idex75),thick=thk,linestyle=2
;          xyouts,25,ang_spec(Nspec/1.8,idex45),'i=45!Uo!N',color=60
;          xyouts,25,ang_spec(Nspec/1.8,idex60),'i=60!Uo!N',color=120
;          xyouts,25,ang_spec(Nspec/1.8,idex75),'i=75!Uo!N'
          xyouts,25,-58,'i=45!Uo!N',color=60
          xyouts,25,-35,'i=60!Uo!N',color=120
          xyouts,25,-2,'i=75!Uo!N'
      endif
  endif

  if ((fix(ifig) eq 5)or(fix(ifig) eq 6)or(fix(ifig) eq 16)) then begin

      !P.Charsize=2
      if (fix(ifig) eq 5) then run_id=143
      if (fix(ifig) eq 6) then run_id=145
      if (fix(ifig) eq 16) then run_id=147
      Nt = 21
      Nspec = 201
      idex75 = 5
      idex60 = 10
      idex45 = 14
      
      spec = dblarr(Nspec,Nt)
      Qspec = dblarr(Nspec,Nt)
      Uspec = dblarr(Nspec,Nt)
      spec_s = dblarr(6,Nspec,Nt)
      Qspec_s = dblarr(6,Nspec,Nt)
      Uspec_s = dblarr(6,Nspec,Nt)
      spec2 = dblarr(Nspec,Nt)
      Qspec2 = dblarr(Nspec,Nt)
      Uspec2 = dblarr(Nspec,Nt)
      spec_dir = dblarr(Nspec,Nt)
      Qspec_dir = dblarr(Nspec,Nt)
      Uspec_dir = dblarr(Nspec,Nt)
      spec_ref = dblarr(Nspec,Nt)
      Qspec_ref = dblarr(Nspec,Nt)
      Uspec_ref = dblarr(Nspec,Nt)
      spec_scat = dblarr(Nspec,Nt)
      Qspec_scat = dblarr(Nspec,Nt)
      Uspec_scat = dblarr(Nspec,Nt)
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
      
      spec_dir(*,*)=spec_s(0,*,*)
      Qspec_dir(*,*)=Qspec_s(0,*,*)
      Uspec_dir(*,*)=Uspec_s(0,*,*)
      spec_ref(*,*)=spec_s(1,*,*)
      Qspec_ref(*,*)=Qspec_s(1,*,*)
      Uspec_ref(*,*)=Uspec_s(1,*,*)
      spec_scat(*,*)=total(spec_s(2:5,*,*),1)
      Qspec_scat(*,*)=total(Qspec_s(2:5,*,*),1)
      Uspec_scat(*,*)=total(Uspec_s(2:5,*,*),1)
      
;LOG ENERGY SPACING
      emin = 0.1
      emax = 1000
      shorte = emin*10.^(findgen(Nspec)/(Nspec-1.)*alog10(emax/emin))
      de = deriv(shorte)

      deg_spec = sqrt(Qspec^2+Uspec^2)/spec
      ang_spec = atan(Uspec,Qspec)/2.*!radeg
      deg_spec_dir = sqrt(Qspec_dir^2+Uspec_dir^2)/spec_dir
      ang_spec_dir = atan(Uspec_dir,Qspec_dir)/2.*!radeg-1
      deg_spec_ref = sqrt(Qspec_ref^2+Uspec_ref^2)/spec_ref
      ang_spec_ref = atan(Uspec_ref,Qspec_ref)/2.*!radeg
      deg_spec_scat = sqrt(Qspec_scat^2+Uspec_scat^2)/spec_scat
      ang_spec_scat = atan(Uspec_scat,Qspec_scat)/2.*!radeg

      if (fix(ifig) eq 5) then min_deg = -80
      if (fix(ifig) eq 6) then min_deg = 10
      if (fix(ifig) eq 16) then min_deg = 10
      neg_ang = where(ang_spec lt min_deg)
      if (neg_ang(0) ge 0) then ang_spec(neg_ang)=ang_spec(neg_ang)+180.
      neg_ang = where(ang_spec_dir lt min_deg)
      if (neg_ang(0) ge 0) then ang_spec_dir(neg_ang)=ang_spec_dir(neg_ang)+180.
      neg_ang = where(ang_spec_ref lt min_deg)
      if (neg_ang(0) ge 0) then ang_spec_ref(neg_ang)=ang_spec_ref(neg_ang)+180.

      if ((ifig ge 5.11)and(ifig le 5.13)) then it = 14
      if ((ifig ge 5.21)and(ifig le 5.23)) then it = 10
      if ((ifig ge 5.31)and(ifig le 5.33)) then it = 5
      if ((ifig ge 6.11)and(ifig le 6.13)) then it = 14
      if ((ifig ge 6.21)and(ifig le 6.23)) then it = 10
      if ((ifig ge 6.31)and(ifig le 6.33)) then it = 5
      if ((ifig ge 16.11)and(ifig le 16.13)) then it = 14
      if ((ifig ge 16.21)and(ifig le 16.23)) then it = 10
      if ((ifig ge 16.31)and(ifig le 16.33)) then it = 5

      spec_ymax = max(spec)
      spec_ymax = 300.
      spec_ymin = 1d-3*spec_ymax
      spec_ymin = 0.01
      if ((fix(100*ifig+0.5) mod 10) eq 1) then begin
          plot_oo,shorte,smooth(spec(*,it),1),$
            xrange=[0.1,20.],xstyle=1,yrange=[spec_ymin,spec_ymax],$
            ystyle=1,xtitle='E!Lobs!N (keV)',ytitle='Intensity',thick=thk
          oplot,shorte,spec_dir(*,it),linestyle=1,thick=thk
          oplot,shorte,spec_ref(*,it),linestyle=2,thick=thk
;          oplot,shorte,spec_scat(*,it),linestyle=3,thick=thk
;          plots,[7,20],[1,1]*(70.),thick=thk,linestyle=1
;          plots,[7,20],[1,1]*(45.),thick=thk,linestyle=2
;          plots,[7,20],[1,1]*(30.),thick=thk
;          xyouts,22,70,'direct'
;          xyouts,22,45,'reflected'
;          xyouts,22,30,'total'
          if (fix(ifig) eq 5) then xyouts,2,100,'a/M = 0'
          if (fix(ifig) eq 6) then xyouts,2,100,'a/M = 0.9'
          if (fix(ifig) eq 16) then xyouts,2,100,'a/M = 0.998'
          if (it eq 5) then xyouts,2,30,'i = 75!Uo!N'
          if (it eq 10) then xyouts,2,30,'i = 60!Uo!N'
          if (it eq 14) then xyouts,2,30,'i = 45!Uo!N'
      endif
      if ((fix(100*ifig+0.5) mod 10) eq 2) then begin
          plot_oo,shorte,smooth(deg_spec(*,it),1)*100.,$
            xrange=[emin,20.],yrange=[1d-1,100],xstyle=1,ystyle=1,$
            xtitle='E!Lobs!N (keV)',ytitle='polarization degree (%)',thick=thk
          oplot,shorte,deg_spec_dir(*,it)*100.,linestyle=1,thick=thk
          oplot,shorte,deg_spec_ref(*,it)*100.,linestyle=2,thick=thk
;          oplot,shorte,deg_spec_scat(*,it)*100.,linestyle=3,thick=thk
;          plots,[0.15,0.3],[1,1]*(70.),thick=thk,linestyle=1
;          plots,[0.15,0.3],[1,1]*(45.),thick=thk,linestyle=2
;          plots,[0.15,0.3],[1,1]*(30.),thick=thk
;          xyouts,0.32,70,'direct'
;          xyouts,0.32,45,'reflected'
;          xyouts,0.32,30,'total'
          if (fix(ifig) eq 5) then xyouts,0.15,50,'a/M = 0'
          if (fix(ifig) eq 6) then xyouts,0.15,50,'a/M = 0.9'
          if (fix(ifig) eq 16) then xyouts,0.15,50,'a/M = 0.998'
          if (it eq 5) then xyouts,0.15,20,'i = 75!Uo!N'
          if (it eq 10) then xyouts,0.15,20,'i = 60!Uo!N'
          if (it eq 14) then xyouts,0.15,20,'i = 45!Uo!N'
      endif
      if ((fix(100*ifig+0.5) mod 10) eq 3) then begin
          plot_oi,shorte,smooth(ang_spec(*,it),3),$
            ;xticks=3,$      ;xtickname=['0','0.5','1.0','1.5','2.0'],$
            xrange=[0.1,20],yrange=[min_deg,min_deg+180],xstyle=1,ystyle=1,$
            xtitle='E!Lobs!N (keV)',ytitle='polarization angle (deg)',thick=thk
          oplot,shorte,smooth(ang_spec_dir(*,it),3),linestyle=1,thick=thk
          oplot,shorte,smooth(ang_spec_ref(*,it),3),linestyle=2,thick=thk
;          oplot,shorte,smooth(ang_spec_scat(*,it),3),linestyle=3,thick=thk
;          plots,[0.15,0.3],[1,1]*(60.),thick=thk,linestyle=1
;          plots,[0.15,0.3],[1,1]*(45.),thick=thk,linestyle=2
;          plots,[0.15,0.3],[1,1]*(30.),thick=thk
;          xyouts,0.32,60,'direct'
;          xyouts,0.32,45,'reflected'
;          xyouts,0.32,30,'total'
          if (fix(ifig) eq 5) then xyouts,0.14,35+min_deg,'a/M = 0'
          if (fix(ifig) eq 6) then xyouts,0.14,35+min_deg,'a/M = 0.9'
          if (fix(ifig) eq 16) then xyouts,0.14,35+min_deg,'a/M = 0.998'
          if (it eq 5) then xyouts,0.14,15+min_deg,'i = 75!Uo!N'
          if (it eq 10) then xyouts,0.14,15+min_deg,'i = 60!Uo!N'
          if (it eq 14) then xyouts,0.14,15+min_deg,'i = 45!Uo!N'
      endif
      !P.Charsize=1.5
;      stop
  endif

  if ((fix(ifig) ge 7)and(fix(ifig) le 9)) then begin
      Nt = 21
      Nspec = 201
      
      if (Nt eq 11) then begin
          inc_lbl = ['i=87!Uo!N','i=82!Uo!N','i=77!Uo!N','i=71!Uo!N','i=66!Uo!N',$
                     'i=60!Uo!N','i=54!Uo!N','i=47!Uo!N','i=39!Uo!N','i=30!Uo!N','i=17!Uo!N']
          idex75 = 2
          idex60 = 5
          idex45 = 7
      endif
      if (Nt eq 21) then begin
          idex75 = 5
          idex60 = 10
          idex45 = 14
          
          inc_lbl = ['i=89!Uo!N','i=86!Uo!N','i=83!Uo!N','i=80!Uo!N','i=78!Uo!N',$
                     'i=75!Uo!N','i=72!Uo!N','i=69!Uo!N','i=66!Uo!N','i=63!Uo!N',$
                     'i=60!Uo!N','i=57!Uo!N','i=53!Uo!N','i=50!Uo!N','i=46!Uo!N',$
                     'i=42!Uo!N','i=38!Uo!N','i=33!Uo!N','i=28!Uo!N','i=22!Uo!N',$
                     'i=12!Uo!N']
      endif

      spec = dblarr(Nspec,Nt)
      Qspec = dblarr(Nspec,Nt)
      Uspec = dblarr(Nspec,Nt)

      emin = 0.1
      emax = 1000
      shorte = emin*10.^(findgen(Nspec)/(Nspec-1.)*alog10(emax/emin))
      de = deriv(shorte)

      if (fix(ifig) eq 7) then run_ids = [143,144,145,146,147]
      if (fix(ifig) eq 8) then run_ids = [148,149,150,151,152]
;      if (fix(ifig) eq 7) then run_ids = [113,114,115,116,117]
;      if (fix(ifig) eq 8) then run_ids = [118,119,120,121,122]
;      if (fix(ifig) eq 9) then run_ids = [143,148,189,189,189]

      for rundex=0,4 do begin
          run_id = run_ids(rundex)
          rdatafile = 'scat_spec.0000.dat'
          dumpstr = string(run_id,format='(I4.4)')
          strput,rdatafile,dumpstr,10
          openr,1,rdatafile
          readf,1,spec,Qspec,Uspec
          close,1
          deg_spec = sqrt(Qspec^2+Uspec^2)/spec
          ang_spec = atan(Uspec,Qspec)/2.*!radeg

          if (fix(ifig) eq 7) then min_deg = -110
          if (fix(ifig) eq 8) then min_deg = -110
          neg_ang = where(ang_spec lt min_deg)
          if (neg_ang(0) ge 0) then ang_spec(neg_ang)=ang_spec(neg_ang)+180.
          pos_ang = where(ang_spec gt min_deg+180)
          if (pos_ang(0) ge 0) then ang_spec(pos_ang)=ang_spec(pos_ang)-180.
          
          if ((ifig eq 7.1)or(ifig eq 8.1)or(ifig eq 9.1)) then begin
              if (rundex eq 0) then begin
                  plot_oo,shorte,deg_spec(*,idex75)*100.,$
                    xrange=[0.1,20.],yrange=[1d-1,100],xstyle=1,ystyle=1,$
                    xtitle='E!Lobs!N (keV)',ytitle='polarization degree (%)',thick=thk
                  plots,[0.15,0.3],[1,1]*10.^(-0.1),thick=thk,color= 60
                  plots,[0.15,0.3],[1,1]*10.^(-0.25),thick=thk,color= 90
                  plots,[0.15,0.3],[1,1]*10.^(-0.4),thick=thk,color= 120
                  plots,[0.15,0.3],[1,1]*10.^(-0.55),thick=thk,color= 150
                  plots,[0.15,0.3],[1,1]*10.^(-0.7),thick=thk,color= 180
                  xyouts,0.32,10.^(-0.12),'a/M=0'
                  xyouts,0.32,10.^(-0.27),'     0.5'
                  xyouts,0.32,10.^(-0.42),'     0.9'
                  xyouts,0.32,10.^(-0.57),'   0.99'
                  xyouts,0.32,10.^(-0.72),' 0.998'
              endif
              oplot,shorte,deg_spec(*,idex75)*100,thick=thk,color=60+rundex*30
          endif
          
          if ((ifig eq 7.2)or(ifig eq 8.2)or(ifig eq 9.2)) then begin
              if (rundex eq 0) then begin
                  plot_oi,shorte,smooth(ang_spec(*,idex75),1),$
                    xrange=[0.1,20.],yrange=[min_deg,min_deg+180],xstyle=1,ystyle=1,$
                    xtitle='E!Lobs!N (keV)',ytitle='polarization angle (deg)'
                  plots,[0.15,0.3],[1,1]*60+min_deg,thick=thk,color=60
                  plots,[0.15,0.3],[1,1]*50+min_deg,thick=thk,color=90
                  plots,[0.15,0.3],[1,1]*40+min_deg,thick=thk,color=120
                  plots,[0.15,0.3],[1,1]*30+min_deg,thick=thk,color=150
                  plots,[0.15,0.3],[1,1]*20+min_deg,thick=thk,color=180
                  xyouts,0.32,58+min_deg,'a/M=0'
                  xyouts,0.32,48+min_deg,'     0.5'
                  xyouts,0.32,38+min_deg,'     0.9'
                  xyouts,0.32,28+min_deg,'   0.99'
                  xyouts,0.32,18+min_deg,' 0.998'
              endif
              oplot,shorte,smooth(ang_spec(*,idex75),1),thick=thk,color=60+rundex*30
          endif
      endfor
  endif

  if ((fix(ifig) ge 10)and(fix(ifig) le 11)) then begin
      !P.charsize=1.3
      Nt = 11
      Nspec = 201
      
      if (Nt eq 11) then begin
          inc_lbl = ['i=87!Uo!N','i=82!Uo!N','i=77!Uo!N','i=71!Uo!N','i=66!Uo!N',$
                     'i=60!Uo!N','i=54!Uo!N','i=47!Uo!N','i=39!Uo!N','i=30!Uo!N','i=17!Uo!N']
          idex75 = 2
          idex60 = 5
          idex45 = 7
      endif
      if (Nt eq 21) then begin
          idex75 = 5
          idex60 = 10
          idex45 = 14
          
          inc_lbl = ['i=89!Uo!N','i=86!Uo!N','i=83!Uo!N','i=80!Uo!N','i=78!Uo!N',$
                     'i=75!Uo!N','i=72!Uo!N','i=69!Uo!N','i=66!Uo!N','i=63!Uo!N',$
                     'i=60!Uo!N','i=57!Uo!N','i=53!Uo!N','i=50!Uo!N','i=46!Uo!N',$
                     'i=42!Uo!N','i=38!Uo!N','i=33!Uo!N','i=28!Uo!N','i=22!Uo!N',$
                     'i=12!Uo!N']
      endif

      spec = dblarr(Nspec,Nt)
      Qspec = dblarr(Nspec,Nt)
      Uspec = dblarr(Nspec,Nt)

      emin = 0.1
      emax = 1000
      shorte = emin*10.^(findgen(Nspec)/(Nspec-1.)*alog10(emax/emin))
      de = deriv(shorte)

      if (fix(ifig) eq 10) then run_ids = [169,170,113,171,172]
      if (fix(ifig) eq 11) then run_ids = [194,195,118,196,197]

      for rundex=0,4 do begin
          run_id = run_ids(rundex)
          rdatafile = 'scat_spec.0000.dat'
          dumpstr = string(run_id,format='(I4.4)')
          strput,rdatafile,dumpstr,10
          openr,1,rdatafile
          readf,1,spec,Qspec,Uspec
          close,1
          deg_spec = sqrt(Qspec^2+Uspec^2)/spec
          ang_spec = atan(Uspec,Qspec)/2.*!radeg

          if (fix(ifig) eq 10) then min_deg = -110
          if (fix(ifig) eq 11) then min_deg = -110
          neg_ang = where(ang_spec lt min_deg)
          if (neg_ang(0) ge 0) then ang_spec(neg_ang)=ang_spec(neg_ang)+180.
          pos_ang = where(ang_spec gt min_deg+180)
          if (pos_ang(0) ge 0) then ang_spec(pos_ang)=ang_spec(pos_ang)-180.
          
          if ((ifig eq 10.1)or(ifig eq 11.1)or(ifig eq 11.1)) then begin
              if (rundex eq 0) then begin
                  plot_oo,shorte,deg_spec(*,idex75)*100.,$
                    xrange=[0.1,20],yrange=[1d-1,100],xstyle=1,ystyle=1,$
                    xtitle='E!Lobs!N (keV)',ytitle='polarization degree (%)',thick=thk
                  plots,[.15,.30],[1,1]*10.^(1.6),thick=thk,color= 60
                  plots,[.15,.30],[1,1]*10.^(1.45),thick=thk,color= 90
                  plots,[.15,.30],[1,1]*10.^(1.3),thick=thk,color= 120
                  plots,[.15,.30],[1,1]*10.^(1.15),thick=thk,color= 150
                  plots,[.15,.30],[1,1]*10.^(1.0),thick=thk,color= 180
                  xyouts,.32,10.^(1.58),'L/L!LEdd!N=0.01'
                  xyouts,.32,10.^(1.43),'           0.03'
                  xyouts,.32,10.^(1.28),'             0.1'
                  xyouts,.32,10.^(1.13),'             0.3'
                  xyouts,.32,10.^(0.98),'             1.0'
                  plots,[0.15,.3],[1,1]*10^(-0.6),thick=thk
                  plots,[0.15,.3],[1,1]*10^(-0.8),thick=thk,linestyle=2
                  xyouts,.32,10^(-0.62),'a/M=0'
                  xyouts,.32,10^(-0.82),'     0.9'
              endif
              oplot,shorte,deg_spec(*,idex75)*100,thick=thk,color=60+rundex*30
          endif
          
          if ((ifig eq 10.2)or(ifig eq 11.2)) then begin
              neg_ang = where(ang_spec lt min_deg)
              if (neg_ang(0) ge 0) then ang_spec(neg_ang)=ang_spec(neg_ang)+180.
              pos_ang = where(ang_spec gt min_deg+180)
              if (pos_ang(0) ge 0) then ang_spec(pos_ang)=ang_spec(pos_ang)-180.
              if (rundex eq 0) then begin
                  plot_oi,shorte,smooth(ang_spec(*,idex75),1),$
                    xrange=[0.1,20.],yrange=[min_deg,min_deg+180.],$
                    xstyle=1,ystyle=1,$
                    xtitle='E!Lobs!N (keV)',ytitle='polarization angle (deg)'
                  max_deg = min_deg+190
                  plots,[0.15,0.3],[1,1]*(-30)+max_deg,thick=thk,color=60
                  plots,[0.15,0.3],[1,1]*(-40)+max_deg,thick=thk,color=90
                  plots,[0.15,0.3],[1,1]*(-50)+max_deg,thick=thk,color=120
                  plots,[0.15,0.3],[1,1]*(-60)+max_deg,thick=thk,color=150
                  plots,[0.15,0.3],[1,1]*(-70)+max_deg,thick=thk,color=180
                  xyouts,0.32,-32+max_deg,'L/L!LEdd!N=0.01'
                  xyouts,0.32,-42+max_deg,'           0.03'
                  xyouts,0.32,-52+max_deg,'             0.1'
                  xyouts,0.32,-62+max_deg,'             0.3'
                  xyouts,0.32,-72+max_deg,'             1.0'
                  max_deg = min_deg+90
                  plots,[.15,.3],[1,1]*(-60)+max_deg,thick=thk
                  plots,[.15,.3],[1,1]*(-70)+max_deg,thick=thk,linestyle=2
                  xyouts,.32,-62+max_deg,'a/M=0'
                  xyouts,.32,-72+max_deg,'     0.9'
              endif
              oplot,shorte,smooth(ang_spec(*,idex75),1),thick=thk,color=60+rundex*30
          endif
      endfor

      if (fix(ifig) eq 10) then run_ids = [190,191,115,192,193]
      if (fix(ifig) eq 11) then run_ids = [198,199,120,200,201]
;      if (fix(ifig) eq 10) then run_ids = [101,125,115,192,193]

      for rundex=0,4 do begin
          run_id = run_ids(rundex)
          rdatafile = 'scat_spec.0000.dat'
          dumpstr = string(run_id,format='(I4.4)')
          strput,rdatafile,dumpstr,10
          openr,1,rdatafile
          readf,1,spec,Qspec,Uspec
          close,1
          deg_spec = sqrt(Qspec^2+Uspec^2)/spec
          ang_spec = atan(Uspec,Qspec)/2.*!radeg
          
          if ((ifig eq 10.1)or(ifig eq 11.1)or(ifig eq 11.1)) then begin
              oplot,shorte,deg_spec(*,idex75)*100,thick=thk,color=60+rundex*30,linestyle=2
          endif
          
          if ((ifig eq 10.2)or(ifig eq 11.2)) then begin
              neg_ang = where(ang_spec lt min_deg)
              if (neg_ang(0) ge 0) then ang_spec(neg_ang)=ang_spec(neg_ang)+180.
              pos_ang = where(ang_spec gt min_deg+180)
              if (pos_ang(0) ge 0) then ang_spec(pos_ang)=ang_spec(pos_ang)-180.
              oplot,shorte,smooth(ang_spec(*,idex75),1),thick=thk,color=60+rundex*30,linestyle=2
          endif
      endfor
  endif

  if ((fix(ifig) eq 100)or(fix(ifig) eq 101)) then begin
      !P.charsize=1.3
      Nt = 11
      Nspec = 201
      
      if (Nt eq 11) then begin
          inc_lbl = ['i=87!Uo!N','i=82!Uo!N','i=77!Uo!N','i=71!Uo!N','i=66!Uo!N',$
                     'i=60!Uo!N','i=54!Uo!N','i=47!Uo!N','i=39!Uo!N','i=30!Uo!N','i=17!Uo!N']
          idex75 = 8
          idex60 = 5
          idex45 = 7
      endif
      if (Nt eq 21) then begin
          idex75 = 5
          idex60 = 10
          idex45 = 14
          
          inc_lbl = ['i=89!Uo!N','i=86!Uo!N','i=83!Uo!N','i=80!Uo!N','i=78!Uo!N',$
                     'i=75!Uo!N','i=72!Uo!N','i=69!Uo!N','i=66!Uo!N','i=63!Uo!N',$
                     'i=60!Uo!N','i=57!Uo!N','i=53!Uo!N','i=50!Uo!N','i=46!Uo!N',$
                     'i=42!Uo!N','i=38!Uo!N','i=33!Uo!N','i=28!Uo!N','i=22!Uo!N',$
                     'i=12!Uo!N']
      endif

      spec = dblarr(Nspec,Nt)
      Qspec = dblarr(Nspec,Nt)
      Uspec = dblarr(Nspec,Nt)

      emin = 0.1
      emin = 0.001
      emax = 1000
      shorte = emin*10.^(findgen(Nspec)/(Nspec-1.)*alog10(emax/emin))
      de = deriv(shorte)

      if (fix(ifig) eq 100) then run_ids = [123,110,202,171,172]
      if (fix(ifig) eq 101) then run_ids = [166,167,101,168,168]
      if (fix(ifig) eq 101) then run_ids = [215,216,218,217,168]

      for rundex=0,3 do begin
          run_id = run_ids(rundex)
          rdatafile = 'scat_spec.0000.dat'
          dumpstr = string(run_id,format='(I4.4)')
          strput,rdatafile,dumpstr,10
          openr,1,rdatafile
          readf,1,spec,Qspec,Uspec
          close,1
          deg_spec = sqrt(Qspec^2+Uspec^2)/spec
          ang_spec = atan(Uspec,Qspec)/2.*!radeg

          if (fix(ifig) eq 100) then min_deg = -120
          if (fix(ifig) eq 101) then min_deg = -100
          neg_ang = where(ang_spec lt min_deg)
          if (neg_ang(0) ge 0) then ang_spec(neg_ang)=ang_spec(neg_ang)+180.
          pos_ang = where(ang_spec gt min_deg+180)
          if (pos_ang(0) ge 0) then ang_spec(pos_ang)=ang_spec(pos_ang)-180.
          
          if ((ifig eq 100.1)) then begin
              if (rundex eq 0) then begin
                  plot_oo,shorte,deg_spec(*,idex75)*100.,$
                    xrange=[0.1,20],yrange=[1d-1,10],xstyle=1,ystyle=1,$
                    xtitle='E!Lobs!N (keV)',ytitle='polarization degree (%)',thick=thk
;                  plots,[15,30],[1,1]*10.^(-0.4),thick=thk,color= 60
;                  plots,[15,30],[1,1]*10.^(-0.5),thick=thk,color= 90
;                  plots,[15,30],[1,1]*10.^(-0.6),thick=thk,color= 120
;                  plots,[15,30],[1,1]*10.^(-0.7),thick=thk,color= 150
;                  plots,[15,30],[1,1]*10.^(-0.8),thick=thk,color= 180
                  plots,[.15,.30],[1,1]*10.^(-0.3),thick=thk,color= 120
                  plots,[.15,.30],[1,1]*10.^(-0.4),thick=thk,color= 90
                  plots,[.15,.30],[1,1]*10.^(-0.5),thick=thk,color= 60
                  xyouts,0.32,10.^(-0.32),'tau!Les!N=2.0, T!Lcor!N=50 keV'
                  xyouts,0.32,10.^(-0.42),'tau!Les!N=1.0, T!Lcor!N=100 keV'
                  xyouts,0.32,10.^(-0.52),'tau!Les!N=0.5, T!Lcor!N=200 keV'
;                  xyouts,32,10.^(-0.42),'L/L!LEdd!N=0.01'
;                  xyouts,32,10.^(-0.52),'           0.03'
;                  xyouts,32,10.^(-0.62),'             0.1'
;                  xyouts,32,10.^(-0.72),'             0.3'
;                  xyouts,32,10.^(-0.82),'             1.0'
                  plots,[.15,.4],[1,1]*10^(-0.7),thick=thk
                  plots,[.15,.4],[1,1]*10^(-0.8),thick=thk,linestyle=2
                  xyouts,.42,10^(-0.72),'a/M=0'
                  xyouts,.42,10^(-0.82),'     0.9'
              endif
              oplot,shorte,deg_spec(*,idex75)*100,thick=thk,color=60+rundex*30
          endif
          
          if ((ifig eq 101.1)) then begin
              if (rundex eq 0) then begin
                  plot_oo,shorte,deg_spec(*,idex75)*100.,$
                    xrange=[0.1,20],yrange=[1d-1,20],xstyle=1,ystyle=1,$
                    xtitle='E!Lobs!N (keV)',ytitle='polarization degree (%)',thick=thk
                  plots,[.15,.30],[1,1]*10.^(-0.4),thick=thk,color= 60
                  plots,[.15,.30],[1,1]*10.^(-0.5),thick=thk,color= 90
                  plots,[.15,.30],[1,1]*10.^(-0.6),thick=thk,color= 120
                  plots,[.15,.30],[1,1]*10.^(-0.7),thick=thk,color= 150
;                  plots,[.15,.30],[1,1]*10.^(-0.8),thick=thk,color= 180
                  xyouts,.32,10.^(-0.42),'corona H/R=0.01'
                  xyouts,.32,10.^(-0.52),'                    0.03'
                  xyouts,.32,10.^(-0.62),'                      0.1'
                  xyouts,.32,10.^(-0.72),'                      0.3'
;                  xyouts,.32,10.^(-0.82),'                       1.0'
              endif
              oplot,shorte,deg_spec(*,idex75)*100,thick=thk,color=60+rundex*30
          endif
          
          if (ifig eq 100.2) then begin
              if (rundex ge 1) then min_deg = -50
              neg_ang = where(ang_spec lt min_deg)
              if (neg_ang(0) ge 0) then ang_spec(neg_ang)=ang_spec(neg_ang)+180.
              pos_ang = where(ang_spec gt min_deg+180)
              if (pos_ang(0) ge 0) then ang_spec(pos_ang)=ang_spec(pos_ang)-180.
              if (rundex eq 0) then begin
                  plot_oi,shorte,smooth(ang_spec(*,idex75),1),$
                    xrange=[0.1,20.],yrange=[min_deg,min_deg+220],xstyle=1,ystyle=1,$
                    xtitle='E!Lobs!N (keV)',ytitle='polarization angle (deg)'
                  plots,[0.15,0.3],[1,1]*(60)+min_deg,thick=thk,color=120
                  plots,[0.15,0.3],[1,1]*(70)+min_deg,thick=thk,color=90
                  plots,[0.15,0.3],[1,1]*(80)+min_deg,thick=thk,color=60
;                  plots,[0.15,0.3],[1,1]*(-50),thick=thk,color=120
;                  plots,[0.15,0.3],[1,1]*(-60),thick=thk,color=150
;                  plots,[0.15,0.3],[1,1]*(-70),thick=thk,color=180
                  xyouts,0.32,58+min_deg,'tau!Les!N=2.0, T!Lcor!N=50 keV'
                  xyouts,0.32,68+min_deg,'tau!Les!N=1.0, T!Lcor!N=100 keV'
                  xyouts,0.32,78+min_deg,'tau!Les!N=0.5, T!Lcor!N=200 keV'
;                  xyouts,0.32,-32,'L/L!LEdd!N=0.01'
;                  xyouts,0.32,-42,'           0.03'
;                  xyouts,0.32,-52,'             0.1'
;                  xyouts,0.32,-62,'             0.3'
;                  xyouts,0.32,-72,'             1.0'
                  plots,[.15,.3],[1,1]*(30+min_deg),thick=thk
                  plots,[.15,.3],[1,1]*(20+min_deg),thick=thk,linestyle=2
                  xyouts,.32,28+min_deg,'a/M=0'
                  xyouts,.32,18+min_deg,'     0.9'
              endif
              oplot,shorte,smooth(ang_spec(*,idex75),1),thick=thk,color=60+rundex*30
          endif
          if (ifig eq 101.2) then begin
              if (rundex ge 2) then min_deg = -10
              neg_ang = where(ang_spec lt min_deg)
              if (neg_ang(0) ge 0) then ang_spec(neg_ang)=ang_spec(neg_ang)+180.
              pos_ang = where(ang_spec gt min_deg+180)
              if (pos_ang(0) ge 0) then ang_spec(pos_ang)=ang_spec(pos_ang)-180.
              if (rundex eq 0) then begin
                  plot_oi,shorte,smooth(ang_spec(*,idex75),1),$
                    xrange=[0.1,20.],yrange=[-100,100],xstyle=1,ystyle=1,$
                    xtitle='E!Lobs!N (keV)',ytitle='polarization angle (deg)'
                  plots,[0.15,0.3],[1,1]*60,thick=thk,color=60
                  plots,[0.15,0.3],[1,1]*50,thick=thk,color=90
                  plots,[0.15,0.3],[1,1]*40,thick=thk,color=120
                  plots,[0.15,0.3],[1,1]*30,thick=thk,color=150
;                  plots,[0.15,0.3],[1,1]*20,thick=thk,color=180
                  xyouts,0.32,58,'corona H/R=0.01'
                  xyouts,0.32,48,'                    0.03'
                  xyouts,0.32,38,'                      0.1'
                  xyouts,0.32,28,'                      0.3'
;                  xyouts,0.32,18,'                       1.0'
              endif
              oplot,shorte,smooth(ang_spec(*,idex75),1),thick=thk,color=60+rundex*30
          endif
      endfor

      if (fix(ifig) eq 100) then run_ids = [125,101,204,192,193]
      if (fix(ifig) eq 11) then run_ids = [198,199,120,200,201]
;      if (fix(ifig) eq 10) then run_ids = [101,125,115,192,193]

      if (fix(ifig) eq 100) then begin
          for rundex=0,2 do begin
              run_id = run_ids(rundex)
              rdatafile = 'scat_spec.0000.dat'
              dumpstr = string(run_id,format='(I4.4)')
              strput,rdatafile,dumpstr,10
              openr,1,rdatafile
              readf,1,spec,Qspec,Uspec
              close,1
              deg_spec = sqrt(Qspec^2+Uspec^2)/spec
              ang_spec = atan(Uspec,Qspec)/2.*!radeg
              
              if ((ifig eq 100.1)or(ifig eq 11.1)or(ifig eq 11.1)) then begin
                  oplot,shorte,deg_spec(*,idex75)*100,thick=thk,color=60+rundex*30,linestyle=2
              endif
              
              if ((ifig eq 100.2)or(ifig eq 11.2)) then begin
                  if (rundex eq 0) then min_deg = -120
                  if (rundex ge 1) then min_deg = -50
                  neg_ang = where(ang_spec lt min_deg)
                  if (neg_ang(0) ge 0) then ang_spec(neg_ang)=ang_spec(neg_ang)+180.
                  pos_ang = where(ang_spec gt min_deg+180)
                  if (pos_ang(0) ge 0) then ang_spec(pos_ang)=ang_spec(pos_ang)-180.
                  oplot,shorte,smooth(ang_spec(*,idex75),1),thick=thk,color=60+rundex*30,linestyle=2
              endif
          endfor
      endif
  endif

  if ((fix(ifig) ge 102)and(fix(ifig) le 102)) then begin
      !P.charsize=1.3
      Nt = 11
      Nspec = 201
      
      if (Nt eq 11) then begin
          inc_lbl = ['i=87!Uo!N','i=82!Uo!N','i=77!Uo!N','i=71!Uo!N','i=66!Uo!N',$
                     'i=60!Uo!N','i=54!Uo!N','i=47!Uo!N','i=39!Uo!N','i=30!Uo!N','i=17!Uo!N']
          idex75 = 2
          idex60 = 5
          idex45 = 7
      endif
      if (Nt eq 21) then begin
          idex75 = 5
          idex60 = 10
          idex45 = 14
          
          inc_lbl = ['i=89!Uo!N','i=86!Uo!N','i=83!Uo!N','i=80!Uo!N','i=78!Uo!N',$
                     'i=75!Uo!N','i=72!Uo!N','i=69!Uo!N','i=66!Uo!N','i=63!Uo!N',$
                     'i=60!Uo!N','i=57!Uo!N','i=53!Uo!N','i=50!Uo!N','i=46!Uo!N',$
                     'i=42!Uo!N','i=38!Uo!N','i=33!Uo!N','i=28!Uo!N','i=22!Uo!N',$
                     'i=12!Uo!N']
      endif

      spec = dblarr(Nspec,Nt)
      Qspec = dblarr(Nspec,Nt)
      Uspec = dblarr(Nspec,Nt)

      emin = 0.001
      emax = 1000
      shorte = emin*10.^(findgen(Nspec)/(Nspec-1.)*alog10(emax/emin))
      de = deriv(shorte)

      run_id = 212
      rdatafile = 'scat_spec.0000.dat'
      dumpstr = string(run_id,format='(I4.4)')
      strput,rdatafile,dumpstr,10
      openr,1,rdatafile
      readf,1,spec,Qspec,Uspec
      close,1
      deg_spec = sqrt(Qspec^2+Uspec^2)/spec
      ang_spec = atan(Uspec,Qspec)/2.*!radeg

      min_deg = 0
      neg_ang = where(ang_spec lt min_deg)
      if (neg_ang(0) ge 0) then ang_spec(neg_ang)=ang_spec(neg_ang)+180.
      pos_ang = where(ang_spec gt min_deg+180)
      if (pos_ang(0) ge 0) then ang_spec(pos_ang)=ang_spec(pos_ang)-180.
      it = idex75
      if (ifig eq 102.1) then begin
          plot_oo,shorte,smooth(deg_spec(*,it),1)*100.,$
            xrange=[0.002,20.],yrange=[2d-1,20],xstyle=1,ystyle=1,$
            xtitle='E!Lobs!N (keV)',ytitle='polarization degree (%)',thick=thk
          oplot,shorte,deg_spec(*,idex60)*100.,thick=thk,color = 120
          oplot,shorte,deg_spec(*,idex45)*100.,thick=thk,color = 60
          xyouts,1.5,deg_spec(Nspec/2,idex75)*120,'i = 75!Uo!N'
          xyouts,1.5,deg_spec(Nspec/2,idex60)*120,'i = 60!Uo!N',color=120
          xyouts,1.5,deg_spec(Nspec/2,idex45)*120,'i = 45!Uo!N',color=60
      endif
      if (ifig eq 102.2) then begin
          plot_oi,shorte,smooth(ang_spec(*,it),3),$
            xrange=[0.002,20],yrange=[min_deg,min_deg+180],xstyle=1,ystyle=1,$
            xtitle='E!Lobs!N (keV)',ytitle='polarization angle (deg)',thick=thk
          oplot,shorte,ang_spec(*,idex60),thick=thk,color = 120
          oplot,shorte,ang_spec(*,idex45),thick=thk,color = 60
          if (it eq 5) then xyouts,0.14,15+min_deg,'i = 75!Uo!N'
          if (it eq 10) then xyouts,0.14,15+min_deg,'i = 60!Uo!N'
          if (it eq 14) then xyouts,0.14,15+min_deg,'i = 45!Uo!N'
      endif
      
  endif

  if (ifig eq 20) then begin
      Nr = 1001.
      Rmin = 1.5
      Rmax = 100.
      rr = Rmin*10.^(2.*findgen(Nr)/(Nr-1.))
      dr = deriv(rr)

      aa = 0.9
      em_model = 2
      Mbh = 10.
      Mstar = Mbh/3.
      Rcm = Mbh*1.45d5
      Mdstar = 1.0
      Mdot = Mdstar*(1.0d17) 
      f_hard = 1.8
      Lstar = 0.1
      
      kapp = 0.4
      ak = 7.566d-15
      sigma = 5.67d-5
      kB = 1.38d-16
      kB_ev = 8.6173d-5
      c = 3.d10
      G = 6.67d-8
      Z1 = 1.+(1-aa*aa)^(1./3.)*((1.+aa)^(1./3.)+(1-aa)^(1./3.))
      Z2 = sqrt(3.*aa*aa+Z1^2.)
      Risco = (3.+Z2-((3.-Z1)*(3+Z1+2.*Z2))^0.5)
      
      if (em_model ge 2) then begin
          
          ntA = 1.+aa^2./rr^2.+2.*aa^2./rr^3.
          ntB = 1.+aa/rr^1.5
          ntC = 1.-3./rr+2.*aa/rr^1.5
          ntD = 1.-2./rr+aa^2./rr^2.
          ntE = 1.+4.*aa^2./rr^2.-4.*aa^2./rr^3.+3.*aa^4./rr^4.
          ntE2 = 2.+5.*aa^2./rr^2.-2.*aa^2./rr^3.+3.*aa^4./rr^4.
          ntF = 1.-2.*aa/rr^1.5+aa^2./rr^2.
          ntJ = dblarr(Nr)
          Nu = Nr
          for i=0,Nr-1 do begin
              if (rr(i) gt Risco) then begin
                  u_star = 1./rr(i)
                  du = u_star/(Nu-1.)
                  subu = findgen(Nu)*du
                  ntJ(i) = exp(1.5*du*total((1.-2.*aa*subu^1.5+aa^2.*subu^2.)/ $
                                            ((1.+aa*subu^1.5)*(1.-3.*subu+2.*aa*subu^1.5))))
              endif
          endfor
          
          L_ms = 2./sqrt(3.)*(3.*sqrt(Risco)-2.*aa)/sqrt(Risco)
          ntL = ntF/sqrt(ntC)-L_ms/sqrt(rr)
          ntQ = dblarr(Nr)
          iisco = where(rr gt Risco)
          iisco = iisco(0)
          for i=iisco,Nr-1 do begin
              Q_int = ntL(iisco:i)*ntF(iisco:i)/(ntB(iisco:i)*ntC(iisco:i)*ntJ(iisco:i)*rr(iisco:i)^1.5)
              ntQ(i) = ntL(i)-3.*ntJ(i)/(2.*sqrt(rr(i)))* $
                (total(Q_int(0:i-iisco)*dr(iisco:i))-0.5*(Q_int(0)*dr(iisco)+Q_int(i-iisco)*dr(i)))
          endfor
          nt_flux = 6.d25*(Mdstar/Mstar^2.)*ntQ/(rr^3.*ntB*sqrt(ntC))
          nt_flux(0:iisco)=0.
          Ltot = 4.*!PI*total(nt_flux*rr*dr)*Rcm^2.
          nt_flux = nt_flux*Lstar/(Ltot/(Mbh*1.3d38))
      endif
      dLdr = deriv(alog10(rr),alog10(nt_flux+0.1))
      nt_flux = nt_flux/max(nt_flux)+0.00001
      plot_oo,rr,nt_flux,xtitle='R (M)',ytitle = 'local flux',$
        yrange=[1d-3,100],ystyle=1,xrange=[1,100],xstyle=1,thick=thk
      alpha = -2
      if ((-alpha) lt dLdr(Nr-2)) then alpha = -dLdr(Nr-2) 
      Rmatch = where(abs(dLdr(iisco:Nr-1)+alpha) eq min(abs(dLdr(iisco:Nr-1)+alpha)))
      pdex = Rmatch(0)+iisco
      Rmatch = rr(pdex)
      pl_flux = nt_flux
      pl_flux(0:pdex)=(rr(0:pdex)/Rmatch)^(-alpha)*pl_flux(pdex)
      oplot,rr,pl_flux,thick=thk,linestyle=2
      alpha = 0.
      if ((-alpha) lt dLdr(Nr-2)) then alpha = -dLdr(Nr-2) 
      Rmatch = where(abs(dLdr(iisco:Nr-1)+alpha) eq min(abs(dLdr(iisco:Nr-1)+alpha)))
      pdex = Rmatch(0)+iisco
      Rmatch = rr(pdex)
      pl_flux = nt_flux
      pl_flux(0:pdex)=(rr(0:pdex)/Rmatch)^(-alpha)*pl_flux(pdex)
      oplot,rr,pl_flux,thick=thk,linestyle=2
      alpha = 2
      if ((-alpha) lt dLdr(Nr-2)) then alpha = -dLdr(Nr-2) 
      Rmatch = where(abs(dLdr(iisco:Nr-1)+alpha) eq min(abs(dLdr(iisco:Nr-1)+alpha)))
      pdex = Rmatch(0)+iisco
      Rmatch = rr(pdex)
      pl_flux = nt_flux
      pl_flux(0:pdex)=(rr(0:pdex)/Rmatch)^(-alpha)*pl_flux(pdex)
      oplot,rr,pl_flux,thick=thk,linestyle=2
      xyouts,2.7,0.1,'NT'
      xyouts,1.2,0.1,'!9a!3=-2'
      xyouts,1.2,0.5,'!9a!3=0'
      xyouts,1.2,10,'!9a!3=2'
      plots,[10,20],[5,5],thick=thk,linestyle=1
      xyouts,22,4,'HARM3D'

      rdata = fltarr(2,192)
      openr,1,'scn_data.dat'
      readf,1,rdata
      close,1
      rdata(1,*)=rdata(1,*)/1.e-5
      oplot,rdata(0,*),rdata(1,*),thick=thk,linestyle=1

  endif

  if ((fix(ifig) eq 21)or(fix(ifig) eq 22)) then begin
      Naa2 = 44.
      Nal2 = 40.
      alpha_max = 3.
      if (fix(ifig) eq 21) then nu = 10.-4.
      if (fix(ifig) eq 22) then nu = 76.-4.
      signu = sqrt(2.*nu)
      new_aa = (findgen(Naa2)/(Naa2-1.)*(-3))#(fltarr(Nal2)+1.)
      new_al = (fltarr(Naa2)+1.)#(findgen(Nal2)/(Nal2-1.))
      new_con = fltarr(Naa2,Nal2)
      iframe = ifig-fix(ifig)
      if ((ifig eq 21.1)or(ifig eq 22.1)) then begin
          if (fix(ifig) eq 21) then openr,1,'polcont_1a.dat'
          if (fix(ifig) eq 22) then openr,1,'polcont_2a.dat'
          atarg = 0
          alphatarg = -30
      endif
      if ((ifig eq 21.2)or(ifig eq 22.2)) then begin
          if (fix(ifig) eq 21) then openr,1,'polcont_1b.dat'
          if (fix(ifig) eq 22) then openr,1,'polcont_2b.dat'
          atarg = 0.97
          alphatarg = 1
      endif
      if ((ifig eq 21.3)or(ifig eq 22.3)) then begin
          if (fix(ifig) eq 21) then openr,1,'polcont_1c.dat'
          if (fix(ifig) eq 22) then openr,1,'polcont_2c.dat'
          atarg = 0.9977
          alphatarg = -30
      endif
      if ((ifig eq 21.4)or(ifig eq 22.4)) then begin
          if (fix(ifig) eq 21) then openr,1,'polcont_1d.dat'
          if (fix(ifig) eq 22) then openr,1,'polcont_2d.dat'
          atarg = 0
          alphatarg = -30
      endif
      if ((ifig eq 21.5)or(ifig eq 22.5)) then begin
          if (fix(ifig) eq 21) then openr,1,'polcont_1e.dat'
          if (fix(ifig) eq 22) then openr,1,'polcont_2e.dat'
          atarg = 0.97
          alphatarg = 1
      endif
      if ((ifig eq 21.6)or(ifig eq 22.6)) then begin
          if (fix(ifig) eq 21) then openr,1,'polcont_1f.dat'
          if (fix(ifig) eq 22) then openr,1,'polcont_2f.dat'
          atarg = 0.9977
          alphatarg = -30
      endif
      if ((ifig eq 21.7)or(ifig eq 22.7)) then begin
          if (fix(ifig) eq 21) then openr,1,'polcont_1g.dat'
          if (fix(ifig) eq 22) then openr,1,'polcont_2g.dat'
          atarg = 0
          alphatarg = -30
      endif
      if ((ifig eq 21.8)or(ifig eq 22.8)) then begin
          if (fix(ifig) eq 21) then openr,1,'polcont_1h.dat'
          if (fix(ifig) eq 22) then openr,1,'polcont_2h.dat'
          atarg = 0.97
          alphatarg = 1
      endif
      if ((ifig eq 21.9)or(ifig eq 22.9)) then begin
          if (fix(ifig) eq 21) then openr,1,'polcont_1i.dat'
          if (fix(ifig) eq 22) then openr,1,'polcont_2i.dat'
          atarg = 0.9977
          alphatarg = -30
      endif
      readf,1,new_con
      if ((ifig eq 22.8)or(ifig eq 22.7)) then begin
          ldata = fltarr(5,(1L)*(22L)*(20L)*(31L)*(5*2+3))
          readf,1,ldata
      endif
      close,1
      contour,new_con,new_aa,new_al,xrange=[0,-3],$
        levels=alog10([1,2,3,4,5]*signu),/fill,c_colors=findgen(5)*40+40.,$
        xtickname=['0','0.9','0.99','0.999'],xtitle='spin a/M',$
        ytickname=['-Inf','-1.8','0.25','1.5','2.3','3.0'],ytitle='!9a!3'

      xtrg = alog10(1.-atarg)
      ytrg = exp(alphatarg/alpha_max-1.)
      xyouts,xtrg,ytrg,'X'
;      stop
  endif
      
device, /close
END

