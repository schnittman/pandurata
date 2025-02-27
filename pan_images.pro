PRO pan_images,ifig,run_id
!P.font = 0
!P.charsize = 1.5
thk = 5

if (ifig eq 1) then begin
   N = 201.
   Nt = 41.
   Ne_i = 11
   Nspec = 101
   run_sort = 0
;   run_id = 1251
   Nend = Nt

   Ixy = dblarr(N,N)
   wght = dblarr(N,N)
   mov = dblarr(N,N,Nt)
   movx = dblarr(N,N,Nt)
   movy = dblarr(N,N,Nt)
   smov = dblarr(Ne_i,N,N,Nt)
   smovx = dblarr(Ne_i,N,N,Nt)
   smovy = dblarr(Ne_i,N,N,Nt)
   mov2 = dblarr(Ne_i,N,N)
   movx2 = dblarr(Ne_i,N,N)
   movy2 = dblarr(Ne_i,N,N)

   inc_lbl = ['i=168!Uo!N',$
              'i=158!Uo!N','i=152!Uo!N','i=147!Uo!N','i=142!Uo!N','i=138!Uo!N',$
              'i=134!Uo!N','i=130!Uo!N','i=127!Uo!N','i=123!Uo!N','i=120!Uo!N',$
              'i=117!Uo!N','i=114!Uo!N','i=111!Uo!N','i=108!Uo!N','i=105!Uo!N',$
              'i=102!Uo!N','i=100!Uo!N','i=97!Uo!N','i=94!Uo!N','i=90!Uo!N',$
              'i=86!Uo!N','i=83!Uo!N','i=80!Uo!N','i=78!Uo!N',$
              'i=75!Uo!N','i=72!Uo!N','i=69!Uo!N','i=66!Uo!N','i=63!Uo!N',$
              'i=60!Uo!N','i=57!Uo!N','i=53!Uo!N','i=50!Uo!N','i=46!Uo!N',$
              'i=42!Uo!N','i=38!Uo!N','i=33!Uo!N','i=28!Uo!N','i=22!Uo!N',$
              'i=12!Uo!N']
   rdatafile = 'data/scat_imag.0000.dat'
   dumpstr = string(run_id,format='(I4.4)')
   strput,rdatafile,dumpstr,15
   openr,1,rdatafile
   readf,1,mov,movx,movy
   close,1
      
   for i=0,Nt-1 do begin
      mov(*,*,i)=transpose(mov(*,*,i))
      movx(*,*,i)=transpose(movx(*,*,i))
      movy(*,*,i)=transpose(movy(*,*,i))
   endfor
   sortmov = mov(sort(mov))
   movmax = max(mov)
   movmax=sortmov(1.*N*N*Nt-100.)
;   movmax=sortmov(1.*N*N*Nt-1.)
   outliers = where(mov gt movmax)
   if (outliers(0) ge 0) then begin
      movx(outliers)=movx(outliers)/mov(outliers)*movmax
      movy(outliers)=movy(outliers)/mov(outliers)*movmax
      mov(outliers)=movmax
   endif
   
   for it=Nend-1,0,-1 do begin
      y_min = 0
      if (N eq 81) then Ixy = dblarr(121,121)
      if (N eq 201) then Ixy = dblarr(301,301)
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
      data = byte(255*(alog10(Ixy/movmax+1e-5)+5.01)/5.)
      N15 = (N-1)*1.5+1
      N2 = (fix(600./N15))*N15
      enlarge = rebin(data,N2,N2)
      
      erase
      tvscl,enlarge,0,y_min,/Device
      dd = 100.
      pstep=10
      maxIxy = max(Ixy)
      for i=N15/6.,N15*(5./6.),pstep do begin
         for j=N15/3.,N15,pstep do begin
            x0 = i/(N15-1.)*12600.
            y0 = j/(N15-1.)*12600.
            dff = 100.*dd*([Xpol(i,j),Ypol(i,j),0])
              if (Ixy(i,j) gt 1d-4*maxIxy) then $
                plots,[x0-dff(0)/2.,x0+dff(0)/2.],[y0-dff(1)/2.,y0+dff(1)/2.],$
                color=255,/Device,thick=4
         endfor
      endfor
      plots,[11000,11000+5.*dd],[12000,12000],color=255,thick=thk,/Device
      xyouts,11000,11500,'deg=5%',color=255,/Device
      
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
      contour,enlarge,Position=[1600,6000,2000,12000],/Device,$
              yticks=1,yminor=1,ytickname=[' ',' '],$
              xticks=1,xminor=1,xtickname=[' ',' '],$
              levels=indgen(255),c_color=indgen(255),/fill,/Noerase
      plot_oi,yscale,xscale,Position=[1600,6000,2000,12000],xticks=1,$
              xtickname=[' ',' '],yticks=5,yminor=1,$
              ytickname=['10!U-5!N','10!U-4!N','10!U-3!N','10!U-2!N','10!U-1!N','1'],$
              ytitle='I/I!Lmax!N',/Noerase,/Device,color=255
      xyouts,11000,10500,inc_lbl(it),/Device,color=255
   endfor
endif

if (fix(ifig) eq 2) then begin
   !P.Charsize=1
   N = 201.
   Nth = 41.
   Nph = 40.
   Ne_i = 11
   Nspec = 101
   run_sort = 0
;   run_id = 1251
   
   if ((Nth eq 41)and(ifig eq 2.1)) then $
      inc_lbl = ['i=168!Uo!N',$
                 'i=158!Uo!N','i=152!Uo!N','i=147!Uo!N','i=142!Uo!N','i=138!Uo!N',$
                 'i=134!Uo!N','i=130!Uo!N','i=127!Uo!N','i=123!Uo!N','i=120!Uo!N',$
                 'i=117!Uo!N','i=114!Uo!N','i=111!Uo!N','i=108!Uo!N','i=105!Uo!N',$
                 'i=102!Uo!N','i=100!Uo!N','i=97!Uo!N','i=94!Uo!N','i=90!Uo!N',$
                 'i=86!Uo!N','i=83!Uo!N','i=80!Uo!N','i=78!Uo!N',$
                 'i=75!Uo!N','i=72!Uo!N','i=69!Uo!N','i=66!Uo!N','i=63!Uo!N',$
                 'i=60!Uo!N','i=57!Uo!N','i=53!Uo!N','i=50!Uo!N','i=46!Uo!N',$
                 'i=42!Uo!N','i=38!Uo!N','i=33!Uo!N','i=28!Uo!N','i=22!Uo!N',$
                 'i=12!Uo!N']
   if ((Nph eq 40)and(ifig eq 2.2)) then $
      inc_lbl = ['f=0!Uo!N','f=9!Uo!N','f=18!Uo!N','f=27!Uo!N','f=36!Uo!N',$
                 'f=45!Uo!N','f=54!Uo!N','f=63!Uo!N','f=72!Uo!N','f=81!Uo!N',$
                 'f=90!Uo!N','f=99!Uo!N','f=108!Uo!N','f=117!Uo!N','f=126!Uo!N',$
                 'f=135!Uo!N','f=144!Uo!N','f=153!Uo!N','f=162!Uo!N','f=171!Uo!N',$
                 'f=180!Uo!N','f=189!Uo!N','f=198!Uo!N','f=207!Uo!N',$
                 'f=216!Uo!N','f=225!Uo!N','f=234!Uo!N','f=243!Uo!N','f=252!Uo!N',$
                 'f=261!Uo!N','f=270!Uo!N','f=279!Uo!N','f=288!Uo!N','f=297!Uo!N',$
                 'f=306!Uo!N','f=315!Uo!N','f=324!Uo!N','f=333!Uo!N','f=342!Uo!N',$
                 'f=351!Uo!N']

   Ixy = dblarr(N,N)
   tpmov = dblarr(N,N,Nph,Nth)
   tpmovx = dblarr(N,N,Nph,Nth)
   tpmovy = dblarr(N,N,Nph,Nth)
      
   rdatafile = 'data/scat_ithp.0000.dat'
   dumpstr = string(run_id,format='(I4.4)')
   strput,rdatafile,dumpstr,15
   openr,1,rdatafile
   readf,1,tpmov,tpmovx,tpmovy
   close,1
      
   for i=0,Nth-1 do begin
      for j=0,Nph-1 do begin
         tpmov(*,*,j,i)=transpose(tpmov(*,*,j,i))
         tpmovx(*,*,j,i)=transpose(tpmovx(*,*,j,i))
         tpmovy(*,*,j,i)=transpose(tpmovy(*,*,j,i))
      endfor
   endfor

   if (ifig eq 2.1) then begin
      jph = 30
      mov = fltarr(N,N,Nth)
      movx = fltarr(N,N,Nth)
      movy = fltarr(N,N,Nth)
      mov(*,*,*)=tpmov(*,*,jph,*)
      movx(*,*,*)=tpmovx(*,*,jph,*)
      movy(*,*,*)=tpmovy(*,*,jph,*)
      Nend = Nth
   endif
   if (ifig eq 2.2) then begin
      jth = 36
      mov = fltarr(N,N,Nph)
      movx = fltarr(N,N,Nph)
      movy = fltarr(N,N,Nph)
      mov(*,*,*)=tpmov(*,*,*,jth)
      movx(*,*,*)=tpmovx(*,*,*,jth)
      movy(*,*,*)=tpmovy(*,*,*,jth)
      Nend = Nph
   endif
   sortmov = mov(sort(mov))
   movmax = max(mov)
   movmax=sortmov(1.*N*N*Nend-100.)
;      movmax=sortmov(1.*N*N*Nend-1.)
   outliers = where(mov gt movmax)
   if (outliers(0) ge 0) then begin
      movx(outliers)=movx(outliers)/mov(outliers)*movmax
      movy(outliers)=movy(outliers)/mov(outliers)*movmax
      mov(outliers)=movmax
   endif
   
   for it=Nend-1,0,-1 do begin
      y_min = 0
      if (N eq 81) then Ixy = dblarr(121,121)
      if (N eq 201) then Ixy = dblarr(301,301)
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
      data = byte(255*(alog10(Ixy/movmax+1e-4)+4.01)/4.)
      N15 = (N-1)*1.5+1
      N2 = (fix(600./N15))*N15
      enlarge = rebin(data,N2,N2)
   
      erase
      tvscl,enlarge,0,y_min,/Device
      dd = 100.
      pstep=10
      maxIxy = max(Ixy)
      for i=N15/6.,N15*(5./6.),pstep do begin
         for j=N15/3.,N15,pstep do begin
            x0 = i/(N15-1.)*12600.
            y0 = j/(N15-1.)*12600.
            dff = 100.*dd*([Xpol(i,j),Ypol(i,j),0])
            if (Ixy(i,j) gt 1d-4*maxIxy) then $
               plots,[x0-dff(0)/2.,x0+dff(0)/2.],[y0-dff(1)/2.,y0+dff(1)/2.],$
                     color=255,/Device,thick=4
         endfor
      endfor
      plots,[11000,11000+5.*dd],[12000,12000],color=255,thick=thk,/Device
      xyouts,11000,11500,'deg=5%',color=255,/Device

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
      contour,enlarge,Position=[1600,6000,2000,12000],/Device,$
              yticks=1,yminor=1,ytickname=[' ',' '],$
              xticks=1,xminor=1,xtickname=[' ',' '],$
              levels=indgen(255),c_color=indgen(255),/fill,/Noerase
      plot_oi,yscale,xscale,Position=[1600,6000,2000,12000],xticks=1,$
              xtickname=[' ',' '],yticks=5,yminor=1,$
              ytickname=['10!U-5!N','10!U-4!N','10!U-3!N','10!U-2!N','10!U-1!N','1'],$
              ytitle='I/I!Lmax!N',/Noerase,/Device,color=255
      xyouts,11000,10500,inc_lbl(it),/Device,color=255
   endfor
endif

;SORT IMAGES BY SPECTRAL BINS
if (fix(ifig) eq 3) then begin
   !P.Charsize=1
   N = 201.
   Nth = 41.
   Ne_i = 11
   
   if (Nth eq 41) then $
      inc_lbl = ['i=168!Uo!N',$
                 'i=158!Uo!N','i=152!Uo!N','i=147!Uo!N','i=142!Uo!N','i=138!Uo!N',$
                 'i=134!Uo!N','i=130!Uo!N','i=127!Uo!N','i=123!Uo!N','i=120!Uo!N',$
                 'i=117!Uo!N','i=114!Uo!N','i=111!Uo!N','i=108!Uo!N','i=105!Uo!N',$
                 'i=102!Uo!N','i=100!Uo!N','i=97!Uo!N','i=94!Uo!N','i=90!Uo!N',$
                 'i=86!Uo!N','i=83!Uo!N','i=80!Uo!N','i=78!Uo!N',$
                 'i=75!Uo!N','i=72!Uo!N','i=69!Uo!N','i=66!Uo!N','i=63!Uo!N',$
                 'i=60!Uo!N','i=57!Uo!N','i=53!Uo!N','i=50!Uo!N','i=46!Uo!N',$
                 'i=42!Uo!N','i=38!Uo!N','i=33!Uo!N','i=28!Uo!N','i=22!Uo!N',$
                 'i=12!Uo!N']
   if (Ne_i eq 11) then $
      e_lbl = ['E = 0.002 keV','E = 0.007 keV','E = 0.02 keV','E = 0.08 keV',$
               'E = 0.3 keV','E = 1.0 keV','E = 3.5 keV','E = 12 keV',$
               'E = 43 keV','E = 150 keV','E = 500 keV']
   Ixy = dblarr(N,N)
   smov = dblarr(Ne_i,N,N,Nth)
   smovx = dblarr(Ne_i,N,N,Nth)
   smovy = dblarr(Ne_i,N,N,Nth)
      
   rdatafile = 'data/scat_ispc.0000.dat'
   dumpstr = string(run_id,format='(I4.4)')
   strput,rdatafile,dumpstr,15
   openr,1,rdatafile
   readf,1,smov,smovx,smovy
   close,1
      
   for i=0,Nth-1 do begin
      for j=0,Ne_i-1 do begin
         smov(j,*,*,i)=transpose(smov(j,*,*,i))
         smovx(j,*,*,i)=transpose(smovx(j,*,*,i))
         smovy(j,*,*,i)=transpose(smovy(j,*,*,i))
      endfor
   endfor

   if (ifig eq 3.1) then begin
      je = 4
      mov = fltarr(N,N,Nth)
      movx = fltarr(N,N,Nth)
      movy = fltarr(N,N,Nth)
      mov(*,*,*)=smov(je,*,*,*)
      movx(*,*,*)=smovx(je,*,*,*)
      movy(*,*,*)=smovy(je,*,*,*)
      sortmov = mov(sort(mov))
      movmax = max(mov)
      movmax=sortmov(1.*N*N*Nth-100.)
;      movmax=sortmov(1.*N*N*Nth-1.)
      Nend = Nth
   endif
   if (ifig eq 3.2) then begin
      jth = 35
      mov = fltarr(N,N,Ne_i)
      movx = fltarr(N,N,Ne_i)
      movy = fltarr(N,N,Ne_i)
      for je=0,Ne_i-1 do begin
         mov(*,*,je)=smov(je,*,*,jth)
         movx(*,*,je)=smovx(je,*,*,jth)
         movy(*,*,je)=smovy(je,*,*,jth)
;         mov(*,*,je) = transpose(mov(*,*,je))
;         movx(*,*,je) = transpose(movx(*,*,je))
;         movy(*,*,je) = transpose(movy(*,*,je))
      endfor
      sortmov = mov(sort(mov))
      movmax = max(mov)
      movmax=sortmov(1.*N*N*Ne_i-100.)
;      movmax=sortmov(1.*N*N*Nth-1.)
      Nend = Ne_i
      inc_lbl = e_lbl
   endif
   outliers = where(mov gt movmax)
   if (outliers(0) ge 0) then begin
      movx(outliers)=movx(outliers)/mov(outliers)*movmax
      movy(outliers)=movy(outliers)/mov(outliers)*movmax
      mov(outliers)=movmax
   endif
   
   for it=Nend-1,0,-1 do begin
      y_min = 0
      if (N eq 81) then Ixy = dblarr(121,121)
      if (N eq 201) then Ixy = dblarr(301,301)
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
      data = byte(255*(alog10(Ixy/movmax+1e-4)+4.01)/4.)
      N15 = (N-1)*1.5+1
      N2 = (fix(600./N15))*N15
      enlarge = rebin(data,N2,N2)
   
      erase
      tvscl,enlarge,0,y_min,/Device
      dd = 100.
      pstep=10
      maxIxy = max(Ixy)
      for i=N15/6.,N15*(5./6.),pstep do begin
         for j=N15/3.,N15,pstep do begin
            x0 = i/(N15-1.)*12600.
            y0 = j/(N15-1.)*12600.
            dff = 100.*dd*([Xpol(i,j),Ypol(i,j),0])
;            if (Ixy(i,j) gt 1d-4*maxIxy) then $
;               plots,[x0-dff(0)/2.,x0+dff(0)/2.],[y0-dff(1)/2.,y0+dff(1)/2.],$
;                     color=255,/Device,thick=4
         endfor
      endfor
      plots,[11000,11000+5.*dd],[12000,12000],color=255,thick=thk,/Device
      xyouts,11000,11500,'deg=5%',color=255,/Device

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
      contour,enlarge,Position=[1600,6000,2000,12000],/Device,$
              yticks=1,yminor=1,ytickname=[' ',' '],$
              xticks=1,xminor=1,xtickname=[' ',' '],$
              levels=indgen(255),c_color=indgen(255),/fill,/Noerase
      plot_oi,yscale,xscale,Position=[1600,6000,2000,12000],xticks=1,$
              xtickname=[' ',' '],yticks=5,yminor=1,$
              ytickname=['10!U-5!N','10!U-4!N','10!U-3!N','10!U-2!N','10!U-1!N','1'],$
              ytitle='I/I!Lmax!N',/Noerase,/Device,color=255
      xyouts,11000,10500,inc_lbl(it),/Device,color=255
   endfor
endif

;SORT IMAGES BY PHOTON HISTORY (ISORT)
if (fix(ifig) eq 4) then begin
   !P.Charsize=1
   N = 201.
   Nth = 41.
   Ns = 6
;   run_id = 1251
   
   if (Nth eq 41) then $
      inc_lbl = ['i=168!Uo!N',$
                 'i=158!Uo!N','i=152!Uo!N','i=147!Uo!N','i=142!Uo!N','i=138!Uo!N',$
                 'i=134!Uo!N','i=130!Uo!N','i=127!Uo!N','i=123!Uo!N','i=120!Uo!N',$
                 'i=117!Uo!N','i=114!Uo!N','i=111!Uo!N','i=108!Uo!N','i=105!Uo!N',$
                 'i=102!Uo!N','i=100!Uo!N','i=97!Uo!N','i=94!Uo!N','i=90!Uo!N',$
                 'i=86!Uo!N','i=83!Uo!N','i=80!Uo!N','i=78!Uo!N',$
                 'i=75!Uo!N','i=72!Uo!N','i=69!Uo!N','i=66!Uo!N','i=63!Uo!N',$
                 'i=60!Uo!N','i=57!Uo!N','i=53!Uo!N','i=50!Uo!N','i=46!Uo!N',$
                 'i=42!Uo!N','i=38!Uo!N','i=33!Uo!N','i=28!Uo!N','i=22!Uo!N',$
                 'i=12!Uo!N']
   Ixy = dblarr(N,N)
   mov = fltarr(N,N,Nth)
   movx = fltarr(N,N,Nth)
   movy = fltarr(N,N,Nth)
   smov = dblarr(Ns,N,N,Nth)
   smovx = dblarr(Ns,N,N,Nth)
   smovy = dblarr(Ns,N,N,Nth)
      
   rdatafile = 'data/scat_isrt.0000.dat'
   dumpstr = string(run_id,format='(I4.4)')
   strput,rdatafile,dumpstr,15
   openr,1,rdatafile
   readf,1,smov,smovx,smovy
   close,1
      
   for i=0,Nth-1 do begin
      for j=0,Ns-1 do begin
         smov(j,*,*,i)=transpose(smov(j,*,*,i))
         smovx(j,*,*,i)=transpose(smovx(j,*,*,i))
         smovy(j,*,*,i)=transpose(smovy(j,*,*,i))
      endfor
   endfor

   if (ifig eq 4.1) then begin
      ji = 0
      mov = fltarr(N,N,Nth)
      movx = fltarr(N,N,Nth)
      movy = fltarr(N,N,Nth)
      mov(*,*,*)=smov(ji,*,*,*)
      movx(*,*,*)=smovx(ji,*,*,*)
      movy(*,*,*)=smovy(ji,*,*,*)
      sortmov = mov(sort(mov))
      movmax = max(mov)
      movmax=sortmov(1.*N*N*Nth-100.)
;      movmax=sortmov(1.*N*N*Nth-1.)
      Nend = Nth
   endif
   if (ifig eq 4.2) then begin
      jth = 40
      mov = fltarr(N,N,Ns)
      movx = fltarr(N,N,Ns)
      movy = fltarr(N,N,Ns)
      for j=0,Ns-1 do begin
         mov(*,*,j)=smov(j,*,*,jth)
         movx(*,*,j)=smovx(j,*,*,jth)
         movy(*,*,j)=smovy(j,*,*,jth)
      endfor
      sortmov = mov(sort(mov))
      movmax = max(mov)
      movmax=sortmov(1.*N*N*Ns-100.)
;      movmax=sortmov(1.*N*N*Nth-1.)
      Nend = Ns
   endif
   outliers = where(mov gt movmax)
   if (outliers(0) ge 0) then begin
      movx(outliers)=movx(outliers)/mov(outliers)*movmax
      movy(outliers)=movy(outliers)/mov(outliers)*movmax
      mov(outliers)=movmax
   endif
   
   for it=Nend-1,0,-1 do begin
      y_min = 0
      if (N eq 81) then Ixy = dblarr(121,121)
      if (N eq 201) then Ixy = dblarr(301,301)
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
      data = byte(255*(alog10(Ixy/movmax+1e-4)+4.01)/4.)
      N15 = (N-1)*1.5+1
      N2 = (fix(600./N15))*N15
      enlarge = rebin(data,N2,N2)
   
      erase
      tvscl,enlarge,0,y_min,/Device
      dd = 100.
      pstep=10
      maxIxy = max(Ixy)
      for i=N15/6.,N15*(5./6.),pstep do begin
         for j=N15/3.,N15,pstep do begin
            x0 = i/(N15-1.)*12600.
            y0 = j/(N15-1.)*12600.
            dff = 100.*dd*([Xpol(i,j),Ypol(i,j),0])
            if (Ixy(i,j) gt 1d-4*maxIxy) then $
               plots,[x0-dff(0)/2.,x0+dff(0)/2.],[y0-dff(1)/2.,y0+dff(1)/2.],$
                     color=255,/Device,thick=4
         endfor
      endfor
      plots,[11000,11000+5.*dd],[12000,12000],color=255,thick=thk,/Device
      xyouts,11000,11500,'deg=5%',color=255,/Device

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
      contour,enlarge,Position=[1600,6000,2000,12000],/Device,$
              yticks=1,yminor=1,ytickname=[' ',' '],$
              xticks=1,xminor=1,xtickname=[' ',' '],$
              levels=indgen(255),c_color=indgen(255),/fill,/Noerase
      plot_oi,yscale,xscale,Position=[1600,6000,2000,12000],xticks=1,$
              xtickname=[' ',' '],yticks=5,yminor=1,$
              ytickname=['10!U-5!N','10!U-4!N','10!U-3!N','10!U-2!N','10!U-1!N','1'],$
              ytitle='I/I!Lmax!N',/Noerase,/Device,color=255
      xyouts,11000,10500,inc_lbl(it),/Device,color=255
   endfor
   stop
endif

  if (fix(ifig) eq 7) then begin
      !P.charsize=1.3
      Nt = 41
      Nspec = 101
;sandwich corona, vary tau and T
      run_ids = [102,101,100,103]
      run_ids = [1250,1250,1250,1250]
      idex75 = 4
      idex60 = 10
      idex45 = 15

      spec = dblarr(Nspec,Nt)
      Qspec = dblarr(Nspec,Nt)
      Uspec = dblarr(Nspec,Nt)

      emin = 0.1
      emax = 1000
      shorte = emin*10.^(findgen(Nspec)/(Nspec-1.)*alog10(emax/emin))
      de = deriv(shorte)

      for rundex=0,0 do begin
          run_id = run_ids(rundex)
          rdatafile = 'data/scat_spec.0000.dat'
          dumpstr = string(run_id,format='(I4.4)')
          strput,rdatafile,dumpstr,15
          openr,1,rdatafile
          readf,1,spec,Qspec,Uspec
          close,1
          deg_spec = sqrt(Qspec^2+Uspec^2)/spec
          ang_spec = atan(Uspec,Qspec)/2.*!radeg

          if (rundex eq 0) then min_deg = -110
          if (rundex eq 1) then min_deg = -110
          if (rundex eq 2) then min_deg = -110
          if (rundex eq 3) then min_deg = -120
          neg_ang = where(ang_spec lt min_deg)
          if (neg_ang(0) ge 0) then ang_spec(neg_ang)=ang_spec(neg_ang)+180.
          pos_ang = where(ang_spec gt min_deg+180)
          if (pos_ang(0) ge 0) then ang_spec(pos_ang)=ang_spec(pos_ang)-180.
          
          if ((ifig eq 7.1)) then begin
              if (rundex eq 0) then begin
                  plot_oi,shorte,deg_spec(*,idex75)*100.,$
                    xrange=[0.1,20],yrange=[0,12],xstyle=1,ystyle=1,$
                    xtitle='E!Lobs!N (keV)',ytitle='polarization degree (%)',$
                    thick=thk
                  plots,[.15,.30],[1,1]*10,thick=thk,color= 60
                  plots,[.15,.30],[1,1]*9,thick=thk,color= 90
                  plots,[.15,.30],[1,1]*8,thick=thk,color= 120
                  plots,[.15,.30],[1,1]*7,thick=thk,color= 150
                  xyouts,0.32,9.9,'H/R = 0.05'
                  xyouts,0.32,8.9,'H/R = 0.1'
                  xyouts,0.32,7.9,'H/R = 0.2'
                  xyouts,0.32,6.9,'H/R = 0.4'
              endif
              oplot,shorte,deg_spec(*,idex75)*100,thick=thk,color=60+rundex*30
          endif
          
          if ((ifig eq 7.11)) then begin
              if (rundex eq 0) then begin
                  plot_oi,shorte,deg_spec(*,idex75)*100.,$
                    xrange=[0.1,20],yrange=[0,12],xstyle=1,ystyle=1,$
                    xtitle='E!Lobs!N (keV)',ytitle='polarization degree (%)',$
                    thick=thk
                  plots,[.15,.30],[1,1]*10,thick=thk,color= 60
                  plots,[.15,.30],[1,1]*9,thick=thk,color= 90
                  plots,[.15,.30],[1,1]*8,thick=thk,color= 120
                  xyouts,0.32,9.9,'!9t!3!Les!N=1.4, T!Le!N=50 keV'
                  xyouts,0.32,8.9,'!9t!3!Les!N=1.0, T!Le!N=100 keV'
                  xyouts,0.32,7.9,'!9t!3!Les!N=0.5, T!Le!N=200 keV'
              endif
              oplot,shorte,deg_spec(*,idex75)*100,thick=thk,color=60+rundex*30
          endif
          
          if (ifig eq 7.2) then begin
;              if (rundex le 1) then min_deg = -50
              neg_ang = where(ang_spec lt min_deg)
              if (neg_ang(0) ge 0) then ang_spec(neg_ang)=ang_spec(neg_ang)+180.
              pos_ang = where(ang_spec gt min_deg+180)
              if (pos_ang(0) ge 0) then ang_spec(pos_ang)=ang_spec(pos_ang)-180.
              if (rundex eq 0) then begin
                  min_deg = -120
                  plot_oi,shorte,smooth(ang_spec(*,idex75),1),$
                    xrange=[0.1,20.],yrange=[min_deg,min_deg+180],xstyle=1,ystyle=1,$
                    xtitle='E!Lobs!N (keV)',ytitle='polarization angle (deg)'
                  plots,[0.12,0.2],[1,1]*(160)+min_deg,thick=thk,color=60
                  plots,[0.12,0.2],[1,1]*(150)+min_deg,thick=thk,color=90
                  plots,[0.12,0.2],[1,1]*(140)+min_deg,thick=thk,color=120
                  plots,[0.12,0.2],[1,1]*(130)+min_deg,thick=thk,color=150
                  xyouts,0.22,158+min_deg,'H/R=0.05'
                  xyouts,0.22,148+min_deg,'H/R=0.1'
                  xyouts,0.22,138+min_deg,'H/R=0.2'
                  xyouts,0.22,128+min_deg,'H/R=0.4'
              endif
              oplot,shorte,smooth(ang_spec(*,idex75),1),thick=thk,color=60+rundex*30
          endif
          if (ifig eq 7.21) then begin
;              if (rundex le 1) then min_deg = -50
              neg_ang = where(ang_spec lt min_deg)
              if (neg_ang(0) ge 0) then ang_spec(neg_ang)=ang_spec(neg_ang)+180.
              pos_ang = where(ang_spec gt min_deg+180)
              if (pos_ang(0) ge 0) then ang_spec(pos_ang)=ang_spec(pos_ang)-180.
              if (rundex eq 0) then begin
                  min_deg = -100
                  plot_oi,shorte,smooth(ang_spec(*,idex75),1),$
                    xrange=[0.1,20.],yrange=[min_deg,min_deg+200],xstyle=1,ystyle=1,$
                    xtitle='E!Lobs!N (keV)',ytitle='polarization angle (deg)'
                  plots,[0.12,0.2],[1,1]*(50)+min_deg,thick=thk,color=60
                  plots,[0.12,0.2],[1,1]*(35)+min_deg,thick=thk,color=90
                  plots,[0.12,0.2],[1,1]*(20)+min_deg,thick=thk,color=120
                  xyouts,0.22,48+min_deg,'!9t!3!Les!N=1.4, T!Lcor!N=50 keV'
                  xyouts,0.22,33+min_deg,'!9t!3!Les!N=1.0, T!Lcor!N=100 keV'
                  xyouts,0.22,18+min_deg,'!9t!3!Les!N=0.5, T!Lcor!N=200 keV'
              endif
              oplot,shorte,smooth(ang_spec(*,idex75),1),thick=thk,color=60+rundex*30
          endif
      endfor
      !P.charsize=1.5
  endif
  if ((fix(ifig) eq 8)) then begin
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

      emin = 0.1
      emax = 1000
      shorte = emin*10.^(findgen(Nspec)/(Nspec-1.)*alog10(emax/emin))
      de = deriv(shorte)

      if (fix(ifig) eq 8) then run_ids = [143,144,145,146,147]
;      if (fix(ifig) eq 8) then run_ids = [148,149,150,151,152]
;      if (fix(ifig) eq 7) then run_ids = [113,114,115,116,117]
;      if (fix(ifig) eq 8) then run_ids = [118,119,120,121,122]
;      if (fix(ifig) eq 9) then run_ids = [143,148,189,189,189]

      for rundex=0,4 do begin
          run_id = run_ids(rundex)
          rdatafile = '~/diskspec/scat_spec.0000.dat'
          dumpstr = string(run_id,format='(I4.4)')
          strput,rdatafile,dumpstr,21
          openr,1,rdatafile
          readf,1,spec,Qspec,Uspec
          close,1
          deg_spec = sqrt(Qspec^2+Uspec^2)/spec
          ang_spec = atan(Uspec,Qspec)/2.*!radeg

          if (fix(ifig) eq 8) then min_deg = -120
          neg_ang = where(ang_spec lt min_deg)
          if (neg_ang(0) ge 0) then ang_spec(neg_ang)=ang_spec(neg_ang)+180.
          pos_ang = where(ang_spec gt min_deg+180)
          if (pos_ang(0) ge 0) then ang_spec(pos_ang)=ang_spec(pos_ang)-180.
          
          if ((ifig eq 7.1)or(ifig eq 8.1)or(ifig eq 9.1)) then begin
              if (rundex eq 0) then begin
                  plot_oi,shorte,deg_spec(*,idex75)*100.,$
                    xrange=[0.1,20.],yrange=[0,12],xstyle=1,ystyle=1,$
                    xtitle='E!Lobs!N (keV)',ytitle='polarization degree (%)',thick=thk
                  plots,[0.15,0.3],[1,1]*11,thick=thk,color= 60
                  plots,[0.15,0.3],[1,1]*10,thick=thk,color= 90
                  plots,[0.15,0.3],[1,1]*9,thick=thk,color= 120
                  plots,[0.15,0.3],[1,1]*8,thick=thk,color= 150
                  plots,[0.15,0.3],[1,1]*7,thick=thk,color= 180
                  xyouts,0.32,10.9,'a/M=0'
                  xyouts,0.32,9.9,'     0.5'
                  xyouts,0.32,8.9,'     0.9'
                  xyouts,0.32,7.9,'   0.99'
                  xyouts,0.32,6.9,' 0.998'
              endif
              oplot,shorte,deg_spec(*,idex75)*100,thick=thk,color=60+rundex*30
          endif
          
          if ((ifig eq 7.2)or(ifig eq 8.2)or(ifig eq 9.2)) then begin
              if (rundex eq 0) then begin
                  plot_oi,shorte,smooth(ang_spec(*,idex75),1),$
                    xrange=[0.1,20.],yrange=[min_deg,min_deg+180],xstyle=1,ystyle=1,$
                    xtitle='E!Lobs!N (keV)',ytitle='polarization angle (deg)'
                  plots,[0.15,0.3],[1,1]*70+min_deg,thick=thk,color=60
                  plots,[0.15,0.3],[1,1]*60+min_deg,thick=thk,color=90
                  plots,[0.15,0.3],[1,1]*50+min_deg,thick=thk,color=120
                  plots,[0.15,0.3],[1,1]*40+min_deg,thick=thk,color=150
                  plots,[0.15,0.3],[1,1]*30+min_deg,thick=thk,color=180
                  xyouts,0.32,68+min_deg,'a/M=0'
                  xyouts,0.32,58+min_deg,'     0.5'
                  xyouts,0.32,48+min_deg,'     0.9'
                  xyouts,0.32,38+min_deg,'   0.99'
                  xyouts,0.32,28+min_deg,' 0.998'
              endif
              oplot,shorte,smooth(ang_spec(*,idex75),1),thick=thk,color=60+rundex*30
          endif
      endfor
  endif

END
