PRO control_iterate2
  it_start = 201
  it_stop = 250
  it_step = 1
  Nframes = 251
  rdata = fltarr(4,Nframes)
  blnk = ' '
  openr,1,'datt/run_status2.dat'
  readf,1,blnk
  readf,1,rdata
  close,1
  runid_ = intarr(Nframes)
  Niter_ = intarr(Nframes)
  Nphot_ = intarr(Nframes)
  rms_ = fltarr(Nframes)
  runid_(*)=fix(rdata(0,*))
  Niter_(*)=fix(rdata(1,*))
  Nphot_(*)=fix(rdata(2,*))
  rms_(*)=rdata(3,*)
  rms_threshold = 0.15

  for it=it_start,it_stop,it_step do begin
      if (Nphot_(it) eq 5) then ilevel = 0
      if (Nphot_(it) eq 10) then ilevel = 1
      if (Nphot_(it) eq 20) then ilevel = 2
      if (Nphot_(it) le 20) then iterate_corona2,runid_(it),runid_(it),ilevel,rmsw
      rms_(it)=rmsw
      if (rmsw lt rms_threshold) then begin
          Niter_(it) = 0
          if (Nphot_(it) eq 20) then Nphot_(it) = 30
          if (Nphot_(it) eq 10) then Nphot_(it) = 20
          if (Nphot_(it) eq 5) then Nphot_(it) = 10
      endif
      if (rmsw ge rms_threshold) then Niter_(it) = Niter_(it)+1
  endfor
  wdata = rdata
  openw,1,'datt/run_status2.dat'
  printf,1,blnk
  for it=0,Nframes-1 do $
    printf,1,runid_(it),Niter_(it),Nphot_(it),rms_(it)
  close,1

END

