pro sky_angles,theta,phi,ecp=ecp,dipole=dipole,deg=deg

  VP= [3.8131929]*!rtod
  VT = [1.9864008]*!rtod

  ;;year = [4,5,6,5,6,7,7]      ;[2004,2007,2009,2009,2013,2013,2013,2013,2013]
  ;;legend=['WMAP1-PA','WMAP5-DP','WMAP9-PA','Planck-DP','Planck-PA','WMAP9-VA']
  WMT=[100,112,117,105,92,VT[*]]*!dtor
  WMP=[237,224,227,227,196,VP[*]]*!dtor

;;CMB dipole direction
;;http://www.sciencedirect.com/science/article/pii/S1387647306001990
     ndipole = [90d0-41.75,263.85]         ;*!dtor
     thph = make_thph_mat(90d0-41.75,263.85,/deg,/outdeg)

     ndipole = [thph[0,0],thph[0,1]]
     sdipole = [thph[1,0],thph[1,1]]

;;Eclipticl poles
     sep=[90+29.81, 276.38]
     nep=[90-29.81, 96.38]

if keyword_set(deg) then begin
   print, 'sep=',sep
   print, 'nep=',nep
   print, 'ndipole',ndipole
   print, 'sdipole',sdipole

   theta = [nep[0],ndipole[0],92.,VT[0]]
   phi = [nep[1],ndipole[1],196.,VP[0]]
endif else begin
   print, 'sep=',sep*!dtor
   print, 'nep=',nep*!dtor
   print, 'ndipole',ndipole*!dtor
   print, 'sdipole',sdipole*!dtor

   theta = [nep[0],ndipole[0],92.,VT[0]]*!dtor
   phi = [nep[1],ndipole[1],196.,VP[0]]*!dtor
endelse

end
