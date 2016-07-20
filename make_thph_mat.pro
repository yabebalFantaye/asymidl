function make_thph_mat,theta,phi,deg=deg,outdeg=outdeg

  thph = dblarr(2,2)

  if keyword_set(deg) then begin
     theta = theta*!dtor
     phi = phi*!dtor
  endif


  thph(0,0) = theta*!rtod
  thph(0,1) = phi*!rtod

  
  thph(1,0) = (!pi - theta)*!rtod
  
  x = !pi + phi
  if x gt 2d0*!pi then x = x - 2*!pi
  thph(1,1) = x*!rtod
  
  
if not keyword_set(outdeg) then  thph = thph*!dtor

  
return,thph

end
