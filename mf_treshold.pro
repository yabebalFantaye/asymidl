pro mf_treshold, phiu, sig=sig,nvec=nvec,bound=bound

  if not keyword_set(sig) then sig=1.
  if not keyword_set(nvec) then nvec=200
  if not keyword_set(bound) then bound=4.
 
  phiu=2*sig*bound*(findgen(nvec)/float(nvec-1) - 0.5)

end
