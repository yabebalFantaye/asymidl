pro simconflevel, mat, y1,y2, cl=cl,pvalue=pvalue,dim=dim,mmat=mmat

  if not keyword_set(pvalue) then pvalue=[0., 0.01, 0.05,0.32, 0.68, 0.95, 0.99, 1.]


mat = reform(mat)

s = size(mat)

npval=n_elements(pvalue)

case s[0] of
   1: begin
      p = percentiles(mat,value=pvalue)
      mmat=mean(mat)
   end
   2: begin
      n = n_elements(mat[*,0])
      p=dblarr(n,npval)
      mmat=dblarr(n)
      for i=0,n-1 do begin
         p[i,*] = percentiles(reform(mat[i,*]),value=pvalue)
         mmat[i]=mean(reform(mat[i,*]))
      endfor
   end
   ELSE: BEGIN
      PRINT, 'morethan two dimensions not implemented'
   END
endcase

cl=p

;help, p
;print,'ind lower ',0,npval/2-1
;print,'ind higher ',npval/2,npval-1

ia=0
ib=npval/2-1
;print, 'npval, ind lower',npval,ia, ib
y1=p[*,ia:ib]

ia=npval/2
ib=npval-1
;print, 'ind higher',npval,ia, ib
y2=p[*,ib:ia:-1]

end
   
