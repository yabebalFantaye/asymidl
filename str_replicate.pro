function str_replicate,str,num

newstr=''
for i=0,num-1 do begin
   newstr = newstr+str
endfor

return,newstr
end
