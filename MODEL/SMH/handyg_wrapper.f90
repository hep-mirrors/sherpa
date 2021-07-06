function gpl1(w1,x)
  use handyG
  complex(kind=prec) :: gpl1, x, w1
  call clearcache
  gpl1 = G([w1], x)
  return
end function gpl1

function gpl2(w1,w2,x)
  use handyG
  complex(kind=prec) :: gpl2, x, w1, w2
  call clearcache
  gpl2 = G([w1,w2], x)
  return
end function gpl2

function gpl3(w1,w2,w3,x)
  use handyG
  complex(kind=prec) :: gpl3, x, w1, w2, w3
  call clearcache
  gpl3 = G([w1,w2,w3], x)
  return
end function gpl3

