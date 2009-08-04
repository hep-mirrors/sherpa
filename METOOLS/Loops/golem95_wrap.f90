module golem95_wrap
use precision
use matrice_s, only: set_ref, s_mat
use form_factor_type, only: form_factor
use form_factor_3p, only: a32, b32
use form_factor_4p, only: a40, a44, b44, c44
use constante, only: s_null

implicit none

contains

subroutine set_s_mat(i, j, v)
   implicit none
   integer, intent(in) :: i, j
   real(ki), intent(in) :: v
   s_mat(i, j) = v
end subroutine set_s_mat

subroutine set_ref3(i, j, k)
   implicit none
   integer, intent(in) :: i, j, k
   set_ref = (/i,j,k/)
end subroutine set_ref3

subroutine set_ref4(i, j, k, l)
   implicit none
   integer, intent(in) :: i, j, k, l
   set_ref = (/i,j,k,l/)
end subroutine set_ref4

subroutine a40w(r1, r2, r3, r4, r5, r6)
  implicit none
  real(ki), intent(out) :: r1, r2, r3, r4, r5, r6
  type(form_factor) :: result
  result =  a40(s_null)
  r1=real(result%a,ki)
  r2=aimag(result%a)
  r3=real(result%b,ki)
  r4=aimag(result%b)
  r5=real(result%c,ki)
  r6=aimag(result%c)
end subroutine a40w

subroutine a32w(l1,l2, r1, r2, r3, r4, r5, r6)
  implicit none
  integer, intent(in) :: l1, l2
  real(ki), intent(out) :: r1, r2, r3, r4, r5, r6
  type(form_factor) :: result
  result =  a32(l1,l2,s_null)
  r1=real(result%a,ki)
  r2=aimag(result%a)
  r3=real(result%b,ki)
  r4=aimag(result%b)
  r5=real(result%c,ki)
  r6=aimag(result%c)
end subroutine a32w

subroutine b32w(r1, r2, r3, r4, r5, r6)
  implicit none
  real(ki), intent(out) :: r1, r2, r3, r4, r5, r6
  type(form_factor) :: result
  result =  b32(s_null)
  r1=real(result%a,ki)
  r2=aimag(result%a)
  r3=real(result%b,ki)
  r4=aimag(result%b)
  r5=real(result%c,ki)
  r6=aimag(result%c)
end subroutine b32w


subroutine a44w(l1, l2, l3, l4, r1, r2, r3, r4, r5, r6)
  implicit none
  integer, intent(in) :: l1, l2, l3, l4
  real(ki), intent(out) :: r1, r2, r3, r4, r5, r6
  type(form_factor) :: result
  result =  a44(l1,l2, l3, l4, s_null)
  r1=real(result%a,ki)
  r2=aimag(result%a)
  r3=real(result%b,ki)
  r4=aimag(result%b)
  r5=real(result%c,ki)
  r6=aimag(result%c)
end subroutine a44w

subroutine b44w(l1,l2, r1, r2, r3, r4, r5, r6)
  implicit none
  integer, intent(in) :: l1, l2
  real(ki), intent(out) :: r1, r2, r3, r4, r5, r6
  type(form_factor) :: result
  result =  b44(l1,l2,s_null)
  r1=real(result%a,ki)
  r2=aimag(result%a)
  r3=real(result%b,ki)
  r4=aimag(result%b)
  r5=real(result%c,ki)
  r6=aimag(result%c)
end subroutine b44w

subroutine c44w(r1, r2, r3, r4, r5, r6)
  implicit none
  real(ki), intent(out) :: r1, r2, r3, r4, r5, r6
  type(form_factor) :: result
  result =  c44(s_null)
  r1=real(result%a,ki)
  r2=aimag(result%a)
  r3=real(result%b,ki)
  r4=aimag(result%b)
  r5=real(result%c,ki)
  r6=aimag(result%c)
end subroutine c44w

end module golem95_wrap
