
!=================================================================
function multi(ndim,frecs, x1)
!extrae al azar un haplotipo de una multinomial de 4 frecs haplotipicas
integer ndim,i
real*8 frecs(ndim), x1, u, unif
integer multi
u = unif(x1)
do i=1,ndim-1
   if (u .lt. sum(frecs(1:i))) then
      multi = i
      return
   end if
enddo
multi = ndim
return
end
!=================================================================
function lunif(x1,m)
real*8 x1,unif
integer lunif, m
x1 = mod (16807.0d0 * x1, 2147483647.0d0)
unif = x1/2147483647.0d0
lunif = unif * m + 1
return
end function lunif
!=================================================================
function unif(x1)
real*8 x1,unif
x1 = mod (16807.0d0 * x1, 2147483647.0d0)
unif = x1/2147483647.0d0
return
end function unif
!=================================================================
subroutine calc_frecs(frecs,p1,p2,d)
! calcula frecuencias haplotipicas en funcion de p1, p2 y d
real*8 frecs(4), p1, p2, d
! frecs (1) es la frecuencia de AB
! frecs (2) es la frecuencia de Ab
! frecs (3) es la frecuencia de aB
! frecs (4) es la frecuencia de ab
frecs (1) = p1 * p2 + d
frecs (2) = p1 * (1.0d0 - p2) - d
frecs (3) = (1.0d0 - p1) * p2 - d
frecs (4) = (1.0d0 - p1) * (1.0d0 - p2) + d
if (abs(1.0d0 - frecs(1) - frecs(2)- frecs(3) - frecs(4)).gt.0.000001)then
   print *, 'error en calc_frecs'
   stop
end if
return
end
!=================================================================
subroutine calc_d(frecs,p1,p2,d)
! calcula frecuencias alelicas y d en funcion de frecuencias haplotipicas
real*8 frecs(4), p1, p2, d
! frecs (1) es la frecuencia de AB
! frecs (2) es la frecuencia de Ab
! frecs (3) es la frecuencia de aB
! frecs (4) es la frecuencia de ab
p1 = frecs(1) + frecs(2)
p2 = frecs(1) + frecs(3)
d = frecs(1) - p1 * p2
if (d .lt. 0.0d0 .or. d .gt. 1.0d0)then
   !print *, 'error en calc_d ', d
   !stop
end if
return
end
!==================================================================================
subroutine multidev(ndim,frecs,n_par,n,x1)
! muestreo de multinomial
real*8 frecs(ndim), x1, unif, bnldev
integer n_par(ndim), n, m
real*8 ft(ndim)

n_par=0
do j=1,n
   i=multi(ndim,frecs,x1)
   n_par(i)=n_par(i)+1
enddo

return

m = n
do i = 1, ndim
  ft(i) = frecs(i)
enddo

do i = 1, ndim -1
  n_par(i) = bnldev(ft(i),m,x1)
  m = m - n_par(i)
  do j = i+1, ndim
    ft(j) = ft(j) / (1.0d0 - ft(i))
  enddo
enddo
n_par(ndim) = m

!n_par(1) = bnldev(ft(1),m,x1)
!do i = 2, 3
!  m = m - n_par(i-1)
!  do j = i, 4
!    ft(j) = ft(j) / (1.0d0 - ft(i-1))
!  enddo
!  n_par(i) = bnldev(ft(i),m,x1)
!enddo
!n_par(4) = n - n_par(1) - n_par(2) - n_par(3)
return
end
!==================================================================================
function bnldev(pp,n,x1)
INTEGER n
real*8 x1, unif
REAL*8 bnldev,pp,PI
INTEGER j,nold
REAL*8 am,em,en,g,oldg,p,pc,pclog,plog,pold,sq,t,y,gammln, alngam
SAVE nold,pold,pc,plog,pclog,en,oldg
DATA nold /-1/, pold /-1./
PI = acos(-1.0d0)
if(pp .le. 0.5d0)then
    p = pp
  else
    p = 1.0d0-pp
endif
am = n * p
if (n .lt. 25)then
  bnldev = 0.0d0
  do j = 1, n
    if (unif(x1) .lt. p) bnldev = bnldev + 1.0d0
  enddo
else if (am .lt. 1.0d0) then
  g = exp(-am)
  t=1.0d0
  do j = 0, n
    t = t * unif(x1)
    if (t.lt.g) go to 1
  enddo
  j = n
1 bnldev = j
else
  if (n.ne.nold) then
    en = n
    oldg = alngam(en+1.0d0)
    nold = n
  endif
  if (p .ne. pold) then
  pc = 1.0d0 - p
  plog = log(p)
  pclog = log(pc)
  pold = p
endif
sq = sqrt(2.0d0 * am * pc)
2 y = tan(PI * unif(x1))
em = sq * y + am
if (em .lt. 0.0d0 .or. em .ge. en+1.0d0) goto 2
em = int(em)
t = 1.2d0 * sq * (1.0d0 + y**2) * exp(oldg - alngam(em + 1.0d0) &
  - alngam (en - em + 1.0d0) + em * plog + (en - em) * pclog)
if (unif(x1) .gt. t) goto 2
  bnldev = em
endif
if (p .ne. pp) bnldev = n - bnldev
return
END

!========================================================
function gammln(xx)
real*8 gammln,xx
integer j
real*8 ser,stp,tmp,x,y,cof(6)
SAVE cof,stp
DATA cof,stp/76.18009172947146d0,-86.50532032941677d0, 24.01409824083091d0, &
             -1.231739572450155d0,.1208650973866179d-2,-.5395239384953d-5,  &
              2.5066282746310005d0/
x = xx
y = x
tmp = x + 5.5d0
tmp = (x + 0.5d0) * log(tmp) - tmp
ser = 1.000000000190015d0
do j=1,6
  y = y + 1.d0
  ser=ser+cof(j)/y
enddo
gammln = tmp + log(stp * ser/x)
return
END




FUNCTION alngam(xvalue) RESULT(fn_val)

!     ALGORITHM AS245  APPL. STATIST. (1989) VOL. 38, NO. 2

!     Calculation of the logarithm of the gamma function

! ELF90-compatible version by Alan Miller
! Latest revision - 29 November 1997

! N.B. Argument IFAULT has been removed

IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(15, 100)
REAL (dp), INTENT(IN) :: xvalue
REAL (dp)             :: fn_val

! Local variables
REAL (dp) :: x, x1, x2, y

!     Coefficients of rational functions

REAL (dp), PARAMETER :: r1(9) = (/ -2.66685511495_dp, -24.4387534237_dp,  &
                                   -21.9698958928_dp,  11.1667541262_dp,  &
                                    3.13060547623_dp,  0.607771387771_dp, &
                                    11.9400905721_dp,  31.4690115749_dp,  &
                                    15.2346874070_dp /)
REAL (dp), PARAMETER :: r2(9) = (/ -78.3359299449_dp, -142.046296688_dp,  &
                                    137.519416416_dp,  78.6994924154_dp,  &
                                    4.16438922228_dp,  47.0668766060_dp,  &
                                    313.399215894_dp,  263.505074721_dp,  &
                                    43.3400022514_dp /)
REAL (dp), PARAMETER :: r3(9) = (/ -2.12159572323E5_dp,  2.30661510616E5_dp,  &
                                    2.74647644705E4_dp, -4.02621119975E4_dp,  &
                                   -2.29660729780E3_dp, -1.16328495004E5_dp,  &
                                   -1.46025937511E5_dp, -2.42357409629E4_dp,  &
                                   -5.70691009324E2_dp /)
REAL (dp), PARAMETER :: r4(5) = (/ 0.279195317918525_dp, 0.4917317610505968_dp, &
                                   0.0692910599291889_dp, 3.350343815022304_dp, &
                                   6.012459259764103_dp /)

!     Fixed constants

REAL (dp), PARAMETER :: alr2pi = 0.918938533204673_dp, four = 4._dp,  &
                        half = 0.5_dp, one = 1._dp, onep5 = 1.5_dp,   &
                        twelve = 12._dp, zero = 0._dp

!     Machine-dependant constants.
!     A table of values is given at the top of page 399 of the paper.
!     These values are for the IEEE double-precision format for which
!     B = 2, t = 53 and U = 1023 in the notation of the paper.

REAL (dp), PARAMETER :: xlge = 5.10E6_dp, xlgst = HUGE(1.0_dp)

x = xvalue
fn_val = zero

!     Test for valid function argument

IF (x >= xlgst) THEN
  WRITE(*, *) 'AS 245: Argument x too large'
  RETURN
END IF
IF (x <= zero) THEN
  WRITE(*, *) 'AS 245: Argument x <= 0'
  RETURN
END IF

!     Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined

IF (x < onep5) THEN
  IF (x < half) THEN
    fn_val = -LOG(x)
    y = x + one
    
!     Test whether X < machine epsilon
    
    IF (y == one) RETURN
  ELSE
    fn_val = zero
    y = x
    x = (x - half) - half
  END IF
  fn_val = fn_val + x * ((((r1(5)*y + r1(4))*y + r1(3))*y + r1(2))*y + r1(1)) / &
                    ((((y + r1(9))*y + r1(8))*y+ r1(7))*y + r1(6))
  RETURN
END IF

!     Calculation for 1.5 <= X < 4.0

IF (x < four) THEN
  y = (x - one) - one
  fn_val = y * ((((r2(5)*x + r2(4))*x + r2(3))*x + r2(2))*x + r2(1)) /  &
               ((((x + r2(9))*x + r2(8))*x + r2(7))*x+ r2(6))
  RETURN
END IF

!     Calculation for 4.0 <= X < 12.0

IF (x < twelve) THEN
  fn_val = ((((r3(5)*x + r3(4))*x + r3(3))*x + r3(2))*x + r3(1)) /  &
           ((((x + r3(9))*x + r3(8))*x + r3(7))*x + r3(6))
  RETURN
END IF

!     Calculation for X >= 12.0

y = LOG(x)
fn_val = x * (y - one) - half * y + alr2pi
IF (x > xlge) RETURN
x1 = one / x
x2 = x1 * x1
fn_val = fn_val + x1 * ((r4(3)*x2 + r4(2))*x2 + r4(1)) /  &
         ((x2 + r4(5))*x2 + r4(4))
RETURN
END FUNCTION alngam
