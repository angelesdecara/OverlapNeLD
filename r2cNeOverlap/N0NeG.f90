!Program to calculate N0 and Ne 
!from a list of pairs of loci (c and r2)

implicit none
! n_ant is the number of parental cohorts 
! c is the recombination frequency
! p is the contribution os each cohort to the newborns
! b is the fertility of individuals of each cohort
! l is proportion of surviving animals at each cohort
! ni is the number of individuals of each cohort
integer n_ant, n_antp,io
parameter (n_antp = 100)
real*8 c, p(0:n_antp), f1, b(0:n_antp)
integer i, j, k, ni(0:n_antp), npairs, long
real*8 r2_mean   !promedio de las ultimas 1000 generaciones
real*8 l(0:n_antp), n0, ne, s_inversos, harm
character*100 inputfile, outputfile, r2file

call getarg(1,inputfile)
open(10,file=inputfile)
read(10,*)outputfile
read(10,*)long
open(15,file=outputfile)
read(10,*)n_ant
l(0) = 1.0d0
do i=1, n_ant
  read(10,*)l(i)
enddo
p(0) = 0.0d0
p(n_ant) = 1.0d0
do i=1, n_ant-1
  read(10,*)p(i)
  p(n_ant) = p(n_ant) - p(i)
enddo
read(10,*)r2file
open(25,file=r2file)
!read(25,*)npairs
s_inversos = 0.0d0
if(long .eq.1) then
  write(15,'(a)')'       pair         c               r2              N0'
end if
io=0;npairs=1
do while(io.eq.0)
!do i = 1, npairs
   read(25,*,iostat=io)c, r2_mean
   if(io.ne.0)exit
   npairs=npairs+1
  call calc_N(n_ant, c, p, l, r2_mean, n0)
  if(long .eq.1) then
    write(15,'(i10,2f16.8,f16.4)')i,c,r2_mean,n0
  end if
  s_inversos = s_inversos + 1.0d0/n0
enddo
npairs=npairs-1
!write(6,*)npairs
!harmonic mean of N0's
harm = npairs / s_inversos
write(15,'(a,19x,f16.4)') 'harmonic mean N0    ', harm
write(6,'(a,19x,f16.4)') 'harmonic mean N0    ', harm
call NefromN0(n_ant, harm,l,p,ne)
stop
end program

subroutine calc_N(n_ant, c, p, l, r2, n0)
  implicit none
  integer n_antp
  parameter (n_antp = 100)
  integer n_ant, i, j,k, iter
  real*8 c, p(0:n_antp),l(0:n_antp),r2, n0
  real*8 tau(0:n_antp), beta (n_antp, n_antp), eps(n_antp)
  real*8 delta(n_antp, n_antp), zeta(n_antp,n_antp)
  real*8 rho(n_antp), alpha(n_antp), xi(n_antp)
  real*8 LL(0:n_antp, 0:n_antp)
  real*8 sum_alphas
  
  delta = 0.0d0
  do i = 1, n_ant
     delta(i,i) = 1.0d0
  enddo
  
  tau = (1.0d0 - c) * p
  beta = 0.0d0
  do i = 1, n_ant
     do k = 1, n_ant
        do j = 1, n_ant
           beta(i,k) = beta(i,k) + delta(abs(i-j),k)*p(j)
        enddo
     enddo  
  enddo

  zeta = 0.0d0
  zeta = (1.0d0 - c) * beta

  alpha = 0.0d0
  do k = 1, n_ant
     do i = 1, n_ant
        do j = 1, n_ant
           alpha(k) = alpha(k) + delta(abs(i-j),k) * p(i) * p(j)
        enddo
     enddo
  enddo

  !starting values for LL (same structure as the H matrix in Felsenstein 1971)
  LL = 0.001d0
  
  do iter = 1, 100
     !First equation
     sum_alphas = 0.0d0
     do i=1, n_ant
        sum_alphas = sum_alphas + alpha(i)* LL(0,i) + p(i)**2*LL(i,i)
     enddo
     n0 = 1.0d0 - (1.0d0 - c)**2*sum_alphas    !n0 es realmente 2*n0
     n0 = n0 / (r2 - (1.0d0 - c)**2 * sum_alphas)
     !Other equations
     do i = 1, n_ant
        LL(0,i) = tau(i)*LL(i,i)
        do j = 1, n_ant
           if (j /= i) LL(0,i) = LL(0,i) + zeta (i,j) * LL(0,j)
        enddo
        LL(0,i)=LL(0,i)/(1-zeta(i,i))
     enddo
     !account for survival
     rho=0.d0
     eps=0.d0
     do i=1,n_ant
        eps(i)= (1.d0 - 1.d0/(n0*l(i)) )/(1.d0- 1.d0/(n0*l(i-1)))        
        rho(i)=1.d0/(n0*l(i)) - eps(i)/(n0*l(i-1))
        if(i.eq.1)then
           LL(i,i)=rho(i)+r2
        else
           LL(i,i)=rho(i)+eps(i)*LL(i-1,i-1)
        endif
     enddo
  enddo
end

subroutine NefromN0(n_ant, n0,l,p, ne)
implicit none
integer n_antp
parameter (n_antp = 100)
integer n_ant, i,j, ii, n
real*8 n0
integer niter
real*8 h(0:n_antp,0:n_antp), h2(0:n_antp,0:n_antp), p(0:n_antp)
real*8 diff, ne, plot(10000,3)
real*8 ni(0:n_antp), temp
real*8 ne_aprox
real*8 l(0:n_antp),s(0:n_antp),d(0:n_antp),v(0:n_antp),t, b(0:n_antp)

niter = 100

!cohort sizes
do i = 0, n_ant
  ni(i) = n0 * l(i)
enddo

!fertilities
do i = 0, n_ant
  b(i) = p(i) / l(i)
enddo

!Formula 10 Felsenstein 1971
t = 0.0d0
do i = 1, n_ant
  !reproductive value
  v(i) = 0.0d0
  do j = i, n_ant
    v(i) = v(i) + l(j)*b(j)
  enddo
  v(i) = v(i) / l(i)
  t = t + (i+1) * l(i)*b(i)
enddo
write(15,'(a,19x,f16.4)')'Generation interval ',t
write(6,'(a,19x,f16.4)')'Generation interval ',t
! t is the generation interval
do i=0, n_ant-1
  s(i) = l(i+1)/l(i)
  d(i) = 1.0d0 - s(i)
enddo

ne = 0.0d0
do i = 0, n_ant-1 
  ne = ne + l(i)*s(i)*d(i)*(v(i+1)**2)
enddo
ne = (ni(0)*t)/(1.0d0+ne)
!write(15,'(a,f16.4)')'Aprox. Ne from Felsenstein (formula 10)', ne
ne_aprox = ne

h  = 1.0d0
!exact calculation with formula 6 of Felsenstein 1971
do ii = 1, niter
  h2(0,0) = 0.0d0
  do i = 1, n_ant
    do j = 1, n_ant
      h2(0,0) = h2(0,0) + p(i)*p(j) * h(i,j) 
    enddo
  enddo
  h2(0,0) = (1.0d0 - 1.0d0/ni(0)) * h2(0,0)
  do i = 1, n_ant
    h2(0,i) = 0.0d0
    do j = 1, n_ant
      h2(0,i) = h2(0,i) + p(j) * h(j,i-1)
    enddo
    h2(i,0) = h2(0,i)
  enddo  
  do i = 1, n_ant
    do j = 1, n_ant
      if (i == j) then
        h2(i,j) = ((1.0d0-1.0d0/ni(i))/(1.0d0-1.0d0/ni(i-1)))*h(i-1,j-1)
      else
        h2(i,j) = h(i-1,j-1)
      end if
    enddo
  enddo
  diff = h2(0,0)/h(0,0)
  ne = 1.0d0/(t*(1.0d0-diff))
  h = h2
  plot(ii,1) = h(0,0)
enddo
write(15,'(a,f16.4)')'Exact  Ne from Felsenstein (formula 6) ', ne
write(6,'(a,f16.4)')'Exact  Ne from Felsenstein (formula 6) ', ne

return
end

