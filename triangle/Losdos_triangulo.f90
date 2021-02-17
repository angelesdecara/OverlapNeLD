
implicit none
integer nn, i, j, ii, iii, i1, n, jj, kk, io
parameter (nn = 3)
integer ni(0:nn)
real*8 c, c1, c2, h(0:nn,0:nn), h2(0:nn,0:nn), p(0:nn), ratio
real*8  Q(0:nn,0:nn), Q2(0:nn,0:nn)
real*8 n1(0:nn), n2(0:nn), Qtot
real*8 b, suma, pepe, beta
integer k
integer ncasos 

ncasos = 20

open (26,file='salida50_caso1.txt')
!choose recombination frequency 
c = 0.31606d0
!c = 0.16484d0
!c = 0.0906d0
!c = 0.0475813d0

c1 = (1.0d0 - c)
c2 = c1*c1

n = 1000
ni(0) = n * 1.0d0
ni(1) = n * 1.0d0 
ni(2) = n * 0.8d0 
ni(3) = n * 0.5d0 

do i = 0, nn
  n1(i) = 1.0d0 / (2.0d0*ni(i))
  n2(i) = 1.0d0 - n1(i)
enddo

p(1) = 0.3d0
p(2) = 0.5d0
p(3) = 0.2d0

do iii = 0, ncasos
do jj = 0, ncasos

kk = ncasos - iii - jj
if (iii+jj > ncasos) goto 3  
p(1) = dble(iii)/ncasos
p(2) = dble(jj)/ncasos
p(3) = 1.0d0 - p(1) - p(2)

!Starting varues for Q
Q  = 0.000d0

do i = 0, nn
  do j = 0, nn
    if (i/=j) then
      h(i,j) = Q(i,j)
    else
      h(i,i) = n1(i) + n2(i) * Q(i,i)
    end if
  enddo
enddo

! for each unit time
do ii = 1, 500
  h2 = 0.0d0
  Q2 = 0.0d0
  !Calculate L00
  do i = 1, nn
    do j = 1, nn
      h2(0,0) = h2(0,0) +  p(i)*p(j) * h(i,j)
    enddo
  enddo
  h2(0,0)=n1(0)+h2(0,0)*c2*n2(0)

  !Calculate L0i
  do i = 1, nn
    h2(0,i)=0.d0
    do j = 1, nn
        h2(0,i) = h2(0,i) + c1* p(j) * h(i,j) 
    enddo
    h2(i,0) = h2(0,i)
  enddo
  !calculate LAA
  do i = nn, 1, -1 
    do j = nn, 1, -1
      if (i /= j) then
        h2(i,j) = h(i-1,j-1)
      else
        h2(i,i) = (h(i-1,i-1)-n1(i-1)) * (n2(i)/n2(i-1)) +n1(i)  
      end if        
    enddo
  enddo 
  !calculate Q
  !Calculate Q00
  do i = 1, nn
    do j = 1, nn
      if (i==j)then
        Q2(0,0) = Q2(0,0) + p(i)*p(i) * (n1(i)+n2(i)*Q(i,i)) 
      else
        Q2(0,0) = Q2(0,0) + p(i)*p(j) * Q(i,j) 
      end if
    enddo
  enddo
  Q2(0,0) = c2 * Q2(0,0)
  !Calculate Q0i
  do i = 1, nn
    do j = 1, nn
      if (i==j) then
        Q2(0,i) = Q2(0,i) + p(i) * (n1(i)+n2(i)*Q(i,i)) 
      else
        Q2(0,i) = Q2(0,i) + p(j) * Q(i,j)
      end if
    enddo
    Q2(0,i) = c1 * Q2(0,i)
    Q2(i,0) = Q2(0,i)
  enddo
  !calculate QAA
  do i =  nn, 1, -1
    do j = nn, 1, -1
      Q2(i,j) = Q(i-1,j-1)
    enddo
  enddo 
  ! writing Q(0,0) y L(0,0)
  write (26,'(i10,4f12.8)')ii,Q(0,0),h(0,0)
  h = h2
  Q = Q2
  !escribir las dos matrices
  if (mod(ii,1)==0) then 
    !alert, using h2 to compare Q with L
    do i = 0, nn
      do j = 0, nn
        if (i/=j) then
          h2(i,j) = Q(i,j)
        else
          h2(i,i) = n1(i) + n2(i) * Q(i,i)
        end if
      enddo
    enddo
    !print *,ii,' -------- This is Q transformed in L -----------'
    do i=0,nn
    !  print '(22f20.10)',(h2(i,j),j=0,nn)
    enddo
  end if
  !pause
 
enddo

print '(3i6,4f12.6)',iii,jj,kk,(p(i1),i1=1,3),h2(0,0)
write (62,'(3i6,4f12.6)')iii,jj,kk,(p(i1),i1=1,3),h2(0,0)


3 continue
enddo
enddo

print *,'  '
do i=0,nn
   print '(22f12.5)',(h2(i,j),j=0,nn)
enddo

print *,h2(0,0)

!checking

pepe = h2(0,0)

do k=1, nn
  do i = 1, nn
    do j = 1, nn
      !if (i /= j) then
        if (abs(i-j) == k) then
          pepe = pepe - n2(0) * c2 * p(i) * p(j) * h2 (k, 0)
        end if
      !end if
    enddo
  enddo
enddo
do k = 1, nn  
  pepe = pepe - n2(0) * c2 * p(k) * p(k) * h2 (k, k)
enddo

print *, ' '
print '(a30,f14.12)', 'first equation left = ',pepe


pepe = n1(0)

print '(a30,f14.12)', 'first equation right = ',pepe


do i = 1, nn
  pepe = 1.0d0*h2(i,0)
  do k = 1, nn
    beta = 0.0d0
    do j = 1, nn
      if (abs(i-j) == k)then ! .and. i /= j) then
        beta = beta + p(j) 
      end if
    enddo
    pepe = pepe - beta * c1 * h2 (k,0)
  enddo
  !termino Lii
  pepe = pepe -  c1 * p(i) * h2 (i, i)
  print '(a30,f14.12)', 'other equations left = ',pepe
  !segundas ecuaciones derechas
  pepe = 0.0d0
  print '(a30,f14.12)', 'other equations right = ',pepe
  print *,'----------'
enddo
print '(a,f5.3)', 'value of Q_eq_00 ',Q(0,0)*1000
print '(a,f5.3)', 'value of L_eq_00 ',h(0,0)*1000
stop
end





