implicit none
integer nn, i,j, ii, n
parameter (nn = 6)
integer niter
real*8 h(nn,nn), h2(nn,nn), p(nn), diff, ne, plot(10000,3)
real*8 ni(nn), temp

!para el calculo aproximado
real*8 l(nn),s(nn),d(nn),v(nn),t, b(nn)

niter = 100
ni = 1418.88 !33.543!51.615d0   !todos los n1 son igual a n de momento

!parametros de sparow de waples y yokota
l(1) = 1.000d0
l(2) = 0.180d0
l(3) = 0.095d0
l(4) = 0.051d0      
l(5) = 0.027d0
l(6) = 0.014d0

b(1) = 0.000d0
b(2) = 2.546d0
b(3) = 2.754d0
b(4) = 2.921d0      
b(5) = 3.130d0
b(6) = 3.339d0

!calculo de los tamaños de cada cohorte
temp = ni(1)
do i = 1, nn
  ni(i) = temp * l(i)
print *, i, ni(i)
enddo

!calculo de las contribuciones de cada cohorte
temp = 0.0d0
do i = 1, nn
  p(i) = l(i) * b(i)
  temp = temp + p(i)
enddo
do i = 1, nn   !cuadrar a uno
  p(i) = p(i) / temp
enddo

!calculo aproximado

t = 0.0d0
do i = 1, nn
  !valor reproductivo
  v(i) = 0.0d0
  do j = i, nn
    v(i) = v(i) + l(j)*b(j)
  enddo
  v(i) = v(i) / l(i)
  t = t + i * l(i)*b(i)
enddo
print *, 't ', t
do i=1,6
  print *,'v ',i,v(i),p(i)
enddo

do i=1, nn-1
  s(i) = l(i+1)/l(i)
  d(i) = 1.0d0 - s(i)
enddo

ne = 0.0d0
do i = 1, nn-1 
  ne = ne + l(i)*s(i)*d(i)*(v(i+1)**2)
enddo
ne = (ni(1)*t)/(1.0d0+ne)
print *, 'aprox ',ne


h  = 1.0d0
h2 = 0.0d0
do i = 1, nn
  h2(i,i) = 1.0d0 !no entiendo porque h2 debe ser la identidad. Luego se machaca
enddo

!solapantes
do ii = 1, niter
  h2(1,1) = 0.0d0
  do i = 1, nn
    do j = 1, nn
      h2(1,1) = h2(1,1) + p(i)*p(j) * h(i,j) 
    enddo
  enddo
  h2(1,1) = (1.0d0 - 1.0d0/ni(1)) * h2(1,1)
  do i = 2, nn
    h2(1,i) = 0.0d0
    do j = 1, nn
      h2(1,i) = h2(1,i) + p(j) * h(j,i-1)
    enddo
    h2(i,1) = h2(1,i)
  enddo  
  !ecuaciones 4 y 5 del paper 
  do i = 2, nn
    do j = 2, nn
      if (i == j) then
        h2(i,j) = ((1.0d0-1.0d0/ni(i))/(1.0d0-1.0d0/ni(i-1)))*h(i-1,j-1)
      else
        h2(i,j) = h(i-1,j-1)
      end if
    enddo
  enddo
  diff = h2(1,1)/h(1,1)
  ne = 1.0d0/(t*(1.0d0-diff))
  if (mod(ii,10)==0)print *,ii, ne
  h = h2
  plot(ii,1) = h(1,1)
enddo
print *,'hola ',ne

!do i=1,niter
!  write(88,'(i10,4f12.5)')i,(plot(i,j),j=1,1)
!enddo

stop
end





