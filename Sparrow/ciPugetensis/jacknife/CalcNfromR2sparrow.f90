!modified to read from a file with nind pops

  !este es una mezcla de Sved_definitivo y FelsconQNew
!se van a poner ambos metodos juntos para ver si se puede
!verificar la equivalencia con la formulita que
!relaciona Q con L

!Ambos metodos encajan bien, pero ahora no encajan 
!porque habia un error en el manejo de los n(i)
!esta cosa se pone mas abajo en un comentario
  !
  ! parameters for the white crowned sparrow from Waples Yokpta
  ! r2 from Lipshutz et al
!  
implicit none
integer nn, i, j, ii, n, icase, npop,nf
parameter (nn = 5,npop=102,nf=3)
integer ni(0:nn)
real*8 c, c1, c2, h(0:nn,0:nn), h2(0:nn,0:nn), p(nn), ratio
real*8  Q(0:nn,0:nn), Q2(0:nn,0:nn), super(0:nn)
! h y h2 es para iterar L, Q y Q2 es para iterar Q
real*8 n1(0:nn), n2(0:nn), Qtot
real*8 b, suma, pepe, beta,n0
real*8 or2(npop,nf),er2(npop,nf),ss(npop,nf) !ss is not sample size but ne as estimated by NeEstimator
real*8 ne5,ne1,ne0,varne5,varne1,varne0
real*8 nest5,nest1,nest0 ! means of nestimator estimates
real*8 varnest5,varnest1,varnest0
real*8 nen0(nf),varnen0(nf),ne
integer k,if
character(len=25)::stuff

nest5=0.d0;nest1=0.d0;nest0=0.d0
varnest5=0.d0;varnest1=0.d0;varnest0=0.d0
open(22,file='r2neNeEst.txt')
do i=1,npop
   read(22,'(a25)',advance='no')stuff
   read(22,*)or2(i,1:nf)
!   write(6,*)stuff,or2(i,1:nf)
   read(22,'(a25)',advance='no')stuff
   read(22,*)er2(i,1:nf)
!   read(22,'(a25,3(f12.7,X))')stuff,er2(i,1:nf)
!   write(6,*)stuff,er2(i,1:nf)
!   pause
   read(22,'(a25)',advance='no')stuff
   read(22,*)ss(i,1:nf)

!   read(22,'(a25,3(f12.7,X))')stuff,ss(i,1:nf)
!   write(6,*)stuff
   if(i.lt.npop)read(22,*)
   nest5=nest5+ss(i,1)
   varnest5=varnest5+(ss(i,1))**2
   nest1=nest1+ss(i,2)
   varnest1=varnest1+(ss(i,2))**2
   nest0=nest0+ss(i,3)
   varnest0=varnest0+(ss(i,3))**2
   write(6,*)ss(i,1:nf)
   !pause
enddo
nest5=nest5/dble(npop)
nest1=nest1/dble(npop)
nest0=nest0/dble(npop)
varnest5=varnest5/dble(npop)-nest5**2
varnest1=varnest1/dble(npop)-nest1**2
varnest0=varnest0/dble(npop)-nest0**2

open (26,file='ciPugetensis.txt')
c = 0.5d0

c1 = (1.0d0 - c)
c2 = c1*c1

n = 1000

p(1)=0.458217d0
p(2)=0.261594d0
p(3)=0.148951d0
p(4)=0.0844984d0
p(5)=1.d0-sum(p(1:4))
super(1)=0.18d0
super(2)=0.095d0
super(3)=0.051d0
super(4)=0.027d0
super(5)=0.014d0
super(0)=1.d0
do i = 0, nn
   n1(i) = 1.0d0 / (2.0d0*ni(i))
   n2(i) = 1.0d0 - n1(i)
enddo
   
!valores iniciales para iterar sobre Q
Q  = 0.000d0

!valores inicales para iterar sobre L
!para partir del mismo sitio se transforma

do i = 0, nn
  do j = 0, nn
    if (i/=j) then
      h(i,j) = Q(i,j)
    else
      h(i,i) = n1(i) + n2(i) * Q(i,i)
    end if
  enddo
enddo

! unidad de tiempo
do ii = 1, 500
  !borrar nueva h y Q
  h2 = 0.0d0
  Q2 = 0.0d0
  !ITERACION DE L
  !Calcular L00 , fichero LinkageDisequilibrium
  do i = 1, nn
    do j = 1, nn
      h2(0,0) = h2(0,0) +  p(i)*p(j) * h(i,j)
    enddo
  enddo
  h2(0,0)=n1(0)+h2(0,0)*c2*n2(0)

  !Calcular L0i
  do i = 1, nn
    h2(0,i)=0.d0
    do j = 1, nn
        h2(0,i) = h2(0,i) + c1* p(j) * h(i,j) 
    enddo
    h2(i,0) = h2(0,i)
  enddo
  !calcular LAA
  do i = nn, 1, -1 
    do j = nn, 1, -1
      if (i /= j) then
        h2(i,j) = h(i-1,j-1)
      else
        h2(i,i) = (h(i-1,i-1)-n1(i-1)) * (n2(i)/n2(i-1)) +n1(i)  
!        h2(i,i) = (h(i-1,i-1)-n1(i)) * (n2(i)/n2(i-1)) +n1(i) 
!        esta es la linea que estaba mal, ponia n1(i) en lugar de
!        n1(i-1) que es lo que debe ser y lo que coincide con Sved
!        Pero ahora no encaja con la simulacion cuando hay n(i) 
!        diferentes en cada cohorte 
      end if        
   enddo
enddo
  !ITERACION DE Q
  !Calcular Q00
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
  !Calcular Q0i
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
  !calcular QAA
  do i =  nn, 1, -1
    do j = nn, 1, -1
      Q2(i,j) = Q(i-1,j-1)
    enddo
  enddo 
  ! escribir en fichero Q(0,0) y L(0,0)
  write (26,'(i10,4f12.8)')ii,Q(0,0),h(0,0)
  h = h2
  Q = Q2
  !escribir las dos matrices
  if (mod(ii,1)==0) then 
!    print *,ii,' -------- Esta es L -----------'
 !   do i=0,nn
 !     print '(22f20.10)',(h(i,j),j=0,nn)
 !   enddo
  !  print *,ii,' -------- Esta es Q -----------'
   ! do i=0,nn
   !   print '(22f20.10)',(Q(i,j),j=0,nn)
   ! enddo
    !ojo, usamos h2 para comparar con Q, luego se borra
    do i = 0, nn
      do j = 0, nn
        if (i/=j) then
          h2(i,j) = Q(i,j)
        else
          h2(i,i) = n1(i) + n2(i) * Q(i,i)
        end if
      enddo
    enddo
!    print *,ii,' -------- Esta es la Q transformada en L -----------'
!    do i=0,nn
!      print '(22f20.10)',(h2(i,j),j=0,nn)
!    enddo
  end if
  !pause
  

enddo


!primera ecuacion izquierda

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
do k = 1, nn  !terminos gamma
  pepe = pepe - n2(0) * c2 * p(k) * p(k) * h2 (k, k)
enddo

!print *, ' '
!print '(a30,f14.12)', 'primera ecuacion izquierda = ',pepe

!primera ecuacion derecha

pepe = n1(0)

!print '(a30,f14.12)', 'primera ecuacion derecha = ',pepe


do i = 1, nn
  !segundas ecuaciones izquierdas
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
!  print '(a30,f14.12)', 'segunda ecuacion izquierda = ',pepe
  !segundas ecuaciones derechas
  pepe = 0.0d0
 ! print '(a30,f14.12)', 'segunda ecuacion derecha = ',pepe
 ! print *,'----------'
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!1
! here replace h(0,0) from observed r2
ne5=0.d0;ne1=0.d0;ne0=0.d0
varne5=0.d0;varne1=0.d0;varne0=0.d0
nen0=0.d0;varnen0=0.d0
do i=1,npop
   
   h(0,0)=or2(i,1)-er2(i,1)
call NfromL00(nn,c,p,super,h(0,0),n0)
!!!!!!!!!!!!!!!!!
!write(6,*)'all WCsparrow',n0,n0/2.d0
ne5=n0/2.d0+ne5
varne5=varne5+(n0/2.d0)**2
!
call NeFromN0(n0/2.d0,ne)
nen0(1)=nen0(1)+ne
varnen0(1)=varnen0(1)+ne**2

!freq 0.1
h(0,0)=or2(i,2)-er2(i,2)
!0 .000766!0.01454654
call NfromL00(nn,c,p,super,h(0,0),n0)
!write(6,*)'ne f=0.1 ',n0,n0/2.d0,h(0,0)
!!!!!!!!!!!!!!!!!
ne1=n0/2.d0+ne1
varne1=varne1+(n0/2.d0)**2

!
call NeFromN0(n0/2.d0,ne)
nen0(2)=nen0(2)+ne
varnen0(2)=varnen0(2)+ne**2

!--- freq 0.d0
h(0,0)=or2(i,3)-er2(i,3)
!0.03213202
call NfromL00(nn,c,p,super,h(0,0),n0)
!!!!!!!!!!!!!!!!!
ne0=n0/2.d0+ne0
varne0=varne0+(n0/2.d0)**2
!write(6,*)'ne f=0.d0 ',n0,n0/2.d0,h(0,0)
!
call NeFromN0(n0/2.d0,ne)
nen0(3)=nen0(3)+ne
varnen0(3)=varnen0(3)+ne**2
write(6,*)ne
enddo

ne5=ne5/dble(npop)
ne1=ne1/dble(npop)
ne0=ne0/dble(npop)
varne5=varne5/dble(npop)-ne5**2
varne1=varne1/dble(npop)-ne1**2
varne0=varne0/dble(npop)-ne0**2
!write(6,*)ne5,ne1,ne0
!write(6,*)varne5,varne1,varne0
!write(6,*)varne5*dble(npop-1)/dble(npop),varne1*dble(npop-1)/dble(npop),varne0*dble(npop-1)/dble(npop)
!
write(26,*)'our estimates of n0'
write(26,*)ne5,ne1,ne0
write(26,*)sqrt(varne5),sqrt(varne1),sqrt(varne0)
write(26,*)sqrt(varne5*dble(npop-1)/dble(npop)),sqrt(varne1*dble(npop-1)/dble(npop)),&
     sqrt(varne0*dble(npop-1)/dble(npop))
write(26,*)'ne estimator'
write(26,*)nest5,nest1,nest0
write(26,*)varnest5,varnest1,varnest0
write(26,*)sqrt(varnest5*dble(npop-1)/dble(npop)),sqrt(varnest1*dble(npop-1)/dble(npop)),&
     sqrt(varnest0*dble(npop-1)/dble(npop))
write(26,*)'our estimates of ne'
write(26,*)nen0(1)/dble(npop),nen0(2)/dble(npop),nen0(3)/dble(npop)
write(26,*)sqrt(varnen0(1)/dble(npop)-(nen0(1)/dble(npop))**2),sqrt(varnen0(2)/dble(npop)-(nen0(2)/dble(npop))**2),&
     sqrt(varnen0(3)/dble(npop)-(nen0(3)/dble(npop))**2)
write(26,*)sqrt(dble(npop-1)/dble(npop))*sqrt(varnen0(1)/dble(npop)-(nen0(1)/dble(npop))**2),&
     sqrt(dble(npop-1)/dble(npop))*sqrt(varnen0(2)/dble(npop)-(nen0(2)/dble(npop))**2),&
     sqrt(dble(npop-1)/dble(npop))*sqrt(varnen0(3)/dble(npop)-(nen0(3)/dble(npop))**2)
stop
end

!subroutine to calculate N from L00
subroutine NfromL00(n_ant,c,p,super,r2,n0)
  implicit none
  integer n_ant, i, j,k, iter
  real*8 c, p(n_ant),super(0:n_ant),r2, n0
  real*8 tau(n_ant), beta (n_ant, n_ant), eps(n_ant)
  real*8 delta(n_ant, n_ant), zeta(n_ant,n_ant)
  real*8 rho(n_ant), alpha(n_ant), xi(n_ant)
  real*8 L(0:n_ant, 0:n_ant)
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

  !valores iniciales para L
  L = 0.001d0
  L(0,0)=r2
  
  do iter = 1, 100
     !PRIMERA ECUACION
     sum_alphas = 0.0d0
     do i=1, n_ant
        sum_alphas = sum_alphas + alpha(i)* L(0,i) + p(i)**2*L(i,i)
     enddo
     n0 = 1.0d0 - (1.0d0 - c)**2*sum_alphas    !n0 es realmente 2*n0
     n0 = n0 / (r2 - (1.0d0 - c)**2 * sum_alphas)
 !    write(6,*)n0,r2,sum_alphas,c
!    pause
     !SEGUNDA ECUACION
     do i = 1, n_ant
        L(0,i) = tau(i)*L(i,i)
        do j = 1, n_ant
           if (j /= i) L(0,i) = L(0,i) + zeta (i,j) * L(0,j)
        enddo
        L(0,i)=L(0,i)/(1-zeta(i,i))
     enddo
     !Tercera eq , L_ii, i!=0
     rho=0.d0
     eps=0.d0
     do i=1,n_ant
        eps(i)= (1.d0 - 1.d0/(n0*super(i)) )/(1.d0- 1.d0/(n0*super(i-1)))        
        rho(i)=1.d0/(n0*super(i)) - eps(i)/(n0*super(i-1))
        if(i.eq.1)then
           L(i,i)=rho(i)+r2
        else
           L(i,i)=rho(i)+eps(i)*L(i-1,i-1)
        endif
     enddo
  enddo

end subroutine NfromL00
!
subroutine NeFromN0(n0,ne)
  implicit none
  integer nn, i,j, ii, n
parameter (nn = 6)
integer niter
real*8 h(nn,nn), h2(nn,nn), p(nn), diff, ne, plot(10000,3)
real*8 ni(nn), temp,n0

!para el calculo aproximado
real*8 l(nn),s(nn),d(nn),v(nn),t, b(nn)

ni=n0
niter = 100
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
!print *, i, ni(i)
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
!print *, 't ', t
!do i=1,6
!  print *,'v ',i,v(i),p(i)
!enddo

do i=1, nn-1
  s(i) = l(i+1)/l(i)
  d(i) = 1.0d0 - s(i)
enddo

ne = 0.0d0
do i = 1, nn-1 
  ne = ne + l(i)*s(i)*d(i)*(v(i+1)**2)
enddo
ne = (ni(1)*t)/(1.0d0+ne)
!print *, 'aprox ',ne


h  = 1.0d0
h2 = 0.0d0
do i = 1, nn
  h2(i,i) = 1.0d0 !no entiendo porque h2 debe ser la identidad. Luego se machaca
enddo

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
!    if (mod(ii,10)==0)print *,ii, ne
  h = h2
  plot(ii,1) = h(1,1)
enddo
!print *,'hola ',ne

!do i=1,niter
!  write(88,'(i10,4f12.5)')i,(plot(i,j),j=1,1)
!enddo

return
end subroutine NeFromN0
