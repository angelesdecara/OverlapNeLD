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
integer nn, i, j, ii, n, icase
parameter (nn = 5)
integer ni(0:nn)
real*8 c, c1, c2, h(0:nn,0:nn), h2(0:nn,0:nn), p(nn), ratio
real*8  Q(0:nn,0:nn), Q2(0:nn,0:nn), super(0:nn)
! h y h2 es para iterar L, Q y Q2 es para iterar Q
real*8 n1(0:nn), n2(0:nn), Qtot
real*8 b, suma, pepe, beta,n0
integer k

open (26,file='salida50_caso1.txt')
c = 0.5d0
open(11,file='Exact_5cM.txt')

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
h(0,0)=0.009026165d0
call NfromL00(nn,c,p,super,h(0,0),n0)
!!!!!!!!!!!!!!!!!
write(6,*)'all WCsparrow',n0,n0/2.d0
h(0,0)=0.000264
!0 .000766!0.01454654
call NfromL00(nn,c,p,super,h(0,0),n0)
!!!!!!!!!!!!!!!!!
write(6,*)'pugetensis',n0,n0/2.d0
h(0,0)=0.000766
!0.03213202
call NfromL00(nn,c,p,super,h(0,0),n0)
!!!!!!!!!!!!!!!!!
write(6,*)'nuttalli',n0,n0/2.d0




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
