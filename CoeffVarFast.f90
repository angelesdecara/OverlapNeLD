!
!returns nloci ncrom coeffVar HarmMean and so on
implicit none
integer nloci,n, ngen, ncoh, ncro,nrep
!nloci, n, ncoh   numero de loci, numero de gametos, numero de cohortes
!numero de cromosomas
parameter (nloci = 1000, n = 2000, ngen = 500, ncoh = 3, ncro = 5,nrep=1000)
real*8 x1, unif, tt, bnldev
integer lunif, i,j,k,l, igen, kk,irep
real*8 locipos(nloci)  !posicion de cada locus
real*8 p(nloci)        !frecuencia alelica inicial
real*8 longitud        !longitud del cromosoma
real*8 recomb(10000)     !puntos de recombinacion
integer genotipo(nloci, n, 0:ncoh) !genotipos de padres e hijos 
integer ancestro(2,2),cro(nloci)   !ancestro (numero cohorte, numero individuo)
real*8 p1(nloci),p2, p12, d, c, r2, temp
real*8 p_coh(ncoh), sup(0:ncoh)      !contribución, supervivencia de cada cohorte
real*8 pco(ncoh)                   ! copia de p_coh sin acumular para la llamada a NfromL00
integer vivos(ncoh)       !numero de vivos de cada cohorte
real*8 n_est              !tamaño estimado para cada pareja de loci
real*8 n_harm             !media armonica de los tamaños
real*8 i_harm,meaninv,coeffvar,invvar
character(len=3)::acrom
character(len=5)::aloci,aind

write(acrom,'(i3)')ncro
acrom=trim(acrom)
write(aloci,'(i5)')nloci
aloci=trim(aloci)
write(aind,'(i5)')n
aind=trim(aind)
open(11,file='n'//trim(adjustl(aind))//'loc'//trim(adjustl(aloci))//'_chr'//trim(adjustl(acrom))//'_ne.dat')

x1 = 0.53648532d+10
longitud = ncro   !ncro cromosomas de longitud 1

!contribuciones de cada cohorte
p_coh(1) = 0.3d0
p_coh(2) = 0.50d0
p_coh(3) = 0.20d0
pco = p_coh
do i = 2, ncoh
  p_coh(i) = p_coh(i) + p_coh(i-1)
enddo
!superviventes de cada cohorte
sup(0) = 1.0d0
sup(1) = 1.0d0
sup(2) = 0.8d0
sup(3) = 0.5d0
!vivos iniciales
do i = 1, ncoh
  vivos(i) = bnldev(sup(i),n,x1)
enddo
print *,'contribuciones ', p_coh
print *,'supervivencia ',  sup
print *,'vivos ',          vivos

do irep=1,nrep
!simular localizaciones de los loci
do i = 1, nloci
  locipos(i) = unif(x1) * longitud
  p(i)       = unif(x1)
  !print *,i,p(i)
enddo

!ordenar localizaciones de los loci y asignar cromosoma
do i = 1, nloci
  do j = 1, nloci-i
    if (locipos(j+1) < locipos(j)) then
      temp = locipos(j+1) 
      locipos(j+1) = locipos(j)
      locipos(j) = temp
    end if
  enddo
enddo
do i = 1, nloci
  cro(i) = int(locipos(i)) + 1
enddo

!simular genotipos iniciales
genotipo = 1
do k = 1, ncoh
  do j = 1, n
    do i = 1, nloci
      if (unif(x1) .gt. p(i)) genotipo (i,j,k) = 0
    enddo
    !print '(2i6,x,200i1)',k,j,genotipo(:,j,k)
  enddo
enddo
!generaciones
do igen = 1, ngen
!   print *, igen
   do i = 1, n
     !Elegir gametos parentales (autofecundacion incluida)
     do k = 1, 2
       temp = unif(x1)
       ancestro(k,1) = 1
       do j = ncoh, 1, -1
         if (p_coh(j) > temp)ancestro(k,1) = j
       enddo
       ancestro(k,2) = lunif(x1,vivos(ancestro(k,1)))
     enddo
     !print *, 'padre ',ancestro(1,1),ancestro(1,2),' madre ',ancestro(2,1),ancestro(2,2)
     !pause
     l = 0
     recomb = 3.0d0
     do j = 1, ncro
       !recombinacion entre cromosomas (j-1) y j
       if (unif(x1) < 0.5d0) then
         l = l + 1
         recomb(l) = dble(j-1)
       end if
       !Elegir puntos de recombinación cromosoma j
       tt = 0.0d0     
       do
         temp = -dlog(unif(x1))
         tt = tt + temp
         !print *, tt, temp
         if (tt > 1.0d0) exit
         l = l + 1
         recomb(l) = tt + dble(j-1)
       enddo
     enddo
     l=l+1
     recomb(l)=ncro+1
     !print '(a,i5,100f10.6)', 'recomb ',l,recomb(1:l)     

     k = 1  !empieza con el ancestro 1
     kk = 1
     do j = 1, nloci
       !print *,locipos(j),recomb(kk)
       do while (locipos(j) > recomb(kk))
         kk = kk + 1   !pasa a la siguiente recombinacion
         k = 3 - k     !cambia de ancestro
         !print *,'salta'
       enddo
       !print *,i,j,k,ancestro(k,1),ancestro(k,2)
       genotipo(j,i,0) = genotipo(j,ancestro(k,2),ancestro(k,1))
     enddo
     !print '(a,i6,x,200i1)','nuevo ',i,genotipo(:,i,0)
   enddo
   !cambiar padres por hijos
   do i=ncoh,1,-1
     genotipo(:,:,i) = genotipo(:,:,i-1)
  enddo
!  if(mod(igen,10).eq.0)then
!     write(6,*)igen,count(genotipo()==1)
!  endif
  !genotipo(:,:,0) = 0
enddo
!   tachan!
!salida
! calc allelic freqs
p1=0.d0
do i=1,nloci
   do j=1,n
      if(genotipo(i,j,0).eq.1)p1(i)=p1(i)+1.d0
      
   enddo
enddo
p1=p1/n
!
n_harm = 0.0d0
i_harm = 0
meaninv=0.d0;invvar=0.d0
do i = 1, nloci
  do j = i+1, nloci
    !frecuencias
    p12 = 0.0d0
    do k = 1, n
      if (genotipo(i,k,0) .eq. 1 .and. genotipo(j,k,0) .eq. 1) p12 = p12 + 1
    enddo
    p12 = p12 / n
    r2 = p12 - p1(i) * p1(j)
    r2 = (p12 - p1(i) * p1(j))*(p12 - p1(i) * p1(j))/(p1(i)*p1(j)*(1.0d0-p1(i))*(1.0d0-p1(j)))
    d = abs(locipos(i)-locipos(j))  !distancia de mapa
      !print *, ' p ',p1,p2,p12
    if (abs(p1(i)-0.5) < 0.4d0 .and. abs(0.5-p1(j)) < 0.4d0) then
      !print *, ' p ',p1,p2,p12
      if (cro(i) /= cro(j)) then
        c = 0.5d0
      else
        c = 0.5d0 * (1.0d0 -dexp(-2.0d0 * d))
      end if 
      call NfromL00(ncoh,c ,pco,sup,r2,n_est)
      n_harm = n_harm + 1.0d0 / n_est
      invvar=invvar+(1.d0/n_est)**2
      i_harm = i_harm + 1
!      write(5,'(2i7,4f12.9,2f15.3)')i, j, locipos(i), locipos(j), r2, c, n_est !,1.0d0/(n_harm/i_harm)
!      if(mod(j,500)==0 .and. mod(i,10)==0)print '(2i7,4f12.9,f12.3)',i, j, locipos(i), locipos(j), r2, c &
!                                                             &,1.0d0/(n_harm/i_harm)
    end if
 enddo
enddo

n_harm=n_harm/dble(i_harm)
invvar=invvar/dble(i_harm)-(n_harm)**2

coeffvar=(invvar)/n_harm**2
write(11,'(2i7,4f14.4)')nloci,ncro,(invvar)/n_harm**2,1.d0/n_harm,2.d0/(coeffvar/i_harm),i_harm

enddo
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

	
include 'multis.f90'
