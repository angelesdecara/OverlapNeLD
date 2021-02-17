
implicit none
! n is the cohort size (as 2N)
! ngen is the numer of generations
! nrep is the numer of replicates
! n_ant is the numer of cohorts contributing 
! c is for recombination frequency
! g1 and g2 are the allelic frequencies p(A) and p(B)
! d is for LD
integer n, ngen, n_ant, nrep, icase
parameter ( ngen = 500, n_ant = 3, nrep = 100000)
real*8 x1, c, p(n_ant), g1, g2, f1
integer irep, nt, i, j, k, ni(0:n_ant)
real*8 frecs(4,0:n_ant,nrep), frecst(4)
real*8 g1_temp(n_ant), g2_temp(n_ant)
integer t(4), tt(4), nrep_vivas, iant, norecomb, recomb, igen
integer haplohijos(n_ant)
real*8 r2_hist(ngen), bnldev, d, dd(ngen)
real*8 r2_mean   !promedio de las ultimas 1000 generaciones
integer hist(0:100,2,3), itemp, imuertos
real*8 super(0:n_ant), n_resul, r2_aqui
real*8 r2perrep(nrep), temp(500), ehtemp, ehtemp0

ehtemp0=0.d0;ehtemp=0.d0
hist = 0
!simulation parameters

!=========================================
!comment lines to choose the current c
!ghange filenames accordingly 
!c = 0.31606d0     !for 50cm
!c = 0.16484d0     !for 20cm
!c = 0.0906d0      !for 10cm
c = 0.0475813d0   !for 5cm
open(33,file='infr2ne_5cm.dat')
open(44,file='n_r2_nest_5cm.dat')
!=========================================

do icase=1,3 !the three scenarios

   super=1.d0
   p=1.d0/3.d0
   if(icase.eq.3)then
      super(0) = 1.0d0
      super(1) = 1.0d0
      super(2) = 0.8d0
      super(3) = 0.5d0
      p(1) = 0.3d0 ;         
      p(2) = 0.5d0 ;         
      p(3) = 1.0d0-p(1)-p(2) 
   elseif(icase.eq.2)then
      p(1) = 0.3d0 ;         
      p(2) = 0.5d0 ;         
      p(3) = 1.0d0-p(1)-p(2) 
   endif

!random number seed
x1 = 0.5489389222d0

!
g1 = 0.5d0     !atarting allelic frequencies A
g2 = 0.5d0     !Id.                          B

do n=100,10000,100
   write(6,*)n,icase

   ni(0) = n

ni = n * super

nt = sum(ni)   ! total size 

!starting haplotipic frequencies for each cohort
frecs(1,0:n_ant,:) = g1*g2
frecs(2,0:n_ant,:) = g1*(1.0d0-g2)
frecs(3,0:n_ant,:) = (1.0d0-g1)*g2
frecs(4,0:n_ant,:) = (1.0d0-g1)*(1.0d0-g2)

open (5, file = 'salidanueva3_5cm_long_caso3.txt')
r2_mean = 0.0d0

do igen = 1, ngen
  nrep_vivas = 0
  r2_hist(igen) = 0.0d0
  do irep = 1, nrep
    ! calculate frequencies
     f1=0.d0
     ehtemp=0.d0
    do iant=1,n_ant
      f1=f1+p(iant)*(frecs(1,iant,irep)+frecs(3,iant,irep))
    enddo
    !calculate newborn gametes
    tt = 0  
    !how many gametes come from each cohort?
    call multidev(n_ant, p(1:n_ant), haplohijos, n, x1)
    do iant = 1, n_ant
      !how many recombine?
      recomb = bnldev(c,haplohijos(iant),x1)
      norecomb = haplohijos(iant) - recomb
      !frequencies for recombinants
      call calc_frecs(frecst, frecs(1,iant,irep) + frecs(2,iant,irep),f1,0.0d0)
      call multidev(4, frecst, t, recomb, x1)
      !print *,' frecs recomb ', t
      tt = tt + t
      !frequencies for non recombinants
      call multidev(4, frecs(:,iant,irep), t, norecomb, x1)
      !print *,' frecs norecomb ', t
      tt = tt + t
    enddo
    !calculate frequencies for newborns
    !print '(6i10,7f12.4)',irep, nrep, tt
    frecs(:,0,irep) = dble(tt) / n
    !calculate disequibrium for newborns
    frecst = frecs(:,0,irep)
    call calc_d(frecst, g1, g2, d)


    if (g1 .gt. 0.001d0 .and. g1 .lt. 0.999d0 .and. g2 .gt. 0.001d0 .and. g2 .lt. 0.999d0)then
      nrep_vivas = nrep_vivas + 1
      r2_hist(igen) = r2_hist(igen) + d*d/(g1*g2*(1.0d0-g1)*(1.0d0-g2))
      r2_aqui = d*d/(g1*g2*(1.0d0-g1)*(1.0d0-g2))
      dd(igen)=d
      !handle fixed replicates
    else
      if(irep.eq.1)then
        frecs(:,:,irep) = frecs(:,:,nrep) 
      else
        frecs(:,:,irep) = frecs(:,:,irep - 1)
      end if
    end if 

   do iant = n_ant, 1, -1
      !apply survival
      imuertos = ni(iant-1) - ni(iant)
      call multidev(4, frecs(:, iant-1, irep), t, imuertos, x1)

      do i = 1, 4
         frecs(i,iant,irep)=(dble(ni(iant-1))*frecs(i,iant-1,irep)-t(i))/dble(ni(iant))
     enddo
   enddo
     if(igen.eq.ngen)then
       call calc_N(n_ant, c, p, super, r2_aqui, n_resul)
       write(22,*)irep,n_resul/2,r2_aqui
       r2perrep(irep)=r2_aqui
    endif
 enddo
  r2_hist(igen) = r2_hist(igen) / nrep_vivas
  if (mod(igen,1).eq.0) then
    write(*,'(2i10,4f14.8)') igen, nrep_vivas, r2_hist(igen), g1, g2 
    write(5,'(2i10,4f14.8)') igen, nrep_vivas, r2_hist(igen) 
    write(15,'(2i10,4f14.8)') igen, nrep_vivas, r2_hist(igen) 
  end if 
  if (igen .gt. ngen-100) r2_mean = r2_mean + r2_hist(igen)
  write(11,*)dd(igen)
  if(sum(r2_hist(igen-9:igen))-sum(r2_hist(igen-10:igen-1)).lt.1.d-5)then
     call calc_N(n_ant, c, p, super, r2_hist(igen), n_resul)
     write(6,*)igen,r2_hist(igen),n_resul
     write(44,*)icase,n,igen,r2_hist(igen),n_resul
     exit
  endif

enddo
    close(5)
    print '(a,f18.10)','r2 medio = ', r2_mean / 100
    call calc_N(n_ant, c, p, super, r2_mean/100, n_resul)
    write(22,*)'#aver ',n_resul/2,r2_mean/100

 enddo
enddo
    close(33)

    stop
end program

subroutine calc_N(n_ant, c, p, super, r2, n0)
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

  !starting values for L
  L = 0.001d0
  
  do iter = 1, 100
     !first equation
     sum_alphas = 0.0d0
     do i=1, n_ant
        sum_alphas = sum_alphas + alpha(i)* L(0,i) + p(i)**2*L(i,i)
     enddo
     n0 = 1.0d0 - (1.0d0 - c)**2*sum_alphas   
     n0 = n0 / (r2 - (1.0d0 - c)**2 * sum_alphas)
     !2nd equation
     do i = 1, n_ant
        L(0,i) = tau(i)*L(i,i)
        do j = 1, n_ant
           if (j /= i) L(0,i) = L(0,i) + zeta (i,j) * L(0,j)
        enddo
        L(0,i)=L(0,i)/(1-zeta(i,i))
     enddo
     !3rd equation , L_ii, i!=0
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
  
end

include 'multis.f90'
