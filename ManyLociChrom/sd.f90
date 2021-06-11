
!prog to calc sd of a harmonic variable
!we've got phi=sigma^2_R/mu^2_R and 1/mu_R

program readstuff
 implicit none
 integer::io,i,j,k,l,irep
 integer,parameter::dp=kind(1.d0)!,nrep=1000
 integer::nloci,ncrom,ncomp,nrep
 real(dp)::nest,mur,varrep,phi,realcomp
 real(dp)::nemean,varne

 nemean=0.d0;varne=0.d0
 open(11,file='n2000loc1000_chr5_ne.dat')
 io=0;irep=0
! do irep=1,nrep
 do while(io.eq.0)
  read(11,*,iostat=io)nloci,ncrom,phi,nest,realcomp,ncomp
  if(io.ne.0)exit
  if(.not.isnan(nest))then
     mur=1.d0/nest
     varrep=phi*mur**2

     write(6,*)sqrt(varrep),sqrt(1.d0/varrep),sqrt(varrep/(ncomp*mur**4))
     nemean=nemean+nest
     varne=varne+nest**2
     irep=irep+1
  endif
 enddo
 nrep=irep

 nemean=nemean/dble(nrep)
 varne=varne/dble(nrep)-nemean**2
 write(6,*)nemean/2.d0,sqrt(varne)/2.d0

end program
