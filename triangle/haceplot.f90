!This program prepare fort.62 written in Losdos...f90
!it writes "esquina" to be plotted in the gnuplot script "plot" 
implicit none
integer n, i, j, k, icicles, ii, jj, ix, iy, iz
real*8 xx,yy,x1,x2,x3
real*8 nn, nnt, tt
integer puntos, io
real*8 ies(0:500,0:500,0:500)  !combinaciones de W's
real*8 c, n2, n2m, t, sw, xmed, xmin, xmax
real*8 Q(0:200,0:200), w(200)
real*8 xnn
real*8 matriz(0:1000,0:1000)
puntos = 20
n = 3  !numero de cohortes
nn = 10000.0d0
c = 1.0d0/nn
ies = 0.0d0
open(77,file='esquina')
matriz = 0.0d0

open(62,file='fort.62')
do while (io==0)
  read (62,'(3i6,4f12.6)',iostat=io)ix,iy,iz,(w(ii),ii=1,3),Q(1,1)
  ies(ix,iy,iz) = Q(1,1)
  print *,ix,iy,iz,(w(ii),ii=1,3),Q(1,1)
enddo
print *,'pasa'
!dibujar eje x
do ix = 0, puntos
do iy = 0, puntos
  do iz = 0, puntos
    if(ies(ix,iy,iz) .gt. 0.000005)then
      x1 = dble(ix)/puntos
      x2 = dble(iy)/puntos
      x3 = dble(iz)/puntos 
      xx = x1 + x3 * 0.5d0
      yy = x3
      print *,ix,iy,iz, ' hola'
      if (ix.eq.puntos.or.iy.eq.puntos.or.iz.eq.puntos) then
        write(77,*)xx,yy,ies(ix,iy,iz),0
        print *, 'hola'
      else
        write(77,*)xx,yy,ies(ix,iy,iz)
        print *, 'hola2'
      end if
    end if
  enddo
enddo
do iy = puntos, 0, -1
  do iz = 0, puntos
    if(ies(ix,iy,iz) .gt. 0.000005)then
      x1 = dble(ix)/puntos
      x2 = dble(iy)/puntos
      x3 = dble(iz)/puntos 
      xx = x1 + x3 * 0.5d0
      yy = x3
      if (ix.eq.puntos.or.iy.eq.puntos.or.iz.eq.puntos) then
        write(77,*)xx,yy,ies(ix,iy,iz),0
      else
        write(77,*)xx,yy,ies(ix,iy,iz)
      end if
    end if
  enddo
enddo
enddo
!dibujar eje z
do iz = 0, puntos
do iy = 0, puntos
  do ix = 0, puntos
    if(ies(ix,iy,iz) .gt. 0.000005)then
      x1 = dble(ix)/puntos
      x2 = dble(iy)/puntos
      x3 = dble(iz)/puntos 
      xx = x1 + x3 * 0.5d0
      yy = x3
      if (ix.eq.puntos.or.iy.eq.puntos.or.iz.eq.puntos) then
        write(77,*)xx,yy,ies(ix,iy,iz),0
      else
        write(77,*)xx,yy,ies(ix,iy,iz)
      end if
    end if
  enddo
enddo
do iy = puntos, 0, -1
  do ix = 0, puntos
    if(ies(ix,iy,iz) .gt. 0.000005)then
      x1 = dble(ix)/puntos
      x2 = dble(iy)/puntos
      x3 = dble(iz)/puntos 
      xx = x1 + x3 * 0.5d0
      yy = x3
      if (ix.eq.puntos.or.iy.eq.puntos.or.iz.eq.puntos) then
        write(77,*)xx,yy,ies(ix,iy,iz),0
      else
        write(77,*)xx,yy,ies(ix,iy,iz)
      end if
    end if
  enddo
enddo
enddo
!dibujar eje y
do iy = 0, puntos
do ix = 0, puntos
  do iz = 0, puntos
    if(ies(ix,iy,iz) .gt. 0.000005)then
      x1 = dble(ix)/puntos
      x2 = dble(iy)/puntos
      x3 = dble(iz)/puntos 
      xx = x1 + x3 * 0.5d0
      yy = x3
      if (ix.eq.puntos.or.iy.eq.puntos.or.iz.eq.puntos) then
        write(77,*)xx,yy,ies(ix,iy,iz),0
      else
        write(77,*)xx,yy,ies(ix,iy,iz)
      end if
    end if
  enddo
enddo
do ix = puntos, 0, -1
  do iz = 0, puntos
    if(ies(ix,iy,iz) .gt. 0.000005)then
      x1 = dble(ix)/puntos
      x2 = dble(iy)/puntos
      x3 = dble(iz)/puntos 
      xx = x1 + x3 * 0.5d0
      yy = x3
      if (ix.eq.puntos.or.iy.eq.puntos.or.iz.eq.puntos) then
        write(77,*)xx,yy,ies(ix,iy,iz),0
      else
        write(77,*)xx,yy,ies(ix,iy,iz)
      end if
    end if
  enddo
enddo
enddo

stop

!---------------------------------------------
! SE conoce N y se estiman todos los Q
!estimar Q(1,0) hasta Q(1,n)
Q = 0.5d0


!ecuacion para Q(0,0)
do icicles = 1, 100000
  xmax = 1.0d0
  xmin = 0.0d0
  call eval_cero(Q,n,nn,c,sw,w,xmin,t)
  if (t .gt. 0) then
    print *,'error en xmin ', t
    stop
  end if
  !print *, xmin, t
  call eval_cero(Q,n,nn,c,sw,w,xmax,t)
  if (t .lt. 0) then
    print *,'error en xmax ', t
    stop 
  end if
  !print *, xmax, t
  do 
    xmed = (xmax + xmin)/2.0d0
    call eval_cero(Q,n,nn,c,sw,w,xmed,t)
    !print '(4f16.13)',xmin,xmax,xmed,t
    if (t .lt. 0)xmin = xmed
    if (t .gt. 0)xmax = xmed
    if (dabs(t) .lt. 1.0d-15) exit
  enddo
  Q(0,0) = xmed
  !print *, 0, Q(0,0)

  !ecuaciones para Q(0,i)
  do i = 1, n
    xmax = 1.0d0
    xmin = 0.0d0
    call eval_i(Q,i,n,nn,c,sw,w,xmin,t)
    if (t .gt. 0) then
      print *,'error en xmin'
      stop
    end if
    !print *, xmin, t
    call eval_i(Q,i,n,nn,c,sw,w,xmax,t)
    if (t .lt. 0) then
      print *,'error en xmax'
      stop
    end if
    !print *, xmax, t
    do 
      xmed = (xmax + xmin)/2.0d0
      call eval_i(Q,i,n,nn,c,sw,w,xmed,t)
      !print '(4f16.13)',xmin,xmax,xmed,t
      if (t .lt. 0)xmin = xmed
      if (t .gt. 0)xmax = xmed
      if (dabs(t) .lt. 1.0d-15) exit
    enddo
    Q(0,i) = xmed
    !print *, i, Q(0,i)
  enddo
  if (mod(icicles,1000).eq.0)print '(i10,6f15.12)',icicles,(Q(0,i),i=0,5)
enddo


stop
end

subroutine eval_n(Q,n,xi,c,sw,w,t)
implicit none
!xi es ahora el numero de bichos por cohorte
integer n, nn, j, k, idif
real*8 c, n2, n2m, sw, xi, t
real*8 Q(0:200,0:200), w(200)

t = Q(0,0) - (1.0d0 - Q(0,0)) * (1.0d0 - c)**2 * sw / (2.0d0 * xi)
do j = 1, n
  do k = 1, n
    t = t - (1.0d0 - c) ** 2 * w(j) * w(k) * Q(0,abs(k-j))
  enddo
enddo
return
end

subroutine eval_cero(Q,n,nn,c,sw,w,xi,t)
implicit none
integer n, j, k, idif
real*8 c, n2, n2m, sw, xi, t, nn
real*8 Q(0:200,0:200), w(200)

t = xi - (1-0d0 - xi) * (1.0d0 - c)**2 * sw / (2.0d0 * nn)
do j = 1, n
  do k = 1, n
    idif = abs(k-j)
    if (idif .eq. 0) then
      t = t - (1.0d0 - c) ** 2 * w(j) * w(k) * xi
    else
      t = t - (1.0d0 - c) ** 2 * w(j) * w(k) * Q(0,abs(k-j))
    end if
  enddo
enddo
return
end

subroutine eval_i(Q,i,n,nn,c,sw,w,xi,t)
implicit none
integer n, i, j, k
real*8 c, n2, n2m, sw, xi, t, nn
real*8 Q(0:200,0:200), w(200)

t = xi - (1-0d0 - Q(0,0)) * w(i) * (1.0d0 - c) / (2.0d0 * nn)
do k = 1, n
j = abs(k-i)
  if (j .eq. i) then
    t = t - (1.0d0 - c) * w(k) * xi
  else
    t = t - (1.0d0 - c) * w(k) * Q(0,abs(k-i))
  end if
enddo

return
end
