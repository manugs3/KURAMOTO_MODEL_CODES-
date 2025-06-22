      program kuramoto_model
      implicit none
!     VECTOR CON LAS FRECUENCIAS
      real*8, allocatable :: omegas(:), thetas(:), f(:), aux(:)
      real*8, allocatable :: k1(:), k2(:), k3(:), k4(:)
      real*8, allocatable :: x_pos(:), y_pos(:) 
      
      integer :: n, i, nciclos
      real*8 :: pi, arg,k, gamma, x, delta_x, t, h, u1, u2, z1, z2
      real*8 :: R_x, R_y
      complex*16 :: R
      
      parameter(pi=4*datan(1.d0))
      parameter(gamma=1.d0)
      real*8 dran_u
      call dran_ini(12371)
      
      print*, "Introduzca numero de osciladores"
      read*, n
      allocate(omegas(n), thetas(n), f(n), aux(n))
      allocate(k1(n), k2(n), k3(n), k4(n), x_pos(n), y_pos(n))
      
!---------------------------------------     
!----INICIALIZO LAS FASES---------------
!---------------------------------------
!      thetas=0.d0
      do i=1, n
         thetas(i)=-pi+dran_u()*2.d0*pi
      end do
      


!---------------------------------------
!---- INICIALIZO LAS FRECUENCIAS--------
!---------------------------------------

      
!     INICIALIZO LAS FRECUENCIAS SEGÚN UNA DISTRIBUCIÓN LORENTZIANA TRUNCADA
      do i = 1, n-1
         do
            omegas(i) = gamma * dtan(pi * (dran_u() - 0.5d0))
            if (abs(omegas(i)) .le. 3.d0) exit
         end do
      end do
      omegas(n)=5.d0

!     INICIALIZO LAS OMEGAS SEGÚN UNA DISTRIBUCIÓN GAUSSIANA (Box-Muller)
!      i=1
!      do while(i.le.n)
!         u1=dran_u()
!         u2=dran_u()
!         if(u1.eq.0.d0) cycle
!         z1=sqrt(-2.d0*log(u1))*cos(2.d0*pi*u2)
!         z2=sqrt(-2.d0*log(u1))*sin(2.d0*pi*u2)
!         omegas(i)=z1
!         if(i+1.le.n) then
!            omegas(i+1)=z2
!         end if
!         i=i+2
!      end do

!     INICIALIZO LAS OMEGAS CON UNA DISTRIBUCIÓN BIMODAL TIPO DIRAC USANDO dran_u()
!      omegas=1.d0
!      do i = 1, n
!         if (dran_u().le.0.5d0) then
!            omegas(i) = -1.d0
!         end if
!      end do

      
!     ABRO LOS FICEROS
!      open(3, file='coordenadas-t.dat')
      open(4, file='omegas-t.dat')
      open(5, file='thetas-t.dat')
!      open(7, file='R-k.dat')
    
!     ALGORITMO DE RUNGE KUTTA
      k=4.d0
!     do while(k.lt.100.d0)
      nciclos=0
      t=0.d0
      h=0.002d0
      do while(t.lt.10.d0)
         call derivadas_mf(omegas, thetas, f, n, k, R)
         do i=1, n
            k1(i)=h*f(i)
            aux(i)=thetas(i)+k1(i)/2.d0
         end do
         call derivadas_mf(omegas, aux, f, n, k, R)
         do i=1, n
            k2(i)=h*f(i)
            aux(i)=thetas(i)+k2(i)/2.d0
         end do
         call derivadas_mf(omegas, aux, f, n, k, R)
         do i=1, n
            k3(i)=h*f(i)
            aux(i)=thetas(i)+k3(i)
         end do
         call derivadas_mf(omegas, aux, f, n, k, R)
         do i=1, n
            k4(i)=h*f(i)
         end do
!     call derivadas_mf(omegas, thetas, f, n, k, R)
         
!     ACTUALIZO VALOR DE LOS ANGULOS y CALCULO POSICIONES
         do i=1, n
            thetas(i)=thetas(i)+(1.d0/6.d0)*(k1(i)+2.d0*k2(i)+2.d0
     &           *k3(i)+k4(i))
            if (thetas(i).gt. pi) then
               thetas(i) = thetas(i) - 2.d0*pi
            else if (thetas(i).lt. -pi) then
               thetas(i) = thetas(i) + 2.d0*pi
            end if
!     x_pos(i)=dcos(thetas(i))
!     y_pos(i)=dsin(thetas(i))
         end do
!     R_x=real(R)
!     R_y=aimag(R)           
            
         
         if(mod(nciclos, 50).eq.0)then
!     write(3,*) t,R_x, R_y, (x_pos(i), y_pos(i), i=1,n)
            write(4,*) t, (f(i), i=100, n, 100)
            write(5,*) t, (thetas(i), i=100, n, 100)
            print*, t
         end if
         
!     ACTUALIZO PASO
         t=t+h
         nciclos=nciclos+1
      end do
      print*, abs(R) 
!     write(7,*) k, abs(R)
      
!     if(k .le. 2.5) then
!     k = k + 0.1
!     elseif(k .gt. 2.5 .and. k .le. 10) then 
!     k = k + 1.0
!     elseif(k .gt. 10 .and. k .le. 100) then
!     k = k + 10.0
!     end if
!     end do
      
!      close(3)
      close(4)
      close(5)
!      close(7)
        
         
      end program
      
!     subroutine derivadas(omegas, thetas, f, n)
!      integer ::  n, i
!      real*8 :: k, suma
!      real*8 :: omegas(n), thetas(n), f(n)
!      k=1.d0
!      suma=0.d0
!      do i=1, n
!         do j=1, n
!            suma=suma+sin(thetas(j)-thetas(i))
!         end do
!         f(i)=omegas(i)+(k/real(n))*suma 
!      end do
!      end subroutine


      subroutine derivadas_mf(omegas, thetas, f, n, k, R)
      integer ::  n, i, l
      real*8 :: k, arg
      real*8 :: omegas(n), thetas(n), f(n)
      complex*16 :: j, R
      parameter (j=complex(0.d0, 1.d0))
      R=(0.d0, 0.d0)
      do l=1, n
         R=R+exp(j*thetas(l))
      end do
      R=R/real(n)
      do i=1, n         
         f(i)=omegas(i)+k*abs(R)*sin(arg(R)-thetas(i))
      end do 
      end subroutine
      
      real*8 function arg(z)
      implicit none
      complex*16 :: z
      arg = atan2(aimag(z), real(z))
      end function arg

      include 'dranxor2_new.f'
 
