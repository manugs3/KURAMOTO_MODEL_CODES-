
      program kuramoto_model2D
      implicit none
      real*8, allocatable :: omegas(:,:), thetas(:,:), f(:,:), aux(:,:)
      real*8, allocatable :: k1(:,:), k2(:,:), k3(:,:), k4(:,:)
      
      integer :: n, i, j, nciclos
      real*8 :: pi, k, t, h,  tol, gamma
      real*8 dran_u
      parameter(pi=4*datan(1.d0))
      call dran_ini(80)
      
      print*, "Introduzca numero de osciladores de la red NxN"
      read*, n
      allocate(omegas(n,n), thetas(n,n), f(n,n), aux(n,n))
      allocate(k1(n,n), k2(n,n), k3(n,n), k4(n,n))
      
      do i=1, n
         do j=1, n
            thetas(i,j) = -pi+dran_u()*2.d0*pi
         end do
      end do
         
      omegas = 0.d0
         
!     ABRO LOS FICHEROS
      open(5, file='thetas-t-2D.dat')
      do i = 1, n
         write(5,*) (thetas(i,j), j=1,n)
      end do
      write(5,*)''
      write(5,*)''
      
      open(3, file='4-tiempos.dat')
      do i = 1, n
         write(3,*) (thetas(i,j), j=1,n)
      end do
      write(3,*)''
      write(3,*)''
      
!     PROGRAMAMOS EL ALGORITMO DE RUNGE KUTTA
      k = 1.d0      
      nciclos = 0
      t = 0.d0
      h = 0.02d0
      tol = 1.d-5  
      do while(t.lt.5001.d0)
         
         call derivadas_2D(omegas, thetas, f, n, k)
         k1 = h * f
         aux = thetas + 0.5d0 * k1
         
         call derivadas_2D(omegas, aux, f, n, k)
         k2 = h * f
         aux = thetas + 0.5d0 * k2
            
         call derivadas_2D(omegas, aux, f, n, k)
         k3 = h * f
         aux = thetas + k3
            
         call derivadas_2D(omegas, aux, f, n, k)
         k4 = h * f
            
         thetas = thetas + (1.d0/6.d0)*(k1 + 2.d0*k2 + 2.d0*k3 + k4)
         do i = 1, n
            do j = 1, n
               if (thetas(i,j) > pi) then
                  thetas(i,j) = thetas(i,j) - 2.d0*pi
               else if (thetas(i,j) <= -pi) then
                  thetas(i,j) = thetas(i,j) + 2.d0*pi
               end if
            end do
         end do
         
         if (mod(nciclos, 50)==0) then
            do i = 1, n
               write(5,*) (thetas(i,j), j=1,n)
            end do
            write(5,*)''
            write(5,*)''
         end if
         
         if (abs(t - 5.d0) < tol .or. abs(t - 200.d0) < tol.or.
     &    abs(t - 5000.d0) < tol) then
            do i = 1, n
               write(3,*) (thetas(i,j), j=1,n)
            end do
            write(3,*)''
            write(3,*)''
            print*, t 
         end if
         
         t = t + h
         nciclos = nciclos + 1
      end do
      
      close(3)
      close(5)
      end program

      subroutine derivadas_2D(omegas, thetas, f, n, k)
      implicit none
      integer :: n, i, j
      real*8 :: thetas(n,n), f(n,n), omegas(n,n), k
      real*8 :: thetas_aux(0:n+1, 0:n+1)
      
!     CREO UNA MATRIZ AUXILIAR PARA METER LAS CONDICIONES DE CONTORNO PERIODICAS
      thetas_aux=0.d0
      do i = 1, n
         do j = 1, n
            thetas_aux(i,j)=thetas(i,j)
         end do
      end do
      
      do i = 1, n
         thetas_aux(0,i) = thetas_aux(n,i)
         thetas_aux(n+1,i) = thetas_aux(1,i)
         thetas_aux(i,0) = thetas_aux(i,n)
         thetas_aux(i,n+1) = thetas_aux(i,1)
      end do     
      
      do i = 1, n
         do j = 1, n
            f(i,j) = omegas(i,j) + k * (
     &           sin(thetas_aux(i+1,j) - thetas_aux(i,j)) +
     &           sin(thetas_aux(i-1,j) - thetas_aux(i,j)) +
     &           sin(thetas_aux(i,j+1) - thetas_aux(i,j)) +
     &           sin(thetas_aux(i,j-1) - thetas_aux(i,j)) )
         end do
      end do
      end subroutine
      
      include 'dranxor2_new.f'
