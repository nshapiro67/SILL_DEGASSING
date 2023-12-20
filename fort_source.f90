!--------------- derivations of basic functions --------------
!   element number num  point coordinate inside the element
! dr(k,i);  k=node number (1,2,3,4);  i=coordinate(x,y,z)
!--------------------------------------------------------

subroutine derivation(dr,n)
  use sizes
  use node_data
  use element_data
  implicit none

  integer::i,n
  real(kind=8):: xn(4),yn(4),zn(4)
  real(kind=8):: dr(4,3)
  real(kind=8):: a,b,c,d,e,vol,v1,v2,v3,v4,onv

  ! coord 4 element nodes
  do i = 1,4
    xn(i) = cord(nop(n,i),1)
    yn(i) = cord(nop(n,i),2)
    zn(i) = cord(nop(n,i),3)
  end do

  ! element volume  (G = 6 * volume )
  !| 1 x1 y1 z1 |
  !| 1 x2 y2 z2 |
  !| 1 x3 y3 z3 |
  !| 1 x4 y4 z4 |

  v1 = xn(2) * (yn(3)*zn(4)-zn(3)*yn(4)) &
     - xn(3) * (yn(2)*zn(4)-zn(2)*yn(4)) &
     + xn(4) * (yn(2)*zn(3)-zn(2)*yn(3))

  v2 = xn(1) * (yn(3)*zn(4)-zn(3)*yn(4)) &
     - xn(3) * (yn(1)*zn(4)-zn(1)*yn(4)) &
     + xn(4) * (yn(1)*zn(3)-zn(1)*yn(3))

  v3 = xn(1) * (yn(2)*zn(4)-zn(2)*yn(4)) &
     - xn(2) * (yn(1)*zn(4)-zn(1)*yn(4)) &
     + xn(4) * (yn(1)*zn(2)-zn(1)*yn(2))

  v4 = xn(1) * (yn(2)*zn(3)-zn(2)*yn(3)) &
     - xn(2) * (yn(1)*zn(3)-zn(1)*yn(3)) &
     + xn(3) * (yn(1)*zn(2)-zn(1)*yn(2))

  vol = v1 - v2 + v3 - v4

  if(vol .ge.0.0_8 )then
    write(6,*)'	Negative volume ',-vol/6,' element ', n
    write(6,*)-v1,v2,-v3,v4
    stop
  end if
  onv= 1.0_8/vol

  !A = | 1 y2 z2 |
  !    | 1 y3 z3 |
  !    | 1 y4 z4 |
  a = yn(3)*zn(4)-zn(3)*yn(4) &
    -(yn(2)*zn(4)-zn(2)*yn(4)) &
    + yn(2)*zn(3)-zn(2)*yn(3)

  !B = | 1 x2 z2 |
  !    | 1 x3 z3 |
  !    | 1 x4 z4 |
  b = xn(3)*zn(4)-zn(3)*xn(4)  &
    -(xn(2)*zn(4)-zn(2)*xn(4)) &
    + xn(2)*zn(3)-zn(2)*xn(3)

  !C = | 1 x2 y2 |
  !    | 1 x3 y3 |
  !    | 1 x4 y4 |
  c = xn(3)*yn(4)-yn(3)*xn(4)  &
    -(xn(2)*yn(4)-yn(2)*xn(4)) &
    + xn(2)*yn(3)-yn(2)*xn(3)

  dr(1,1) = -a * onv
  dr(1,2) =  b * onv
  dr(1,3) = -c * onv

  !A = | 1 y1 z1 |
  !    | 1 y3 z3 |
  !    | 1 y4 z4 |
  a = yn(3)*zn(4)-zn(3)*yn(4)  &
    -(yn(1)*zn(4)-zn(1)*yn(4)) &
    + yn(1)*zn(3)-zn(1)*yn(3)

  !B = | 1 x1 z1 |
  !    | 1 x3 z3 |
  !    | 1 x4 z4 |
  b = xn(3)*zn(4)-zn(3)*xn(4)  &
    -(xn(1)*zn(4)-zn(1)*xn(4)) &
    + xn(1)*zn(3)-zn(1)*xn(3)

  !C = | 1 x1 y1 |
  !    | 1 x3 y3 |
  !    | 1 x4 y4 |
  c = xn(3)*yn(4)-yn(3)*xn(4)  &
    -(xn(1)*yn(4)-yn(1)*xn(4)) &
    + xn(1)*yn(3)-yn(1)*xn(3)

  dr(2,1) =  a * onv
  dr(2,2) = -b * onv
  dr(2,3) =  c * onv

  !A = | 1 y1 z1 |
  !    | 1 y2 z2 |
  !    | 1 y4 z4 |
  a = yn(2)*zn(4)-zn(2)*yn(4)  &
    -(yn(1)*zn(4)-zn(1)*yn(4)) &
    + yn(1)*zn(2)-zn(1)*yn(2)

  !B = | 1 x1 z1 |
  !    | 1 x2 z2 |
  !    | 1 x4 z4 |
  b = xn(2)*zn(4)-zn(2)*xn(4)  &
    -(xn(1)*zn(4)-zn(1)*xn(4)) &
    + xn(1)*zn(2)-zn(1)*xn(2)

  !C = | 1 x1 y1 |
  !    | 1 x2 y2 |
  !    | 1 x4 y4 |
  c = xn(2)*yn(4)-yn(2)*xn(4)  &
    -(xn(1)*yn(4)-yn(1)*xn(4)) &
    + xn(1)*yn(2)-yn(1)*xn(2)

  dr(3,1) = -a * onv
  dr(3,2) =  b * onv
  dr(3,3) = -c * onv

  !A = | 1 y1 z1 |
  !    | 1 y2 z2 |
  !    | 1 y3 z3 |
  a = yn(2)*zn(3)-zn(2)*yn(3)  &
    -(yn(1)*zn(3)-zn(1)*yn(3)) &
    + yn(1)*zn(2)-zn(1)*yn(2)

  !B = | 1 x1 z1 |
  !    | 1 x2 z2 |
  !    | 1 x3 z3 |
  b = xn(2)*zn(3)-zn(2)*xn(3)  &
    -(xn(1)*zn(3)-zn(1)*xn(3)) &
    + xn(1)*zn(2)-zn(1)*xn(2)

  !C = | 1 x1 y1 |
  !    | 1 x2 y2 |
  !    | 1 x3 y3 |
  c = xn(2)*yn(3)-yn(2)*xn(3)  &
    -(xn(1)*yn(3)-yn(1)*xn(3)) &
    + xn(1)*yn(2)-yn(1)*xn(2)

  dr(4,1) =  a * onv
  dr(4,2) = -b * onv
  dr(4,3) =  c * onv

  return
end subroutine derivation
! Elastic stress-strain for damage-breakage material
! returns stress in element and eigenvalue

subroutine elastic(n,se)
use sizes
use element_data
implicit none

integer:: n,j
real(kind=8):: rl,rm,se(6)
real(kind=8):: i1,b

! stress calculation
! elastic moduli
    rl  = lambda(n)
    rm  = mu(n)
 
! stress for solid
      i1 = str_e(1,n) + str_e(2,n) + str_e(3,n)
      b  = rl*i1

      do j=1,3
        se(j)   = b + 2.0_8*rm*str_e(j,n)
        se(j+3) =     2.0_8*rm*str_e(j+3,n)
      end do

  return 

end subroutine elastic

subroutine flac(dt,vvv)
  use sizes
  use node_data
  use element_data
  implicit none

  integer::n,j,i,noff,nmax,nn
  real(kind=8)::dt,vvv,se(6)
  
  ! Initialization 
  force=0.0_8
  balance=0.0_8
  
! loop through all elements 
  !$OMP PARALLEL
  !$OMP DO PRIVATE(n,se,j)
  do n=1,ne
!calculate elastic strains
    do j=1,6
      str_e(j,n) = strain(j,n) - strainp(j,n)
    end do

    call elastic(n,se)

    do j=1,6
      stress(j,n)=se(j)
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  

  do n=1,ne    
! calculate and transfer element forces onto nodes
! output: node forces and balance=sum(abs(force))
    call force_n_balance(n)
  end do

! include boundary conditions - force and velocity 
! calculate  velocity wih damping 
  call velocity_n_balance(dt)

!  do n=1,np    
!  print *,'  Force ',n,force(n,1),force(n,2),force(n,3)
!  end do  


! Estimation max balance-off in system (boff)
  call boff_calc(dt,vvv)
  
  return
end subroutine flac

!*****************************************************
subroutine force_n_balance(n)
  use sizes
  use node_data
  use element_data
  implicit none

  integer::n,nside,nn,j,ii,n1,n2,n3,baltoo
  real(kind=8)::x1,x2,y1,y2,z1,z2,xnorm,ynorm,znorm
  real(kind=8)::fs(4,3),f(4,3)

! Calculation forces (4 sides, one element)
! On the sides (fs(i,j), i-side, j-dimension) fi = sij*nj
  do nside = 1,4
! Three nodes of the side
    n1 = nop(n,side(1,nside))
    n2 = nop(n,side(2,nside))
    n3 = nop(n,side(3,nside))
! Two vectors of the side
    x1 = cord(n2,1) - cord(n1,1)
    y1 = cord(n2,2) - cord(n1,2)
    z1 = cord(n2,3) - cord(n1,3)

    x2 = cord(n3,1) - cord(n1,1)
    y2 = cord(n3,2) - cord(n1,2)
    z2 = cord(n3,3) - cord(n1,3)

! Normal vector to the side (NORM = X1 X X2 / 2.)
    xnorm = 0.5_8*(y1*z2 - z1*y2)
    ynorm = 0.5_8*(z1*x2 - x1*z2)
    znorm = 0.5_8*(x1*y2 - y1*x2)

    fs(nside,1) = stress(1,n)*xnorm + stress(4,n)*ynorm + stress(5,n)*znorm
    fs(nside,2) = stress(4,n)*xnorm + stress(2,n)*ynorm + stress(6,n)*znorm
    fs(nside,3) = stress(5,n)*xnorm + stress(6,n)*ynorm + stress(3,n)*znorm
 end do


! In the nodes (f(i,j), i-node,j-dimension) 
! each node belongs to three sides (see "side" above)
  do j=1,3
    f(1,j) = (fs(1,j)+fs(2,j)+fs(4,j))/3.0_8
    f(2,j) = (fs(1,j)+fs(3,j)+fs(4,j))/3.0_8
    f(3,j) = (fs(2,j)+fs(3,j)+fs(4,j))/3.0_8
    f(4,j) = (fs(1,j)+fs(2,j)+fs(3,j))/3.0_8
  end do

! Insert forces in whole system:
  do ii=1,4
    nn = nop(n,ii)
    do j=1,3
      force(nn,j) = force(nn,j) - f(ii,j)
      balance(nn,j) = balance(nn,j) + abs(f(ii,j))
!  if(nn.eq.2) print *,nn,n,j,force(nn,j),f(ii,j)
    end do
  end do
   
  return
end subroutine force_n_balance

!*****************************************************
subroutine velocity_n_balance(dt)
  use sizes
  use node_data
  use boundary_node_data
  use element_data
  implicit none

  integer::j,nb,i,ii,i3,nrd,node,nn,nd
  real(kind=8)::fs(4,3),f(4,3),fv,vvs,fff,bbb,vvv,dt
  
  ! Including volumetric force
  !$OMP PARALLEL
  !$OMP DO PRIVATE(i,i3,fv)
  do i = 1, np
    do i3 = 1,3
      fv = mass(i)*g(i3)
      force(i,i3)   = force(i,i3) + fv	! g(3) negative
      balance(i,i3) = balance(i,i3) + abs(fv)
    end do
  end do
  !$OMP END DO

  ! BOUNDARY CONDITION IN STRESSES
  !$OMP DO PRIVATE(j,nb,nn)
  do j=1,3
    do nb = 1,nbp
        nn = numbn( nb )
        bforce(nb,j)=force(nn,j)
      if( force_code(nb,j) )then
        force(nn,j)   = force(nn,j)   + value(nb,j)
        balance(nn,j) = balance(nn,j) + abs (value(nb,j))
      end if
    end do
  end do
  !$OMP END DO


  !$OMP DO PRIVATE(i,j,vvs,fff)
  ! Damping
  do i=1,np
    do j=1,3
      vvs = vel(i,j) + dt*force(i,j)/(mass(i)*den_scale)
      if ( abs(vvs) .gt. 1.0e-13_8) then
        fff = abs(force(i,j))
        force(i,j) = force(i,j) - demf(i)*sign(fff,vvs)
        vel(i,j)=vel(i,j) + dt*force(i,j)/(mass(i)*den_scale)
      end if
    end do
  end do
  !$OMP END DO

    
  !
  ! Boundary condition in velocities:
  !$OMP DO PRIVATE(i,j,node)
  do i=1,nbp 
    do j=1,3 
      if( vel_code(i,j) )then
        node = numbn(i)
        vel(node,j)= value(i,j)
        force(node,j)= 0.0_8
        balance(node,j)= 0.0_8
      end if
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  
  return
end subroutine velocity_n_balance

!*****************************************************
! Estimation max balance-off in system (boff)
subroutine boff_calc(dt,vvv)
  use sizes
  use node_data
  implicit none
  
  integer::n,j,i,noff,nmax,nn
  real(kind=8)::dt
  real(kind=8)::bbb,vvv,vmax
  
  boff = 0.0_8
  noff = 0
  bbb = 1.0e-6_8	!0.0_8
  
  do j =1,3
    do i =1,np
      bbb = max(bbb,balance(i,j))
    end do
  end do


  do j =1,3
    do i =1,np
!   print *,i,j,force(i,j),vel(i,j)
      balance(i,j)=abs (force(i,j)/bbb)
      if(balance(i,j).gt.boff ) then
        boff = balance(i,j)
        noff = i
      end if
    end do
  end do

  ! Maximum vel
  vvv = vmax (np,vel,nmax)

  ! Screen output
  write(6,200)dt,vvv,boff,noff
200 format(' dt,vmax,boff:', 3g11.3, i10)
    
  return
end subroutine boff_calc
	program fract 	!  test vol fraction

implicit none

integer, dimension(4) ::  nop

!	Working variables
integer :: i, ne, mat
real*8 :: tot,melt, frc

melt = 0.
open(111,file='elements.dat')

read(111,*) ne
tot = ne

do i = 1,ne
read(111,* )nop,mat

if ( mat .eq. 2 ) melt = melt + 1.

end do

close(111)

frc = melt/tot * 100.

print *,tot,melt,frc

stop
end program fract
! **************************************************
! minmax
 subroutine minmax(a,n,amin,amax)
  implicit none
  integer:: n,i
  real(kind=8):: a(n),amax,amin

  amin = a(1)
  amax = a(1)
  do i =2,n
    if( a(i) .lt. amin ) amin = a(i) 
    if( a(i) .gt. amax ) amax = a(i)
  end do

  return
 end subroutine minmax
 
 ! **************************************************
 ! Determination abs max velocity
 function  vmax(np,r1,nmax) 
  implicit none
  integer:: np,nmax,nnp
  real(kind=8)::r1(np,3),vm,rr,vmax
    
  vm=0.0_8
  do nnp =1,np
    rr = r1(nnp,1)**2 + r1(nnp,2)**2 + r1(nnp,3)**2
    if(rr.gt.vm) then
      vm = rr
      nmax = nnp
    end if
  end do
  vmax = sqrt ( vm )
  
  return
 end function vmax
 
subroutine input_mgm
  use sizes
  use element_data
  use node_data
  use boundary_node_data
  implicit none
  
  integer:: nn,nomat,i,j,k,nbm,ind,nlm,nmat
  integer:: nbr
  
  logical:: there
  
  real(kind=8):: rl0,rm0,vp,vs,dist,dtp,dst
  real(kind=8):: rlmax,rmumax,densmin,rl,rm
  real(kind=8):: xx,yy,zz,rr

  real(kind=8),dimension(:,:),allocatable::prop
  

 ! read the dimensions of the grid
    inquire(file='grid_size.dat',exist=there)
    if(there)then
      open(1,file="grid_size.dat")
      read(1,*) ne
      read(1,*) np
      read(1,*) nbp
      read(1,*) ddmcm
      close(1)
    else
      write(6,*)' can not find file "grid_size.dat".  Stopping '
      stop
    end if

!-----------------
        write(6,*)' Model size: '
        write(6,*)'np=',np,'   ne=',ne
        write(6,*)'min node to node distance : ',ddmcm ,' [m]'
        
        dist = ddmcm*ddmcm

!-----------------

 ! allocate arrays to the correct sizes
 ! Tetra structure 
      allocate(nop(ne,4))

 ! Boundary conditions
  allocate(numbn(nbp),vel_code(nbp,3),force_code(nbp,3),  &
  	 value(nbp,3),bforce(nbp,3))

 ! Attributes of elements
  allocate(lambda(ne),mu(ne),dens(ne),ductile(ne),src_code(ne,3),src_tm(ne), &
           stress(6,ne),strain(6,ne),strainp(6,ne),str_e(6,ne),strain0(6,ne),el_vol(ne) )
           
 ! Attributes of nodes
  allocate(mass(np),demf(np),force(np,3),             &
           balance(np,3),cord(np,3),disp(np,3),vel(np,3) )

!---------- Pressure source --------------------  
	src_code = .false.
	src_tm   = 0.0

!------	READ node coordinates & boundary condition--------
    inquire(file='coordinates.dat',exist=there)
    if(there)then
      open(1,file='coordinates.dat')
      read(1,*) nn

      if ( nn .ne. np ) then
        write(6,*)' Number of points is not correct ',nn,np
        stop
      end if
      
      		nbm = 0
      do i = 1,np
        read(1,*) cord(i,1),cord(i,2),cord(i,3),ind
 
  		vel(i,1) = 0.0_8
        vel(i,2) = 0.0_8
        vel(i,3) = 0.0_8
 
 		disp(i,1) = 0.0_8
        disp(i,2) = 0.0_8
        disp(i,3) = 0.0_8
                
                demf(i)   = 0.3_8

		if ( ind .ne. 0 ) then
		   nbm = nbm + 1	
	                numbn(nbm) = i

!---- type of boundary condition
                do k = 1,3
        vel_code(nbm,k)   = .true.
        force_code(nbm,k) = .false.

		value(nbm,k) = 0.0
   				end do
		end if

      end do
      close(1)
      write(6,*) 'end READ node coordinates'  
    else
      write(6,*)'can not find file "coordinates.dat". Stopping'
      stop
    end if

        if ( nbm .ne. nbp ) then
        print *,'  Number of boundary points is not correct ',nbm,nbp
        stop
        endif

!-----  search for minimum distance between nodes -----------
!	put value instead
!	dist = 100000.
!	do i=1,np
!		do j=1,np
!	if ( i .ne. j ) then
!	dst =       ( cord(i,1)-cord(j,1) )**2
!	dst = dst + ( cord(i,2)-cord(j,2) )**2
!	dst = dst + ( cord(i,3)-cord(j,3) )**2
!	if( dst .le. dist ) dist = dst
!	end if
!		end do
!	end do
!	dist = ddmcm*ddmcm
!
!	ddmcm = sqrt(dist)
!
!        write(6,*)'calculated node to node distance : ',ddmcm ,' [m]'

! READ material properties ------------------
!	Units:
!	density - g/m^3
!	stress  - MPa
!	time scale = sqrt(1000.) second for MPa - km
!	time scale = sqrt(1.e-3) second for MPa - m
  	
  	tsc = dsqrt(1.0e-3_8)
  write(6,*)' Numerical time scale: ',tsc
 
  inquire(file='material.dat',exist=there)
  if(there)then
    open(1,file='material.dat')
    read(1,*) nomat
    write(6,*)' number of different materials found =',nomat
  else
    write(6,*)'can not find file "material.dat". Stopping'
    stop
  end if
  
  allocate(prop(nomat,4))

  rlmax   = 0.0_8
  rmumax  = 0.0_8
  densmin = 999999.0_8
  dt_vis  = 999999.0_8

  do i = 1,nomat
	read(1,*) (prop(i,k), k=1,4)

!	Bulk, G  ---> Lambda, Mu (G - rigidity)
	rl0 = prop(i,2) - 2.0/3.0 * prop(i,3)
	rm0 = prop(i,3) 

	prop(i,2) = rl0
	prop(i,3) = rm0

!---------------------------------------------------------------------

  write(6,*) ' -------------------------------------------------'
  write(6,*)'	Material number ',i
  write(6,*)'		Density  = ',prop(i,1)
  write(6,*)'		Lambda   = ',prop(i,2)
  write(6,*)'		Mu       = ',prop(i,3)
  write(6,*)'		Fluidity = ',prop(i,4)
  write(6,*) ' -------------------------------------------------'
    
!	scale fluidity to time scale
	prop(i,4) = prop(i,4) * tsc

	if ( prop(i, 4) .ge. 1.0e-9_8 ) then
 	dtp = 0.1 / (prop(i,3) * prop(i,4) )          
	  if(dtp .le. dt_vis ) dt_vis = dtp
	end if
  write(6,*)'		Visc dt = ',dt_vis*tsc*10.0

    if ( rlmax .le. rl0  ) rlmax = rl0
	if ( rmumax .le. rm0 ) rmumax = rm0 
	if ( densmin .ge. prop(i,1) ) densmin = prop(i,1)

      vp = sqrt((rl0+2.*rm0)/prop(i,1))/tsc
      vs = sqrt(        rm0 /prop(i,1))/tsc
  write(6,*)'	Vp [m/s]      = ',vp
  write(6,*)'	Vs [m/s]      = ',vs

  end do

  close (1)

!       calculate wave-related time scale (vp2=time_scale**2)
        vp2 = ( rlmax + 2.0_8* rmumax ) / densmin / dist

!-----  end READ material properties ------------------
        write(6,*)' end read material properties '

!-----  READ element connections and put properties ----
  inquire(file='elements.dat',exist=there)
  if(there)then
    open(1,file='elements.dat')
    read(1,*) nlm
    if ( nlm .ne. ne ) then
      write(6,*)'  Number of elements is not correct ',nlm,ne
      stop
    end if
  else
    write(6,*)'can not find file "elements.dat". Stopping'
    stop
  end if
    
  do i = 1,ne
   read(1,*) nop(i,1),nop(i,2),nop(i,3),nop(i,4),nmat

    if( nmat .gt. nomat )then
      write(6,*)' Number of material types is not correct ',nmat,nomat
      stop
    end if
    
     if(nmat .eq. 3) then            ! sponge
			nmat = 1
        demf( nop(i,1) ) = 0.8_8
        demf( nop(i,2) ) = 0.8_8
        demf( nop(i,3) ) = 0.8_8
        demf( nop(i,4) ) = 0.8_8
     end if
	
	if (nmat .eq. 2) then 		! magma
	xx = (cord(nop(i,1),1)+cord(nop(i,2),1)+cord(nop(i,3),1)+cord(nop(i,4),1))/4.0
	yy = (cord(nop(i,1),2)+cord(nop(i,2),2)+cord(nop(i,3),2)+cord(nop(i,4),2))/4.0
	zz = (cord(nop(i,1),3)+cord(nop(i,2),3)+cord(nop(i,3),3)+cord(nop(i,4),3))/4.0
	
	rr = sqrt( xx*xx+yy*yy)
	
		if ( rr .le.  5.0 .and. abs(zz) .le. 0.5) src_code(i,2) = .true.

		if ( rr .le. 50.0 .and. abs(zz) .le. 0.5) then
			src_code(i,3) = .true.
		else
			nmat = 1
		end if
	end if

!  Test for correct volume (det<0)
   call volume_correct_sign(i)

!	Put material properties

        dens  (i) = prop(nmat,1)
        lambda(i) = prop(nmat,2)
        mu    (i) = prop(nmat,3)

!------------------ fluidity ----------------------
        ductile(i) = prop(nmat,4)
        
!------  Initial strain distribution ---------------- 
      
        strain(1,i) = -(787.0-700.0) / (3.*lambda(i) + 2.*mu(i)) 
        strain(2,i) = strain(1,i)
        strain(3,i) = strain(1,i)
        strain(4,i) = 0.0_8
        strain(5,i) = 0.0_8
        strain(6,i) = 0.0_8

	  strainp(1,i) = 0.0_8
	  strainp(2,i) = 0.0_8
	  strainp(3,i) = 0.0_8
      strainp(4,i) = 0.0_8
	  strainp(5,i) = 0.0_8
	  strainp(6,i) = 0.0_8

  end do
!-----  end READ element connections and put properties ----
 close(1)
  write(6,*)' end READ element connections and put properties '

  return
 end subroutine input_mgm

!==================================================================
 
subroutine volume_correct_sign(n)
  use sizes
  use element_data
  use node_data
  implicit none
  
  integer:: j,n,itemp
  
  real(kind=8):: x(4),y(4),z(4),det

          do j = 1,4
               x(j) = cord ( abs(nop(n,j)),1 )
               z(j)=  cord ( abs(nop(n,j)),3 )
               y(j) = cord ( abs(nop(n,j)),2 )
	  end do

          det=0
         det=det+(x(2)-x(1))*(y(3)-y(1))*(z(4)-z(1))
         det=det+(y(2)-y(1))*(z(3)-z(1))*(x(4)-x(1))
         det=det+(z(2)-z(1))*(x(3)-x(1))*(y(4)-y(1))

         det=det-(z(2)-z(1))*(y(3)-y(1))*(x(4)-x(1))
         det=det-(y(2)-y(1))*(x(3)-x(1))*(z(4)-z(1))
         det=det-(x(2)-x(1))*(z(3)-z(1))*(y(4)-y(1))


         if( det > 0.0_8 ) then
!---      the tetra vertices must be reordered for a NEGATIVE volume
          itemp=nop(n,1)
          nop(n,1)=nop(n,2)
          nop(n,2)=itemp
         end if

 return
end subroutine volume_correct_sign
 
!       Input files:
!               corodinates.dat - node coordinates
!               elememt.dat - element connections and properties
!               material.dat - material properties
!		grid_size.dat - size of the mesh
!-------------------------------------------------------------
!       ne   - number of elements
!       np   - number of nodes
!       nbp  - number of boundary nodes
!       nx, ny, nz size of the grid ( nx by  ny by nz  nodes)
!       nslip    number of pairs of nodes above-bellow the fault
!       ddmcm   min node-to-node distance in mcm
!-------------------------------------------------------------
program main_mgm
  use sizes
  use element_data
  use node_data
  use boundary_node_data
  implicit none
  
  integer:: nloop,nout,n
  real(kind=8)::dt,dtp,adp,s_mean
  real(kind=8):: t_start,t_stop,t_delta,tout,boff_min,boff_max,vmax

  ! Set general model parameters        
  
	time    = 0.0
	t_start = 0.0	
	t_stop  = 1.0
	t_delta = 0.005
  
  	adp = 0.01_8	!0.05_8
  
  boff_min = 1.e-3
  boff_max = 1.e-2

        g(1) = 0.
        g(2) = 0.
        g(3) = 0.

  
!  input  model geometry, material properties and initial conditions 

  call input_mgm


! time step for elasticity
      den_scale = 1.0e+4_8
  
  dt = 0.25_8 / sqrt(vp2 / den_scale)
  
  write(6,*)' Wave-related time ',sqrt(vp2)
  write(6,*)' Initial time step dt = ',dt,' density scale = ',den_scale
  
 write(6,*)' Time: ',t_start,t_stop,tout
  call node_mass
  
  
!      -------------- Loop of time ---------------------------
                time = 0.
                nloop = 0
  				nout=1
  boff = 1.

       do while ( boff .ge. 1.e-5 .or. nloop .le. 10)  		!(time .le. t_stop)
!       do while (nloop .le. 1000)

	nloop = nloop + 1
    write(6,600) nloop,time
600	format(I8,'''s ',g12.5,$)
  	time = time + dt*tsc

    call flac(dt,vmax)
    
    call move_grid(dt)

! New plastic strain
!    call plastic(dt)
	
    end do

 	  call output_mgm(nout)

!      -------------- Loop of time ---------------------------
                time = 0.
				tout = t_delta
				nout = 2
                nloop = 0
		den_scale = 1.0
        dt = 0.25_8 * sqrt(den_scale/vp2)
        print *,' Time step elastic ',dt*tsc,' Viscouse ',dt_vis*tsc
       	if( dt .ge. dt_vis ) dt = dt_vis

      open(120,file='potency.dat')

	print *,' ================================================= '
	print *,' Include source, dt = ',dt*tsc 
		src_code(:,1) = src_code(:,2)
		strain0 = strain

       do while (time .le. t_stop)

	nloop = nloop + 1
    write(6,600) nloop,time
  	time = time + dt*tsc
  	
! Nucleation
  do n = 1,ne
    if (src_code(n,3) .and. .not. src_code(n,1) ) then
!  Mean stress
      s_mean = -(stress(1,n) + stress(2,n) + stress(3,n))/3.0_8
      
      	if ( s_mean .le. 86.99 ) then
      	  src_code(n,1) = .true.
      	  src_tm(n) = time
!      	print *,' Nuc ',n,src_code(n,1),s_mean
      	end if
     end if
   end do

    call flac(dt,vmax)
    
    call move_grid(dt)

! New plastic strain
    call plastic(dt)

    call output_src(nloop)

! change time step ------------------------
!  if( boff .ge. boff_max ) then
!            den_scale = den_scale*(1.0_8-adp)
!	if( den_scale .lt. 1.0e+4_8 ) den_scale = 1.0e+4_8
!      dt = 0.25_8 * sqrt(den_scale / vp2 )
!  else if ( boff .le. boff_min ) then
!      den_scale = den_scale*(1.0_8+adp)
!       if(den_scale .ge. 1.0e+6_8) den_scale = 1.0e+6_8
!       if(den_scale .ge. 1.0e+10_8) den_scale = 1.0e+10_8
!      dt = 0.25_8 * sqrt(den_scale/vp2)
!    end if

! OUTPUT
    if(time .ge. tout) then
 	  call output_mgm(nout)
	nout = nout + 1	
 	tout = tout + t_delta
    end if
    
    end do
  
      close(120)

  stop
end program main_mgm
module sizes
! contains dimesnions and commonly used scaling factors 
  integer:: ne,np,nbp
  integer, dimension(3,4):: side / 2,1,4, &
                                   1,3,4, &
                                   3,2,4, &
                                   2,3,1 /
  real(kind=8) g(3),den_scale,boff,vp2,tsc,ddmcm,dt_vis,time
  
end module sizes
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

module element_data
  integer,dimension(:,:),allocatable:: nop
  real(kind=8),dimension(:),allocatable::lambda,mu,dens,ductile,src_tm,el_vol
  real(kind=8),dimension(:,:),allocatable::stress,strain,strainp,str_e,strain0
  logical,dimension(:,:),allocatable::src_code
    
end module element_data

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

module node_data
  real(kind=8),dimension(:),allocatable:: mass,demf
  real(kind=8),dimension(:,:),allocatable:: balance,cord,disp,vel
  real(kind=8),dimension(:,:),allocatable::force  
  
end module node_data
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
module boundary_node_data
  logical,dimension(:,:),allocatable::vel_code,force_code
  integer,dimension(:),allocatable:: numbn
  real(kind=8),dimension(:,:),allocatable::value,bforce
  
end module boundary_node_data
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!Move grid along veocity field (Lagrangian) and new total strain tensor  
subroutine move_grid(dt)
  
  use sizes
  use node_data
  use element_data

  implicit none
  
  integer:: ii,j,i,n,nn
  real(kind=8)::dt,de(6),dr(4,3)
  
  ! moving of mesh
  !$OMP PARALLEL
  !$OMP DO PRIVATE(i,j)
  do i=1,np
    do j=1,3
      cord(i,j) = cord(i,j) + vel(i,j)*dt
      disp(i,j) = disp(i,j) + vel(i,j)*dt
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  ! New elastic strain tensor
  !$OMP PARALLEL
  !$OMP DO PRIVATE(dr,n,ii,de,j,nn)
  do n = 1,ne
 ! derivation of basic functions
    call derivation(dr,n)
    
 ! Strain rate
    de = 0.0_8
    do ii=1,4
      nn = nop(n,ii)
      de(1)=de(1) + vel(nn,1)*dr(ii,1)
      de(2)=de(2) + vel(nn,2)*dr(ii,2)
      de(3)=de(3) + vel(nn,3)*dr(ii,3)
      de(4)=de(4) + 0.5_8*(vel(nn,1)*dr(ii,2)+vel(nn,2)*dr(ii,1))
      de(5)=de(5) + 0.5_8*(vel(nn,1)*dr(ii,3)+vel(nn,3)*dr(ii,1))
      de(6)=de(6) + 0.5_8*(vel(nn,2)*dr(ii,3)+vel(nn,3)*dr(ii,2))
    end do

 ! New total strains
    do j=1,6
      strain(j,n) = strain(j,n) + de(j)*dt
      !write(6,*)n,j,de(j),dt
    end do
    
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  
return
end subroutine move_grid
subroutine node_mass
  ! calculate element volumes and masses
  use sizes
  use node_data
  use element_data
  implicit none

  integer::num,i
  real(kind=8)::xn(4),yn(4),zn(4)
  real(kind=8)::vol,vola,volmin,volmax,elmass
  real(kind=8)::v1,v2,v3,v4

  mass=0.0_8 !zero mass array

  vola = 0.0_8
  volmin = 999999999999.0_8
  volmax = 0.0_8

  ! loop for all elements
  do num = 1,ne
    ! coord 4 element nodes
    do i = 1,4
      xn(i) = cord(nop(num,i),1)
      yn(i) = cord(nop(num,i),2)
      zn(i) = cord(nop(num,i),3)
    end do

    !element volume  (G = 6 * volume )
    !| 1 x1 y1 z1 |
    !| 1 x2 y2 z2 |
    !| 1 x3 y3 z3 |
    !| 1 x4 y4 z4 |
    v1= xn(2) * (yn(3)*zn(4)-zn(3)*yn(4)) &
      - xn(3) * (yn(2)*zn(4)-zn(2)*yn(4)) &
      + xn(4) * (yn(2)*zn(3)-zn(2)*yn(3)) 

    v2= xn(1) * (yn(3)*zn(4)-zn(3)*yn(4)) &
      - xn(3) * (yn(1)*zn(4)-zn(1)*yn(4)) &
      + xn(4) * (yn(1)*zn(3)-zn(1)*yn(3))

    v3= xn(1) * (yn(2)*zn(4)-zn(2)*yn(4)) &
      - xn(2) * (yn(1)*zn(4)-zn(1)*yn(4)) &
      + xn(4) * (yn(1)*zn(2)-zn(1)*yn(2))

    v4= xn(1) * (yn(2)*zn(3)-zn(2)*yn(3)) &
      - xn(2) * (yn(1)*zn(3)-zn(1)*yn(3)) &
      + xn(3) * (yn(1)*zn(2)-zn(1)*yn(2))

    vol = v1 - v2 + v3 - v4
    vol = abs(vol)/6.0_8

    vola = vola + vol

    if ( volmin .ge. vol ) volmin = vol
    if ( volmax .le. vol ) volmax = vol

    el_vol(num) = vol
    elmass = dens(num) * vol
    do i = 1,4
      mass(nop(num,i)) = mass(nop(num,i)) + 0.25_8*elmass
    end do
  end do

  vola = vola / ne
  write(6,*)' Average element volume is ',vola
  write(6,*)' Volume ranges from ',volmin,'  to ',volmax

  return
end subroutine node_mass
!----------- Output data for VIEWER-----------------
subroutine output_mgm(nstep)
  use sizes
  use element_data
  use node_data
  implicit none
  
  integer:: nstep,i
  
  character(len=14):: file1=('    .nodes.dat')
  character(len=14):: file2=('    .tetra.dat')
  character(len=4):: step1,step2
  character(len=1):: st1(4),st2(4)
  
  equivalence (step1,file1),(step1,st1(1))
  equivalence (step2,file2),(step2,st2(1))
  
  real(kind=8)::eps_v

  ! generate file-names for output
  write(step1,5) nstep
  write(step2,5) nstep
5 format(i4)

  if (st1(1) .eq. ' ' ) st1(1) = '0'
  if (st1(2) .eq. ' ' ) st1(2) = '0'
  if (st1(3) .eq. ' ' ) st1(3) = '0'

  if (st2(1) .eq. ' ' ) st2(1) = '0'
  if (st2(2) .eq. ' ' ) st2(2) = '0'
  if (st2(3) .eq. ' ' ) st2(3) = '0'

        open(79,file=file1)

       do i=1,np
	    
        write(79,*) i,				&
      		cord(i,1),cord(i,2),cord(i,3),  &
     		disp(i,1),disp(i,2),disp(i,3)
 
         end do

        close(79)

        open(79,file=file2)

       do i=1,ne

  eps_v = (strain (1,i) + strain (2,i) + strain (3,i)) - 	&
  	      (strain0(1,i) + strain0(2,i) + strain0(3,i))
	    
        write(79,*) stress(1,i),stress(2,i),stress(3,i), &
                    stress(4,i),stress(5,i),stress(6,i),eps_v, &
                    src_code(i,1),src_tm(i)

!        write(79,*) strain(1,i),strain(2,i),strain(3,i), &
!                    strain(4,i),strain(5,i),strain(6,i)
 
         end do

        close(79)

  return
end subroutine output_mgm
!----------- Output data for every step -----------------
subroutine output_src(nstep)
  use sizes
  use element_data
  use node_data
  implicit none
  
  integer:: nstep,i,k
  real(kind=8):: pot(6)


!	potency

	pot = 0.0
	do i=1,ne
		if ( src_code(i,3) ) then
			do k = 1,6
			pot(k) = pot(k) + (strain(k,i)-strain0(k,i)) * el_vol(i)
		    end do
		end if
	end do

	write (120,'(7(1x,e12.5))') time,pot

  return
end subroutine output_src

! Plastic strain accumulation
subroutine plastic(dt)
  use sizes
  use element_data
  implicit none
  
  integer n,i,j
  real(kind=8) dt,t(6),e(6),s_mean

  !$OMP PARALLEL
  !$OMP DO PRIVATE(n,s_mean,e,t,j)
  do n = 1,ne
  
!  Mean stress
      s_mean = (stress(1,n) + stress(2,n) + stress(3,n))/3.0_8
      
      	if ( src_code(n,1) ) call source(n,s_mean)
        
! Deviatoric stress
      do j=1,3
        t(j)   = stress(j,n)  - s_mean
        t(j+3) = stress(j+3,n)
      end do

      do j=1,6
        e(j) = t(j)*ductile(n)
      end do
      
! New plastic strain
      do j=1,6
        strainp(j,n) = strainp(j,n) + e(j)*dt
      end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  
  return
end subroutine plastic
subroutine source(n,s_mean)
  use sizes
  use element_data
  implicit none

  integer::n,j
  real(kind=8)::tmm,s_mean,p_source,delta,power,b
  real(kind=8)::eps_v,sm,c1,p1,p2,ro_m,P0,p_eq

		tmm = time - src_tm(n)
		sm = -s_mean + 700.0_8
		
! Final equilibrium pressure

  c1   =  	1.44E-05
  p1   = 	0.5826
  p2   =   	419
  ro_m =	2800
  P0   =	825.5
  
  eps_v = (strain (1,n) + strain (2,n) + strain (3,n)) - 	&
  	      (strain0(1,n) + strain0(2,n) + strain0(3,n))

	p_eq = (P0*ro_m*c1-eps_v*p2) / (eps_v*c1+ro_m*c1)

! Source function
		b = -187.5
		power = 2.4
	p_source = p_eq - (p_eq-787.0)*exp(b*tmm**power)  
	
	if(n .eq. 1437590) print *,' sr ',n,tmm,sm,p_source,p_eq,eps_v
	
	delta = (p_source - sm) / (3.0*lambda(n) + 2.0*mu(n))
	
	  do j=1,3
        strainp(j,n) = strainp(j,n) + delta
      end do	
	
  return
end subroutine source
