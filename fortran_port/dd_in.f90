/******************************************************************************************
* Copyright (C) 2014 Dr. Ati Moncef <a_moncef@live.fr> and contributors.
* All rights reserved.
* 
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
******************************************************************************************/

program main

  implicit none

  real ( kind = 8 ), allocatable :: phi(:,:),ddphi(:,:),dphi(:,:)
  real ( kind = 8 ), allocatable :: box(:)
  real ( kind = 8 ) ctime
  real ( kind = 8 ) ctime1
  real ( kind = 8 ) ctime2
  real ( kind = 8 ), parameter :: dt = 0.0002D+00
  real ( kind = 8 ) EneTot
  real ( kind = 8 ), allocatable :: Moment_F(:,:)
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) id
  real ( kind = 8 ) kinetic
  integer ( kind = 4 ) last
  real ( kind = 8 ), parameter :: Inertie = 1.0D+00,pi = 3.141592653589793D+00
  integer ( kind = 4 ), parameter :: nd = 1
  integer ( kind = 4 ) Nspin
  real ( kind = 8 ) potential,Temperature
  integer ( kind = 8 ) step,temps
  integer ( kind = 4 ) i,j
  integer ( kind = 4 ) step_print
  integer ( kind = 4 ) step_print_index
  integer ( kind = 4 ) step_print_num
  real ( kind = 8 ) wtime

 
OPEN(UNIT=1,FILE="in.txt",FORM="formatted",ACTION="read",STATUS="old")
open(unit=2,file="out.txt",status='replace',action='write')
!  Set initial Angles, moment angulaire, et accelerations.
! Lire l'angle initiale de chaque Spin
temps = 2000
read(1,*) Nspin

! Allocate memory.

  allocate ( ddphi(nd,Nspin))
  allocate ( box(nd) )
  allocate ( Moment_F(nd,Nspin))
  allocate ( phi(nd,Nspin))
  allocate ( dphi(nd,Nspin))


! Dimensions de Boite de simulation.

box(1:nd) = 10.0D+00

write ( *, '(a)' ) '  Initialisation des angles, moment angulaire, et accelerations.'

 call initialize ( Nspin, nd, box, phi, dphi, ddphi )

  write ( *, '(a)' ) ' Calcul de Moment de force et le Potentiel d''interaction.'

  call compute ( Nspin, nd, phi, dphi, Inertie, Moment_F, potential, kinetic, eneTot, Temperature )

  call cpu_time ( ctime1 )
  write (*, '(a)' ) '     ------    ----------  -----------------    ------------      ------------'
  write (*, '(a)' ) '#    Temps     Potentiel   Energie Cénetique    Energie Total      Température '
  write (*, '(a)' ) '     ------    ----------  -----------------    ------------      ------------'
  

do step = 0,temps

    call compute ( Nspin, nd, phi, dphi, Inertie, Moment_F, potential, kinetic, eneTot, Temperature )
  
write ( 2 , '(1x,i8,2x,g14.6,2x,g14.6,2x,g17.6,2x,g17.6)' ) step, potential, kinetic,EneTot, Temperature
write ( * , '(1x,i8,2x,g14.6,2x,g14.6,2x,g17.6,2x,g17.6)' ) step, potential, kinetic,EneTot, Temperature

    call update ( Nspin, nd, phi, dphi, Moment_F, ddphi, Inertie, dt )

  end do

  call cpu_time ( ctime2 )

  ctime = ctime2 - ctime1
 write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Temps du CPU pour cette opération de calcul est :' 
  write ( *, '(2x,g14.6,a)' ) ctime, ' seconds'

!  Free memory.
!
  deallocate ( ddphi )
  deallocate ( box )
  deallocate ( Moment_F )
  deallocate ( phi )
  deallocate ( dphi )
 
  stop
end
subroutine compute ( Nspin, nd, phi, dphi, Inertie, F, pot, kin,EneT, T )

!*****************************************************************************80

  implicit none

  integer ( kind = 4 ) Nspin,nd,i,j
  real ( kind = 8 ) theta, Inertie,kin,M_F(nd,Nspin),phi(nd,Nspin),pot,thetaij(nd,Nspin),dphi(nd,Nspin),rij,EneT,T,F_fr(nd,Nspin)
  real ( kind = 8 ) F(nd,Nspin),mu,mu0,R2,R3,ddphi(nd,Nspin)
  real*8 , parameter :: r = 12e-9,kb = 1, nu = 1,muB = 9.27400968 * 1E-24; !Bohr magneton
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00,rop = 5240 
  real ( kind = 8 ), parameter :: M = 4.5e5 ! magnitisation saturation A/S
  real ( kind = 8 ), parameter :: d = 12e-9,I0=1,dt = 0.0002D+00
mu = (M*pi*d**3)/6
mu0 = 4*pi*1e-7
pot = 0.0D+00
  do i = 1, Nspin-1

!  Calcul de l'energie potentiel et le moment de force.

    M_F(1:nd,i) = 0.0D+00

    do j = i+1, Nspin

!Angle entre deux spins

        thetaij(:,j) = phi(:,i) - phi(:,j)

! Distance entre deux spins

        rij = abs(r*(j-i))
R2=rij**2
R3=rij**3
! potentiel d'interaction

        pot  = pot + (cos(thetaij(1,j))-(3*cos(phi(1,i))*cos(phi(1,j))/R2))*(mu0*mu*mu/(4*pi*R3))

!Force de frotement

F_fr(:,j)=-(8*pi*rop*(0.02)**3)*dphi(1,j)

! Moments de Force

M_F(1,j) =M_F(1,j)-(mu0*mu*mu/(4*pi*R3))*(sin(thetaij(1,j)-3*cos(phi(1,j))*sin(phi(1,i))))

! Force Total 

F(:,j)=M_F(:,j)+F_fr(:,j)

! Caclude la vitesse angulaire
dphi(:,j) = M_F(1,j) * dt / 0.2 + (1 - M_F(1,j) / 0.2 ) * (1 - exp(- 0.2 * dt / I0)) * I0 / 0.2

!dphi(:,j) = (1/pi*(d**3)*0.2)*F(:,j)*cos(phi(1,j))*phi(1,j)

! Caclul de l'acceleration

ddphi(:,j) = (F(:,j))/I0
    
	end do
end do

!  Calcul de l'energie cénetique.

  kin = sum ( dphi(1:nd,1:Nspin)**2 )/(2*I0)
  T = kin/3*Nspin*kb
  EneT = pot + kin  
  return
end

subroutine initialize ( Nspin, nd, box, phi, dphi, ddphi )

!*****************************************************************************80

! Inisialisation des Angles, moment angulaire, et accelerations.

  implicit none

  integer ( kind = 4 ) Nspin,i,j,nd
  real ( kind = 8 ) phi(nd,Nspin),dphi(nd,Nspin),ddphi(nd,Nspin),box(nd)


do i = 1,Nspin
read(1,*) phi(:,i)
read ( 1, * ) dphi(:,i)
enddo

!Le moment angulaire et accelerations sont initialement a 0.

 return
end
!******************************************************************************


subroutine update ( Nspin, nd, phi, dphi, F, ddphi, Inertie, dt )

!*****************************************************************************80
!! UPDATE updates Angles, moment angulaire et accelerations.
!
!  Discussion:
!
!    The time integration is fully parallel.
!d
!    A moment angulaire Verlet algorithm is used for the updating.
!
!    x(t+dt) = x(t) + v(t) * dt + 0.5 * a(t) * dt * dt
!    v(t+dt) = v(t) + 0.5 * ( a(t) + a(t+dt) ) * dt
!    a(t+dt) = f(t) / m
!
  implicit none

  integer ( kind = 4 ) Nspin,nd,i,j
  real ( kind = 8 ) ddphi(nd,Nspin),phi(nd,Nspin),dphi(nd,Nspin),Inertie,F(nd,Nspin),dt

  do j = 1, Nspin
    do i = 1, nd 
      phi(:,j) = phi(1,j) + dphi(1,j) * dt + 0.5D+00 * ddphi(1,j) * dt * dt
      
     dphi(:,j) = dphi(1,j) + 0.5D+00 * dt * ddphi(1,j)
    
 end do
  end do

  return
end
