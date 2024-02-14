program PlasmaPastAPotSphere
use kernel
implicit none

type(wld) :: world
type(vec3d) :: x0,xl
type(vec3i) :: nn

type(sphere) :: sph1

double precision, dimension(:,:,:), allocatable :: phi,rho,den_ions,node_vol
type(vec3d), dimension(:,:,:), allocatable :: ef


integer :: t

call RANDOM_SEED()

x0=vec3d(-0.1d0,-0.1d0,0.0d0)
xl=vec3d(0.1d0,0.1d0,0.4d0)
nn=vec3i(21,21,41)

world%x0=x0
world%xl=xl
world%nn=nn

world%dd=vec3d( (xl%x-x0%x)/(nn%i-1) , &
		(xl%y-x0%y)/(nn%j-1) , &
		(xl%z-x0%z)/(nn%k-1) )


world%dt=1.0d-7
world%Nmax=400

world%tol=1.0d-4
world%max_it=20000

world%n0=1.0d10
world%T0=1.5d0
world%phi0=0.0d0

world%nbody=0
world%ntotal=400000
ALLOCATE(world%p1(1:world%ntotal),&
	 world%p2(1:world%ntotal))
	 
!-------------------------------------------------
sph1=sphere(vec3d(0.d0,0.0d0,0.15d0),0.05d0,-100.0d0)

call world%initFlagObs()
call world%initScalar(phi)
call world%initScalar(rho)
call world%initScalar(den_ions)
call world%initScalar(node_vol)
call world%initVecField(ef)

world%flagObs(:,:,1)=1
call world%addPotentialSphere(sph1,phi)

call world%nodeVolCalc( node_vol )

do t=1,world%Nmax
	rho=0.0d0
	den_ions=0.0d0
	
	call world%coldBeam(7.0d3,1.0d10,5.0d1)
	call world%scatter(den_ions)
	den_ions=den_ions/node_vol
	rho=rho+den_ions*QE
	
	call world%solvePot(phi,rho)
	call world%computeEF(ef,phi)
	call world%gather( world%p1%efx , ef%x )
	call world%gather( world%p1%efy , ef%y )
	call world%gather( world%p1%efz , ef%z )
	
	
	
	call animation(t,1,world,phi,rho,den_ions,ef)
	
	call world%integrateVelocity()
	call world%integratePosition()
	call world%removeParticles()
	
	print*, '# part,',world%nbody,t
	call animationParticles(world,t,1,'particles_')
		
	world%p1=world%p2
enddo

!--------------------------------------------------
call world%destroyFlagObs()
call world%destroyScalar(phi)
call world%destroyScalar(rho)
call world%destroyScalar(den_ions)
call world%destroyScalar(node_vol)
call world%destroyVecField(ef)

DEALLOCATE(world%p1,world%p2)
end program
