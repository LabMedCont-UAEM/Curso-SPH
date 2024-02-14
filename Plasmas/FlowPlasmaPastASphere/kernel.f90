module kernel

double precision, parameter :: EPS0=8.8541878e-12	!C/V/m, Permitividad del vacio
double precision, parameter :: QE=1.602176565e-19	!C, carga electrica del electron
double precision, parameter :: AMU=1.660538921e-27	!Kg, unidad de masa atomica
double precision, parameter :: ME=9.10938215e-31	!Kg, masa del electron
double precision, parameter :: K0=1.380648E-23		!J/K, constante de Boltzmann
double precision, parameter :: PI=3.141592653		!pi
double precision, parameter :: EvToK=QE/K0		!1 eV in K - 11604

type :: vec3d
	double precision :: x,y,z
end type

type :: vec3i
	integer :: i,j,k
end type

type, extends(vec3d) :: particle
	double precision :: u,v,w,mpw,efx,efy,efz
end type

type :: sphere
	type(vec3d) :: center
	double precision :: r,pot
end type

type :: wld
	type(vec3d) :: x0,xl,dd
	type(vec3i) :: nn
	double precision :: dt,tol,n0,T0,phi0
	integer :: Nmax,max_it,nbody,ntotal
	
	integer, dimension(:,:,:), allocatable :: flagObs
	type(particle), dimension(:), allocatable :: p1,p2 
	
	contains
	
	procedure :: initFlagObs
	procedure :: destroyFlagObs
	!-------------------------
	procedure :: initScalar
	procedure :: destroyScalar
	procedure :: initVecField
	procedure :: destroyVecField
	procedure :: solvePot
	procedure :: addPotentialSphere
	procedure :: coldBeam
	procedure :: integrateVelocity
	procedure :: integratePosition
	procedure :: removeParticles
	procedure :: scatter
	procedure :: gather
	procedure :: nodeVolCalc
	procedure :: computeEF
end type

contains

subroutine initFlagObs(this)
class(wld) :: this
allocate(this%flagObs( this%nn%i,this%nn%j,this%nn%k))
flagObs=0
end subroutine

subroutine destroyFlagObs(this)
class(wld) :: this
deallocate(this%flagObs)
end subroutine

subroutine initScalar(this,f)
	implicit none
	class(wld) :: this
	double precision, dimension(:,:,:), allocatable :: f

	ALLOCATE(f(this%nn%i,this%nn%j,this%nn%k))

	f=0.0d0

	return
end subroutine

subroutine destroyScalar(this,f)
	implicit none
	class(wld) :: this
	double precision, dimension(:,:,:), allocatable :: f

	DEALLOCATE(f)
	
end subroutine

subroutine initVecField(this,f)
	implicit none
	class(wld) :: this
	type(vec3d), dimension(:,:,:), allocatable :: f

	ALLOCATE(f(this%nn%i,this%nn%j,this%nn%k))

	f=vec3d(0.0d0,0.0d0,0.0d0)
	
	return
end subroutine

subroutine destroyVecField(this,f)
	implicit none
	class(wld) :: this
	type(vec3d), dimension(:,:,:), allocatable :: f

	DEALLOCATE(f)
	
end subroutine

subroutine solvePot(this,phi,rho)
implicit none
class(wld) :: this
double precision, dimension(:,:,:), allocatable :: phi,rho

integer :: i,j,k,n
double precision :: idx2,idy2,idz2,phi_new,summ,R,L2,ne

	idx2=1.0/(this%dd%x*this%dd%x)
	idy2=1.0/(this%dd%y*this%dd%y)
	idz2=1.0/(this%dd%z*this%dd%z)
	
	do n=1,this%max_it
		
		
		
		do i=1,this%nn%i
		do j=1,this%nn%j
		do k=1,this%nn%k
		if( this%flagObs(i,j,k) == 0) then
		
		if(i==1) then
			phi(i,j,k)=phi(i+1,j,k)
		else if(i==this%nn%i) then
			phi(i,j,k)=phi(i-1,j,k)
		else if(j==1) then
			phi(i,j,k)=phi(i,j+1,k)
		else if(j==this%nn%j) then
			phi(i,j,k)=phi(i,j-1,k)
		else if(k==1) then
			phi(i,j,k)=phi(i,j,k+1)
		else if(k==this%nn%k) then
			phi(i,j,k)=phi(i,j,k-1)
		else
			ne=this%n0*exp( (phi(i,j,k)-this%phi0)/this%T0)	
			phi_new=( (rho(i,j,k) - QE*ne)/EPS0 &
				+ idx2*( phi(i-1,j,k) + phi(i+1,j,k) ) &
				+ idy2*( phi(i,j-1,k) + phi(i,j+1,k) ) &
				+ idz2*( phi(i,j,k-1) + phi(i,j,k+1) ) )/&
				( 2.0*idx2 + 2.0*idy2 + 2.0*idz2 )
				
			phi(i,j,k)=phi(i,j,k)+1.4*(phi_new - phi(i,j,k))
		end if
		
		end if
		enddo
		enddo
		enddo
	
		if( mod(n,25) == 0 ) then
	
		summ=0.0d0
		
		do i=2,this%nn%i-1
		do j=2,this%nn%j-1
		do k=2,this%nn%k-1
		if( this%flagObs(i,j,k) == 0) then
		if(i==1) then
			R=phi(i,j,k)-phi(i+1,j,k)
		else if(i==this%nn%i) then
			R=phi(i,j,k)-phi(i-1,j,k)
		else if(j==1) then
			R=phi(i,j,k)-phi(i,j+1,k)
		else if(j==this%nn%j) then
			R=phi(i,j,k)-phi(i,j-1,k)
		else if(k==1) then
			R=phi(i,j,k)-phi(i,j,k+1)
		else if(k==this%nn%k) then
			R=phi(i,j,k)-phi(i,j,k-1)
		else
		ne=this%n0*exp( (phi(i,j,k)-this%phi0)/this%T0)		
		R=-phi(i,j,k)*( 2.0*idx2 + 2.0*idy2 + 2.0*idz2  ) &
		+ (rho(i,j,k)-QE*ne)/EPS0 &
		+ idx2*( phi(i-1,j,k) + phi(i+1,j,k) ) &
		+ idy2*( phi(i,j-1,k) + phi(i,j+1,k) ) &
		+ idz2*( phi(i,j,k-1) + phi(i,j,k+1) ) 
		summ=summ+R*R
		end if
		end if
		enddo
		enddo
		enddo
		
		L2=sqrt(summ/(this%nn%i*this%nn%j*this%nn%k))
		
		if(L2 < this%tol) then
			write(*,*) '** Converged reached **',n
			exit
		end if
		
		end if
	enddo
	
end subroutine

subroutine addPotentialSphere(this,sph1,phi)
	implicit none
	class(wld) :: this
	type(sphere) :: sph1
	integer :: i,j,k
	double precision :: xcp,ycp,zcp,dist
	double precision, dimension(:,:,:), allocatable :: phi
	
	do i=2,this%nn%i-1
	do j=2,this%nn%j-1
	do k=2,this%nn%k-1
	xcp=(i-1)*this%dd%x+this%x0%x
	ycp=(j-1)*this%dd%y+this%x0%y
	zcp=(k-1)*this%dd%z+this%x0%z
	
	dist=sqrt( (sph1%center%x-xcp)**2 + (sph1%center%y-ycp)**2 + (sph1%center%z-zcp)**2 )
	
	if( dist <= sph1%r) then
		this%flagObs(i,j,k)=2
		phi(i,j,k)=sph1%pot
	end if
	
	enddo
	enddo
	enddo
	
	return
end subroutine

subroutine scatter(this,f)
	implicit none
	class(wld) :: this
	double precision, dimension(:,:,:), allocatable :: f
	double precision :: lx,ly,lz,di,dj,dk,val
	integer :: i,j,k,n

	do n=1,this%nbody
	lx=(this%p1(n)%x-this%x0%x)/this%dd%x 
	ly=(this%p1(n)%y-this%x0%y)/this%dd%y
	lz=(this%p1(n)%z-this%x0%z)/this%dd%z

	i=int(lx)
	j=int(ly)
	k=int(lz)

	di=lx-i
	dj=ly-j
	dk=lz-k

	i=i+1
	j=j+1
	k=k+1

	val=this%p1(n)%mpw
	
	f(i,j,k)=f(i,j,k)+val*(1-di)*(1-dj)*(1-dk)
	f(i+1,j,k)=f(i+1,j,k)+val*(di)*(1-dj)*(1-dk)
	f(i+1,j+1,k)=f(i+1,j+1,k)+val*(di)*(dj)*(1-dk)
	f(i,j+1,k)=f(i,j+1,k)+val*(1-di)*(dj)*(1-dk)
	f(i,j,k+1)=f(i,j,k+1)+val*(1-di)*(1-dj)*(dk)
	f(i+1,j,k+1)=f(i+1,j,k+1)+val*(di)*(1-dj)*(dk)
	f(i+1,j+1,k+1)=f(i+1,j+1,k+1)+val*(di)*(dj)*(dk)
	f(i,j+1,k+1)=f(i,j+1,k+1)+val*(1-di)*(dj)*(dk)
	enddo
	
	return
end subroutine

subroutine gather(this,Tval,f)
	implicit none
	class(wld) :: this
	double precision, dimension(:,:,:) :: f
	double precision, dimension(:) :: Tval
	double precision :: lx,ly,lz,di,dj,dk
	integer :: i,j,k,n

	do n=1,this%nbody
	lx=(this%p1(n)%x-this%x0%x)/this%dd%x 
	ly=(this%p1(n)%y-this%x0%y)/this%dd%y
	lz=(this%p1(n)%z-this%x0%z)/this%dd%z

	i=int(lx)
	j=int(ly)
	k=int(lz)

	di=lx-i
	dj=ly-j
	dk=lz-k

	i=i+1
	j=j+1
	k=k+1

	
	Tval(n) = f(i,j,k)*(1-di)*(1-dj)*(1-dk) &
		+f(i+1,j,k)*(di)*(1-dj)*(1-dk) &
		+f(i+1,j+1,k)*(di)*(dj)*(1-dk) &
		+f(i,j+1,k)*(1-di)*(dj)*(1-dk) &
		+f(i,j,k+1)*(1-di)*(1-dj)*(dk) &
		+f(i+1,j,k+1)*(di)*(1-dj)*(dk) &
		+f(i+1,j+1,k+1)*(di)*(dj)*(dk) &
		+f(i,j+1,k+1)*(1-di)*(dj)*(dk)
	enddo
	
	return
end subroutine

!----------------------------------------------------------------

subroutine coldBeam(this,speed,density,mpw0)
	implicit none
	class(wld) :: this
	real :: rdn1,rdn2
	integer :: m,current_nbody,i,num_sim
	double precision :: speed,num_real,density,mpw0
	
	
	
	num_real=(density)*( (this%xl%x-this%x0%x)*(this%xl%y-this%x0%y) )&
		*(this%dt)*( speed )
	num_sim=nint( num_real/mpw0 )
	
	this%nbody=this%nbody+num_sim
	
	if(this%nbody > this%ntotal) then
	DEALLOCATE(this%p2)
	
	this%ntotal=int( this%ntotal*2 )
	print*,'* Allocating new array with...*',this%ntotal
	ALLOCATE(this%p2(1:this%ntotal))
	
	this%p2(1:size(this%p1))=this%p1(1:size( this%p1))
	DEALLOCATE(this%p1)
	ALLOCATE(this%p1(1:this%ntotal))
	this%p1=this%p2
	end if
	
	do i=this%nbody-num_sim + 1,this%nbody
	call RANDOM_NUMBER(rdn1)
	call RANDOM_NUMBER(rdn2)
	this%p1(i)%x=this%x0%x+rdn1*(this%xl%x-this%x0%x)
	this%p1(i)%y=this%x0%y+rdn2*(this%xl%y-this%x0%y)
	this%p1(i)%z=this%x0%z
	
	this%p1(i)%w=speed
	this%p2(i)%w=this%p1(i)%w
	
	this%p1(i)%mpw=mpw0
	
	this%p2(i)=this%p1(i)
	enddo
	

	
	return
end subroutine

subroutine removeParticles(this)
implicit none

	class(wld) :: this
	integer :: i,m,ii,jj,kk
	
	i=0
	do while(i <= this%nbody)
	i=i+1
	ii=int( (this%p2(i)%x-this%x0%x)/this%dd%x ) + 1
	jj=int( (this%p2(i)%y-this%x0%y)/this%dd%y ) + 1
	kk=int( (this%p2(i)%z-this%x0%z)/this%dd%z ) + 1
	
	if( (this%p2(i)%z > this%xl%z) .or.( this%flagObs(ii,jj,kk)==2 ) ) then
	this%p2(i)=this%p2(this%nbody)
	this%p2(this%nbody)=particle(0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
					0.0,0.0,0.0)
	i=i-1
	this%nbody=this%nbody-1
	end if
	enddo
return
end subroutine

subroutine integratePosition(this)
implicit none
	class(wld) :: this
	integer :: i
	
	do i=1,this%nbody
	
	this%p2(i)%x = this%p1(i)%x + this%dt*( this%p1(i)%u )
	this%p2(i)%y = this%p1(i)%y + this%dt*( this%p1(i)%v )
	this%p2(i)%z = this%p1(i)%z + this%dt*( this%p1(i)%w )
	
	if( (this%p2(i)%x < this%x0%x ) ) then 
		this%p2(i)%x=2.0*this%x0%x-this%p2(i)%x
		this%p2(i)%u=-1.0*this%p2(i)%u
	else if( (this%p2(i)%x > this%xl%x) ) then
		this%p2(i)%x=2.0*this%xl%x-this%p2(i)%x
		this%p2(i)%u=-1.0*this%p2(i)%u
	end if
	
	if( (this%p2(i)%y < this%x0%y ) ) then 
		this%p2(i)%y=2.0*this%x0%y-this%p2(i)%y
		this%p2(i)%v=-1.0*this%p2(i)%v
	else if( (this%p2(i)%y > this%xl%y) ) then
		this%p2(i)%y=2.0*this%xl%y-this%p2(i)%y
		this%p2(i)%v=-1.0*this%p2(i)%v
	end if
		
	enddo
return
end subroutine

subroutine integrateVelocity(this)
	implicit none
	class(wld) :: this
	integer :: i
	
	do i=1,this%nbody
	this%p2(i)%u=this%p1(i)%u + this%dt*(QE/(16.0*AMU))*this%p1(i)%efx
	this%p2(i)%v=this%p1(i)%v + this%dt*(QE/(16.0*AMU))*this%p1(i)%efy
	this%p2(i)%w=this%p1(i)%w + this%dt*(QE/(16.0*AMU))*this%p1(i)%efz
	enddo
	
	return
end subroutine

subroutine nodeVolCalc( this, f )
	implicit none
	class(wld) :: this
	double precision, dimension(:,:,:), allocatable :: f
	double precision :: V
	integer :: i,j,k

	do i=1,this%nn%i
	do j=1,this%nn%j
	do k=1,this%nn%k
		V=this%dd%x*this%dd%y*this%dd%z
		if( (i==1).or.(i==this%nn%i) ) V=V*0.5
		if( (j==1).or.(j==this%nn%j) ) V=V*0.5
		if( (k==1).or.(k==this%nn%k) ) V=V*0.5
		f(i,j,k)=V
	enddo
	enddo
	enddo
	return
end subroutine

subroutine computeEF(this,ef,phi)
	implicit none
	class(wld) :: this
	type(vec3d), dimension(:,:,:), allocatable :: ef
	double precision, dimension(:,:,:), allocatable :: phi
	double precision :: dx,dy,dz
	integer :: i,j,k,ni,nj,nk
	
	dx=this%dd%x
	dy=this%dd%y
	dz=this%dd%z
	
	ni=this%nn%i
	nj=this%nn%j
	nk=this%nn%k
	
	do i=1,ni
	do j=1,nj
	do k=1,nk
	
	
		if(i==1) then
			ef(i,j,k)%x=-( -3.0*phi(i,j,k)+4.0*phi(i+1,j,k) -phi(i+2,j,k))/(2.0*dx)
		else if(i==ni) then
			ef(i,j,k)%x=-(phi(i-2,j,k)-4.0*phi(i-1,j,k)+3.0*phi(i,j,k))/(2.0*dx)
		else
			ef(i,j,k)%x=-(phi(i+1,j,k)-phi(i-1,j,k))/(2.0*dx)
		end if
		
		if(j==1) then
			ef(i,j,k)%y=-( -3.0*phi(i,j,k)+4.0*phi(i,j+1,k) -phi(i,j+2,k))/(2.0*dy)
		else if(j==nj) then
			ef(i,j,k)%y=-(phi(i,j-2,k)-4.0*phi(i,j-1,k)+3.0*phi(i,j,k))/(2.0*dy)
		else
			ef(i,j,k)%y=-(phi(i,j+1,k)-phi(i,j-1,k))/(2.0*dy)
		end if
		
		if(k==1) then
			ef(i,j,k)%z=-( -3.0*phi(i,j,k)+4.0*phi(i,j,k+1) -phi(i,j,k+2))/(2.0*dz)
		else if(k==nk) then
			ef(i,j,k)%z=-(phi(i,j,k-2)-4.0*phi(i,j,k-1)+3.0*phi(i,j,k))/(2.0*dz)
		else
			ef(i,j,k)%z=-(phi(i,j,k+1)-phi(i,j,k-1))/(2.0*dz)
		end if

	enddo
	enddo
	enddo
	
	return
end subroutine
!----------------------------------------------------------------


subroutine animationParticles(world,t,nprint,NAMEE)
implicit none
type(wld) :: world
integer :: t,nprint,i,NCOUNT


character :: nline*1,EXT*4,DESTINY*64,FNAME*64,OUTFILE*128,NCTEXT*4,NPOINTSTEXT*7

character(*) :: NAMEE
integer(8) :: offset(2),ii

EXT='.vtk'
DESTINY='anim/'
FNAME=trim(NAMEE)
NCOUNT=t/nprint
WRITE(NCTEXT,'(I4.4)') NCOUNT 
OUTFILE=trim(DESTINY)//trim(FNAME)//trim(NCTEXT)//EXT

WRITE(NPOINTSTEXT,'(I0)') (world%nbody)

offset(1)=0
offset(2)=(world%nbody)

	OPEN(10, file=OUTFILE, form="unformatted", access='stream', status='replace',convert="big_endian")

	WRITE(10) '# vtk DataFile Version 5.1',new_line(nline)
	WRITE(10) 'Really cool data',new_line(nline)
	WRITE(10) 'BINARY',new_line(nline)
	WRITE(10) 'DATASET POLYDATA',new_line(nline)
	WRITE(10) 'POINTS '//trim(NPOINTSTEXT)//' double',new_line(nline)
	
	do i=1,world%nbody
	WRITE(10) world%p1(i)%X,world%p1(i)%Y,world%p1(i)%Z
	enddo

	WRITE(10) 'VERTICES 2 '//trim(NPOINTSTEXT),new_line(nline)
	WRITE(10) 'OFFSETS vtktypeint64',new_line(nline)
	WRITE(10) offset
	
	WRITE(10) 'CONNECTIVITY vtktypeint64',new_line(nline)
	

	do ii=0,(world%nbody-1)
	WRITE(10) ii
	enddo
	
	

close(10)
	
return
end subroutine


subroutine animation(t,nplot,&
			this,phi,rho,den_ions,ef) 
implicit none
type(wld) :: this
double precision, dimension(:,:,:), allocatable :: phi,rho,den_ions
type(vec3d), dimension(:,:,:), allocatable :: ef
integer :: i,j,k,ncount,t,nplot
character ::  NCTEXT*4,FNAME*16,EXT*4,OUTFILE*256

EXT='.dat'
FNAME='TIME_'
NCOUNT=t/Nplot
WRITE(NCTEXT,'(I4.4)') NCOUNT

OUTFILE='anim/'//trim(FNAME)//NCTEXT//EXT

	100 format(9(E14.6))
	
	open(1,file=OUTFILE)
		
		write(1,*) 'TITLE="PLASMA3D"'
		write(1,*) 'VARIABLES="X","Y","Z","PHI","RHO","dn.O+","efx","efy","efz"'
		write(1,*) 'ZONE T="1", I=',this%nn%i,", J=",this%nn%j,", K=",this%nn%k
		
		do k=1,this%nn%k
		do j=1,this%nn%j
		do i=1,this%nn%i
			write(1,100) this%x0%x + (i-1)*this%dd%x, this%x0%y + (j-1)*this%dd%y, this%x0%z + (k-1)*this%dd%z,&
				phi(i,j,k),rho(i,j,k),den_ions(i,j,k),&
				ef(i,j,k)%x,ef(i,j,k)%y,ef(i,j,k)%z
		enddo
		enddo
		enddo
		
	close(1)
return
end subroutine


end module
