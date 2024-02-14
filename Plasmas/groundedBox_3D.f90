program plasma2d
implicit none

double precision, parameter :: EPS0=8.8541878e-12	!C/V/m, Permitividad del vacio
double precision, parameter :: QE=1.602176565e-19	!C, carga electrica del electron
double precision, parameter :: AMU=1.660538921e-27	!Kg, unidad de masa atomica
double precision, parameter :: ME=9.10938215e-31	!Kg, masa del electron
double precision, parameter :: K0=1.380648E-23		!J/K, constante de Boltzmann
double precision, parameter :: PI=3.141592653		!pi
double precision, parameter :: EvToK=QE/K0		!1 eV in K - 11604

real :: xrd,yrd,zrd
double precision :: x0,y0,z0,xm,ym,zm,x0e,y0e,z0e,xme,yme,zme,&
	dx,dy,dz,tol,v,&
	dt,Eke,Eki,PE
double precision, dimension(:), allocatable :: xi,yi,zi,ui,vi,wi,mpwi,xe,ye,ze,ue,ve,we,mpwe,&
						efxe,efye,efze,efxi,efyi,efzi
double precision, dimension(:,:,:), allocatable :: phi,rho,efx,efy,efz,den_eles,den_ions,node_vol

integer :: i,j,k,nd_ions,nd_eles,&
	ni,nj,nk,max_it,t

!-------------------------------------------------------------------| Particles info |
CALL RANDOM_SEED()

dt=1.0d-10

x0=-0.1d0
y0=-0.1d0
z0=0.0d0
xm=0.1d0
ym=0.1d0
zm=0.2d0

x0e=-0.1d0
y0e=-0.1d0
z0e=0.0
xme=0.0d0
yme=0.0d0
zme=0.1d0


nd_ions=100000
nd_eles=50000

allocate( efxi(nd_ions),efyi(nd_ions),efzi(nd_ions),&
	xi(nd_ions),yi(nd_ions),zi(nd_ions),ui(nd_ions),vi(nd_ions),wi(nd_ions),mpwi(nd_ions),&
	efxe(nd_eles),efye(nd_eles),efze(nd_eles),&
	xe(nd_eles),ye(nd_eles),ze(nd_eles),ue(nd_eles),ve(nd_eles),we(nd_eles),mpwe(nd_eles) )
	

do i=1,nd_eles
	call RANDOM_NUMBER(xrd)
	xe(i)=x0e+xrd*(xme-x0e)
	
	call RANDOM_NUMBER(yrd)
	ye(i)=y0e+yrd*(yme-y0e)
	
	call RANDOM_NUMBER(zrd)
	ze(i)=z0e+zrd*(zme-z0e)
enddo

do i=1,nd_ions
	call RANDOM_NUMBER(xrd)
	xi(i)=x0+xrd*(xm-x0)
	
	call RANDOM_NUMBER(yrd)
	yi(i)=y0+yrd*(ym-y0)
	
	call RANDOM_NUMBER(zrd)
	zi(i)=z0+zrd*(zm-z0)
enddo

ui=0.0d0
vi=0.0d0
wi=0.0d0
ue=0.0d0
ve=0.0d0
we=0.0d0
!-----------------------------------------------------------------------------------------

ni=30
nj=30
nk=30

dx=(xm-x0)/(ni-1)
dy=(ym-y0)/(nj-1)
dz=(zm-z0)/(nk-1)

tol=1.0d-4
max_it=10000

allocate(phi(ni,nj,nk),rho(ni,nj,nk),efx(ni,nj,nk),efy(ni,nj,nk),efz(ni,nj,nk),&
	den_eles(ni,nj,nk),den_ions(ni,nj,nk),node_vol(ni,nj,nk))

!phi(1,:)=1.0d0
!phi(:,1)=2.0d0


	
mpwe=(xme-x0e)*(yme-y0e)*(zme-z0e)*(1.0d11)/(nd_eles)
mpwi=(xm-x0)*(ym-y0)*(zm-z0)*(1.0d11)/(nd_ions)


do i=1,ni
do j=1,nj
do k=1,nk
	V=dx*dy*dz
	if( (i==1).or.(i==ni) ) V=V*0.5
	if( (j==1).or.(j==nj) ) V=V*0.5
	if( (k==1).or.(k==nk) ) V=V*0.5
	node_vol(i,j,k)=V
enddo
enddo
enddo

open(100,file='diagnostic.dat')

do t=1,20000

	rho=0.0d0
	den_ions=0.0d0
	den_eles=0.0d0
	
	do i=1,nd_ions
		call scatter(den_ions,mpwi(i),&
			xi(i),x0,dx,yi(i),y0,dy,zi(i),z0,dz,&
			ni,nj,nk)
	enddo

	den_ions=den_ions/node_vol

	do i=1,nd_eles
		call scatter(den_eles,mpwe(i),&
				xe(i),x0e,dx,ye(i),y0e,dy,ze(i),z0e,dz,&
				ni,nj,nk)
	enddo

	den_eles=den_eles/node_vol

	rho=rho+den_ions*QE
	rho=rho+den_eles*(-1.0*QE)
		

	call solve_pot(dx,dy,dz,max_it,tol,ni,nj,nk,phi,rho,&
			EPS0)	

	call computeEF(ni,nj,nk,dx,dy,dz,efx,efy,efz,phi)

	do i=1,nd_eles

		!-------------------------------------| Advance for eles |
		call gather(efxe(i),efx,&
			xe(i),x0,dx,ye(i),y0,dy,ze(i),z0,dz,&
				ni,nj,nk)
				
		call gather(efye(i),efy,&
			xe(i),x0,dx,ye(i),y0,dy,ze(i),z0,dz,&
				ni,nj,nk)
				
		call gather(efze(i),efz,&
			xe(i),x0,dx,ye(i),y0,dy,ze(i),z0,dz,&
				ni,nj,nk)
				
		ue(i)=ue(i)+(-1.0*QE/ME)*efxe(i)*dt
		ve(i)=ve(i)+(-1.0*QE/ME)*efye(i)*dt
		we(i)=we(i)+(-1.0*QE/ME)*efze(i)*dt
		
		xe(i)=xe(i)+ue(i)*dt
		ye(i)=ye(i)+ve(i)*dt
		ze(i)=ze(i)+we(i)*dt
		
		if( (xe(i) < x0) ) then 
			xe(i)=2.0*x0-xe(i)
			ue(i)=-1.0*ue(i)
		else if( (xe(i) > xm) ) then
			xe(i)=2.0*xm-xe(i)
			ue(i)=-1.0*ue(i)
		end if
		
		if( (ye(i) < y0) ) then
			ye(i)=2.0*y0-ye(i)
			ve(i)=-1.0*ve(i)
		else if ( ( ye(i) > ym) ) then
			ye(i)=2.0*ym-ye(i)
			ve(i)=-1.0*ve(i)
		end if
		
		if( (ze(i) < z0) ) then
			ze(i)=2.0*z0-ze(i)
			we(i)=-1.0*we(i)
		else if ( ( ze(i) > zm) ) then
			ze(i)=2.0*zm-ze(i)
			we(i)=-1.0*we(i)
		end if
	enddo
	
	do i=1,nd_ions	
		!-------------------------------------| Advance for ions |
		call gather(efxi(i),efx,&
			xi(i),x0,dx,yi(i),y0,dy,zi(i),z0,dz,&
				ni,nj,nk)
				
		call gather(efyi(i),efy,&
			xi(i),x0,dx,yi(i),y0,dy,zi(i),z0,dz,&
				ni,nj,nk)	
		
		call gather(efzi(i),efz,&
			xi(i),x0,dx,yi(i),y0,dy,zi(i),z0,dz,&
				ni,nj,nk)
						
		
		ui(i)=ui(i)+(QE/(16.0*AMU))*efxi(i)*dt
		vi(i)=vi(i)+(QE/(16.0*AMU))*efyi(i)*dt
		wi(i)=wi(i)+(QE/(16.0*AMU))*efzi(i)*dt
		
		xi(i)=xi(i)+ui(i)*dt
		yi(i)=yi(i)+vi(i)*dt
		zi(i)=zi(i)+wi(i)*dt
		
		if( (xi(i) < x0) ) then 
			xi(i)=2.0*x0-xi(i)
			ui(i)=-1.0*ui(i)
		else if( (xi(i) > xm) ) then
			xi(i)=2.0*xm-xi(i)
			ui(i)=-1.0*ui(i)
		end if
		
		if( (yi(i) < y0) ) then
			yi(i)=2.0*y0-yi(i)
			vi(i)=-1.0*vi(i)
		else if ( ( yi(i) > ym) ) then
			yi(i)=2.0*ym-yi(i)
			vi(i)=-1.0*vi(i)
		end if
		
		if( (zi(i) < z0) ) then
			zi(i)=2.0*z0-zi(i)
			wi(i)=-1.0*wi(i)
		else if ( ( zi(i) > zm) ) then
			zi(i)=2.0*zm-zi(i)
			wi(i)=-1.0*wi(i)
		end if
	enddo


	if(mod(t,100).eq.0.0) then
	
		PE=0.0d0
		Eki=0.0d0
		Eke=0.0d0
		
		call diagnostic(Eke,Eki,PE,EPS0,AMU,ME,&
				nd_eles,nd_ions,mpwi,mpwe,ue,ve,we,ui,vi,wi,&
				ni,nj,nk,node_vol,efx,efy,efz)
		
		print*, t,PE+Eke+Eki
		
		10 format(I0,4(E14.6))
		write(100,10) t,PE,Eke,Eki,PE+Eke+Eki
		
		call animation(t,100,&
			x0,y0,z0,dx,dy,dz,&
			phi,den_ions,den_eles,rho,&
			ni,nj,nk)
	end if
enddo

close(100)

deallocate(phi,rho,efx,efy,&
	den_eles,den_ions,node_vol)

deallocate( efxi,efyi,efzi,&
	xi,yi,zi,ui,vi,wi,mpwi,&
	efxe,efye,efze,&
	xe,ye,ze,ue,ve,we,mpwe )
	
end program

subroutine solve_pot(dx,dy,dz,max_it,tol,ni,nj,nk,phi,rho,&
		EPS0)
implicit none

integer :: i,j,k,n,max_it,ni,nj,nk
double precision :: dx,dy,dz,idx2,idy2,idz2,phi_new,summ,R,EPS0,L2,tol
double precision :: phi(ni,nj,nk),rho(ni,nj,nk)

	idx2=1.0/(dx*dx)
	idy2=1.0/(dy*dy)
	idz2=1.0/(dz*dz)

	do n=1,max_it
		
		do i=2,ni-1
			do j=2,nj-1
				do k=2,nk-1
					phi_new=( rho(i,j,k)/EPS0 &
						+ idx2*( phi(i-1,j,k) + phi(i+1,j,k) ) &
						+ idy2*( phi(i,j-1,k) + phi(i,j+1,k) ) &
						+ idz2*( phi(i,j,k-1) + phi(i,j,k+1) ) )/&
						( 2.0*idx2 + 2.0*idy2 + 2.0*idz2 )
						
					phi(i,j,k)=phi(i,j,k)+1.4*(phi_new - phi(i,j,k))
				enddo
			enddo
		enddo
	
		if( mod(n,25) == 0 ) then
		
			summ=0.0d0
			
			
			do i=2,ni-1
				do j=2,nj-1
					do k=2,nk-1
					
						R=-phi(i,j,k)*( 2.0*idx2 + 2.0*idy2 + 2.0*idz2  ) &
						+ rho(i,j,k)/EPS0 &
						+ idx2*( phi(i-1,j,k) + phi(i+1,j,k) ) &
						+ idy2*( phi(i,j-1,k) + phi(i,j+1,k) ) &
						+ idz2*( phi(i,j,k-1) + phi(i,j,k+1) ) 
						summ=summ+R*R
					enddo
				enddo
			enddo
			
			L2=sqrt(summ/(ni*nj))
			
			if(L2 < tol) then
			!	write(*,*) '** Converged reached **',n
				exit
			end if
		
		end if
	enddo

return	
end subroutine

subroutine computeEF(ni,nj,nk,dx,dy,dz,efx,efy,efz,phi)
implicit none

integer :: ni,nj,nk,i,j,k
double precision :: dx,dy,dz,efx(ni,nj,nk),efy(ni,nj,nk),efz(ni,nj,nk),phi(ni,nj,nk)

do i=1,ni
	do j=1,nj
		do k=1,nk
	
	
		if(i==1) then
			efx(i,j,k)=-( -3.0*phi(i,j,k)+4.0*phi(i+1,j,k) -phi(i+2,j,k))/(2.0*dx)
		else if(i==ni) then
			efx(i,j,k)=-(phi(i-2,j,k)-4.0*phi(i-1,j,k)+3.0*phi(i,j,k))/(2.0*dx)
		else
			efx(i,j,k)=-(phi(i+1,j,k)-phi(i-1,j,k))/(2.0*dx)
		end if
		
		if(j==1) then
			efy(i,j,k)=-( -3.0*phi(i,j,k)+4.0*phi(i,j+1,k) -phi(i,j+2,k))/(2.0*dy)
		else if(j==nj) then
			efy(i,j,k)=-(phi(i,j-2,k)-4.0*phi(i,j-1,k)+3.0*phi(i,j,k))/(2.0*dy)
		else
			efy(i,j,k)=-(phi(i,j+1,k)-phi(i,j-1,k))/(2.0*dy)
		end if
		
		if(k==1) then
			efz(i,j,k)=-( -3.0*phi(i,j,k)+4.0*phi(i,j,k+1) -phi(i,j,k+2))/(2.0*dz)
		else if(k==nk) then
			efz(i,j,k)=-(phi(i,j,k-2)-4.0*phi(i,j,k-1)+3.0*phi(i,j,k))/(2.0*dz)
		else
			efz(i,j,k)=-(phi(i,j,k+1)-phi(i,j,k-1))/(2.0*dz)
		end if

	enddo
	enddo
	enddo

return
end subroutine

subroutine scatter(f,val,&
			x,x0,dx,y,y0,dy,z,z0,dz,&
			ni,nj,nk)			!Particle to Mesh
implicit none

integer :: ni,nj,nk
double precision :: f(ni,nj,nk)
double precision :: lx,ly,lz,di,dj,dk,x,x0,dx,y,y0,dy,val,z,z0,dz
integer :: i,j,k

lx=(x-x0)/dx
ly=(y-y0)/dy
lz=(z-z0)/dz

i=int(lx)
j=int(ly)
k=int(lz)

di=lx-i
dj=ly-j
dk=lz-k

i=i+1
j=j+1
k=k+1

	f(i,j,k)=f(i,j,k)+val*(1-di)*(1-dj)*(1-dk)
	f(i+1,j,k)=f(i+1,j,k)+val*(di)*(1-dj)*(1-dk)
	f(i+1,j+1,k)=f(i+1,j+1,k)+val*(di)*(dj)*(1-dk)
	f(i,j+1,k)=f(i,j+1,k)+val*(1-di)*(dj)*(1-dk)
	f(i,j,k+1)=f(i,j,k+1)+val*(1-di)*(1-dj)*(dk)
	f(i+1,j,k+1)=f(i+1,j,k+1)+val*(di)*(1-dj)*(dk)
	f(i+1,j+1,k+1)=f(i+1,j+1,k+1)+val*(di)*(dj)*(dk)
	f(i,j+1,k+1)=f(i,j+1,k+1)+val*(1-di)*(dj)*(dk)

return
end subroutine

subroutine gather(Tval,f,&
		x,x0,dx,y,y0,dy,z,z0,dz,&
			ni,nj,nk)
implicit none

integer :: ni,nj,nk,i,j,k
double precision :: Tval,f(ni,nj,nk),di,dj,dk,lx,ly,lz,x,x0,y,y0,z,z0,dx,dy,dz

lx=(x-x0)/dx
ly=(y-y0)/dy
lz=(z-z0)/dz

i=int(lx)
j=int(ly)
k=int(lz)

di=lx-i
dj=ly-j
dk=lz-k

i=i+1
j=j+1
k=k+1

Tval = f(i,j,k)*(1-di)*(1-dj)*(1-dk) &
		+f(i+1,j,k)*(di)*(1-dj)*(1-dk) &
		+f(i+1,j+1,k)*(di)*(dj)*(1-dk) &
		+f(i,j+1,k)*(1-di)*(dj)*(1-dk) &
		+f(i,j,k+1)*(1-di)*(1-dj)*(dk) &
		+f(i+1,j,k+1)*(di)*(1-dj)*(dk) &
		+f(i+1,j+1,k+1)*(di)*(dj)*(dk) &
		+f(i,j+1,k+1)*(1-di)*(dj)*(dk)

return		
end subroutine

subroutine animation(t,nplot,&
			x0,y0,z0,dx,dy,dz,phi,den_ions,den_eles,rho,ni,nj,nk)
implicit none

integer :: i,j,k,ncount,t,nplot,ni,nj,nk
double precision :: x0,y0,z0,dx,dy,dz,phi(ni,nj,nk),den_ions(ni,nj,nk),den_eles(ni,nj,nk),rho(ni,nj,nk)
character ::  NCTEXT*4,FNAME*16,EXT*4,OUTFILE*256

EXT='.dat'
FNAME='TIME_'
NCOUNT=t/Nplot
WRITE(NCTEXT,'(I4.4)') NCOUNT

OUTFILE='anim/'//trim(FNAME)//NCTEXT//EXT

	100 format(3(E14.6))
	
	open(1,file=OUTFILE)
		
		write(1,*) 'TITLE="PLASMA2D"'
		write(1,*) 'VARIABLES="X","Y","Z","PHI","dn.e-","dn.O+","RHO"'
		write(1,*) 'ZONE T="1", I=',ni,", J=",nj,", K=",nk
		
		do k=1,nk
		do j=1,nj
		do i=1,ni
			write(1,100) x0 + (i-1)*dx, y0 + (j-1)*dy, z0 + (k-1)*dz,&
				phi(i,j,k),den_eles(i,j,k),den_ions(i,j,k),rho(i,j,k)
		enddo
		enddo
		enddo
		
	close(1)
return
end subroutine

subroutine diagnostic(Eke,Eki,PE,EPS0,AMU,ME,&
		nd_eles,nd_ions,mpwi,mpwe,ue,ve,we,ui,vi,wi,&
		ni,nj,nk,node_vol,efx,efy,efz)
implicit none

integer :: i,j,k,nd_eles,nd_ions,ni,nj,nk
double precision :: mpwi(nd_ions),mpwe(nd_eles),Eke,Eki,PE,EPS0,AMU,ME,&
		ue(nd_eles),ve(nd_eles),we(nd_eles),&
		ui(nd_ions),vi(nd_ions),wi(nd_ions),&
		efx(ni,nj,nk),efy(ni,nj,nk),efz(ni,nj,nk),&
		node_vol(ni,nj,nk)

do i=1,nd_eles
EKe=Eke+0.5*(ME)*mpwe(i)*( ue(i)*ue(i) + ve(i)*ve(i) + we(i)*we(i) )
enddo

do i=1,nd_ions
EKi=Eki+0.5*(16.0*AMU)*mpwi(i)*( ui(i)*ui(i) + vi(i)*vi(i) + wi(i)*wi(i) )
enddo

do k=1,nk
do j=1,nj
do i=1,ni
	PE=PE+0.5*EPS0*node_vol(i,j,k)*( efx(i,j,k)*efx(i,j,k) + efy(i,j,k)*efy(i,j,k) + efz(i,j,k)*efz(i,j,k) )
enddo
enddo
enddo

return
end subroutine
