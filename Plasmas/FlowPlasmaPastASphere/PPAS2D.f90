program plasma2d
implicit none

double precision, parameter :: EPS0=8.8541878e-12	!C/V/m, Permitividad del vacio
double precision, parameter :: QE=1.602176565e-19	!C, carga electrica del electron
double precision, parameter :: AMU=1.660538921e-27	!Kg, unidad de masa atomica
double precision, parameter :: ME=9.10938215e-31	!Kg, masa del electron
double precision, parameter :: K0=1.380648E-23		!J/K, constante de Boltzmann
double precision, parameter :: PI=3.141592653		!pi
double precision, parameter :: EvToK=QE/K0		!1 eV in K - 11604

real :: xrd,yrd
double precision :: x0,y0,xm,ym,x0e,y0e,xme,yme,&
	dx,dy,tol,v,&
	dt
double precision, dimension(:), allocatable :: xi,yi,ui,vi,mpwi,xe,ye,ue,ve,mpwe,&
						efxe,efye,efxi,efyi
double precision, dimension(:,:), allocatable :: phi,rho,efx,efy,den_eles,den_ions,node_vol

integer, dimension(:,:), allocatable ::	flagPot

integer :: i,j,nd_ions,nd_eles,&
	ni,nj,max_it,t

!-------------------------------------------------------------------| Particles info |
CALL RANDOM_SEED()

dt=1.0d-10

x0=0.0d0
y0=0.0d0
xm=0.2d0
ym=0.2d0

x0e=0.0d0
y0e=0.0d0
xme=0.1d0
yme=0.1d0



nd_ions=80000
nd_eles=20000

allocate( efxi(nd_ions),efyi(nd_ions),&
	xi(nd_ions),yi(nd_ions),ui(nd_ions),vi(nd_ions),mpwi(nd_ions),&
	efxe(nd_eles),efye(nd_eles),&
	xe(nd_eles),ye(nd_eles),ue(nd_eles),ve(nd_eles),mpwe(nd_eles) )
	

do i=1,nd_eles
	call RANDOM_NUMBER(xrd)
	xe(i)=x0e+xrd*(xme-x0e)
	
	call RANDOM_NUMBER(yrd)
	ye(i)=y0e+yrd*(yme-y0e)
enddo

do i=1,nd_ions
	call RANDOM_NUMBER(xrd)
	xi(i)=x0+xrd*(xm-x0)
	
	call RANDOM_NUMBER(yrd)
	yi(i)=y0+yrd*(ym-y0)
enddo

ui=0.0
vi=0.0
ue=0.0
ve=0.0
!-----------------------------------------------------------------------------------------

ni=50
nj=50

dx=(xm-x0)/(ni-1)
dy=(ym-y0)/(nj-1)

tol=1.0d-4
max_it=5000

allocate(phi(ni,nj),rho(ni,nj),efx(ni,nj),efy(ni,nj),&
	den_eles(ni,nj),den_ions(ni,nj),node_vol(ni,nj),&
	flagPot(ni,nj))

!phi(1,:)=1.0d0
!phi(:,1)=2.0d0

flagPot=0
call defineSphere(flagPot,phi,ni,nj,dx,dy)
	
mpwe=(xme-x0e)*(yme-y0e)*(1.0d11)/(nd_eles)
mpwi=(xm-x0)*(ym-y0)*(1.0d11)/(nd_ions)


do i=1,ni
do j=1,nj
	V=dx*dy
	if( (i==1).or.(i==ni) ) V=V*0.5
	if( (j==1).or.(j==nj) ) V=V*0.5
	node_vol(i,j)=V
enddo
enddo

do t=1,max_it

	rho=0.0
	den_ions=0.0
	den_eles=0.0
	
	do i=1,nd_ions
		call scatter(den_ions,mpwi(i),&
			xi(i),x0,dx,yi(i),y0,dy,&
			ni,nj)
	enddo

	den_ions=den_ions/node_vol

	do i=1,nd_eles
		call scatter(den_eles,mpwe(i),&
				xe(i),x0e,dx,ye(i),y0e,dy,&
				ni,nj)
	enddo

	den_eles=den_eles/node_vol

	rho=rho+den_ions*QE
	rho=rho+den_eles*(-1.0*QE)
		

	call solve_pot(dx,dy,max_it,tol,ni,nj,phi,rho,&
			EPS0,flagPot)	

	call computeEF(ni,nj,dx,dy,efx,efy,phi,flagPot)

	do i=1,nd_eles

		!-------------------------------------| Advance for eles |
		call gather(efxe(i),efx,&
			xe(i),x0,dx,ye(i),y0,dy,&
				ni,nj)
				
		call gather(efye(i),efy,&
			xe(i),x0,dx,ye(i),y0,dy,&
				ni,nj)	
				
		
		ue(i)=ue(i)+(-1.0*QE/ME)*efxe(i)*dt
		ve(i)=ve(i)+(-1.0*QE/ME)*efye(i)*dt
		
		xe(i)=xe(i)+ue(i)*dt
		ye(i)=ye(i)+ve(i)*dt
		
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
	enddo
	
	do i=1,nd_ions	
		!-------------------------------------| Advance for ions |
		call gather(efxi(i),efx,&
			xi(i),x0,dx,yi(i),y0,dy,&
				ni,nj)
				
		call gather(efyi(i),efy,&
			xi(i),x0,dx,yi(i),y0,dy,&
				ni,nj)	
				
		
		ui(i)=ui(i)+(QE/(16.0*AMU))*efxi(i)*dt
		vi(i)=vi(i)+(QE/(16.0*AMU))*efyi(i)*dt
		
		xi(i)=xi(i)+ui(i)*dt
		yi(i)=yi(i)+vi(i)*dt
		
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
	enddo


	if(mod(t,100).eq.0.0) then
		print*, t
		
		call animation(t,100,&
			flagPot,&
			x0,y0,dx,dy,&
			phi,den_ions,den_eles,rho,&
			ni,nj)
	end if
enddo



deallocate(phi,rho,efx,efy,&
	den_eles,den_ions,node_vol,flagPot)

deallocate( efxi,efyi,&
	xi,yi,ui,vi,mpwi,&
	efxe,efye,&
	xe,ye,ue,ve,mpwe )
	
end program

subroutine solve_pot(dx,dy,max_it,tol,ni,nj,phi,rho,&
		EPS0,flagPot)
implicit none

integer :: i,j,n,max_it,ni,nj,flagPot(ni,nj)
double precision :: dx,dy,idx2,idy2,phi_new,summ,R,EPS0,L2,tol
double precision :: phi(ni,nj),rho(ni,nj)

	idx2=1.0/(dx*dx)
	idy2=1.0/(dy*dy)


	do n=1,max_it
		
		do i=2,ni-1
			do j=2,nj-1
				if(flagPot(i,j)==0) then
					phi_new=( rho(i,j)/EPS0 &
						+ idx2*( phi(i-1,j) + phi(i+1,j) ) &
						+ idy2*( phi(i,j-1) + phi(i,j+1) ) )/&
						( 2.0*idx2 + 2.0*idy2 )
						
					phi(i,j)=phi(i,j)+1.4*(phi_new - phi(i,j))
				end if
			enddo
		enddo
	
		if( mod(n,25) == 0 ) then
		
			summ=0.0d0
			
			
			do i=2,ni-1
				do j=2,nj-1
			
					if(flagPot(i,j)==0) then
						R=-phi(i,j)*( 2.0*idx2 + 2.0*idy2  ) &
						+ rho(i,j)/EPS0 &
						+ idx2*( phi(i-1,j) + phi(i+1,j) ) &
						+ idy2*( phi(i,j-1) + phi(i,j+1) )
						
						summ=summ+R*R
					end if
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

subroutine computeEF(ni,nj,dx,dy,efx,efy,phi,flagPot)
implicit none

integer :: ni,nj,i,j,flagPot(ni,nj)
double precision :: dx,dy,efx(ni,nj),efy(ni,nj),phi(ni,nj)

do i=1,ni
	do j=1,nj

	if(flagPot(i,j)==0) then
	
		if(i==1) then
			efx(i,j)=-( -3.0*phi(i,j)+4.0*phi(i+1,j) -phi(i+2,j))/(2.0*dx)
		else if(i==ni) then
			efx(i,j)=-(phi(i-2,j)-4.0*phi(i-1,j)+3.0*phi(i,j))/(2.0*dx)
		else
			efx(i,j)=-(phi(i+1,j)-phi(i-1,j))/(2.0*dx)
		end if
		
		if(j==1) then
			efy(i,j)=-( -3.0*phi(i,j)+4.0*phi(i,j+1) -phi(i,j+2))/(2.0*dy)
		else if(j==nj) then
			efy(i,j)=-(phi(i,j-2)-4.0*phi(i,j-1)+3.0*phi(i,j))/(2.0*dy)
		else
			efy(i,j)=-(phi(i,j+1)-phi(i,j-1))/(2.0*dy)
		end if

	end if
			
	enddo
	enddo

return
end subroutine

subroutine scatter(f,val,&
			x,x0,dx,y,y0,dy,&
			ni,nj)			!Particle to Mesh
implicit none

integer :: ni,nj
double precision :: f(ni,nj)
double precision :: lx,ly,di,dj,x,x0,dx,y,y0,dy,val
integer :: i,j

lx=(x-x0)/dx
ly=(y-y0)/dy

i=int(lx)
j=int(ly)

di=lx-i
dj=ly-j

i=i+1
j=j+1

f(i,j)=f(i,j)+val*(1-di)*(1-dj)
f(i+1,j)=f(i+1,j)+val*(di)*(1-dj)
f(i+1,j+1)=f(i+1,j+1)+val*(di)*(dj)
f(i,j+1)=f(i,j+1)+val*(1-di)*(dj)

return
end subroutine

subroutine gather(Tval,f,&
		x,x0,dx,y,y0,dy,&
			ni,nj)
implicit none

integer :: ni,nj,i,j
double precision :: Tval,f(ni,nj),di,dj,lx,ly,x,x0,y,y0,dx,dy

lx=(x-x0)/dx
ly=(y-y0)/dy

i=int(lx)
j=int(ly)

di=lx-i
dj=ly-j

i=i+1
j=j+1

Tval = f(i,j)*(1-di)*(1-dj) &
		+f(i+1,j)*(di)*(1-dj) &
		+f(i+1,j+1)*(di)*(dj) &
		+f(i,j+1)*(1-di)*(dj)

return		
end subroutine


subroutine defineSphere(flagPot,phi,ni,nj,dx,dy)
implicit none

integer :: i,j,ni,nj,flagPot(ni,nj)
double precision :: phi(ni,nj),dx,dy,xcenter,ycenter,xcellpos,ycellpos,dist,r

xcenter=0.1d0
ycenter=0.1d0
r=0.02d0

do j=1,nj
do i=1,ni
xcellpos=(i-1)*dx
ycellpos=(j-1)*dy

dist=sqrt( (xcenter-xcellpos)**2 + (ycenter-ycellpos)**2 )


if( dist <= r) then
	flagPot(i,j)=1
	phi(i,j)=-1.0d0
end if
enddo
enddo

end subroutine

subroutine animation(t,nplot,flagpot,&
			x0,y0,dx,dy,phi,den_ions,den_eles,rho,ni,nj)
implicit none

integer :: i,j,ncount,t,nplot,ni,nj,flagPOT(ni,nj)
double precision :: x0,y0,dx,dy,phi(ni,nj),den_ions(ni,nj),den_eles(ni,nj),rho(ni,nj)
character ::  NCTEXT*4,FNAME*16,EXT*4,OUTFILE*256

EXT='.dat'
FNAME='TIME_'
NCOUNT=t/Nplot
WRITE(NCTEXT,'(I4.4)') NCOUNT

OUTFILE='anim/'//trim(FNAME)//NCTEXT//EXT

	100 format(4(E14.6))
	
	open(1,file=OUTFILE)
		
		write(1,*) 'TITLE="PLASMA2D"'
		write(1,*) 'VARIABLES="X","Y","PHI","dn.e-","dn.O+","RHO","FlagPOT"'
		write(1,*) 'ZONE T="1", I=',ni,", J=",nj
		
		do j=1,nj
		do i=1,ni
			write(1,100) x0 + (i-1)*dx, y0 + (j-1)*dy, phi(i,j),den_eles(i,j),den_ions(i,j),rho(i,j),flagPot(i,j)
		enddo
		enddo
	close(1)
return
end subroutine
