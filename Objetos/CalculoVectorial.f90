module kernel
	type :: vec3
		real :: x,y,z
		contains
		procedure :: getMag
		procedure :: getNormal
	end type

	public operator(+)
	public operator(-)
	public operator(*)
	
	interface operator(+)
		procedure :: plus_vec
	end interface
	
	interface operator(-)
		procedure :: minus_vec
	end interface
	
	interface operator(*)
		procedure :: ppunto_vec
	end interface
	
	contains
	
	function getMag(this) result(m)
	class(vec3) :: this
	real :: m 
	m=sqrt(this%x*this%x+this%y*this%y+this%z*this%z)
	end function
	
	function getNormal(this) result(m)
	class(vec3) :: this
	type(vec3) :: m
	
	m%x=this%x/this%getMag()
	m%y=this%y/this%getMag()
	m%z=this%z/this%getMag()
	
	end function
	
	function plus_vec(u,v) result(w)
	type(vec3), intent(in) :: u,v
	type(vec3) :: w
	
	w%x=u%x+v%x
	w%y=u%y+v%y
	w%z=u%z+v%z
	
	end function
	
	function minus_vec(u,v) result(w)
	type(vec3), intent(in) :: u,v
	type(vec3) :: w
	
	w%x=u%x-v%x
	w%y=u%y-v%y
	w%z=u%z-v%z
	
	end function
	
	function ppunto_vec(u,v) result(w)
	type(vec3), intent(in) :: u,v
	real :: w
	
	w=u%x*v%x + u%y*v%y + u%z*v%z
	end function
	
end module

program calculo
use kernel

type(vec3) :: f,nf,a,b

a=vec3(1.0,0.0,-2.0)
b=vec3(-1.0,0.0,2.0)

f=vec3(1.0,5.0,-2.0)

nf=f%getNormal()

print*, f%getMag()
print*, nf

print*, nf%getMag()

print*, a+b
print*, a-b
print*, a*b
print*, nf*f
end program
