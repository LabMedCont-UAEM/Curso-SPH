module kernel

	type :: cliente
		character*36 :: nombre
		integer :: edad
		real :: subtotal,total
		character*8 :: codigoProducto
		logical :: premium
	end type
	
	type :: producto
		character*64 :: nombre
		character*4 :: codigo
		real :: precio
		logical :: descuento
	end type
	
end module

program clientes
use kernel

type(cliente) :: usr
type(producto), dimension(5) :: pdc 
integer :: i,j
real :: total=0.0,subtotal=0.0
character*4, dimension(3) :: listaCompras

pdc(1)=producto('Lavadora MABE','12HN',12890.0,.false.)
pdc(2)=producto('Horno microondas LG','ax23',7800.0,.true.)
pdc(3)=producto('Celular Samsung Galaxy','SM3R',6500.0,.true.)
pdc(4)=producto('Comedor Deluxe','DLE2',26500.0,.true.)
pdc(5)=producto('Pantalla LG','LG56',19990.0,.false.)

usr=cliente('Alexis',29,0.0,0.0,'ABC-WXYZ',.false.)

listaCompras=(/'12HN','SM3R','LG56'/)

	write(*,'("¡Hola ",A,"! Tu edad es de ",I0," años")') trim(usr%nombre),usr%edad
	write(*,'("Aqui tienes tu lista de compras:")')
	
do i=1,size(listaCompras)
	do j=1,size(pdc)
	if(listaCompras(i)==pdc(j)%codigo) then
		write(*,'("	* ",A," Precio $",F12.2)') trim(pdc(j)%nombre),pdc(j)%precio
		usr%subtotal=usr%subtotal+pdc(j)%precio
		
	end if
	enddo
enddo
	print*, '	-----------------------------'
	write(*,'("	Subtotal: $",F12.2)') usr%subtotal
	
	if(usr%premium) then
	write(*,'(" 	Descuento premium sí - 20% descuento")')
	usr%total=usr%subtotal*0.8
	write(*,'("	Total a pagar: $",F12.2)') usr%total
	else
	write(*,'(" 	Descuento premium No")')
	usr%total=usr%subtotal
	write(*,'("	Total a pagar: $",F12.2)') usr%total
	end if
	
	
end program
