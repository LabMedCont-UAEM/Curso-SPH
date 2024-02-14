all: a.out

kernel.o: kernel.f90
	ifx kernel.f90 -c
	
a.out: main.f90 kernel.o
	ifx main.f90 kernel.o -o a.out
