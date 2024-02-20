####  for windows 
	gg=gfortran
	makehtml:=make html
	makelatex:=make latex
mqc: main.o  Global.o Int1e.o Int2e.o IO.o HartreeFock.o NumInt.o 
	$(gg) -o mqc main.o  Global.o Int1e.o Int2e.o IO.o HartreeFock.o\
		NumInt.o Init.o  -L .\lib\  -llapack -lblas  -lcint -fdefault-integer-8
main.o: main.f90 IO.mod Int1e.mod Int2e.mod hf.mod Global.mod  Init.mod
		$(gg) -c main.f90   
HartreeFock.o hf.mod: Global.o HartreeFock.f90
		$(gg) -c   HartreeFock.f90   
Global.o  Global.mod: Global.f90
		$(gg) -c Global.f90  
Int1e.o Int1e.mod: Int1e.f90  gausshermite.mod
		$(gg) -c Int1e.f90  
Int2e.o Int2e.mod: Int2e.f90 
		$(gg) -c   Int2e.f90  
gausshermite.mod NumInt.o: NumInt.f90 
		$(gg) -c  NumInt.f90 		
IO.o IO.mod: IO.f90 Global.o Init.mod
		$(gg) -c IO.f90 
Init.o Init.mod: Init.f90 
		$(gg) -c Init.f90 		
clean:
	del  *.o  *.mod
#### for Linux		
mqclinux: main.o  Global.o Int1e.o Int2e.o IO.o HartreeFock.o NumInt.o 
	$(gg) -o mqc main.o  Global.o Int1e.o Int2e.o IO.o HartreeFock.o NumInt.o  Init.o -L ./lib/  -llapack -lblas -fdefault-integer-8
cleanlinux:
	 rm  *.o  *.mod
html:
	cd Doc  && $(makehtml)
latex:
	cd Doc  && $(makelatex)
cleandoc:
	cd Doc  && make clean