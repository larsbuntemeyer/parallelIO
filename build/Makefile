


include Makefile.h


MODULES = mo_mpi.o grid.o io.o dbase.o mo_parorg.o mo_decomp.o mo_org.o

main: main.o ${MODULES} ${OBJECTS}
	$(F90) main.o -o main ${MODULES} $(LFLAGS)

matrix: matrix3D.o
	$(F90) matrix3D.o -o matrix

main.o: ${OBJECTS} ${MODULES} 

dbase.o: grid.o mo_decomp.o
io.o : dbase.o
grid.o: mo_mpi.o
mo_decomp.o: mo_parorg.o mo_org.o grid.o

%.o: %.f90
	$(F90) -c $(CFLAGS) $<

%.o: %.f
	$(F90) -c $(CFLAGS) $<

clean:
	rm *.o *.mod main
