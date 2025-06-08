include $(AMBERHOME)/config.h

SUFFIXES:
SUFFIXES: .c .o

test:  test_mask test_xmin test_md test_lmod test_nmode

test_xmin:  txmin
	time ./txmin xfin.parm7 xfin.rst7 foo.rst7
	/bin/rm -f foo.rst7

test_rismxmin:  trismxmin
	time ./trismxmin xfin.parm7 xfin.rst7 

test_lmod:  tlmod
	./tlmod trpcage.parm7 trpcage.rst7
	/bin/rm -f conflib.dat

test_md:  tmd
	time ./tmd xfin.parm7 xfin.md1.x
	/bin/rm -f xfin.md2.x

test_md.MPI:  tmd.MPI
	time mpirun -np 4 ./tmd.MPI xfin.parm7 xfin.md1.x

test_nmode:  tnmode
	./tnmode crambin.parm7 crambin.rst7 > nmode.out
	/bin/rm -f vecs

test_mask:  tmask
	./tmask xfin.parm7 xfin.md1.x ":1-3"

test_hcp1:  thcp1
	./thcp1

test_hcp2:  thcp2
	./thcp2

test_nabc:  tnabc
	./tnabc 2.pdb 2.rst7

test_netcdf:  tnetcdf
	./tnetcdf

#---------------------------------------------------------------------------
txmin:  txmin.o
	$(CC) -o txmin txmin.o -L$(AMBERHOME)/lib  $(FLIBS) $(LM)

trismxmin:  trismxmin.o
	$(CC) -o trismxmin trismxmin.o -L$(AMBERHOME)/lib  $(FLIBS) $(LM)

tlmod:  tlmod.o
	$(CC) -o tlmod tlmod.o -L$(AMBERHOME)/lib  $(FLIBS) $(LM)

tmd:  tmd.o
	$(CC) -o tmd tmd.o -L$(AMBERHOME)/lib  $(FLIBS) $(LM)

tmd.MPI:  tmd.MPI.o
	mpicc -o tmd.MPI tmd.MPI.o -L$(AMBERHOME)/lib  $(FLIBS) $(LM)

tnmode:  tnmode.o
	$(CC) -o tnmode tnmode.o -L$(AMBERHOME)/lib  $(FLIBS) $(LM)

tmask:  tmask.o
	$(CC) -o tmask tmask.o -L$(AMBERHOME)/lib  $(FLIBS) $(LM)

thcp1:  thcp1.o
	$(CC) -o thcp1 thcp1.o -L$(AMBERHOME)/lib  $(FLIBS) $(LM)

thcp2:  thcp2.o
	$(CC) -o thcp2 thcp2.o -L$(AMBERHOME)/lib  $(FLIBS) $(LM)

tnabc:  tnabc.o
	$(CC) -o tnabc tnabc.o -L$(AMBERHOME)/lib  -lnabc -lcifparse $(FLIBS) $(LM)

tnetcdf:  tnetcdf.o
	$(CC) -o tnetcdf tnetcdf.o -L$(AMBERHOME)/lib $(FLIBS) $(LM)

#---------------------------------------------------------------------------
.c.o:
	$(CC) -c $(COPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) -I$(INCDIR) -o $@ $<

tmd.MPI.o: tmd.c
	mpicc -DMPI -c $(COPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) -I$(INCDIR) -o $@ $<

clean:
	/bin/rm -f *.o txmin tmd tmd.MPI tnabc tlmod tnmode tmask thcp1 thcp2 \
       trpcageTest.nc
