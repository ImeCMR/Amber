include ../config.h

.c.o:
	@echo "[XTALUTIL]  CC $<"
	$(VB)$(CC) -c $(COPTFLAGS) $(CFLAGS) $(WARNFLAGS) -o $@ $< 

OBJS = crdmanip.o grid.o matrix.o pdbRead.o rand2.o topRead.o vector.o

install: AddToBox$(SFX) ChBox$(SFX) PropPDB$(SFX) UnitCell$(SFX) 
	mv AddToBox$(SFX) ChBox$(SFX) PropPDB$(SFX) UnitCell$(SFX) ${BINDIR}

AddToBox$(SFX): AddToBox.o $(OBJS)
	@echo "[XTALUTIL]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -lm

ChBox$(SFX): ChBox.o $(OBJS)
	@echo "[XTALUTIL]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -lm

CheckMD$(SFX): CheckMD.o $(OBJS)
	@echo "[XTALUTIL]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -lm

PropPDB$(SFX): PropPDB.o $(OBJS)
	@echo "[XTALUTIL]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -lm

UnitCell$(SFX): UnitCell.o $(OBJS)
	@echo "[XTALUTIL]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -lm

uninstall:
	/bin/rm -f $(BINDIR)/AddToBox$(SFX) $(BINDIR)/ChBox$(SFX) \
       $(BINDIR)/PropPDB$(SFX) $(BINDIR)/UnitCell$(SFX)

clean:
	/bin/rm -f *.o

