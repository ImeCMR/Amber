# Makefile for antechamber tests

fast:
	-cd acdoctor && ./Run.cf
	-cd bondtype && ./Run.bondtype
#	-cd charmm && ./Run.charmm
	-cd fpph && ./Run.fpph
	-cd fpph && ./Run.fpph_resp
	-cd g09 && ./Run.g09
	-cd parmchk2 && ./Run.parmchk2
	-cd residuegen && ./Run.ala
	-cd sustiva && ./Run.sustiva
	-cd sustiva && ./Run.sustiva_g2
	-cd tp && ./Run.tp
	-cd multirespfix && ./Run.multirespfix

test: fast
	-cd c60 && ./Run.c60
#	-cd ccif && ./Run.ccif
	-cd charge_equlibration && ./Run.equ_charge
	-cd chemokine && ./Run.1b2t
	-cd tp && ./Run.tp2
#	-cd translate && ./Run.2rh1
	-cd ash && ./Run.ash
	-cd fluorescein && ./Run.fluorescein
	-cd guanine_amber && ./Run.guanine_amber
#	-cd acdoctor && ./Run.acdoctor

test.antechamber: test

test.notreadyforprimetime:
	-cd database && ./Run.database
	-cd gaff2 && ./Run.trimethylammonium
	-cd guanine && ./Run.guanine
	-cd match_atomname && ./Run.match_atomname
	-cd match && ./Run.match
	-cd top2mol2 && ./Run.top2mol2
	-cd tp && ./Run.tp_resp

test.snail:
	-cd nad && ./Run.nad
	-cd nad && ./Run.nad_bcc
	-cd nad && ./Run.nad_resp

test.all: test test.notreadyforprimetime

clean:
	@echo "Cleaning all antechamber test subdirectories:"
	find . -name "*.dif" -print0 | xargs -0 /bin/rm -f
	/bin/rm -f \
		ash/antechamber.out \
		ash/sqm.in \
		ash/sqm.out \
		ash/sqm.pdb \
		bondtype/antechamber.out \
		c60/antechamber.out \
		charmm/antechamber.out \
		fluorescein/antechamber.out \
		fluorescein/sqm.in \
		fluorescein/sqm.out \
		fpph/antechamber.out \
		guanine_amber/antechamber.out \
		guanine_amber/sqm.in \
		guanine_amber/sqm.out \
		guanine_amber/sqm.pdb \
		residuegen/ala_conf1.esp \
		residuegen/ala_conf2.esp \
		residuegen/ANTECHAMBER_AC.AC \
		residuegen/ANTECHAMBER_AC.AC0 \
		residuegen/ANTECHAMBER_BOND_TYPE.AC \
		residuegen/ANTECHAMBER_BOND_TYPE.AC0 \
		residuegen/antechamber.out \
		residuegen/ATOMTYPE.INF \
		residuegen/esout \
		residuegen/espgen.out \
		residuegen/NEWPDB.PDB \
		residuegen/PREP.INF \
		residuegen/punch \
		residuegen/QIN \
		residuegen/qout \
		residuegen/RESIDUE_GEN_MAINCHAIN.DAT \
		residuegen/RESIDUE_GEN_RESP.INPUT1 \
		residuegen/RESIDUE_GEN_RESP.INPUT2 \
		residuegen/RESIDUE_GEN_RESP.OUTPUT1 \
		residuegen/RESIDUE_GEN_RESP.OUTPUT2 \
		residuegen/RESIDUE_GEN_RESPGEN.DAT \
		residuegen/RESIDUE_GEN.AC \
		residuegen/residuegen.out \
		sustiva/antechamber.out \
		sustiva/leap.log \
		sustiva/prmtop \
		sustiva/sqm.in \
		sustiva/sqm.out \
		tp/antechamber.out \
		tp/sqm.out

check:
	(find . -name "*.dif" ! -size 0c -print | \
	while read dif ;\
	do \
	    echo $$dif ;\
	    cat $$dif ;\
	    echo ;\
	done ;\
	)

