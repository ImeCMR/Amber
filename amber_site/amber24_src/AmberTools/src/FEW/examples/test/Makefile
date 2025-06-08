include ../../../config.h

test:
	# Use AMBER_PREFIX for AMBERHOME to use bin/ folder
	# AMBER_SOURCE is also used since FEW run test in source folder
	AMBERHOME=$(AMBER_PREFIX) AMBER_SOURCE=$(AMBER_SOURCE) csh ./Run.few.test

clean:
	AMBERHOME=$(AMBER_SOURCE) csh ./Run.few.clean
