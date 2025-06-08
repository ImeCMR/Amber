include $(AMBERHOME)/config.h

fast: test.pytraj

test.pytraj:
	@(if [ "$(SKIP_PYTHON)" = 'yes' ] ; then \
	    echo "Skipping testing of pytraj because python was skipped." ;\
	elif [ "$(SKIP_PYTHON)" = 'no' ] ; then \
		echo "Testing serial pytraj" ;\
		echo "$(PYTHON) test.py" ;\
		$(PYTHON) test.py || echo "$@:  Program error" ;\
	    $(PYTHON) -m pytest -v $(AMBER_SOURCE)/AmberTools/src/pytraj/tests/test_energy/test_sander_energies.py; \
	else \
	    echo "test.pytraj:  Program error" ;\
	    echo "              Undefined value for config.h variable SKIP_PYTHON." ;\
	    exit 2 ;\
	fi ;\
	$(PYTHON) clean.py \
	)

clean:
	$(PYTHON) clean.py
