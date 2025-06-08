dacdif	= $(AMBERHOME)/test/dacdif
exe	= elsize
input	= trxox.2.4ns.pqr
output	= elsize_out.trxox
save	= $(output).save

test:
	$(exe) $(input) > $(output)
	$(dacdif) $(output) $(save)
	$(MAKE) clean

.PHONY:	clean

clean:
	-rm -rf $(output)