FFLAGS= -lm -g -extend_source -C -check_bounds -ieee

OBJS= SOPHIA20.o eventgen.o sampling.o inpoutput.o\
      jetset74dp.o rndm.o

SLIB= SOPHIA

$(SLIB): $(OBJS)
	f77 $(FFLAGS) $(OBJS) -o $(SLIB)

SOPHIAr:$(SLIB)
	$(SLIB)
