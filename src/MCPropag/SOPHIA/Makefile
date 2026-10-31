.SUFFIXES: .o .f

FFLAGS= -ffixed-line-length-200 -lm -ggdb -extend_source -C 

OBJS= SOPHIA20.o eventgen.o sampling.o inpoutput.o\
      jetset74dp.o rndm.o

OBJSCS= SOPHIA20crosssec.o eventgen.o sampling.o inpoutput.o\
      jetset74dp.o rndm.o

SLIB= SOPHIA
SLIBCS= SOPHIA_CS

all: $(SLIB) $(SLIBCS)

clean:
	rm *.o $(SLIB) $(SLIBCS)

$(SLIB): $(OBJS)
	gfortran $(FFLAGS) $(OBJS) -o $(SLIB) ; cp SOPHIA ../bin/

$(SLIBCS): $(OBJSCS)
	gfortran $(FFLAGS) $(OBJSCS) -o $(SLIBCS)

.f.o:
	gfortran -ffixed-line-length-200 -std=legacy -ggdb -c $<
