# program makefile.
# valid targets are:
#       all stellar polytr lstmod clean realclean install
.SUFFIXES:
.SUFFIXES: .out .o .c .e .r .F .f .y .l .s .p .vmsf .for .vmsc .com
#FOR=f95
#FOR=f90
FOR=g77
#FOR=f77
FFLAGS=-O
CFLAGS=
#       Flags for generating code for transport to a VMS system
VMSFLAGS=-Uunix -DVMS -P
.SUFFIXES: .F
.F.o: ; $(FOR) $(FFLAGS) -c $*.F
.F.for: ; /lib/cpp $(VMSFLAGS) $*.F > $*.for
.f.for: ; cp $*.f $*.for
RM=rm -f
MV=mv -f
#       Yow!
PROGS=stellar polytr lstmod

all: $(PROGS)
  
install: all
	if [ ! -d $(HOME)/bin ]; then mkdir $(HOME)/bin; \
	else true; fi
	for i in $(PROGS);\
	do \
	$(MV) $$i $(HOME)/bin;\
	done

stellar: stellar.o opacity.o atmos.o gi.o henyey.o nucrat.o addsub.o invstate.o
	$(FOR) $(FFLAGS) stellar.o opacity.o atmos.o gi.o henyey.o nucrat.o addsub.o invstate.o -o stellar
polytr: polytr.o
	$(FOR) $(FFLAGS) polytr.o -o polytr
lstmod: lstmod.o invstate.o
	$(FOR) $(FFLAGS) lstmod.o invstate.o -o lstmod

clean:
	$(RM) *.o core

realclean: clean
	$(RM) $(PROGS) datefile

