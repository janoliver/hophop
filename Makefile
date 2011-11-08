#  Created by Jan Oliver Oelerich

COMP=/usr/bin/gcc
CFLAGS=-Wall -O3 -ffast-math -march=native 
LFLAGS=-lgsl -lgslcblas -lm
#GGO=gengetopt
GGO=/import/vtthomes/oelerich/software/bin/gengetopt

all: hop clean

hop: hop.c hop.h cmdline.o init_mc.o hopping_mc.o output.o analyze_mc.o
	$(COMP) -o hop hop.c cmdline.o mc_hopping.o output.o mc_analyze.o mc_init.o $(CFLAGS) $(LFLAGS)

init_mc.o: mc_init.c hop.h
	$(COMP) -c -o mc_init.o mc_init.c $(CFLAGS)

hopping_mc.o: mc_hopping.c hop.h
	$(COMP) -c -o mc_hopping.o mc_hopping.c $(CFLAGS)

analyze_mc.o: mc_analyze.c hop.h
	$(COMP) -c -o mc_analyze.o mc_analyze.c $(CFLAGS)

output.o: output.c hop.h
	$(COMP) -c -o output.o output.c $(CFLAGS)

cmdline.c: options.ggo 
	$(GGO) --conf-parser -i options.ggo 

cmdline.o: cmdline.c
	$(COMP) cmdline.c -c -o cmdline.o $(CFLAGS)

clean: 
	rm -rf *.o
	echo "Done."

.SILENT: clean
