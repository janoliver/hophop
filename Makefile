#  Created by Jan Oliver Oelerich

COMP=/usr/bin/gcc
CFLAGS=-Wall -O3 -ffast-math
LFLAGS=-lgsl -lgslcblas -lm
#GGO=gengetopt
GGO=/import/vtthomes/oelerich/bin/gengetopt

all: hop clean

hop: hop.c hop.h cmdline.o init_mc.o hopping_mc.o output.o analyze_mc.o
	$(COMP) -o hop hop.c cmdline.o mc/hopping_mc.o output.o mc/analyze_mc.o mc/init_mc.o $(CFLAGS) $(LFLAGS)

init_mc.o: mc/init_mc.c hop.h
	$(COMP) -c -o mc/init_mc.o mc/init_mc.c $(CFLAGS)

hopping_mc.o: mc/hopping_mc.c hop.h
	$(COMP) -c -o mc/hopping_mc.o mc/hopping_mc.c $(CFLAGS)

analyze_mc.o: mc/analyze_mc.c hop.h
	$(COMP) -c -o mc/analyze_mc.o mc/analyze_mc.c $(CFLAGS)

output.o: output.c hop.h
	$(COMP) -c -o output.o output.c $(CFLAGS)

cmdline.c: options.ggo 
	$(GGO) --conf-parser -i options.ggo 

cmdline.o: cmdline.c
	$(COMP) cmdline.c -c -o cmdline.o $(CFLAGS)

clean: 
	rm -rf output.o mc/*.o
	echo "Done."

.SILENT: clean
