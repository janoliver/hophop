#  Created by Jan Oliver Oelerich

COMP=/usr/bin/gcc
CFLAGS=-Wall -O3 -ffast-math -march=native 
LFLAGS=-lgsl -lgslcblas -lm
GGO=gengetopt
#GGO=/import/vtthomes/oelerich/software/bin/gengetopt
#CFLAGS=-Wall -O3 -ffast-math -I/import/vtthomes/oelerich/Workspace/cluster/bin/include
#LFLAGS=-L/import/vtthomes/oelerich/Workspace/cluster/bin/lib -lgsl -lgslcblas -lm


all: hop clean

hop: hop.c hop.h cmdline.o mc_init.o mc_hopping.o output.o mc_analyze.o
	$(COMP) -o hop hop.c cmdline.o mc_hopping.o output.o mc_analyze.o mc_init.o $(CFLAGS) $(LFLAGS)

mc_init.o: mc_init.c hop.h
	$(COMP) -c -o mc_init.o mc_init.c $(CFLAGS)

mc_hopping.o: mc_hopping.c hop.h
	$(COMP) -c -o mc_hopping.o mc_hopping.c $(CFLAGS)

mc_analyze.o: mc_analyze.c hop.h
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
