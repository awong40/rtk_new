# makefile for rnx2rtkp

BINDIR  = /usr/local/bin
SRC     = /home/ec2-user/rtk_new

OPTS    = -DTRACE -DENAGLO -DENAQZS -DENAGAL -DNFREQ=3

# for no lapack
CFLAGS  = -Wall -O3 -ansi -pedantic -Wno-unused-but-set-variable -I$(SRC) $(OPTS) -g
LDLIBS  = -lm -lrt

#CFLAGS  = -Wall -O3 -ansi -pedantic -Wno-unused-but-set-variable -I$(SRC) -DLAPACK $(OPTS)
#LDLIBS  = -lm -lrt -llapack -lblas

# for gprof
#CFLAGS  = -Wall -O3 -ansi -pedantic -Wno-unused-but-set-variable -I$(SRC) -DLAPACK $(OPTS) -pg
#LDLIBS  = -lm -lrt -llapack -lblas -pg

# for mkl
#MKLDIR  = /opt/intel/mkl
#CFLAGS  = -O3 -ansi -pedantic -Wno-unused-but-set-variable -I$(SRC) $(OPTS) -DMKL
#LDLIBS  = -L$(MKLDIR)/lib/intel64 -lm -lrt -lmkl_core -lmkl_intel_lp64 -lmkl_gnu_thread -liomp5 -lpthread

all        : rnx2rtkp
rnx2rtkp   : rnx2rtkp.o rtkcmn.o pntpos.o


rnx2rtkp.o : ../rnx2rtkp.c
	$(CC) -c $(CFLAGS) ../rnx2rtkp.c
rtkcmn.o   : $(SRC)/rtkcmn.c
	$(CC) -c $(CFLAGS) $(SRC)/rtkcmn.c
pntpos.o   : $(SRC)/pntpos.c
	$(CC) -c $(CFLAGS) $(SRC)/pntpos.c


rnx2rtkp.o : $(SRC)/rtklib.h
rtkcmn.o   : $(SRC)/rtklib.h
pntpos.o   : $(SRC)/rtklib.h

CMD1    = ./rnx2rtkp
INPUT11 = ../../../test/data/rinex/07590920.05o ../../../test/data/rinex/30400920.05n
INPUT12 = ../../../test/data/rinex/30400920.05o
OPTS1   = -r -3978241.958 3382840.234 3649900.853

test : test1 test2 test3 test4 test5 test6 test7 test8 test9 test10
test : test11 test12 test13 test14 test15 test16 test17 test18 test19 test20
test : test21 test22 test23 test24

test1 :
	$(CMD1) $(INPUT11) -x 5 -o test1.pos
test2 :
	$(CMD1) -t -e $(OPTS1) $(INPUT11) > test2.pos
test3 :
	$(CMD1) -t -p 1 -e $(OPTS1) $(INPUT11) $(INPUT12) > test3.pos
test4 :
	$(CMD1) -t -p 3 -e $(OPTS1) $(INPUT11) $(INPUT12) > test4.pos
test5 :
	$(CMD1) -t -m 15 -e $(OPTS1) $(INPUT11) $(INPUT12) > test5.pos
test6 :
	$(CMD1) -t -f 1 -e $(OPTS1) $(INPUT11) $(INPUT12) > test6.pos
test7 :
	$(CMD1) -t -v 5 -e $(OPTS1) $(INPUT11) $(INPUT12) > test7.pos
test8 :
	$(CMD1) -t -i -e $(OPTS1) $(INPUT11) $(INPUT12) > test8.pos
test9 :
	$(CMD1) -t -p 0 $(OPTS1) $(INPUT11) > test9.pos
test10 :
	$(CMD1) -t -p 0 $(OPTS1) $(INPUT11) -o test10.pos
test11 :
	$(CMD1) -t -p 0 -n $(OPTS1) $(INPUT11) > test11.pos
test12 :
	$(CMD1) -t -p 0 -g $(OPTS1) $(INPUT11) > test12.pos
test13 :
	$(CMD1) -t -p 0 $(OPTS1) $(INPUT11) > test13.pos
test14 :
	$(CMD1) -t -p 0 -u $(OPTS1) $(INPUT11) > test14.pos
test15 :
	$(CMD1) -t -p 0 -d 9 $(OPTS1) $(INPUT11) > test15.pos
test16 :
	$(CMD1) -t -p 0 -s , $(OPTS1) $(INPUT11) > test16.pos
test17 :
	$(CMD1) -t -b -e $(OPTS1) $(INPUT11) $(INPUT12) > test17.pos
test18 :
	$(CMD1) -t -c -e $(OPTS1) $(INPUT11) $(INPUT12) > test18.pos
test19 :
	$(CMD1) -t -h -e $(OPTS1) $(INPUT11) $(INPUT12) > test19.pos
test20 :
	$(CMD1) -t -p 4 -a $(OPTS1) $(INPUT11) $(INPUT12) > test20.pos
test21 :
	$(CMD1) $(INPUT11) $(INPUT12) > test21.pos
test22 :
	$(CMD1) -k opts1.conf $(INPUT11) $(INPUT12) > test22.pos
test23 :
	$(CMD1) -k opts2.conf $(INPUT11) $(INPUT12) > test23.pos
test24 :
	$(CMD1) -k opts3.conf $(INPUT11) $(INPUT12) -y 2 -o test24.pos
test25 :
	$(CMD1) -k opts4.conf $(INPUT11) $(INPUT12) -y 2 -o test25.pos

clean :
	rm -f rnx2rtkp rnx2rtkp.exe *.o *.pos *.trace

prof :
	gprof rnx2rtkp > prof.txt

install :
	cp rnx2rtkp $(BINDIR)

