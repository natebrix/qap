CPP = g++
BNB = bnbdrv

#CF = -pg -O4
#CF = -g -O4
CF = -O3 

CFLAGS = $(CF) -DBNB_STANDALONE

MESCHACH = ../src/MESCHACH
MESCHACHLIB = ../src/MESCHACH/meschach.a

INCLUDES = -I$(MESCHACH) -I$(EVB3) -I./
 
LDLIBS = -lm $(MESCHACHLIB)

OBJS = qpbnd.o qpbnd2.o lap.o dlap.o qap.o assign.o link1.o \
       solution.o branch.o bnb.o util.o linalg.o fw.o estimate.o svd.o \
       sa_qap.o stats.o console.o lgb_train.o

all:    bnb

clean:
	rm -f *.o *~ a.out

.C.o:	
	$(CPP) $(CFLAGS) $(INCLUDES) -c $?

bnb:	$(OBJS) driver.o
	$(CPP) $(CFLAGS) $(OBJS) driver.o -o $(BNB) $(LDLIBS)

