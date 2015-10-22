

all: getPCA getMult

clean: 
	rm getPCA getPCA.o getMult getMult.o

getPCA: getPCA.o
	cc getPCA.o -llapack -lblas -lm -o getPCA
 
getPCA.o: getPCA.c
	cc -c getPCA.c -o getPCA.o

getMult: getMult.o
	cc getMult.o -o getMult -lm

getMult.o: getMult.c
	cc -c getMult.c -o getMult.o
