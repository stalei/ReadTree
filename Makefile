EXECS=ReadTree
CC=gcc
CXX=gcc
OPTIMIZE = -Wall -g
#OBJS  =ReadTree.o core_mymalloc.o
#all:${EXECS}

ReadTree: ReadTree.o core_mymalloc.o
	${CC} ${OPTIMIZE} -o ReadTree ReadTree.o core_mymalloc.o

ReadTree.o: ReadTree.c
	${CC} ${OPTIMIZE} -c ReadTree.c
core_mymalloc.o: core_mymalloc.c
	${CC} ${OPTIMIZE} -c core_mymalloc.c




#all: ${EXECS}
#	${CC} ${OPTIMIZE} ReadTree ReadTree.c
Clean:
	rm -f ${EXECS}
