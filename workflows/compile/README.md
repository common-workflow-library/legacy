# Compile with GNU GCC


## compile1.cwl

compile1.cwl is the CWL equivalent of the following Makefile.

```make
.PHONY:all
all:a.out
a.out: source1.o source2.o
	gcc -o a.out source1.o source2.o
source1.o : source1.c source1.h
	gcc -Wall -c -o source1.o source1.c
source2.o : source2.c
	gcc -Wall -c -o source2.o source2.c
```

it compiles two sources and link them.
