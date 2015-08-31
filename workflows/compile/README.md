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


execute:

```bash
/cwl-runner --basedir ${PWD} compile1.cwl

cwl-runner 1.0.20150728161219
Must provide input in the form of a json file or command line parameters.
[job 179982700] exec gcc -c -Wall -o source2.o cwl-workflows/workflows/compile/source2.c
[job 179983628] exec gcc -c -Wall -o source1.o cwl-workflows/workflows/compile/source1.c
[job 179982412] exec gcc -o a.out /tmp/tmpjzWWeE/source1.o /tmp/tmpEhY9J5/source2.o
[workflow 180010220] outdir is /cwl-workflows/workflows/compile
Final process status is success
{
    "output": {
        "path": "cwl-workflows/workflows/compile/a.out", 
        "checksum": "sha1$d2ab8ec752fec7bd6e124b6f11da4c47cdd52a05", 
        "class": "File", 
        "size": 7176
    }

```

