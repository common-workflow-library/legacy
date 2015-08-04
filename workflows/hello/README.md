A simple CWL workflow printing "Hello World" into a File.

See : https://github.com/common-workflow-language/workflows/pull/2

## echo 'Hello'

Execute:

```bash
$ cwl-runner  hello.cwl 
cwl-runner 1.0.20150728161219
Must provide input in the form of a json file or command line parameters.
[job 177019500] exec echo 'Hello World' > /tmp/tmpz9BXNP/messageout.txt
[workflow 177019180] outdir is /path/to/workflows/workflows/hello
Final process status is success
{
    "output": {
        "path": "//path/to/workflows/workflows/hello/messageout.txt", 
        "checksum": "sha1$648a6a6ffffdaa0badb23b8baf90b6168dd16b3a", 
        "class": "File", 
        "size": 12
    }
}


$ cat /path/to/workflows/workflows/hello/messageout.txt
Hello World
```

## echo 'Hello' with parameters

```bash
$ cwl-runner  hello-param.cwl params.json
cwl-runner 1.0.20150728161219
[job 153175084] exec echo -n -e 'Hello, CWL !
Hello World !' > /tmp/tmp1er8QF/useroutput.txt
[workflow 153177196] outdir is workflows/hello
Final process status is success
{
    "output": {
        "path": "workflows/hello/useroutput.txt", 
        "checksum": "sha1$e8bb28df025c10299db8e73281fbf96d402a1bc0", 
        "class": "File", 
        "size": 26
    }
}

$ cat useroutput.txt 
Hello, CWL !
Hello World !
```

with parameters specified on the command line:

```bash
$ cwl-runner  hello-param.cwl --usermessage "Yellow submarine" --useroutput song.txt
cwl-runner 1.0.20150728161219
[job 154935692] exec echo -n -e 'Yellow submarine' > /tmp/tmpljSIiy/song.txt
[workflow 154938540] outdir is workflows/hello
Final process status is success
{
    "output": {
        "path": "workflows/hello/song.txt", 
        "checksum": "sha1$72fd47e534c95ce205f8fa1d2f78ab144e6a04ff", 
        "class": "File", 
        "size": 16
    }
}lindenb@hardyweinberg:~/src/cwl-workflows/workflows/hello$ cat song.txt 
Yellow submarine
```
