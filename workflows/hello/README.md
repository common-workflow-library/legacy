A simple CWL workflow printing "Hello World" into a File.

See : https://github.com/common-workflow-language/workflows/pull/2


Execute:

```bash
$ cwl-runner  hello.cwl 
/home/lindenb/.local/bin/cwl-runner 1.0.20150728161219
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

