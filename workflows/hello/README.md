A simple CWL workflow printing "Hello World" into a File.

See : https://github.com/common-workflow-language/workflows/pull/2

## echo 'Hello'

Execute:

```bash
$ cwl-runner  hello.cwl
cwl-runner 1.0.20161007181528
[job step0] /tmp/tmp5t6hw9$ echo \
    'Hello World' > /tmp/tmp5t6hw9/response.txt
[step step0] completion status is success
[workflow hello.cwl] outdir is /tmp/tmp7BOFU6
Final process status is success
{
    "response": {
        "checksum": "sha1$648a6a6ffffdaa0badb23b8baf90b6168dd16b3a", 
        "basename": "response.txt",
        "location": "file:///path/to/workflows/workflows/hello/response.txt",
        "path": "/path/to/workflows/workflows/hello/response.txt",
        "class": "File", 
        "size": 12
    }
}


$ cat response.txt
Hello World
```

## echo 'Hello' with parameters

```bash
$ cwl-runner  hello-param.cwl params.yaml
cwl-runner 1.0.20161007181528
[job step0] /tmp/tmp3iO0b_$ echo \
    -n \
    -e \
    'Hello, CWL !
Hello World !' > /tmp/tmp3iO0b_/response.txt
[step step0] completion status is success
[workflow hello-param.cwl] outdir is /tmp/tmpX33xib
Final process status is success
{
    "response": {
        "checksum": "sha1$e8bb28df025c10299db8e73281fbf96d402a1bc0", 
        "basename": "response.txt",
        "location": "file:///path/to/workflows/workflows/hello/response.txt",
        "path": "/path/to/workflows/workflows/hello/response.txt",
        "class": "File", 
        "size": 26
    }
}

$ cat response.txt 
Hello, CWL !
Hello World !
```

with parameters specified on the command line:

```bash
$ cwl-runner  hello-param.cwl --usermessage "Yellow submarine\n"
cwl-runner 1.0.20161007181528
[job step0] /tmp/tmpEC7Ibv$ echo \
    -n \
    -e \
    'Yellow submarine\n' > /tmp/tmpEC7Ibv/response.txt
[step step0] completion status is success
[workflow hello-param.cwl] outdir is /tmp/tmpzcZggz
Final process status is success
{
    "response": {
        "checksum": "sha1$63cce799f2d6cf47c6c364dff126c2e7dcbabd9d", 
        "basename": "response.txt",
        "location": "file:///path/to/workflows/workflows/hello/response.txt",
        "path": "/path/to/workflows/workflows/hello/response.txt",
        "class": "File", 
        "size": 17
    }
}

$ cat response.txt 
Yellow submarine
```
