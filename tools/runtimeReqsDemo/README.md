``` bash
$ cwl-runner runtime-demo.cwl 
cwl-runner 1.0.20160811184335
Final process status is success
{
    "outdirSize": 1024, 
    "ram": 1024, 
    "tmpdirSize": 1024, 
    "cores": 1, 
    "tmpdir": "/tmp/tmpKAeAd6", 
    "outdir": "/tmp/tmp89v6zt"
}
$ cwl-runner runtime-demo.cwl runtime-demo.input.yml 
cwl-runner 1.0.20160811184335
Final process status is success
{
    "outdirSize": 88888, 
    "ram": 4242, 
    "tmpdirSize": 100000000, 
    "cores": 4, 
    "tmpdir": "/tmp/tmpoafYIJ", 
    "outdir": "/tmp/tmpgjfU29"
}
```
