lobSTR is a tool for profiling Short Tandem Repeats (STRs) from high throughput sequencing data.

This lobSTR workflow is based on http://melissagymrek.com/lobstr-code/usage.html

To run locally, you may install the reference implementation of CWL:

1) Install the reference implementation of cwl-runner with  "pip install cwl-runner".

2)  Download reference data:

$ wget http://files.teamerlich.org/lobstr/v3/ref/lobSTR_v3.0.2_hg19_resource_bundle.tar.gz
$ tar xvzf lobSTR_v3.0.2_hg19_resource_bundle.tar.gz

3) Run the demo:

$ ./lobSTR-workflow.cwl lobSTR-demo.json


Alternately, you may run lobSTR on Curoverse cloud.  This does not require
downloading the reference data.

1) Install Arvados cwl-runner using "pip install arvados-cwl-runner".

2) Sign up for a Curoverse Cloud account at https://cloud.curoverse.com

3) Navigate to https://cloud.curoverse.com/current_token and follow the
instructions to set ARVADOS_API_TOKEN and ARVADOS_API_HOST into your shell
session (these are your Curoverse cloud credentials).

4) You can now run the lobSTR workflow and submit jobs to run on Curoverse
cloud:

$ ./lobSTR-workflow.cwl lobSTR-arvados-demo.json

Note that the output data will also be stored on Curoverse cloud, and may be
downloaded through the Workbench web interface or using "arv-get".
