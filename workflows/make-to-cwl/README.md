"I'd like to know if I can easily convert a Makefile based workflow to CWL. For
example, for the simple Makefile, want would be equivalent CWL file please?"

To run for yourself after `pip install cwl-runner`

`./rna.json`
or
`cwl-runner rna.json`

See: Makefile, dna.cwl

Notes:

* Each tool (echo, tr, cat) is wrapped as an individual reusable component.
* The workflow is parameterized (unlike the supplied Makefile).
* The "get_sequences" and "translate_sequences" steps are scatter steps, so they can be processed in parallel.

General advantages of CWL over plain Makefiles:

* You can run components in Docker containers without the details of getting data into/out of the container
* You can push the workflow out to run on a cluster
* Your workflow can run unmodified on multiple cluster platforms
