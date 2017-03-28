## Reverse / Complement in Common Workflow Language

Here are two tools and two workflows, all written to
demonstrate features of the [Common Workflow Language](http://www.commonwl.org/).

In `reverse/` and `complement/` there are simply Python scripts that implement
the reverse and [complement](https://en.wikipedia.org/wiki/Complementary_DNA)
operations on a file containing a DNA string. There are also corresponding
Docker containers and Dockerfiles for building them.

### reverse: a Common Workflow Language tool description

Each tool has a Common Workflow Language tool description file. The
file from `reverse.cwl` is reproduced below:

    cwlVersion: v1.0
    class: CommandLineTool

    hints:
      DockerRequirement:
        dockerPull: pvanheus/reverse:latest

    baseCommand: reverse.py

    inputs:
      dnafile:
        type: File
        inputBinding:
          position: 1

    stdout: $(inputs.dnafile.nameroot)\_reversed$(inputs.dnafile.nameext)

    outputs:
      rev_dnafile:
        type: stdout

For an introduction to the CWL, see the [User Guide](http://www.commonwl.org/v1.0/UserGuide.html).
A few items bear mentioning, however: the Docker image for the tool is associated
with the tool using a hint. I.e. this tool can either be run with or without
(Docker) container support. The `baseCommand` refers to the command line
tool (`reverse.py`) being run, and a single command line argument is specified. Given this command line, `reverse.py` writes output to `stdout`. This is captured and redirected to a file.

Instead of giving a simple string as the output
filename, parameters of the input are used as part of a filename pattern. In
CWL, the `$()` refers to a variable lookup. The `dnafile` is a `File` type,
you can read about its attributes in the corresponding section of the
[CWL specification](http://www.commonwl.org/v1.0/CommandLineTool.html#File).
The effect of this pattern is that if the input is `dna.txt` the output
will be `dna_reversed.txt`.

Finally, a single output is captured from the tool, given the name `rev_dnafile` which is associated with the previously-captured `stdout`.

This tool description allows the tool to be run, either by command line
invocation of something like the `cwltool` script, or via other mechanisms
provided by CWL-supporting workflow management systems. Effectively this
tool description translates between a description of inputs and outputs
and the command line parameters and output files associated with a tool.

### revcomp: a Common Workflow Language Workflow

The reverse and complement tools can be chained together using a workflow. This
is what such a workflow could look like (`revcomp.cwl`):

    cwlVersion: v1.0
    class: Workflow

    requirements:
      - class: InlineJavascriptRequirement

    inputs:
      infile:
        type: File
        inputBinding:
          position: 1

    outputs:
      revcomp_dnafile:
        type: File
        outputSource: complement/comp_dnafile

    steps:
      reverse:
        run: reverse.cwl
        in:
          dnafile: infile
        out: [rev_dnafile]
      complement:
        run: complement.cwl
        in:
          dnafile: reverse/rev_dnafile
        out: [comp_dnafile]

Note that this is `class: Workflow` instead of `class: CommandLineTool`. It has its own `inputs` and `outputs`, just like a command line tool. Instead of linking to a command (with `baseCommand`), a workflow has `steps`, which in this case each refer to a CWL command line tool description file. In each step the inputs provided and outputs used are mentioned using `in` and `out` respectively. The `infile` input from `inputs` is bound to the `dnafile` parameter of `reverse.cwl`. The `rev_dnafile` output of this tool is captured and then bound to the `dnafile` parameter of `complement.cwl`. Finally `comp_dnafile` from the complement step is bound (using `outputSource`) to the `revcomp_dnafile` output of the workflow as a whole.

At this (workflow) level the implementation details of the steps are hidden.
The workflow has no need to know about the structure of the command line or
Docker containers. It only needs to know steps and how steps relate to each other.

One effect of how the steps are written is that the final output file will
inherit patterns from each step. So `dna.txt` becomes `dna_reversed_complement.txt`.

The final example `revcomp_with_rename.cwl` illustrates some of the more
advanced features of CWL. Firstly this workflow takes an optional string
parameter (`outfile_name`) which is the name to give the output file of the workflow.
This parameter is optional because its type is `string?`. Like in
regular expressions, the `?` denotes that this paramter does not have to be
provided. If the parameter is not provided, the output filename would be the same
as for `revcomp.cwl`.

The second difference from `revcomp.cwl` is the step
that allows output renaming. Instead of being provided as a command line tool,
it runs a Javascript expression. This require a `InlineJavascriptRequirement`
to be specified and then a step of type `ExpressionTool`:

    rename:
      run:
        class: ExpressionTool
        inputs:
          infile:
            type: File
          outfile_name:
            type: string?
        outputs:
          outfile: File
        expression: >
          ${
          var outfile = inputs.infile;
          if (inputs.outfile_name) {
            outfile.basename = inputs.outfile_name;
          }
          return { "outfile": outfile }; }
      in:
        infile: complement/comp_dnafile
        outfile_name: outfile_name
      out: [outfile]

This also has an optional `outfile_name` parameter. The Javascript tests if this is
set and if so, uses it to modify the `outfile.basename`, resulting in the
output file having the name specified by the user. The Javascript expression
returns a dictionary whose keys correspond to the names of outputs specified
for the step.

When using `cwltool` the `InlineJavascriptRequirement` requires the use of NodeJS (the `node`) command, either through having it installed as a packer or via the use of a Docker container that can run `node`. The Javascript code is sandboxed and cannot access anything outside of its limited execution context.

### Conclusion

These examples show just a little of what can be done with the Common
Workflow Language. While they don't do much useful, they should illustrate
how the Common Workflow Language allows you to map command line realities to
a higher level tool and workflow specification. These specifications
can then be shared and improved on, paving the way towards collaborative
and reproductible science!
