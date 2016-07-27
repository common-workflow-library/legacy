#!/usr/bin/env cwl-runner
cwlVersion: "cwl:draft-3"

class: SchemaDefRequirement
types:
  - name: "#depts"
    type: "record"
    fields:
      - name: "#dept"
        type:
          type: enum
          name: "JustDepts"
          symbols: ["-bg","-bga","-d"]
        description: |
          -bg Report depth in BedGraph format. For details, see
          genome.ucsc.edu/goldenPath/help/bedgraph.html
          -bga Report depth in BedGraph format, as above (-bg).
          However with this option, regions with zero
          coverage are also reported. This allows one to
          quickly extract all regions of a genome with 0
          coverage by applying: "grep -w 0$" to the output.
          -d Report the depth at each genome position (with one-based coordinates).
          Default behavior is to report a histogram.
        inputBinding:
          position: 4
      - name: "#scale"
        type: ["null",float ]
        description: |
          Scale the coverage by a constant factor.
          Each coverage value is multiplied by this factor before being reported.
          Useful for normalizing coverage by, e.g., reads per million (RPM).
          - Default is 1.0; i.e., unscaled.
          - (FLOAT)
        inputBinding:
          position: 4
          prefix: -scale
