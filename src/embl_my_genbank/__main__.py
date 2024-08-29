#!/usr/bin/env python3

"""
`__main__.py` provides a redundant interface with the embl_my_genbank functions
that can be accessed when built into a cli tool:

```
usage: emb_my_gbk [-h] --gb_path GB_PATH [--metadata METADATA] --species SPECIES [--out_fmt OUT_FMT] [--view_intermediate] [--multi_output]

options:
  -h, --help            show this help message and exit
  --gb_path GB_PATH, -g GB_PATH
                        Genbank file to be converted.
  --metadata METADATA, -m METADATA
                        Metadata file with, at minimum, a column of allele names and a column of representative animals.
  --species SPECIES, -s SPECIES
                        Scientific name for the species under examination.
  --out_fmt OUT_FMT, -o OUT_FMT
                        Format to convert to. Can convert to EMBL, IPD_EMBL, and FASTA.
  --view_intermediate, -v
                        Whether to write out the intermediate cleaned Genbank file for inspection.
  --multi_output        Whether to split each output record into its own embl file.
```
"""

from .embl_my_genbank import main

if __name__ == "__main__":
    main()
