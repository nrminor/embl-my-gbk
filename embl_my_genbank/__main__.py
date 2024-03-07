#!/usr/bin/env python3

"""
`__main__.py` provides a redundant interface with the embl_my_genbank functions
that can be accessed when built into a cli tool:

```
usage: emb_my_gbk [-h] --gb_path GB_PATH --species SPECIES [--out_fmt OUT_FMT] [--view_intermediate VIEW_INTERMEDIATE]

options:
  -h, --help            show this help message and exit
  --gb_path GB_PATH, -g GB_PATH
                        Genbank file to be converted.
  --species SPECIES, -s SPECIES
                        Scientific name for the species under examination.
  --out_fmt OUT_FMT, -o OUT_FMT
                        Format to convert to. Can convert to EMBL, IPD_EMBL, and FASTA.
  --view_intermediate VIEW_INTERMEDIATE, -v VIEW_INTERMEDIATE
                        Boolean; whether to write out the intermediate cleaned Genbank file for inspection.
```
"""

import os

from Bio import SeqIO

from .embl_my_genbank import clean_genbank, parse_command_line_args, write_output


def main() -> None:
    """
    Coordinate the flow of data through the embl_my_genbank functions.
    """

    # parse command line arguments
    gb_path, species, out_format, view_intermediate = parse_command_line_args()

    assert os.path.isfile(gb_path), "Provided file path does not exist"

    # parse out a name for the output file based on the input file
    basename = os.path.basename(gb_path).replace(".gb", "")
    out_path = f"{basename}.emb"

    # revise and clean all records in the Genbank file
    records = [
        clean_genbank(record, species) for record in SeqIO.parse(gb_path, "genbank")
    ]

    # write output in the desired format
    write_output(records, out_path, out_format, view_intermediate)


if __name__ == "__main__":
    main()
