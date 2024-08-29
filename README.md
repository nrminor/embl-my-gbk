# EMBL-my-Genbank: Convert a complicated Genbank file to a barebones EMBL file

### Usage

This repo contains a minimal-dependency, Ruff-formatted, pure Python module that can be accessed directly through a script, e.g.,

```bash
python3 src/embl_my_genbank/embl_my_genbank.py -g file.gb -s "Homo sapiens"
```

Or, after running `uv venv`, `source .venv/bin/activate`, and `uv sync`, as a command:

```bash
emb_my_gbk -g file.gb -s "Homo sapiens"
```

Or simplest of all, with:

```bash
uv run src/embl_my_genbank/embl_my_genbank.py -g file.gb -s "Homo sapiens"
```

The repo comes with a `uv`-generated `pyproject.toml` to reproduce our dev environment really, the only non-stdlib dependencies are [BioPython](https://biopython.org/), [Polars](https://pola.rs/), and [Loguru](https://github.com/Delgan/loguru).

For complete documentation of the API, see the HTML in `docs/`. Recommended usage looks like this:

```bash
usage: emb_my_gbk [-h] --gb_path GB_PATH --species SPECIES [--out_fmt OUT_FMT] [--view_intermediate VIEW_INTERMEDIATE]

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
```

### Purpose

With some regularity, the [Genomic Services Unit](https://primate.wisc.edu/research-services/genomics-services/) at the Wisconsin National Primate Research Center run an allele discovery pipeline on primate MHC amplicons, the goal of which is to submit previously undocumented MHC alleles to the [Immuno-Polymorphism Database (IPD)](https://www.ebi.ac.uk/ipd/). One of the primary pain points in this process is converting between Genbank format, which we use internally to annotate and review allele candidates alongside exon annotations, and the equally complicated EMBL format, which IPD requires. Additionally, IPD imposes non-standard EMBL format requirements, which we have to adhere to in as many as hundreds of new allele candidates. This python module implements many cleaning and conversion steps to take our internal review Genbank files and convert them into IPD-compatible EMBL files--a surprisingly tricky process!
