import marimo

__generated_with = "0.2.7"
app = marimo.App()


@app.cell
def __():
    import os
    import marimo as mo
    import Bio
    from enum import Enum, auto
    from itertools import filterfalse
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from typing import List, Optional
    return (
        Bio,
        Enum,
        FeatureLocation,
        List,
        Optional,
        Seq,
        SeqFeature,
        SeqIO,
        SeqRecord,
        auto,
        filterfalse,
        mo,
        os,
    )


@app.cell
def __():
    gb_path = "/Users/nickminor/Documents/dholk_experiments/29841/mafa_input/29766_PacBio125_Mafa_annotated.gb"
    view_intermediate = True
    species = "Macaca fascicularis"
    format = "ipd_embl"
    return format, gb_path, species, view_intermediate


@app.cell
def __(gb_path, os):
    basename = os.path.basename(gb_path).replace(".gb", "")
    out_path = f"{basename}.emb"
    return basename, out_path


@app.cell
def __(Enum, auto):
    class OutType(Enum):
        EMBL = auto()
        IPD_EMBL = auto()
        GENBANK = auto()
        FASTA = auto()
    return OutType,


@app.cell
def __(OutType):
    def choose_out_type(format_choice: str) -> str:
        """ """
        try:
            choice_str = OutType.__members__[format_choice.upper()]
        except KeyError:
            raise ValueError(f"Invalid file type selected: {format_choice}")

        out_type: str
        match choice_str:
            case OutType.EMBL:
                out_type = "embl"
            case OutType.IPD_EMBL:
                out_type = "embl"
            case OutType.GENBANK:
                out_type = "genbank"
            case OutType.FASTA:
                out_type = "fasta"
            case _:
                raise ValueError(
                    f"Invalid output file type selected: {format_choice}"
                )
        return out_type
    return choose_out_type,


@app.cell
def __(SeqRecord):
    def handle_allele_id_placement(_record: SeqRecord) -> SeqRecord:
        """ """

        # parse out the working name of the allele
        allele = _record.name

        # parse out the gene name from the allele name
        gene = allele.split("_")[1]

        # move the allele name into the id position, replacing a geneious annotation
        _record.id = allele

        # rename the accession to be the allele name
        _record.annotations["accessions"] = [allele]

        # check that the allele name was properly propogated
        assert (
            _record.annotations["accessions"][0] == _record.name
        ), "Record name could not be properly assigned to accession annotation"

        return _record
    return handle_allele_id_placement,


@app.cell
def __(SeqRecord):
    def assign_species(_record: SeqRecord, species: str) -> SeqRecord:
        """ """
        # assign the correct species to the organism annotation
        _record.annotations["organism"] = species

        return _record
    return assign_species,


@app.cell
def __(SeqRecord, filterfalse):
    def remove_geneious_annotations(_record: SeqRecord) -> SeqRecord:
        """ """

        # get rid of the misc features Geneious adds in for where parts of annotations
        # were manually deleted
        _record.features = list(
            filterfalse(lambda f: f.type == "misc_feature", _record.features)
        )

        # delete geneious notes
        for _feature in _record.features:
            # Check if the feature has a 'note' qualifier and remove it
            if "note" not in _feature.qualifiers:
                continue
            del _feature.qualifiers["note"]

        return _record
    return remove_geneious_annotations,


@app.cell
def __(SeqRecord, filterfalse):
    def clean_record_features(_record: SeqRecord, species: str) -> SeqRecord:
        """ """

        # filter out genes
        _record.features = list(
            filterfalse(lambda f: f.type == "gene", _record.features)
        )

        for _feature in _record.features:
            # Remove standard_name qualifier from CDS features
            if _feature.type == "CDS" and "standard_name" in _feature.qualifiers:
                del _feature.qualifiers["standard_name"]

            # Update source feature with organism and mol_type
            if _feature.type == "source":
                _feature.qualifiers["organism"] = [species]
                _feature.qualifiers["molecule_type"] = ["genomic DNA"]

        return _record
    return clean_record_features,


@app.cell
def __(
    SeqRecord,
    assign_species,
    clean_record_features,
    handle_allele_id_placement,
    remove_geneious_annotations,
):
    def clean_genbank(_record: SeqRecord, species: str) -> SeqRecord:
        """ """

        # shift around the allele ID into different slots than
        # Geneious places it by default (i.e., Accession and )
        _record = handle_allele_id_placement(_record)

        # assign the user provided scientific name for annotation
        # onto the record
        _record = assign_species(_record, species)

        # clean off geneious annotations
        _record = remove_geneious_annotations(_record)

        # clean up the many features appended to each Genbank record,
        # which BioPython represents as a dictionary
        _record = clean_record_features(_record, species)

        return _record
    return clean_genbank,


@app.cell
def __(SeqIO, records):
    def id_line_replacement(int_embl: str, final_out: str, out_type: str) -> None:
        """ """

        NEW_ID_LINE: str = "ID   XXX; XXX; linear; XXX; XXX; XXX;\n"

        # make an intermediate embl file first
        with open(int_embl, "w", encoding="utf8") as int_handle:
            for _record in records:
                SeqIO.write(records, int_handle, out_type)

        # Open the input file for reading and the output file for writing
        with open(int_embl, "r") as infile, open(final_out, "w") as outfile:
            for line in infile:
                # Check if the line starts with "ID"
                if line.startswith("ID"):
                    # Write the replacement text to the output file
                    outfile.write(NEW_ID_LINE)
                    continue
                # Write the original line to the output file
                outfile.write(line)
    return id_line_replacement,


@app.cell
def __(List, SeqIO, SeqRecord, choose_out_type, id_line_replacement):
    def write_output(
        records: List[SeqRecord],
        final_out: str,
        format: str,
        view_intermediate: bool,
    ) -> None:
        """ """

        out_type = choose_out_type(format)

        if view_intermediate:
            SeqIO.write(records, "intermediate.gbk", "genbank")

        # check if an IPC intermediate embl is needed
        if "IPD" in format.upper():
            int_embl = "intermediate.emb"
            id_line_replacement(int_embl, final_out, out_type)
            return

        with open(final_out, "w", encoding="utf8") as out_file:
            for _record in records:
                SeqIO.write(records, out_file, out_type)
    return write_output,


@app.cell
def __(SeqIO, clean_genbank, gb_path, species):
    records = [
        clean_genbank(_record, species)
        for _record in SeqIO.parse(gb_path, "genbank")
    ]
    return records,


@app.cell
def __(format, out_path, records, view_intermediate, write_output):
    write_output(records, out_path, format, view_intermediate)
    return


if __name__ == "__main__":
    app.run()
