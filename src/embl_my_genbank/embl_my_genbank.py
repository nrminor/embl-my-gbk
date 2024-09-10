#!/usr/bin/env python3

"""
This module can be used as a library for the command line interface
`emb_my_gbk`, or as an independent python script. It takes an arbitrarily
complicated Genbank flat file (we tested on a Geneious Prime-exported file
with many added annotations and feature qualifiers) and converts it to a
clean, barebones EMBL file. It can optionally convert this EMBL file with
additional restrictions to the Immuno-Polymorphism Database flavor of the
EMBL format.

```
usage: embl_my_genbank.py [-h] --gb_path GB_PATH [--metapath METAPATH] --species SPECIES [--out_fmt OUT_FMT] [--view_intermediate] [--divide]

options:
  -h, --help            show this help message and exit
  --gb_path GB_PATH, -g GB_PATH
                        Genbank file to be converted.
  --metapath METAPATH, -m METAPATH
                        Metadata file with, at minimum, a column of allele names and a column of representative animals.
  --species SPECIES, -s SPECIES
                        Scientific name for the species under examination.
  --out_fmt OUT_FMT, -o OUT_FMT
                        Format to convert to. Can convert to EMBL, IPD_EMBL, and FASTA.
  --view_intermediate, -i
                        Whether to write out the intermediate cleaned Genbank file for inspection.
  --divide, -d          Whether to divide the batch into one file per record.
  --verbose, -v
                        Increase verbosity level (-v for WARNING, -vv for INFO, -vvv for DEBUG)
```
"""

from __future__ import annotations

import argparse
import sys
from enum import Enum, auto
from itertools import filterfalse
from pathlib import Path
from typing import TYPE_CHECKING

import polars as pl
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature
from loguru import logger

if TYPE_CHECKING:
    from Bio.SeqRecord import SeqRecord


def parse_command_line_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gb_path",
        "-g",
        type=Path,
        required=True,
        help="Genbank file to be converted.",
    )
    parser.add_argument(
        "--meta_path",
        "-m",
        type=Path,
        required=False,
        default=None,
        help="Metadata file with, at minimum, a column of allele names and a column of representative animals.",
    )
    parser.add_argument(
        "--species",
        "-s",
        type=str,
        required=True,
        help="Scientific name for the species under examination.",
    )
    parser.add_argument(
        "--out_fmt",
        "-o",
        type=str,
        required=False,
        default="IPD_EMBL",
        help="Format to convert to. Can convert to EMBL, IPD_EMBL, and FASTA.",
    )
    parser.add_argument(
        "--view_intermediate",
        "-i",
        action="store_true",
        help="Whether to write out the intermediate cleaned Genbank file for inspection.",
    )
    parser.add_argument(
        "--divide",
        "-d",
        action="store_true",
        help="Whether to divide the batch into one file per record.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase verbosity level (-v for WARNING, -vv for INFO, -vvv for DEBUG)",
    )

    return parser.parse_args()


class OutType(Enum):
    """
    Enum for preventing unusuable format requests
    from being represented.
    """

    EMBL = auto()
    IPD_EMBL = auto()
    GENBANK = auto()
    FASTA = auto()


def choose_out_type(format_choice: str) -> str:
    """
    Determines the output file format based on the user's choice.

    This function maps a user-provided string indicating the desired
    output file format to the corresponding format string used by Biopython
    for file writing. It validates the input against a predefined list of
    acceptable formats and raises an error for invalid inputs.

    Args:
        format_choice (str): A string representing the user's choice of output
        file format. Expected values correspond to the keys of the `OutType` enum
        (e.g., 'EMBL', 'IPD_EMBL', 'GENBANK', 'FASTA').

    Returns:
        str: The string identifier for the output file format as recognized by
        Biopython (e.g., 'embl', 'genbank', 'fasta').

    Raises:
        ValueError: If the provided `format_choice` does not correspond to a valid
        output file format, indicating either a typo or an unsupported format choice.

    Example:
        >>> choose_out_type('fasta')
        'fasta'
        >>> choose_out_type('genbank')
        'genbank'
        >>> choose_out_type('invalid')  # This will raise a ValueError

    Note:
        The function is case-insensitive for `format_choice` but expects it to match
        one of the predefined format keys exactly.
    """
    try:
        choice_str = OutType.__members__[format_choice.upper()]
    except KeyError as exc:
        msg = f"Invalid file type selected: {format_choice}"
        raise ValueError(msg) from exc

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
            msg = f"Invalid output file type selected: {format_choice}"
            raise ValueError(msg)
    return out_type


def handle_allele_id_placement(record: SeqRecord) -> SeqRecord:
    """
    Corrects the placement of allele identifiers in sequence records exported
    from Geneious Prime.

    Geneious Prime sometimes places allele names in a location that does not
    conform to expected standards for certain analyses or database submissions.
    This function addresses the issue by relocating the allele name from its original
    position to both the record's ID and accession fields. It ensures that the allele
    name is consistently represented across the sequence record, facilitating proper
    identification and reference.

    Args:
        record (SeqRecord): The sequence record to be corrected. This `SeqRecord` object,
        as defined by the Biopython library, represents a single sequence and its associated
        data. The record is expected to have an allele name as part of its naming convention,
        which will be utilized to update the ID and accession fields.

    Returns:
        SeqRecord: The modified `SeqRecord` object with the allele name correctly positioned
        within the ID and accession annotations. This adjustment ensures that the allele name
        is prominently and correctly displayed in the record, aligning with standard practices
        for sequence identification.

    Raises:
        AssertionError: If the allele name is not properly propagated to the accession annotation
        after the function's execution, an AssertionError is raised. This serves as a check to
        ensure that the intended modifications have been successfully applied to the sequence
        record.

    Example:
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> original_record = SeqRecord(Seq("AGTC"), id="geneious_id", name="Allele_A1")
        >>> corrected_record = handle_allele_id_placement(original_record)
        >>> corrected_record.id
        'Allele_A1'
        >>> corrected_record.annotations["accessions"]
        ['Allele_A1']

    Note:
        This function modifies the input `SeqRecord` object in place by adjusting its ID and
        annotations. It is designed to ensure that allele names are correctly positioned within
        the sequence record, reflecting best practices for sequence data management and submission.
    """

    # parse out the working name of the allele
    allele = record.name
    logger.debug(f"Processing allele with working name {allele}")

    # move the allele name into the id position, replacing a geneious annotation
    logger.debug(
        f"Replacing the original id field, {record.id}, with the allele name {allele}"
    )
    record.id = allele

    # rename the accession to be the allele name
    logger.debug(
        f"Attempting to set the accessions annotation, {record.annotations['accessions']}, to the allele name."
    )
    record.annotations["accessions"] = [allele]

    # check that the allele name was properly propogated
    assert (
        record.annotations["accessions"][0] == record.name
    ), "Record name could not be properly assigned to accession annotation"

    return record


def assign_species(record: SeqRecord, species: str) -> SeqRecord:
    """
    Assigns a specified species name to the sequence record.

    This function updates the 'organism' annotation within the `SeqRecord` object
    to reflect the specified species' scientific name. This is particularly useful
    for ensuring consistency in species annotation across multiple sequence records,
    especially when preparing data for analysis or database submission.

    Args:
        record (SeqRecord): The sequence record to be updated. This is a `SeqRecord`
        object as defined by the Biopython library which represents a single sequence
        and its associated data, including annotations. species (str): The scientific
        name of the species to be assigned to the sequence record. This name is used
        to update the 'organism' field within the record's annotations, ensuring that
        the record accurately reflects the specified species.

    Returns:
        SeqRecord: The modified `SeqRecord` object with the 'organism' annotation updated
        to the specified species' scientific name. This ensures that the sequence record
        is correctly annotated with the appropriate species information for downstream
        applications.

    Example:
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> record = SeqRecord(Seq("ACTG"), id="example_id", description="An example sequence.")
        >>> updated_record = assign_species(record, "Homo sapiens")
        >>> updated_record.annotations["organism"]
        'Homo sapiens'

    Note:
        The function modifies the input `SeqRecord` object in place by updating its
        'organism' annotation. It also returns the modified `SeqRecord` object for
        convenience, allowing further manipulation or analysis if desired.
    """
    # assign the correct species to the organism annotation
    logger.debug(
        f"Replacing the previous organism entry, {record.annotations['organism']}, with the provided species '{species}'"
    )
    record.annotations["organism"] = species

    return record


def remove_geneious_annotations(record: SeqRecord) -> SeqRecord:
    """
    Removes Geneious-specific annotations from a sequence record.

    Geneious software is known to add certain annotations to sequence records that may not
    be desired in all contexts. This function specifically targets and removes:

    1. "misc_feature" features that Geneious adds, often related to manual editing history.
    2. "/note" qualifiers within any feature, which can include comments or metadata added
    by Geneious.

    The removal of these annotations can help in simplifying the sequence record for downstream
    analysis or for submission to databases where such annotations are not required or could
    cause confusion.

    Args:
        record (SeqRecord): The sequence record from which Geneious-specific annotations are
            to be removed. The record should be a SeqRecord object as defined by the Biopython
            library, potentially containing features and qualifiers added by Geneious software.

    Returns:
        SeqRecord: The same SeqRecord object passed as input, but modified to remove Geneious-
            specific annotations. Specifically, all "misc_feature" features are removed, and any
            "note" qualifiers within other features are deleted. This ensures the record is cleaner
            and more compliant with general sequence record standards.

    Example:
        >>> from Bio.SeqRecord import SeqRecord
        >>> from Bio.SeqFeature import SeqFeature, FeatureLocation
        >>> record = SeqRecord(seq="ACTG", id="example")
        >>> feature = SeqFeature(FeatureLocation(start=0, end=4), type="misc_feature", qualifiers={"note": "Geneious annotation"})
        >>> record.features.append(feature)
        >>> cleaned_record = remove_geneious_annotations(record)
        >>> len(cleaned_record.features)
        0

    Note:
        The function modifies the input SeqRecord object in place, but also returns it for convenience.
        This allows for the modified record to be reassigned to the same or a different variable if
        desired.
    """

    # get rid of the misc features Geneious adds in for where parts of annotations
    # were manually deleted
    logger.debug(
        f"Attempting to take {len(record.features)} features for {record.id} and remove all Geneious fields labelled 'misc_feature'"
    )
    record.features = list(
        filterfalse(lambda f: f.type == "misc_feature", record.features),
    )
    logger.debug(f"Keeping {len(record.features)} features.")

    # delete geneious notes
    for feature in record.features:
        # Check if the feature has a 'note' qualifier and remove it
        if "note" not in feature.qualifiers:
            continue
        logger.debug(
            f"Deleting the note feature qualifier, {feature.qualifiers['note']}"
        )
        del feature.qualifiers["note"]

    return record


def replace_special_chars(allele_name: str) -> str:
    """
    Replaces special characters in an allele name with underscores.

    This function takes an allele name as input and replaces asterisks (*),
    hyphens (-), and colons (:) with underscores (_). This is often necessary
    to ensure compatibility with various file systems and software that may
    have restrictions on special characters in file names or identifiers.

    Args:
        allele_name (str): The original allele name that may contain special characters.

    Returns:
        str: The modified allele name with special characters replaced by underscores.

    Example:
        >>> replace_special_chars("HLA-A*01:01")
        'HLA_A_01_01'
    """
    logger.debug(f"Replacing special characters in {allele_name}")
    no_special_chars = allele_name.replace("*", "_").replace("-", "_").replace(":", "_")
    logger.debug(
        f"Gene name will be split out with `.split('_')[1]` from {no_special_chars}"
    )
    return no_special_chars


def clean_record_features(
    record: SeqRecord, species: str, isolate_dict: None | dict[str, str] = None
) -> SeqRecord:
    """
        Geneious also exports a number of features and feature qualifiers that are
        unnecessary for our purposes, namely "standard_name" and "gene". This function does
        away with those as well, and also adds species and DNA information to the source
        feature, if available.

    Args:
        record (SeqRecord): The sequence record to be cleaned and updated. A SeqRecord object
            from Biopython that contains features which may include "gene", "CDS", and "source"
            types among others.
        species (str): The species name to be used for updating the "organism" qualifier in the
            "source" feature of the record. This string should represent the scientific name of
            the organism, e.g., "Homo sapiens".

    Returns:
        SeqRecord: The modified SeqRecord object with cleaned features according to the
            specifications. The returned record will have "gene" features removed, "standard_name"
            qualifiers deleted from "CDS" features where applicable, and the "source" feature
            updated with the provided species name and molecule type set to "genomic DNA".

    Note:
        The function directly modifies the input SeqRecord object and also returns this modified
        object. Therefore, the changes will be reflected in the original SeqRecord passed to the
        function.
    """

    # separate out non-species-specific gene name
    gene = replace_special_chars(record.name).split("_")[1]
    logger.debug(f"Parsed out gene name for {record.id}: {gene}.")

    # separate out record sequence length
    seq_length = len(record.seq)

    # filter out genes
    record.features = list(filterfalse(lambda f: f.type == "gene", record.features))

    for i, feature in enumerate(record.features):
        # skips current iteration immediately
        if feature.type not in ("CDS", "source"):
            continue

        # Remove standard_name qualifier from CDS features
        if feature.type == "CDS":
            if "standard_name" in feature.qualifiers:
                logger.debug(
                    f"Deleting the 'standard_name' field for {record.id}: {feature.qualifiers['standard_name']}"
                )
                del feature.qualifiers["standard_name"]

            # in-place mutate the allele name in the CDS, if present
            if "allele" in feature.qualifiers:
                logger.debug(
                    f"Setting the allele feature qualifier to {record.name} for {record.id}"
                )
                feature.qualifiers["allele"] = [record.name]

            # in-place mutate the allele name in the CDS
            feature.qualifiers["gene"] = gene

            # bring in an isolate ID, if provided
            if isolate_dict:
                logger.debug(
                    f"Isolate animal IDs provided. Setting '/isolate' to {isolate_dict[record.name]} for {record.id}."
                )
                feature.qualifiers["isolate"] = isolate_dict[record.name]

        # Update source feature with organism and mol_type, while making
        # sure it has the correct coordinates
        if feature.type == "source":
            feature.qualifiers["organism"] = [species]
            feature.qualifiers["mol_type"] = ["genomic DNA"]
            if (
                int(feature.location.end) != seq_length
            ):  # .end method returns a int-like type but I think the int() annotation handles errors well just in case
                logger.debug(
                    f"Mismatch in sequence length and recorded coordinates for {record.id}. Attempting to fix."
                )
                source_start = feature.location.start
                source_end = seq_length
            else:
                source_start = feature.location.start  # immutable
                source_end = feature.location.end

            source_location = FeatureLocation(source_start, source_end)
            new_source = SeqFeature(location=source_location, type="source")
            new_source.qualifiers = record.features[i].qualifiers
            record.features[i] = new_source

    return record


def read_metadata(meta_path: Path) -> dict[str, str]:
    """
    Reads metadata from an Excel file and converts it to a dictionary.

    This function reads an Excel file containing metadata, specifically looking for
    columns named 'Local designation' and 'Animal ID'. It then constructs a dictionary
    where 'Local designation' values serve as keys and 'Animal ID' values as corresponding
    values.

    Args:
        meta_path (Path): A Path object pointing to the Excel file containing metadata.

    Returns:
        Dict[str, str]: A dictionary mapping 'Local designation' to 'Animal ID'.

    Raises:
        AssertionError: If either 'Local designation' or 'Animal ID' columns are missing
        from the Excel file.

    Example:
        >>> read_metadata(Path('metadata.xlsx'))
        {'Designation1': 'ID1', 'Designation2': 'ID2', ...}

    Note:
        The function assumes the Excel file is structured with the specified columns. It
        raises an assertion error if the expected columns are not found, ensuring the
        presence of necessary data for further processing.
    """

    # read the provided table
    meta_df = pl.read_excel(source=meta_path)

    # run some checks to make sure the expected column headers are present
    assert (
        "Local designation" in meta_df.columns
    ), "Expected column named 'Local designation' is missing"
    assert (
        "Animal ID" in meta_df.columns
    ), "Expected column named 'Animal ID' is missing"

    # select the two columns
    meta_dicts = (
        meta_df.select("Animal ID", "Local designation")
        .with_columns(
            [
                pl.col("Animal ID").str.strip_chars().alias("Animal ID"),
                pl.col("Local designation")
                .str.strip_chars()
                .alias("Local designation"),
            ],
        )
        .to_dicts()
    )

    return {
        mapping["Local designation"]: mapping["Animal ID"] for mapping in meta_dicts
    }


def construct_description(
    record: SeqRecord,
    species: str,
    isolate_dict: dict[str, str],
) -> SeqRecord:
    """
    Constructs a descriptive string for a SeqRecord based on provided metadata.

    This function updates the description field of a SeqRecord object based on allele
    information, species name, and a dictionary mapping alleles to animal isolates. The
    constructed description follows a specific format incorporating these elements.

    Args:
        record (SeqRecord): The SeqRecord object to be updated with a new description.
        species (str): The scientific name of the species.
        isolate_dict (Dict[str, str]): A dictionary mapping allele names (keys) to
            animal isolate identifiers (values).

    Returns:
        SeqRecord: The updated SeqRecord object with a new description.

    Example:
        >>> from Bio.Seq import Seq
        >>> record = SeqRecord(Seq("ATGC"), id="allele_1", name="allele_1")
        >>> construct_description(record, "Homo sapiens", {"allele_1": "isolate_1"})
        SeqRecord(description="Homo sapiens gene for MHC class I antigen (allele_1 gene),
        isolate isolate_1, allele allele_1")

    Note:
        If the allele from the SeqRecord is not found in the provided isolate_dict, the function
        prints a warning message and returns the SeqRecord unchanged. This ensures the integrity
        of the SeqRecord's description in cases where metadata may be incomplete or incorrect.
    """

    # parse out the working name of the allele
    allele = record.name

    # check that the allele is actually represented in the provided dict
    if allele not in isolate_dict:
        logger.warning(f"Allele {allele} is missing in the provided table of metadata.")
        return record

    # parse out a gene name with the locus
    gene = "-".join(allele.split("_")[:2])

    # pull out the representative animal isolate for this allele
    isolate = isolate_dict.get(allele)

    # construct the description
    description = f"{species} gene for MHC class I antigen ({gene} gene), isolate {isolate}, allele {allele}"

    # add the description to the SeqRecord object
    record.description = description

    return record


def clean_genbank(
    record: SeqRecord,
    species: str,
    meta_path: Path | None,
) -> SeqRecord:
    """
    Cleans a GenBank record by updating allele IDs, assigning species names, and
    removing Geneious annotations.

    This function performs a series of cleaning and updating steps on a given GenBank
    sequence record. These steps include correcting the placement of allele IDs, assigning
    a scientific name to the species annotation, removing specific Geneious annotations,
    and cleaning various feature annotations added by Geneious software.

    Args:
        record (SeqRecord): The sequence record to be cleaned and updated.
        species (str): The scientific name of the species to be assigned to the record.

    Returns:
        SeqRecord: The cleaned and updated sequence record.

    The function is a wrapper that sequentially applies several specific cleaning and
    updating operations, making it convenient to perform multiple adjustments through a
    single function call.
    """

    # shift around the allele ID into different slots than
    # Geneious places it by default (i.e., Accession and )
    record = handle_allele_id_placement(record)

    # assign the user provided scientific name for annotation
    # onto the record
    record = assign_species(record, species)

    # clean off geneious annotations
    record = remove_geneious_annotations(record)

    # if no metadata path was provided, clean record features and finish here
    if meta_path is None:
        return clean_record_features(record, species)

    # if isolate metadata is available, use it to construct record
    # descriptions
    isolate_dict = read_metadata(meta_path)

    # clean up the many features appended to each Genbank record,
    # which BioPython represents as a dictionary
    record = clean_record_features(record, species, isolate_dict)

    return construct_description(record, species, isolate_dict)


def handle_multi_output(records: list[SeqRecord], out_type: str) -> None:
    """
    TODO
    """
    for record in records:
        allele_embl_name = f"{record.name}.tmp.{out_type}"
        with Path(allele_embl_name).open("w", encoding="utf8") as per_record:
            SeqIO.write(record, per_record, out_type)


def id_line_replacement(
    records: list[SeqRecord],
    int_embl: str,
    final_out: str,
    out_type: str,
    divide: bool,
) -> None:
    """
    Replaces the ID line in EMBL-formatted files and writes the modified
    records to a new file.

    This function first creates an intermediate file from the given sequence
    records in the specified output format. It then reads this intermediate file,
    replaces lines starting with "ID" with a predetermined new ID line, and writes
    the results to a final output file.

    Args:
        records (List[SeqRecord]): A list of SeqRecord objects to be processed.
        int_embl (str): The path for the intermediate file to be created.
        final_out (str): The path for the final output file with modified ID lines.
        out_type (str): The format of the output file (e.g., 'embl', 'genbank').

    The function is useful for adjusting file formats that require specific header
    formats not directly supported by the original exporting tool.
    """

    _new_id_line: str = "ID   XXX; XXX; linear; genomic DNA; XXX; XXX;\n"

    # make an intermediate embl file first
    with Path(int_embl).open("w", encoding="utf8") as int_handle:
        for record in records:
            record.annotations["keywords"] = ["."]
            SeqIO.write(record, int_handle, out_type)

    # split up the intermediate file if multi output mode is on
    if divide:
        with Path(int_embl).open("r", encoding="utf8") as multi_input:
            records = list(SeqIO.parse(multi_input, out_type))
            handle_multi_output(records, out_type)

    # collect the intermediate files to be corrected if ipd files are called for
    intermediate_files = [path.name for path in Path.cwd().glob(f"*tmp*.{out_type}")]

    # go through the intermediate files
    for int_file in intermediate_files:
        final_handle = (
            final_out if len(intermediate_files) == 1 else int_file.replace(".tmp", "")
        )
        # Open the input file for reading and the output file for writing
        with Path(int_file).open(encoding="utf8") as infile, Path(final_handle).open(
            "w",
            encoding="utf8",
        ) as outfile:
            for line in infile:
                # Check if the line starts with "ID"
                if line.startswith("ID"):
                    # Write the replacement text to the output file
                    new_id_line = line.replace("; DNA", "; genomic DNA")
                    outfile.write(new_id_line)
                    continue
                # Write the original line to the output file
                outfile.write(line)

        if int_file != int_embl:
            Path(int_file).unlink()
    Path(int_embl).unlink()


def write_output(
    records: list[SeqRecord],
    final_out: str,
    requested_format: str,
    view_intermediate: bool,
    divide: bool,
) -> None:
    """
    Writes sequence records to an output file in a specified format, with
    options for intermediate viewing and format adjustments.

    This function supports writing the sequence records to a file in various
    formats. It includes options for creating an intermediate file for review
    and handling special format requirements such as ID line replacements in EMBL
    files for IPC entries.

    Args:
        records (List[SeqRecord]): A list of SeqRecord objects to be written to the file.
        final_out (str): The path for the output file.
        requested_format (str): The desired format for the output file (e.g., 'EMBL',
        'IPD_EMBL', 'GENBANK', 'FASTA').
        view_intermediate (bool): Whether to create and view an intermediate file before
        final output. Useful for debugging or intermediate analysis.

    The function adapts to the requested format, providing flexibility in output file
    generation and format-specific handling.
    """

    out_type = choose_out_type(requested_format)

    if view_intermediate:
        logger.info(
            "Writing intermediate cleaned genbank file called 'intermediate.gbk' for manual inspection.",
        )
        SeqIO.write(records, "intermediate.gbk", "genbank")

    # check if an IPC intermediate embl is needed
    if "IPD" in requested_format.upper():
        int_embl = "tmp.embl"
        id_line_replacement(records, int_embl, final_out, out_type, divide)
        return

    # if multi output mode is on, write each EMBL record into its own file
    if divide:
        logger.info("Writing outputs with one file per record.")
        handle_multi_output(records, out_type)
        return

    # if multi output mode is not on, write all the records into a single file
    with Path(final_out).open("w", encoding="utf8") as out_file:
        for record in records:
            SeqIO.write(record, out_file, out_type)


def main() -> None:
    """
    Coordinate the flow of data through the embl_my_genbank functions.
    """

    # parse command line arguments
    args = parse_command_line_args()

    # set up logger based on the requested verbosity
    match args.verbose:
        case 0:
            level = "WARNING"
        case 1:
            level = "SUCCESS"
        case 2:
            level = "INFO"
        case 3:
            level = "DEBUG"
        case _:
            level = "WARNING"
    logger.remove()
    logger.add(sys.stderr, colorize=True, level=level)

    # run some checks to catch if the specified files don't exist
    assert Path(args.gb_path).is_file(), "Provided Genbank file path does not exist"
    if args.meta_path is not None:
        assert Path(
            args.gb_path,
        ).is_file(), "Provided metadata file path does not exist"

    # parse out a name for the output file based on the input file
    basename = Path(args.gb_path).name.replace(".gb", "").replace(".gbk", "")
    out_path = f"{basename}.embl"
    logger.info(f"Parsing records from {args.gb_path}.")

    # revise and clean all records in the Genbank file
    records = [
        clean_genbank(record, args.species, args.meta_path)
        for record in SeqIO.parse(args.gb_path, "genbank")
    ]
    logger.success(f"Successfully cleaned {len(records)} records.")

    # write output in the desired format
    logger.info(f"Writing cleaned and formatted records to {out_path}")
    write_output(records, out_path, args.out_fmt, args.view_intermediate, args.divide)


if __name__ == "__main__":
    main()
