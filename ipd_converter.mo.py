import marimo

__generated_with = "0.2.13"
app = marimo.App()


@app.cell
def __():
    import marimo as mo
    import Bio
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    import xmltodict
    return (
        Bio,
        FeatureLocation,
        Seq,
        SeqFeature,
        SeqIO,
        SeqRecord,
        mo,
        xmltodict,
    )


@app.cell
def __():
    gb_path = "mafa_input/29766_PacBio125_Mafa_annotated.gb"
    out_path = ""
    species = "Macaca fascicularis"
    tmp_gbk = "tmp1.gb"
    return gb_path, out_path, species, tmp_gbk


@app.cell
def __(FeatureLocation, species):
    def quarterback_tweaks(_record):

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

        # assign the correct species to the organism annotation
        _record.annotations["organism"] = species

        # Check if it's a 'misc_feature' without location
        for _feature in _record.features:
            if _feature.type == "misc_feature" and not _feature.location:
                # Assign default coordinates (1..1) to misc_feature
                _feature.location = FeatureLocation(0, 1, strand=1)

        # get rid of the misc features Geneious adds in for where parts of annotations
        # were manually deleted
        _record.features = [
            _feature
            for _feature in _record.features
            if _feature.type != "misc_feature"
        ]

        # get rid of gene features
        _record.features = [
            _feature
            for _feature in _record.features
            if _feature.type != "gene"
        ]

        # Check if it's a 'CDS' or 'gene' and get rid of standard names
        for _feature in _record.features:
            if _feature.type == "CDS" or _feature.type == "gene":
                if "standard_name" in _feature.qualifiers:
                    del _feature.qualifiers["standard_name"]

        # delete geneious notes
        for _feature in _record.features:
            # Check if the feature has a 'note' qualifier and remove it
            if "note" in _feature.qualifiers:
                del _feature.qualifiers["note"]

        for _feature in _record.features:
            if _feature.type == "source":
                # Update organism and mol_type qualifiers
                _feature.qualifiers["organism"] = [species]
                _feature.qualifiers["mol_type"] = ["genomic DNA"]

        return _record
    return quarterback_tweaks,


@app.cell
def __(SeqIO, gb_path, quarterback_tweaks):
    records = [
        quarterback_tweaks(_record) for _record in SeqIO.parse(gb_path, "genbank")
    ]
    return records,


@app.cell
def __(SeqIO, records, tmp_gbk):
    # write the modified genbank
    SeqIO.write(records, tmp_gbk, "genbank")
    return


@app.cell
def __(SeqIO, tmp_gbk):
    with open("tmp1.emb", "w") as tmp_embl1:
        for _record in SeqIO.parse(tmp_gbk, "genbank"):
            SeqIO.write(_record, tmp_embl1, "embl")
    return tmp_embl1,


@app.cell
def __():
    return


if __name__ == "__main__":
    app.run()
