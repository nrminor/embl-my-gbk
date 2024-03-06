#!/bin/bash

pip install pip install git+https://github.com/NBISweden/EMBLmyGFF3.git --upgrade

EMBLmyGFF3 --expose_translations

# modify the json as needed

EMBLmyGFF3 \
mafa_input/29766_PacBio125_Mafa_annotated_corrected.gff \
mafa_input/29766_PacBio125_Mafa_annotated.fasta \
--topology linear \
--molecule_type "genomic DNA" \
--transl_table 1  \
--species 'Macaca fascicularis' \
--locus_tag MHC \
--project_id PRJXXXXXXX \
| grep -vE \
'/note|/locus_tag|/transl_table|RT   |RL   |OC   |PR   |RN   |RP   |RG   |AC   XXX|DT   ' \
| sed 's/\* _/  /g' \
> result.embl
