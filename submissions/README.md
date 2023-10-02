# Submissions
This directory contains the scripts and data to prepare assembled consensus genomes for submission to GISAID, GenBank, WA DOH, and Oregon HA.

## Setup
1. Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) if you don't already have [Conda](https://docs.conda.io/en/latest/) installed
1. Clone this repo:
    ```
    git clone https://github.com/seattleflu/assembly.git
    ```
1. Clone these repos:
    ```
    git clone https://github.com/seattleflu/hcov19-sequence-identifiers.git
    git clone https://github.com/seattleflu/rsv-flu-sequence-identifiers.git
    ```
    Ask the `@dev-team` if permission is restricted
1. Create conda environment for submissions:
    ```
    conda env create -f ./envs/submissions.yaml
    ```
1. Follow [instructions](../docs/SFS-assembly-steps.md#globus) to set up Globus account and connection.
    - _Note_: make sure to **not** install with "High Assurance" option since the NWGC endpoint does not support it
1. Install [Docker](https://docs.docker.com/get-docker/)
1. Follow [instructions](https://github.com/seattleflu/documentation/wiki/Linelists#connect-to-the-production-id3c-database) to set up connection to ID3C.
    - Ask `@dev-team` to grant your user the `assembly-exporter` role within ID3C
1. Ask `@dev-team` to grant and provide instructions for access to LIMS API.
1. Create an account for [GISAID](https://www.gisaid.org/)
1. Create an account for [NCBI](https://www.ncbi.nlm.nih.gov/) (You can create an account linked to your UW net id)
1. Ask `@dev-team` to be added to the NCBI SFS submissions group.
1. Ask `@dev-team` to reach out to WA DOH to set up new account for SFT.


## Data

Files stored in the `source-data/` directory are:

* `authors/` subdirectory contains separate files that lists all the authors for each sample source.
These files are named with their respective sources:
    * `altius.txt`
    * `cascadia.txt`
    * `scan.txt`
    * `sfs.txt`
    * `wa-doh.txt`

* `lab_addresses.tsv` contains the lab name and addresses for external labs and SFS.
These are required for GISAID submissions.

* `manual_location_annotations.tsv` contains manual curations of location data from external labs.
External labs only give the sample's county, but not the state.

* `us_state_abbreviations.csv` contains US states and their abbreviations.
This is used to create the appropriate strain name based on the sample's original state.

* `variants_of_concern.tsv` contains mapping of variants of concerns for the clade name from Nextclade, WHO, and Pangolin.

* `washington_counties.txt` contains a list of all counties in the state of Washington.

* `oregon_counties.txt` contains a list of all counties in the state of Oregon.


## Submission Process

1. Transfer the assembly results from S3 to you local machine:
    Consensus genomes:
    - [s3-folder]/assembly/consensus_genomes/sars-cov-2.all-samples.masked_consensus.fasta
    - [s3-folder]/assembly/consensus_genomes/rsv-a.all-samples.masked_consensus.fasta
    - [s3-folder]/assembly/consensus_genomes/rsv-b.all-samples.masked_consensus.fasta
    - [s3-folder]/assembly/consensus_genomes/flu-a-h1n1.all-samples.masked_consensus.fasta
    - [s3-folder]/assembly/consensus_genomes/flu-a-h3n2.all-samples.masked_consensus.fasta
    - [s3-folder]/assembly/consensus_genomes/flu-b.all-samples.masked_consensus.fasta

    Picard metrics:
    - [s3-folder]/assembly/summary/picard/sars-cov-2/*
    - [s3-folder]/assembly/summary/picard/rsv-a/*
    - [s3-folder]/assembly/summary/picard/rsv-b/*
    - [s3-folder]/assembly/summary/picard/flu-a-h1n1/*
    - [s3-folder]/assembly/summary/picard/flu-a-h3n2/*
    - [s3-folder]/assembly/summary/picard/flu-b/*

1. Run file prep script
    ```
    python3 ./submissions/scripts/prep_files_mip.py \
        --fasta [path to FASTA file]
        --metrics [Picard metrics files]
        --output-dir [path to output directory]
        --batch-date YYYYMMDD
        --metadata-file [path to external metadata XLSX file]
        --pathogen [pathogen]
        --subtype [subtype]
    ```
    - pathogen values: `sars-cov-2`, `rsv-a`, `rsv-b`, `flu-a`, or `flu-b`
    - subtype values (only applicable to `flu-a`, otherwise omit): `h1n1`, `h3n2`

1. Follow prompts to include/exclude specific sets of records or individual file preparation steps as needed (use defaults for standard use case).

1. Review outputs of file prep script. Investigate missing or erroneous metadata.

1. Check that LIMS and ID3C metadata files include collection dates for all samples (open the file in a NON-Excel text editor).
    - collection date may be missing if metadata has not been ingested into ID3C yet
    - ping `@dev-team` in `#informatics` with SFS sample barcodes if missing collection dates
    - submission of sequences missing collection date will have to be delayed

1. Download the latest TSV of sequences from Github:

     [SARS-CoV-2](https://github.com/seattleflu/hcov19-sequence-identifiers/blob/master/hcov19-sequence-identifiers.tsv): save as `[output-folder]/sars-cov-2-previous-submissions.tsv`

     or

     [RSV and flu](https://github.com/seattleflu/rsv-flu-sequence-identifiers/blob/master/rsv-flu-sequence-identifiers.tsv): save as `[output-folder]/rsv-flu-previous-submissions.tsv`

1. Create submission files by running:

    SARS-CoV-2:
    ```
    python3 ./submissions/scripts/create_submissions_sars_cov_2.py \
        --batch-name [YYYYMMDD] \
        --metadata [batch-dir]/external-metadata.xlsx \
        --id3c-metadata [batch-dir]/id3c-metadata-with-county.csv \
        --metrics [batch-dir]/[pathogen]-metrics.tsv \
        --nextclade [batch-dir]/nextclade-[pathogen].tsv \
        --previous-submissions [batch-dir]/sars-cov-2-previous-submissions.tsv \
        --previous-submissions-rsv-flu [batch-dir]/rsv-flu-previous-submissions.tsv \
        --strain-id #####
        --fasta [batch-dir]/[FASTA-filename].fa \
        --gisaid-username [your username]
        --output-dir [batch-dir]/ \
        --excluded-vocs [batch-dir]/excluded-vocs.txt \
        --test-name MIPsSEQ \
        --vadr-dir [batch-dir]/genbank-[pathogen]
    ```

    RSV:
    ```
    python3 ./submissions/scripts/create_submissions_rsv.py \
        --batch-name [YYYYMMDD] \
        --metadata [batch-dir]/external-metadata.xlsx \
        --id3c-metadata [batch-dir]/id3c-metadata-with-county.csv \
        --lims-metadata [batch-dir]/lims-metadata-with-county.csv \
        --metrics [batch-dir]/[pathogen]-metrics.tsv \
        --nextclade [batch-dir]/nextclade-[pathogen].tsv \
        --previous-submissions [batch-dir]/rsv-flu-previous-submissions.tsv \
        --previous-submissions-sars-cov-2 [batch-dir]/sars-cov-2-previous-submissions.tsv \
        --strain-id #####
        --fasta [batch-dir]/[FASTA-filename].fa \
        --gisaid-username [your username]
        --output-dir [batch-dir]/ \
        --test-name MIPsSEQ \
        --pathogen [pathogen] \
        --vadr-dir [batch-dir]/genbank-[pathogen]
    ```
    - pathogen values: `rsv-a`, `rsv-b`

    Influenza:
    ```
    python3 ./submissions/scripts/create_submissions_flu.py \
        --batch-name [YYYYMMDD] \
        --metadata [batch-dir]/external-metadata.xlsx \
        --id3c-metadata [batch-dir]/id3c-metadata-with-county.csv \
        --lims-metadata [batch-dir]/lims-metadata-with-county.csv \
        --metrics [batch-dir]/metrics.tsv \
        --nextclade [batch-dir]/nextclade.tsv \
        --previous-submissions [batch-dir]/rsv-flu-previous-submissions.tsv \
        --previous-submissions-sars-cov-2 [batch-dir]/sars-cov-2-previous-submissions.tsv \
        --strain-id #####
        --fasta [batch-dir]/AAAKJKHM5.fa \
        --gisaid-username [your username]
        --output-dir [batch-dir]/ \
        --test-name MIPsSEQ \
        --pathogen [pathogen] \
        --subtype [subtype]
    ```
    - pathogen values: `flu-a`, `flu-b`
    - subtype values (only applicable to `flu-a`): `h1n1`, `h3n2`

    __NOTE__: Replace the batch-dir, ,GISAID username, the strain ID, pathogen and subtype (where applicable).
        - The strain ID should be the next integer after the last strain name in the previous submissions, and there should be no overlap in strain names beween the two previous submission files.

1. Share reports of SFS SARS-CoV-2 VoCs to study point people following the SARS-CoV-2 VoC Reporting SOP pinned in `#sequencing`
    - reports are CSV files with the study name in the filename, e.g. `YYYYMMDD_SCAN_vocs.csv`

1. Post total SARS-CoV-2 VoC counts in `YYYYMMDD_total_vocs.csv` to `#sequencing`

1. Share count of >10% Ns, control, failed, and submitted and attach files, e.g `YYYYMMDD_sample_status.csv` to `#assembly`

1. For SARS-CoV-2, submit sequences to [GISAID](https://www.gisaid.org/).
    - Navigate to the EpiCoV:tm: tab.
    - Click on "Upload" and select "Batch Upload" in the pop-up.
    - Upload `SFS_YYYYMMDD_EpiCoV_BulkUpload.csv` and `SFS_YYYYMMDD_EpiCoV_BulkUpload.fasta` files.
    - Select "Notify me only about NOT PREVIOUSLY REPORTED FRAMESHIFTS in this submission for reconfirmation of affected sequences"

1. Submit sample data to [BioSample](https://submit.ncbi.nlm.nih.gov/subs/biosample/)
    - Follow instructions in [NCBI protocol](https://www.protocols.io/view/sars-cov-2-ncbi-submission-protocol-sra-biosample-14egn8ydmg5d/v5) (ignore instructions for SRA submissions, we are currently not submitting to SRA)
    - Upload `YYYYMMDD_biosample.tsv` to the "Attributes" tab
1. If there were SCH samples in this batch of sequences, email SCH the sequencing results:
    - Email `Batch_YYYYMMDD_SCH_sequencing_results.xlsx` to SCH as an encrypted email.
    - SCH will fill in the original lab accession id and upload to WA DOH SFT
1. Submit sequencing results to [WA DOH SFT](https://sft.wa.gov/)
    - Upload `Batch_YYYYMMDD_sequencing_results.xlsx` to the `NW_Genomics` folder
1. BioSample will send an email when accessions are available for download.
    - Go to the ["My Submissions" page](https://submit.ncbi.nlm.nih.gov/subs/)
    - Download the "attributes file with BioSample accessions" and save as `~/Documents/ivar-releases/Batch-YYYYMMDD/biosample_accessions.tsv`
1. Link BioSample accessions to the GenBank submissions metadata by running:
    ```
    python3 ./submissions/scripts/add_biosample_accessions.py \
        --biosample [batch-dir]/biosample_accessions.tsv \
        --genbank [batch-dir]/YYYYMMDD_*_genbank_metadata.tsv
    ```
1. Submit sequences to [GenBank](https://submit.ncbi.nlm.nih.gov/subs/genbank/)
    - Follow instructions in [NCBI protocol](https://www.protocols.io/view/sars-cov-2-ncbi-consensus-submission-protocol-genb-n92ldy1w9l5b/v3)
    - Create a separate submission for each submission group (SFS, SCAN, WA DOH, Cascadia) since they have different authors
    - There should be a separate TSV and FASTA file for each submission group, e.g. `YYYYMMDD_scan_genbank_metadata_with_biosample.tsv` and `YYYYMMDD_scan_genbank.fasta`
1. GISAID will send an email with accession numbers when the sequences have been published.
Search for these accession numbers to get mapping between strain name and GISAID accessions:
    - Navigate to the EpiCoV:tm: tab.
    - Click on "Search"
    - Click on "Select" in the lower right hand corner
    - Copy and paste list of accession numbers into the pop-up box
    - Click "OK" to select these sequences
    - Click "Download" in the lower right hand corner
    - Select "Patient status metadata" and click "Download" in the pop-up box
    - Save file as `[batch-dir]/gisaid_accessions.tsv`
1. GenBank will send an email when sequences have been published.
    - Go to the GenBank submission portal
    - Find the submissions for this batch and download the "AccessionReport.tsv" for each submission
    - Save the files as `[batch-dir]/genbank_accessions_{1|2|3}.tsv`
1. Merge accessions with identifiers by running:
    ```
    python3 ./submissions/scripts/merge_identifiers.py \
        --identifiers [batch-dir]/identifiers.tsv \
        --gisaid-accessions [batch-dir]/gisaid_accessions.tsv \
        --genbank-accessions [batch-dir]/genbank_accessions_*.tsv
        --output [batch-dir]/identifiers_with_accessions.tsv
    ```
1. Copy the identifiers from `[batch-dir]/identifiers_with_accessions.tsv` and append them to the TSV in the
[hcov19 identifiers repo](https://github.com/seattleflu/hcov19-sequence-identifiers) or [rsv-flu identifiers repo](https://github.com/seattleflu/rsv-flu-sequence-identifiers)
1. Push the new identifiers up the master branch of the repo.
This will email Krisandra from WA DOH that new sequence identifiers have been added.


## Submission Errors

### GISAID
GISAID excludes sequences that have frameshifts in the bulk submission.
They usually email a list of sequences that have been excluded in the submissions.
1. Inspect the variants TSV generated by NWGC (e.g. `/fh/fast/bedford_t/seattleflu/ivar-releases/20210701_fastq/variants/793380.merged.trimmed.sorted.markeddups.tsv`) to confirm the frameshifts of the listed sequences.
1. Contact NWGC to investigate sequences that have frameshifts not listed in the variants TSV or if
the deletion listed in the TSV does not match the GISAID's reported frameshift.
1. Save the list of sequence strain names in a text file,
e.g. `~/Documents/ivar-releases/Batch-20210701/gisaid-errors/20210701-resubmit-subset.txt`
1. Pull out the metadata and sequences that have been confirmed from previous submission files by running:
    ```
    python3 ./submissions/scripts/extract_submission_subset.py \
        --subset ~/Documents/ivar-releases/Batch-20210701/gisaid-errors/resubmit-subset.txt \
        --metadata ~/Documents/ivar-releases/Batch-20210701/SFS_20210701_EpiCoV_BulkUpload.csv \
        --metadata-format GISAID \
        --fasta ~/Documents/ivar-releases/Batch-20210701/SFS_20210701_EpiCoV_BulkUpload.fasta \
        --output-dir ~/Documents/ivar-releases/Batch-20210701/gisaid-errors/
    ```
1. Resubmit sequences to GISAID by uploading the subset files
`SFS_20210701_EpiCoV_BulkUpload-subset.csv` and `SFS_20210701_EpiCoV_BulkUpload-subset.fasta`
1. Reply to GISAID's email to confirm that the frameshifts have been verified and the sequences have been resubmitted.
1. Update the [hcov19 identifiers repo](https://github.com/seattleflu/hcov19-sequence-identifiers) to mark the rejected sequences as "pending" under `status`.
This will make it easier to track with samples are still be investigated for frameshift or assembly anomalies.
### GenBank
- Before using the VADR program, GenBank would often flag sequences that have errors.
- The same `extract_submission_subset.py` script can be used to pull out subset
of GenBank metadata and sequences for resubmissions by passing the option `--metadata-format GenBank`.
- This is no longer required since VADR failed sequences are excluded from submissions.
