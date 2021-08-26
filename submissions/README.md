# Submssions
This directory contains the scripts and data to prepare assembled consensus genomes for submission to GISAID, GenBank, and WA DOH.

## Setup
1. Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) if you don't already have [Conda](https://docs.conda.io/en/latest/) installed
1. Clone this repo:
    ```
    git clone https://github.com/seattleflu/assembly.git
    ```
1. Create conda environment for submissions:
    ```
    conda env create -f ./envs/submissions.yaml
    ```
1. Follow [instructions](../docs/SFS-assembly-steps.md#globus) to set up Globus account and connection.
    - _Note_: make sure to **not** install with "High Assurance" option since the NWGC endpoint does not support it
1. Install [Docker](https://docs.docker.com/get-docker/)
1. Follow [instructions](https://github.com/seattleflu/documentation/wiki/Linelists#connect-to-the-production-id3c-database) to set up connection to ID3C.
    - Ask `@dev-team` to grant your user the `assembly-exporter` role within ID3C
1. Create an account for [GISAID](https://www.gisaid.org/)
1. Create an account for [NCBI](https://www.ncbi.nlm.nih.gov/) (You can create an account linked to your UW net id)
1. Ping Jover to get added to the NCBI SFS submissions group.
1. Ping Jover to reach out to WA DOH to set up new account for SFT.


## Data

Files stored in the `source-data/` directory are:

* `authors/` subdirectory contains separate files that lists all the authors for each sample source.
These files are named with their respective sources:
    * `altius.txt`
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
This will be used to flag samples with counties that are from outside of Washington.


## Submission Process

Example commands and filenames are based on sequence flow cell `AAAKJKHM5` released on `2021-07-01`

1. A member of NWGC will ping `#data-transfer-nwgc` channel in SFS Slack when a new sequence batch is available on Globus.
1. Transfer assembly results from [NWGC Globus endpoint](https://app.globus.org/file-manager?origin_id=178d2980-769b-11e9-8e59-029d279f7e24&origin_path=%2Fseattle_flu_project%2Fivar_releases%2F) to the Fred Hutch rhino cluster
1. Create directory for unzipping the .tar.gz file
    ```
    mkdir /fh/fast/bedford_t/seattleflu/ivar-releases/20210701_fastq
    ```
1. Extract .tar.gz file:
    ```
    tar -xvzf AAAKJKHM5.tar.gz -C /fh/fast/bedford_t/seattleflu/ivar-releases/20210701_fastq
    ```
    _Skip the next two steps if you are not working from the FH rhino cluster_
1. Create a local directory to hold all files related to this batch of sequences:
    ```
    mkdir ~/Documents/ivar-releases/Batch-20210701
    ```
1. Copy the FASTA and metrics TSV to local directory
    ```
    scp joverlee@rhino.fhcrc.org:/fh/fast/bedford_t/seattleflu/ivar-releases/20210701_fastq/AAAKJKHM5.fa ~/Documents/ivar-releases/Batch-20210701/
    scp joverlee@rhino.fhcrc.org:/fh/fast/bedford_t/seattleflu/ivar-releases/20210701_fastq/AAAKJKHM5.metrics.tsv ~/Documents/ivar-releases/Batch-20210701/
    ```
1. Drag and drop the FASTA file to [NextClade](https://clades.nextstrain.org/)
1. Do a quick scan of the sequences on NextClade to note any abnormal sequences:
    - sequences assigned to early 19A/19B clades should be twist positive controls or sequences that have been masked with majority Ns
    - sequences that are VoCs should have major defining S gene mutations listed in [CoVariants](https://covariants.org/)
        - if major defining mutations are masked with Ns, add the NWGC id to `~/Documents/ivar-releases/Batch-20210701/excluded-vocs.txt`
1. Download the NextClade TSV and save as `~/Documents/ivar-releases/Batch-20210701/nextclade.tsv`.
1. Pull down the [latest Docker image for VADR](https://hub.docker.com/r/staphb/vadr) maintained by StaPH-B
    ```
    docker pull staphb/vadr
    ```
1. Run the sequences through the NCBI VADR program by running:
    ```
    ./submissions/scripts/run_vadr ~/Documents/ivar-releases/Batch-20210701/ AAAKJKHM5.fa
    ```
    - If you run into permission denied error, make the script executable by running:
        ```
        chmod +x ./submissions/scripts/run_vadr
        ```
    - This will create a sub-directory `~/Documents/ivar-releases/Batch-20210701/genbank/` with all output files from VADR
1. Drag and drop the FASTA file to [Pangolin](https://pangolin.cog-uk.io/) and click "Start Analysis"
1. When the analysis is complete, download the Pangolin lineages and save as `~/Documents/ivar-releases/Batch-20210701/pangolin.csv`
1. Download the metadata Excel file from SFS Slack and save as `~/Documents/ivar-releases/Batch-20210701/external-metadata.xlsx`
    - The metadata file is usually in the original thread in the `#data-transfer-nwgc` channel
    - Ping Machiko for the metadata file if it has not already been posted on Slack.
1. Check the metadata in the `Metadata` or first sheet:
    - All external samples should have collection dates
    - A majority of external samples should have county.
    (It's fine if some are missing county, they will be assumed to be from Washington state)
    - There should _not_ be any dates earlier than February 2020
    - The dates should _not_ all be identical
1. Report missing or erroneous metadata to Machiko to verify and correct.
1. If there are SFS/SCAN samples listed in the metadata sheet, create a CSV of SFS sample barcodes by running:
    ```
    python3 ./submissions/scripts/extract_sfs_identifiers.py \
        --metadata ~/Documents/ivar-releases/Batch-20210701/external-metadata.xlsx \
        --output ~/Documents/ivar-releases/Batch-20210701/sfs-sample-barcodes.csv
    ```
1. If there are incorrectly formatted barcodes are printed to stdout:
    - Try to find the correct barcode by looking up the associated NWGC ID in the Metabase [NWGC ID lookup query](https://backoffice.seattleflu.org/metabase/question/641).
    - Try to find the correct barcode by looking up the barcode in the Metabase [Unknown barcode query](https://backoffice.seattleflu.org/metabase/question/439)
    - Ping Machiko to confirm and correct barcodes in the Excel metadata file.
    - Re-run the previous step after correcting Excel metadata file
1. Export SFS/SCAN metadata from ID3C by running:
    ```
    PGSERVICE=seattleflu-production ./submissions/scripts/export_id3c_metadata \
        ~/Documents/ivar-releases/Batch-20210701/sfs-sample-barcodes.csv > ~/Documents/ivar-releases/Batch-20210701/id3c-metadata-with-county.csv
    ```
1. Check the ID3C metadata include collection dates for all samples.
    - collection date may be missing if metadata has not been ingested into ID3C yet
    - ping `@dev-team` in `#informatics` with SFS sample barcodes if missing collection dates
    - submission of sequences missing collection date will have to be delayed
1. Download the latest TSV of sequences from [GitHub](https://github.com/seattleflu/hcov19-sequence-identifiers/blob/master/hcov19-sequence-identifiers.tsv) and save as `~/Documents/ivar-releases/previous-submissions.tsv`
1. Create submission files by running:
    ```
    python3 ./submissions/scripts/create_submissions.py \
        --batch-name 20210701 \
        --metadata ~/Documents/ivar-releases/Batch-20210701/external-metadata.xlsx \
        --id3c-metadata ~/Documents/ivar-releases/Batch-20210701/id3c-metadata-with-county.csv \
        --metrics ~/Documents/ivar-releases/Batch-20210701/AAAKJKHM5.metrics.tsv \
        --nextclade ~/Documents/ivar-releases/Batch-20210701/nextclade.tsv \
        --pangolin ~/Documents/ivar-releases/Batch-20210701/pangolin.csv \
        --previous-submissions ~/Documents/ivar-releases/previous-submissions.tsv \
        --strain-id 9161
        --fasta ~/Documents/ivar-releases/Batch-20210701/AAAKJKHM5.fa \
        --gisaid-username joverlee
        --output-dir ~/Documents/ivar-releases/Batch-20210701/ \
        --excluded-vocs ~/Documents/ivar-releases/Batch-20210701/excluded-vocs.txt \
        --vadr-dir ~/Documents/ivar-releases/Batch-20210701/genbank
    ```
    __NOTE__: Replace both the GISAID username and the strain ID.
        - The strain ID should be the next integer after the last strain name in the previous submissions.
1. Share reports of SFS VoCs to study point people following the VoC Reporting SOP pinned in `#sequencing`
    - reports are CSV files with the study name in the filename, e.g. `20210701_SCAN_vocs.csv`
1. Post total VoC counts in `20210701_total_vocs.csv` to `#sequencing`
1. Submit sequences to [GISAID](https://www.gisaid.org/).
    - Navigate to the EpiCoV:tm: tab.
    - Click on "Upload" and select "Batch Upload" in the pop-up.
    - Upload `SFS_20210701_EpiCoV_BulkUpload.csv` and `SFS_20210701_EpiCoV_BulkUpload.fasta` files.
1. Submit sample data to [BioSample](https://submit.ncbi.nlm.nih.gov/subs/biosample/)
    - Follow instructions in [NCBI protocol](https://www.protocols.io/view/sars-cov-2-ncbi-submission-protocol-sra-biosample-bui7nuhn) (ignore instructions for SRA submissions, we are currently not submitting to SRA)
    - Upload `20210701_biosample.tsv` to the "Attributes" tab
1. If there were SCH samples in this batch of sequences, email SCH the sequencing results:
    - Email `Batch_20210701_SCH_sequencing_results.xlsx` to SCH as an encrypted email.
    - SCH will fill in the original lab accession id and upload to WA DOH SFT
1. Submit sequencing results to [WA DOH SFT](https://sft.wa.gov/)
    - Upload `Batch_20210701_sequencing_results.xlsx` to the `NW_Genomics` folder
1. BioSample will send an email when accessions are available for download.
    - Go to the ["My Submissions" page](https://submit.ncbi.nlm.nih.gov/subs/)
    - Download the "attributes file with BioSample accessions" and save as `~/Documents/ivar-releases/Batch-20210701/biosample_accessions.tsv`
1. Link BioSample accessions to the GenBank submissions metadata by running:
    ```
    python3 ./submissions/scripts/add_biosample_accessions.py \
        --biosample ~/Documents/ivar-releases/Batch-20210701/biosample_accessions.tsv \
        --genbank ~/Documents/ivar-releases/Batch-20210701/20210701_*_genbank_metadata.tsv
    ```
1. Submit sequences to [GenBank](https://submit.ncbi.nlm.nih.gov/subs/genbank/)
    - Follow instructions in [NCBI protocol](https://www.protocols.io/view/sars-cov-2-ncbi-consensus-submission-protocol-genb-bid7ka9n?step=2)
    - Create a separate submission for each submission group (SFS, SCAN, WA DOH) since they have different authors
    - There should be a separate TSV and FASTA file for each submission group, e.g. `20210701_scan_genbank_metadata_with_biosample.tsv` and `20210701_scan_genbank.fasta`
1. GISAID will send an email with accession numbers when the sequences have been published.
Search for these accession numbers to get mapping between strain name and GISAID accessions:
    - Navigate to the EpiCoV:tm: tab.
    - Click on "Search"
    - Click on "Select" in the lower right hand corner
    - Copy and paste list of accession numbers into the pop-up box
    - Click "OK" to select these sequences
    - Click "Download" in the lower right hand corner
    - Select "Patient status metadata" and click "Download" in the pop-up box
    - Save file as `~/Documents/ivar-releases/Batch-20210701/gisaid_accessions.tsv`
1. GenBank will send an email when sequences have been published.
    - Go to the GenBank submission portal
    - Find the submissions for this batch and download the "AccessionReport.tsv" for each submission
    - Save the files as `~/Documents/ivar-releases/Batch-20210701/genbank_accessions_{1|2|3}.tsv`
1. Merge accessions with identifiers by running:
    ```
    python3 ./submissions/scripts/merge_identifiers.py \
        --identifiers ~/Documents/ivar-releases/Batch-20210701identifiers.tsv \
        --gisaid-accessions ~/Documents/ivar-releases/Batch-20210701/gisaid_accessions.tsv \
        --genbank-accessions ~/Documents/ivar-releases/Batch-20210701/genbank_accessions_1.tsv ~/Documents/ivar-releases/Batch-20210701/genbank_accessions_2.tsv ~/Documents/ivar-releases/Batch-20210701/genbank_accessions_3.tsv
        --output ~/Documents/ivar-releases/Batch-20210701/identifiers_with_accessions.tsv
    ```
1. Copy the identifiers from `~/Documents/ivar-releases/Batch-20210701/identifiers_with_accessions.tsv` and append them to the TSV in the
[hcov19 identifiers repo](https://github.com/seattleflu/hcov19-sequence-identifiers)
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
