# Submssions
This directory contains the scripts and data to prepare assembled consensus genomes for submission to GISAID, GenBank, and WA DOH.

## Data

Files stored in the `source-data/` directory are:

* `authors/` subdirectory contains separate files that lists all the authors for each sample source.
These files are named with their respective sources:
    * `altius.txt`
    * `scan.txt`
    * `sfs.txt`
    * `wa-doh.txt`

* `manual_location_annotations.tsv` contains manual curations of location data from external labs.
External labs only give the sample's county, but not the state.

* `us_state_abbreviations.csv` contains US states and their abbreviations.
This is used to create the appropriate strain name based on the sample's original state.

* `variants_of_concern.tsv` contains mapping of variants of concerns for the clade name from Nextclade, WHO, and Pangolin.

* `washington_counties.txt` contains a list of all counties in the state of Washington.
This will be used to flag samples with counties that are from outside of Washington.

* `washington_pumas_to_county.csv` contains a mapping of Washington PUMAs to their counties.
This file only contains PUMAs that map to **one** specific county within Washington.
This is taken from the PUMA maps listed at
https://www.census.gov/geographies/reference-maps/2010/geo/2010-pumas/washington.html

## Submission Process

Example commands and filenames are based on sequence flow cell `AAAKJKHM5` released on `2021-07-01`

1. A member of NWGC will ping `#data-transfer-nwgc` channel in SFS Slack when a new sequence batch is available on Globus.
1. Transfer assembly results from [NWGC Globus endpoint](https://app.globus.org/file-manager?origin_id=178d2980-769b-11e9-8e59-029d279f7e24&origin_path=%2Fseattle_flu_project%2Fivar_releases%2F) to the Fred Hutch rhino cluster
    - See [docs](../docs/SFS-assembly-steps.md#globus) on how to set up Globus account and connection.
1. Extract .tar.gz file on FH rhino cluster:
    ```
    tar -xvzf AAAKJKHM5.tar.gz -C /fh/fast/bedford_t/seattleflu/ivar-releases/20210701_fastq
    ```
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
1. Install [Docker](https://docs.docker.com/get-docker/)
1. Pull down the [latest Docker image for VADR](https://hub.docker.com/r/staphb/vadr) maintained by StaPH-B
    ```
    docker pull staphb/vadr
    ```
1. Run the sequences through the NCBI VADR program by running:
    ```
    ./submissions/scripts/run_vadr ~/Documents/ivar-releases/Batch-20210701/ AAAKJKHM5.fa
    ```
    - This will create a sub-directory `~/Documents/ivar-releases/Batch-20210701/genbank/` with all output files from VADR
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
1. If this is the first time exporting SFS metadata from ID3C, follow [instructions](https://github.com/seattleflu/documentation/wiki/Linelists#connect-to-the-production-id3c-database) to set up connection to ID3C.
    - Ask `@dev-team` to grant your user the `assembly-exporter` role within ID3C
1. Export SFS/SCAN metadata from ID3C by running:
    ```
    PGSERVICE=seattleflu-production ./submissions/scripts/export_id3c_metadata \
        ~/Documents/ivar-releases/Batch-20210701/sfs-sample-barcodes.csv > ~/Documents/ivar-releases/Batch-20210701/id3c-export.csv
    ```
1. Check the ID3C metadata include collection dates for all samples.
    - collection date may be missing if metadata has not been ingested into ID3C yet
    - ping `@dev-team` in `#informatics` with SFS sample barcodes if missing collection dates
    - submission of sequences missing collection date will have to be delayed
1. Convert PUMAs in ID3C metadata to counties by running:
    ```
    python3 ./submissions/scripts/match_puma_to_county.py \
        --metadata ~/Documents/ivar-releases/Batch-20210701/id3c-export.csv \
        --output-metadata ~/Documents/ivar-releases/Batch-20210701/id3c-metadata-with-county.csv
    ```
    - records missing county will be output to stdout - these can be ignored as long as it's not __all__ records missing county
1. Download the latest TSV of sequences from [GitHub](https://github.com/seattleflu/hcov19-sequence-identifiers/blob/master/hcov19-sequence-identifiers.tsv) and save as `~/Documents/ivar-releases/previous-submissions.tsv`
1. Create submission files by running:
    ```
    python3 ./submissions/scripts/create_submissions.py \
        --batch-name 20210701 \
        --metadata ~/Documents/ivar-releases/Batch-20210701/external-metadata.xlsx \
        --id3c-metadata ~/Documents/ivar-releases/Batch-20210701/id3c-metadata-with-county.csv \
        --metrics ~/Documents/ivar-releases/Batch-20210701/AAAKJKHM5.metrics.tsv \
        --nextclade ~/Documents/ivar-releases/Batch-20210701/nextclade.tsv \
        --previous-submissions ~/Documents/ivar-releases/previous-submissions.tsv \
        --strain-id 9161
        --fasta ~/Documents/ivar-releases/Batch-20210701/AAAKJKHM5.fa \
        --gisaid-username joverlee
        --output-dir ~/Documents/ivar-releases/Batch-20210701/ \
        --excluded-vocs ~/Documents/ivar-releases/Batch-20210701/excluded-vocs.txt \
        --vadr-dir ~/Documents/ivar-releases/Batch-20210701/genbank
    ```
1. Share reports of SFS VoCs to study point people following the VoC Reporting SOP pinned in `#sequencing`
    - reports are CSV files with the study name in the filename, e.g. `20210701_SCAN_vocs.csv`
1. Post total VoC counts in `20210701_total_vocs.csv` to `#sequencing`
1. Submit sequences to [GISAID](https://www.gisaid.org/).
    - Navigate to the EpiCoV:tm: tab.
    - Click on "Upload" and select "Batch Upload" in the pop-up.
    - Upload `SFS_20210701_EpiCoV_BulkUpload.csv` and `SFS_20210701_EpiCoV_BulkUpload.fasta` files.
1. Submit sequences to [GenBank](https://submit.ncbi.nlm.nih.gov/subs/genbank/)
    - Follow instructions in [NCBI protocol](https://www.protocols.io/view/sars-cov-2-ncbi-consensus-submission-protocol-genb-bid7ka9n?step=2)
    - Create a separate submission for each submission group (SFS, SCAN, WA DOH) since they have different authors
    - The `create_submissions` script should have created a separate TSV and FASTA file for each submission group, e.g. `20210701_scan_genbank_metadata.tsv` and `20210701_scan_genbank.fasta`
1. Submit sequencing results to [WA DOH SFT](https://sft.wa.gov/)
    - Upload `Batch_20210701_sequencing_results.xlsx` to the `NW_Genomics` folder
