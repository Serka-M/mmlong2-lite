## Description of mmlong2-lite custom dataframes

### Column names for <output_name>_bins_prok.tsv

| Category | Description |
| --- | --- |
| bin | Genome bin (MAG) identifier |
| completeness_checkm1 | Genome bin completeness estimate (%), reported by CheckM |
| contamination_checkm1 | Genome bin contamination estimate (%), reported by CheckM |
| strain_heterogeneity_checkm1 | Genome bin strain heterogeneity index, reported by CheckM |
| completeness_checkm2 | Genome bin completeness estimate (%), reported by CheckM2 |
| contamination_checkm2 | Genome bin contamination estimate (%), reported by CheckM2 |
| model_checkm2 | Model used for the completeness estimate (general or specific), reported by CheckM2 |
| coding_density_checkm2 | Genome bin coding density (fraction), reported by CheckM2 |
| contigs | Number of contigs, reported by Quast |
| longest_contig | Longest contig length (bp), reported by Quast |
| genome_size | Genome bin size (bp), reported by Quast |
| gc | Genome bin guanine-cytosine content (%), reported by Quast |
| contig_n50 | Genome bin N50 (bp), reported by Quast |
| contig_n90 | Genome bin N90 (bp), reported by Quast |
| [aun](http://lh3.github.io/2020/04/08/a-new-metric-on-assembly-contiguity) | Nx area under the curve, reported by Quast |
| n_per_100kb | Rate of Ns in a genome bin per 100 kbp, reported by Quast |
| cov | Genome bin mean coverage for the main read dataset, reported by CoverM |
| r_abund | Genome bin relative abundance (%) for the main read dataset, reported by CoverM |
| wf_name | Workflow output name |
| wf_read_mode | Workflow long-read mode |
| wf_bin_mode | Workflow binning mode |
| wf_v | Workflow version |
| wf_date | Date of workflow completion |

<br/>

### Column names for <output_name>_bins_euk.tsv

| Category | Description |
| --- | --- |
| bin | Eukaryotic genome bin (MAG) identifier |
| completeness_eukcc | Genome bin completeness estimate (%), reported by EukCC |
| contamination_eukcc | Genome bin contamination estimate (%), reported by EukCC |
| lineage_eukcc | Lineage of the marker gene set used for the genome bin as dash-separated NCBI taxonomy IDs, reported by EukCC |
| lineage_qc | Status of the marker gene set quality check, reported by EukCC (PASS or FAIL) |
| lineage_name | Name of the major eukaryotic group of the genome bin (e.g. Stramenopiles), derived from the EukCC lineage |
| contigs | Number of contigs, reported by Quast |
| longest_contig | Longest contig length (bp), reported by Quast |
| genome_size | Genome bin size (bp), reported by Quast |
| gc | Genome bin guanine-cytosine content (%), reported by Quast |
| contig_n50 | Genome bin N50 (bp), reported by Quast |
| contig_n90 | Genome bin N90 (bp), reported by Quast |
| [aun](http://lh3.github.io/2020/04/08/a-new-metric-on-assembly-contiguity) | Nx area under the curve, reported by Quast |
| n_per_100kb | Rate of Ns in a genome bin per 100 kbp, reported by Quast |
| cov | Genome bin mean coverage for the main read dataset, reported by CoverM |
| r_abund | Genome bin relative abundance (%) for the main read dataset, reported by CoverM |
| wf_name | Workflow output name |
| wf_read_mode | Workflow long-read mode |
| wf_bin_mode | Workflow binning mode |
| wf_v | Workflow version |
| wf_date | Date of workflow completion |

<br/>

### Column names for <output_name>_usage.tsv

| Category | Description |
| --- | --- |
| stage | Workflow stage (sequentially sorted) |
| step | Workflow step (alphabetically sorted) |
| s | Running time in seconds |
| h:m:s | Running time in hours:minutes:seconds format |
| max_rss | Maximum Resident Set Size (MB) - peak amount of non-swapped physical memory used by the process |
| max_vms | Maximum Virtual Memory Size (MB) - total amount of virtual memory used by the process |
| max_uss | Maximum Unique Set Size (MB) - peak memory unique to the process |
| max_pss | Maximum Proportional Set Size (MB) - peak unique memory plus shared memory, divided equally among processes sharing it |
| io_in | The number of MB read (cumulative) |
| io_out | The number of MB written (cumulative) |
| mean_load | CPU usage (%) over time, divided by the total running time |
| cpu_time | Total CPU time (user + system), in seconds |
| wf_name | Workflow output name |
| wf_read_mode | Workflow long-read mode |
| wf_bin_mode | Workflow binning mode |
| wf_v | Workflow version |
| wf_date | Date of workflow completion |

<br/>

### Row names for column `stage` in <output_name>_usage.tsv
| stage | Description |
| --- | --- |
| assembly | Metagenomic assembly |
| polishing | Optional polishing of the assembled metagenome |
| curation | Assembly curation |
| filtering | Metagenome filtering based on contig length and domain-level classification |
| singletons | Recovery of single-contig genomes |
| coverage | Read mapping and contig coverage calculation |
| binning | Recovery of multi-contig prokaryotic and eukaryotic metagenomic bins |
| summary | Summary statistics for the recovered genomes |
| cleanup | Removal of temporary files |

<br/>

### Row names for column `step` in <output_name>_usage.tsv
| stage | step | Description |
| --- | --- | --- |
| assembly | myloasm | Metagenomic assembly with myloasm |
| assembly | metaflye | Metagenomic assembly with metaFlye |
| assembly | metamdbg | Metagenomic assembly with metaMDBG |
| assembly | custom | Import of custom assembly |
| polishing | prep | Preparation for polishing metagenome with Medaka |
| polishing | consensus | Generation of polishing consensus (per-batch) |
| polishing | stitch | Generation of the full consensus sequence |
| curation | map | Read multi-mapping to the assembly |
| curation | screening | Misassembly screening |
| curation | aggregate | Selection of contigs to curate |
| curation | selection | Generation of curated assembly |
| filtering | length | Metagenome filtering based on contig length |
| filtering | tiara | Detection of eukaryotic contigs with Tiara |
| filtering | whokaryote | Detection of eukaryotic contigs with Whokaryote |
| filtering | eukaryotes | Separation of eukaryotic contigs from the metagenome |
| singletons | circ | Extraction of circular contigs as separate genomes |
| singletons | lin | Extraction of linear contigs as separate genomes |
| singletons | qc | Quality control of single-contig genomes |
| coverage | prep | Preparation for read mapping to metagenome |
| coverage | map | Read mapping and contig coverage calculation (per-dataset) |
| coverage | aggregate | Aggregation of contig coverage profiles |
| binning | prep | Preparation for contig binning (per-round) |
| binning | metabat2 | Metagenomic binning with MetaBAT2 (per-round) |
| binning | semibin2 | Metagenomic binning with SemiBin2 (per-round) |
| binning | vamb | Metagenomic binning with VAMB (per-round) |
| binning | lorbin | Metagenomic binning with LorBin (per-round) |
| binning | quickbin | Metagenomic binning with QuickBin (per-round) |
| binning | binette | Ensemble bin selection with Binette (per-round) |
| binning | prep-comebin-innit | Preparation for contig binning with COMEBin (per-round) |
| binning | prep-comebin | Read mapping preparation for contig binning with COMEBin (per-dataset, per-round) |
| binning | comebin | Metagenomic binning with COMEBin (per-round) |
| binning | qc | Quality control of the recovered bins (per-round) |
| binning | prep-euk | Preparation for binning of eukaryotic contigs |
| binning | remag | Binning of eukaryotic contigs with REMAG |
| binning | eukcc_prep | Preparation for quality control of eukaryotic bins |
| binning | eukcc | Quality control of eukaryotic bins with EukCC |
| binning | eukclust | Clustering of the eukaryotic bins recovered by different binners |
| binning | eukselect | Selection of the final eukaryotic bins |
| binning | eukagr | Aggregation of the final eukaryotic bins |
| binning | aggregate | Aggregation of the recovered bins |
| binning | qc2 | Secondary quality control of the recovered bins |
| summary | cov | Genome coverage and relative abundance calculation with CoverM |
| summary | stats | Genome summary statistic calculation with Quast |
| cleanup | prebin | Removal of temporary files generated before binning |
| cleanup | postbin | Removal of temporary files generated after binning |

<br/>

### Suffixes for column `step` in <output_name>_usage.tsv
| Suffix | Example | Description |
| --- | --- | --- |
| -\<batch\> | consensus-01 | Polishing batch number (per-batch steps) |
| -\<dataset\>-\<type\> | map-1-NP | Dataset number and read type (per-dataset steps): NP - Nanopore, PB - PacBio, IL - Illumina |
| -R\<round\> | metabat2-R1 | Binning round number (per-round steps) |
| -Reuk | remag-Reuk | Binning round for eukaryotic contigs |

[//]: # (Written by Mantas Sereika)
