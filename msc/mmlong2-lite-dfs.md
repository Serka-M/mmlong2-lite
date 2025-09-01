## Description of mmlong2-lite custom dataframes

### Column names for <output_name>_bins.tsv

| Category | Description |
| --- | --- |
| bin | Genome bin (MAG) identifier |
| completeness_checkm1 | Genome bin completeness estimate, reported by CheckM |
| contamination_checkm1 | Genome bin contamination estimate, reported by CheckM |
| strain_heterogeneity_checkm1 | Genome bin strain heterogeneity index, reported by CheckM |
| completeness_checkm2 | Genome bin completeness estimate, reported by CheckM2 |
| contamination_checkm2 | Genome bin contamination estimate, reported by CheckM2 |
| contigs | Number of contigs, reported by Quast |
| longest_contig | Longest contig length (bp), reported by Quast |
| genome_size | Genome bin size (bp), reported by Quast |
| gc | Genome bin guanine-cytosine content (%), reported by Quast |
| contig_n50 | Genome bin N50 (bp), reported by Quast |
| contig_n90 | Genome bin N90 (bp), reported by Quast |
| [aun](http://lh3.github.io/2020/04/08/a-new-metric-on-assembly-contiguity) | Nx area under the curve, reported by Quast |
| n_per_100kb | Rate of Ns in a genome bin per 100 kbp, reported by Quast |
| cov | Genome bin coverage, reported by CoverM |
| r_abund | Genome bin relative abundance (%), reported by CoverM |
| wf_name | Workflow output name |
| wf_read_mode | Workflow long-read mode |
| wf_binning_mode | Workflow binning mode |
| wf_v | Workflow version |
| wf_date | Date of workflow completion |

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
| max_pss | Maximum Proportional Set Size (MB) - peak shared memory, divided equally among processes sharing it |
| io_in | The number of MB read (cumulative) |
| io_out | The number of MB written (cumulative) |
| mean_load | CPU usage over time, divided by the total running time |
| cpu_time | Total CPU time (user + system), in seconds |

### Row names for column `stage` in <output_name>_usage.tsv
| stage | Description |
| assembly | Metagenomic assembly |
| polishing | Optional polishing of the assembled metagenome |
| filtering | Metagenome filtering based on contig length and domain-level classification |
| singletons | Recovery of single-contig genomes |
| coverage | Read mapping and contig coverage calculation |
| binning | Recovery of multi-contig metagenomic bins |
| summary | Summary statistics for the recovered genomes |

### Row names for column `step` in <output_name>_usage.tsv
| stage | step | Description |
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
| filtering | length| Metagenome filtering based on contig length |
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
| binning | binette | Ensemble bin selection with Binette (per-round) |
| binning | prep-comebin-innit | Preparation for contig binning with COMEBin (per-round) |
| binning | prep-comebin | Read mapping preparation for contig binning with COMEBin (per-dataset, per-round) |
| binning | comebin | Metagenomic binning with COMEBin (per-round) |
| binning | qc | Quality control of the recovered bins (per-round) |
| binning | aggregate | Aggregation of the recovered bins |
| binning | qc2 | Secondary quality control of the recovered bins |
| summary | cov | Genome coverage and relative abundance calculation with CoverM |
| summary | stats | Genome summary statistic calculation with QUAST |

[//]: # (Written by Mantas Sereika)
