## Description of mmlong2-lite custom dataframes

### bins_[output_name].tsv

| Category | Description |
| --- | --- |
| bin | Genome bin ID |
| Completeness | Genome bin completeness estimate, reported by CheckM2 |
| Contamination | Genome bin contamination estimate, reported by CheckM2 |
| Completeness_Model_Used | Model used by CheckM2 |
| Coding_Density | Genome bin gene coding density, reported by CheckM2 |
| Contig_N50 | Genome bin N50 (bp), reported by CheckM2 |
| Average_Gene_Length | Average gene length, reported by CheckM2 |
| Genome_Size | Genome bin size (bp), reported by CheckM2 |
| GC_Content | Genome bin guanine-cytosine content (fraction), reported by CheckM2 |
| Total_Coding_Sequences | Number of coding sequences, reported by CheckM2 |
| contigs | Number of contigs, reported by Quast |
| longest_contig | Length of the longest contig (bp), reported by Quast |
| N90 | Genome bin N90 (bp), reported by Quast |
| [auN](http://lh3.github.io/2020/04/08/a-new-metric-on-assembly-contiguity) | Nx area under the curve, reported by Quast |
| N_per_100kb | Rate of Ns in a genome bin per 100 kb, reported by Quast |
| cov | Average genome bin coverage, reported by CoverM |
| r_abund | Average genome bin relative abundance (%), reported by CoverM  |

[//]: # (Written by Mantas Sereika)
