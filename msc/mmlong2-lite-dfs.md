## Description of mmlong2-lite custom dataframes

### bins_[output_name].tsv

| Category | Description |
| --- | --- |
| bin | Genome bin ID |
| completeness_checkm1 | Genome bin completeness estimate, reported by CheckM |
| contamination_checkm1 | Genome bin contamination estimate, reported by CheckM |
| strain_heterogeneity_checkm1 | Genome bin strain heterogeneity index, reported by CheckM |
| completeness_checkm2 | Genome bin completeness estimate, reported by CheckM2 |
| contamination_checkm2 | Genome bin contamination estimate, reported by CheckM2 |
| coding_density | Genome bin gene coding density, reported by CheckM2 |
| genome_size | Genome bin size (bp), reported by CheckM2 |
| contigs | Number of contigs, reported by CheckM2 |
| gc |  Genome bin guanine-cytosine content (%), reported by Quast |
| contig_n50 | Genome bin N50 (bp), reported by Quast |
| contig_n90 | Genome bin N90 (bp), reported by Quast |
| [aun](http://lh3.github.io/2020/04/08/a-new-metric-on-assembly-contiguity) | Nx area under the curve, reported by Quast |
| n_per_100kb | Rate of Ns in a genome bin per 100 kb, reported by Quast |
| cov | Genome bin coverage, reported by CoverM |
| r_abund | Genome bin relative abundance (%), reported by CoverM  |
| wf_name | Workflow output name |
| wf_mode | Workflow mode |
| wf_v | Workflow version |
| wf_date | Date of workflow completion |

[//]: # (Written by Mantas Sereika)
