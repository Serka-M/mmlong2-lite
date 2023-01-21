<p align="center">
<img align="center" width="200" height="200" src="msc/mmlong2-lite-logo.png" alt="logo" style="zoom:100%;" />
</p>

Lightweight workflow for microbial genome recovery using either Nanopore or PacBio HiFi reads. <br/>
mmlong2-lite is the microbial genome production part of the [mmlong2](https://github.com/Serka-M/mmlong2) pipeline. <br/>
<br/>
<br/>
**Core features:**
* [Snakemake](https://snakemake.readthedocs.io) workflow running dependencies from a [Singularity](https://docs.sylabs.io/guides/latest/user-guide/) container for enhanced reproducibility
* Bioinformatics tool and parameter optimizations for high complexity metagenomics samples
* Circular microbial genome extraction as separate genome bins
* Eukaryotic contig removal for reduced microbial genome contamination
* Differential coverage support for improved microbial genome recovery
* Iterative binning strategy for improved microbial genome recovery
<br/>
