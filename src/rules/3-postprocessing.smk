# DESCRIPTION: Snakemake rules for the postprocessing sections (bin summary statistics, intermediate file cleanup, usage statistics, dataframe generation) of the mmlong2 pipeline
# AUTHOR: Mantas Sereika (mase@bio.aau.dk)
# LICENSE: GNU General Public License

rule Summary_coverage:
    conda: config["env_1"]
    input:
        os.path.join(loc, sample, "tmp/binning/contig_bin_prok.tsv") if "TRUE" in config["prok_status"] else [],
        os.path.join(loc, sample, "tmp/binning/contig_bin_euk.tsv") if "TRUE" in config["euk_status"] else [],
    output:
        cov=os.path.join(loc, sample, "tmp/binning/bin_cov.tsv"),
        abund=os.path.join(loc, sample, "tmp/binning/bin_abund.tsv"),
    params:
        cov_method=config["cov_method"],
        ident=get_map_ident(get_read_type(mode),config["np_map_ident"],config["pb_map_ident"],config["il_map_ident"]),
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_summary_cov.tsv")
    threads: proc
    shell:
        """
        coverm genome -t {threads} -b {loc}/{sample}/tmp/binning/mapping/1-*.bam -f {loc}/{sample}/tmp/binning/bins_*/*.fa -o {output.cov} -m {params.cov_method} --min-read-percent-identity {params.ident}
        coverm genome -t {threads} -b {loc}/{sample}/tmp/binning/mapping/1-*.bam -f {loc}/{sample}/tmp/binning/bins_*/*.fa -o {output.abund} -m relative_abundance --min-read-percent-identity {params.ident}
        if ! grep -q "coverm" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "coverm " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi  
        """

rule Summary_stats:
    conda: config["env_1"]
    input:
        os.path.join(loc, sample, "tmp/binning/bin_cov.tsv"),
    output:
        os.path.join(loc, sample, "tmp/binning/quast.tsv"),
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_summary_stats.tsv")
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/quast" ]; then rm -r {loc}/{sample}/tmp/binning/quast; fi
        cat {loc}/{sample}/tmp/binning/contig_bin_*.tsv > {loc}/{sample}/tmp/binning/contig_bin.tsv
        quast.py {loc}/{sample}/tmp/binning/bins_*/*.fa -o {loc}/{sample}/tmp/binning/quast -t {threads} --no-plots --no-html --no-icarus
        cut -f1,14,15,16,17,18,19,20,23 {loc}/{sample}/tmp/binning/quast/transposed_report.tsv > {output}
        if ! grep -q "quast" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "quast " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi  
        """

rule Cleanup_prebin:
    conda: config["env_1"]
    input:
        os.path.join(loc, sample, "tmp/binning/quast.tsv"),
    output:
        os.path.join(loc, sample, "tmp/logs/cleanup_prebin.txt")
    params:
        cleanup=config["cleanup_status"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_cleanup_prebin.tsv")
    resources: usage=50
    threads: proc / 2
    shell:
        """
        if [ "{params.cleanup}" == "TRUE" ]; then

        size_pre=$(du -sb --exclude='{loc}/{sample}/tmp/binning' {loc}/{sample}/tmp | awk '{{printf "%.1f", $1/1024/1024/1024}}')
        files_pre=$(find {loc}/{sample}/tmp -path '{loc}/{sample}/tmp/binning' -prune -o -type f -printf '.' | wc -c)
        echo "Status before cleanup (prebinning): $files_pre files and $size_pre GB of storage" > {output}

        if [ -d "{loc}/{sample}/tmp/assembly/00-assembly" ]; then rm -r {loc}/{sample}/tmp/assembly/00-assembly; fi
        if [ -d "{loc}/{sample}/tmp/assembly/10-consensus" ]; then rm -r {loc}/{sample}/tmp/assembly/10-consensus; fi
        if [ -d "{loc}/{sample}/tmp/assembly/20-repeat" ]; then rm -r {loc}/{sample}/tmp/assembly/20-repeat; fi
        if [ -d "{loc}/{sample}/tmp/assembly/30-contigger" ]; then rm -r {loc}/{sample}/tmp/assembly/30-contigger; fi
        if [ -d "{loc}/{sample}/tmp/assembly/40-polishing" ]; then rm -r {loc}/{sample}/tmp/assembly/40-polishing; fi
        if [ -f "{loc}/{sample}/tmp/assembly/assembly.fasta.map-ont.mmi" ]; then rm {loc}/{sample}/tmp/assembly/assembly.fasta.map-ont.mmi; fi
        if [ -f "{loc}/{sample}/tmp/assembly/assembly.fasta" ]; then pigz {loc}/{sample}/tmp/assembly/assembly.fasta; fi
        if [ -f "{loc}/{sample}/tmp/assembly/assembly_graph.gfa" ]; then pigz {loc}/{sample}/tmp/assembly/assembly_graph.gfa; fi

        if [ -d "{loc}/{sample}/tmp/assembly/0-cleaning_and_unitigs" ]; then rm -r {loc}/{sample}/tmp/assembly/0-cleaning_and_unitigs; fi
        if [ -d "{loc}/{sample}/tmp/assembly/1-light_resolve" ]; then rm -r {loc}/{sample}/tmp/assembly/1-light_resolve; fi
        if [ -d "{loc}/{sample}/tmp/assembly/2-heavy_path_resolve" ]; then rm -r {loc}/{sample}/tmp/assembly/2-heavy_path_resolve; fi
        if [ -d "{loc}/{sample}/tmp/assembly/3-mapping" ]; then rm -r {loc}/{sample}/tmp/assembly/3-mapping; fi
        if [ -d "{loc}/{sample}/tmp/assembly/binary_temp" ]; then rm -r {loc}/{sample}/tmp/assembly/binary_temp; fi
        if [ -d "{loc}/{sample}/tmp/assembly/alternate_assemblies" ]; then rm -r {loc}/{sample}/tmp/assembly/alternate_assemblies; fi
        if [ -d "{loc}/{sample}/tmp/assembly/assembly_graphs" ]; then rm -r {loc}/{sample}/tmp/assembly/assembly_graphs; fi
        if [ -f "{loc}/{sample}/tmp/assembly/assembly.fa.map-ont.mmi" ]; then rm {loc}/{sample}/tmp/assembly/assembly.fa.map-ont.mmi; fi
        if [ -f "{loc}/{sample}/tmp/assembly/assembly.fa" ]; then pigz {loc}/{sample}/tmp/assembly/assembly.fa; fi
        if [ -f "{loc}/{sample}/tmp/assembly/assembly_primary.fa" ]; then pigz {loc}/{sample}/tmp/assembly/assembly_primary.fa; fi

        if [ -d "{loc}/{sample}/tmp/assembly/tmp" ]; then rm -r {loc}/{sample}/tmp/assembly/tmp; fi
        if [ -f "{loc}/{sample}/tmp/assembly/contigs.fasta.map-ont.mmi" ]; then rm {loc}/{sample}/tmp/assembly/contigs.fasta.map-ont.mmi; fi
        if [ -f "{loc}/{sample}/tmp/assembly/contigs.fasta" ]; then pigz {loc}/{sample}/tmp/assembly/contigs.fasta; fi

        if [ -f "{loc}/{sample}/tmp/assembly/assembly_custom.fa.map-ont.mmi" ]; then rm {loc}/{sample}/tmp/assembly/assembly_custom.fa.map-ont.mmi; fi
        if [ -f "{loc}/{sample}/tmp/assembly/assembly_custom.fa" ]; then pigz {loc}/{sample}/tmp/assembly/assembly_custom.fa; fi

        if [ -f "{loc}/{sample}/tmp/polishing/ids.txt" ]; then rm {loc}/{sample}/tmp/polishing/*.txt; fi
        if [ -f "{loc}/{sample}/tmp/polishing/lin_contigs_id_01.hdf" ]; then rm {loc}/{sample}/tmp/polishing/lin_contigs_id_*.hdf; fi
        if [ -f "{loc}/{sample}/tmp/polishing/calls_to_draft.bam" ]; then rm {loc}/{sample}/tmp/polishing/calls_to_draft.bam*; fi
        if [ -f "{loc}/{sample}/tmp/polishing/asm_pol.fasta" ]; then pigz {loc}/{sample}/tmp/polishing/asm_pol.fasta; fi

        if [ -f "{loc}/{sample}/tmp/curation/asm_curated.fasta" ]; then pigz {loc}/{sample}/tmp/curation/asm_curated.fasta; fi
        if [ -f "{loc}/{sample}/tmp/curation/mapping.bam" ]; then rm {loc}/{sample}/tmp/curation/mapping.bam*; fi

        if [ -f "{loc}/{sample}/tmp/filtering/asm_filt_len.fasta" ]; then pigz {loc}/{sample}/tmp/filtering/asm_filt_len.fasta; fi
        if [ -f "{loc}/{sample}/tmp/filtering/asm_filt_prok.fasta" ]; then pigz {loc}/{sample}/tmp/filtering/asm_filt_prok.fasta; fi
        if [ -f "{loc}/{sample}/tmp/filtering/asm_filt_euk.fasta" ]; then pigz {loc}/{sample}/tmp/filtering/asm_filt_euk.fasta; fi
        if [ -d "{loc}/{sample}/tmp/filtering/whokaryote" ]; then rm -r {loc}/{sample}/tmp/filtering/whokaryote; fi
 
        size_post=$(du -sb --exclude='{loc}/{sample}/tmp/binning' {loc}/{sample}/tmp | awk '{{printf "%.1f", $1/1024/1024/1024}}')
        files_post=$(find {loc}/{sample}/tmp -path '{loc}/{sample}/tmp/binning' -prune -o -type f -printf '.' | wc -c)
        echo "Status after cleanup (prebinning): $files_post files and $size_post GB of storage" >> {output}

        else echo "Cleanup skipped for prebinning" > {output}; fi
        """

rule Cleanup_postbin:
    conda: config["env_1"]
    input:
        os.path.join(loc, sample, "tmp/binning/quast.tsv"),
    output:
        os.path.join(loc, sample, "tmp/logs/cleanup_postbin.txt")
    params:
        cleanup=config["cleanup_status"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_cleanup_postbin.tsv")
    resources: usage=50
    threads: proc / 2
    shell:
        """
        if [ "{params.cleanup}" == "TRUE" ]; then

        size_pre=$(du -sb {loc}/{sample}/tmp/binning | awk '{{printf "%.1f", $1/1024/1024/1024}}')
        files_pre=$(find {loc}/{sample}/tmp/binning -type f -printf '.' | wc -c)
        echo "Status before cleanup (postbinning): $files_pre files and $size_pre GB of storage" > {output}

        if [ -f "{loc}/{sample}/tmp/binning/round_2/contigs.fasta" ]; then pigz {loc}/{sample}/tmp/binning/round_2/contigs.fasta; fi
        if [ -f "{loc}/{sample}/tmp/binning/round_3/contigs.fasta" ]; then pigz {loc}/{sample}/tmp/binning/round_3/contigs.fasta; fi
        if [ -f "{loc}/{sample}/tmp/binning/round_4/contigs.fasta" ]; then pigz {loc}/{sample}/tmp/binning/round_4/contigs.fasta; fi
        if [ -f "{loc}/{sample}/tmp/binning/round_euk/contigs.fasta" ]; then pigz {loc}/{sample}/tmp/binning/round_euk/contigs.fasta; fi

        if [ -d "{loc}/{sample}/tmp/binning/singl/innit" ]; then rm -r {loc}/{sample}/tmp/binning/singl/innit; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_1/bins_innit" ]; then rm -r {loc}/{sample}/tmp/binning/round_1/bins_innit; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_2/bins_innit" ]; then rm -r {loc}/{sample}/tmp/binning/round_2/bins_innit; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_3/bins_innit" ]; then rm -r {loc}/{sample}/tmp/binning/round_3/bins_innit; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_4/bins_innit" ]; then rm -r {loc}/{sample}/tmp/binning/round_4/bins_innit; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_euk/bins_innit" ]; then rm -r {loc}/{sample}/tmp/binning/round_euk/bins_innit; fi
        if [ -d "{loc}/{sample}/tmp/binning/prok_innit" ]; then rm -r {loc}/{sample}/tmp/binning/prok_innit; fi

        if [ -d "{loc}/{sample}/tmp/binning/singl/checkm2" ]; then rm -r {loc}/{sample}/tmp/binning/singl/checkm2; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_1/checkm2" ]; then rm -r {loc}/{sample}/tmp/binning/round_1/checkm2; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_2/checkm2" ]; then rm -r {loc}/{sample}/tmp/binning/round_2/checkm2; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_4/checkm2" ]; then rm -r {loc}/{sample}/tmp/binning/round_4/checkm2; fi
        if [ -d "{loc}/{sample}/tmp/binning/checkm1" ]; then rm -r {loc}/{sample}/tmp/binning/checkm1; fi

        if [ -d "{loc}/{sample}/tmp/binning/round_euk/eukcc/refine_workdir/mergers" ]; then rm -r {loc}/{sample}/tmp/binning/round_euk/eukcc/refine_workdir/mergers; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_euk/eukcc/refine_workdir/metaeuk" ]; then rm -r {loc}/{sample}/tmp/binning/round_euk/eukcc/refine_workdir/metaeuk; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_euk/eukcc/refine_workdir/refine" ]; then rm -r {loc}/{sample}/tmp/binning/round_euk/eukcc/refine_workdir/refine; fi
        if [ -f "{loc}/{sample}/tmp/binning/round_euk/eukcc/refine_workdir/contigs.fasta" ]; then rm {loc}/{sample}/tmp/binning/round_euk/eukcc/refine_workdir/contigs.fasta; fi
        if ls {loc}/{sample}/tmp/binning/round_euk/eukcc/refine_workdir/faa/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_euk/eukcc/refine_workdir/faa/*.fa; fi

        if [ -d "{loc}/{sample}/tmp/binning/round_1/binette/temporary_files" ]; then rm -r {loc}/{sample}/tmp/binning/round_1/binette/temporary_files; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_3/binette/temporary_files" ]; then rm -r {loc}/{sample}/tmp/binning/round_3/binette/temporary_files; fi
        if ls {loc}/{sample}/tmp/binning/round_3/binette/final_bins/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_3/binette/final_bins/*.fa; fi

        if [ -d "{loc}/{sample}/tmp/binning/round_2/mapping" ]; then rm -r {loc}/{sample}/tmp/binning/round_2/mapping; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_2/comebin/data_augmentation" ]; then rm -r {loc}/{sample}/tmp/binning/round_2/comebin/data_augmentation; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_2/comebin/comebin_res/cluster_res" ]; then rm -r {loc}/{sample}/tmp/binning/round_2/comebin/comebin_res/cluster_res; fi
        if ls {loc}/{sample}/tmp/binning/round_2/comebin/comebin_res/*.pth.tar >/dev/null 2>&1; then rm {loc}/{sample}/tmp/binning/round_2/comebin/comebin_res/*.pth.tar; fi
        if [ -f "{loc}/{sample}/tmp/binning/round_2/contigs.fasta.frag.ffn" ]; then rm {loc}/{sample}/tmp/binning/round_2/contigs.fasta.frag.*; fi
        if [ -f "{loc}/{sample}/tmp/binning/round_2/contigs.fasta.bacar_marker.2quarter_lencutoff_1001.seed" ]; then rm {loc}/{sample}/tmp/binning/round_2/contigs.fasta.bacar_marker.2quarter_lencutoff_1001.seed*; fi
        if ls {loc}/{sample}/tmp/binning/round_2/comebin/comebin_res/comebin_res_bins/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_2/comebin/comebin_res/comebin_res_bins/*.fa; fi
        if ls {loc}/{sample}/tmp/binning/round_2/lorbin/output_bins/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_2/lorbin/output_bins/*.fa; fi

        if ls {loc}/{sample}/tmp/binning/round_1/metabat2/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_1/metabat2/*.fa; fi
        if ls {loc}/{sample}/tmp/binning/round_3/metabat2/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_3/metabat2/*.fa; fi
        if ls {loc}/{sample}/tmp/binning/round_euk/metabat2/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_euk/metabat2/*.fa; fi

        if ls {loc}/{sample}/tmp/binning/round_euk/remag/bins/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_euk/remag/bins/*.fa; fi
        if ls {loc}/{sample}/tmp/binning/round_euk/quickbin/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_euk/quickbin/*.fa; fi
        if ls {loc}/{sample}/tmp/binning/round_4/quickbin/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_4/quickbin/*.fa; fi

        if ls {loc}/{sample}/tmp/binning/round_1/vamb/bins/*.fna >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_1/vamb/bins/*.fna; fi
        if ls {loc}/{sample}/tmp/binning/round_3/vamb/bins/*.fna >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_3/vamb/bins/*.fna; fi
        if ls {loc}/{sample}/tmp/binning/round_euk/vamb/bins/*.fna >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_euk/vamb/bins/*.fna; fi

        if ls {loc}/{sample}/tmp/binning/bins_euk/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/bins_euk/*.fa; fi
        if ls {loc}/{sample}/tmp/binning/bins_prok/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/bins_prok/*.fa; fi
        if ls {loc}/{sample}/tmp/binning/round_2/bins/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_2/bins/*.fa; fi
        if ls {loc}/{sample}/tmp/binning/round_3/bins/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_3/bins/*.fa; fi
        if ls {loc}/{sample}/tmp/binning/round_4/bins/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_4/bins/*.fa; fi

        if ls {loc}/{sample}/tmp/binning/round_3/semibin/output_bins/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_3/semibin/output_bins/*.fa; fi
        if ls {loc}/{sample}/tmp/binning/round_euk/semibin/output_bins/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_euk/semibin/output_bins/*.fa; fi
        if [ "{binmode}" != "fast" ]; then if ls {loc}/{sample}/tmp/binning/round_1/semibin/output_bins/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_1/semibin/output_bins/*.fa; fi; fi

        if [ "{binmode}" == "extended" ]; then
        if [ -f "{loc}/{sample}/tmp/binning/round_1/contigs.fasta" ]; then pigz {loc}/{sample}/tmp/binning/round_1/contigs.fasta; fi
        if ls {loc}/{sample}/tmp/binning/singl/bins/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/singl/bins/*.fa; fi
        if ls {loc}/{sample}/tmp/binning/round_1/binette/final_bins/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_1/binette/final_bins/*.fa; fi
        if ls {loc}/{sample}/tmp/binning/round_1/bins/*.fa >/dev/null 2>&1; then pigz {loc}/{sample}/tmp/binning/round_1/bins/*.fa; fi; fi

        size_post=$(du -sb {loc}/{sample}/tmp/binning | awk '{{printf "%.1f", $1/1024/1024/1024}}')
        files_post=$(find {loc}/{sample}/tmp/binning -type f -printf '.' | wc -c)
        echo "Status after cleanup (postbinning): $files_post files and $size_post GB of storage" >> {output}

        else echo "Cleanup skipped for postbinning" > {output}; fi
        """

rule Summary_usage:
    conda: config["env_1"]
    input:
        os.path.join(loc, sample, "tmp/logs/cleanup_prebin.txt"),
        os.path.join(loc, sample, "tmp/logs/cleanup_postbin.txt"),
    output:
        os.path.join(loc, sample, "tmp/logs/summary_mmlong2-lite.tsv")
    params:
        cleanup=config["cleanup_status"],
    threads: proc
    shell:
        """
        head -n1 {loc}/{sample}/tmp/logs/usage_summary_stats.tsv | sed 's/^/stage\tstep\t/' > {output}.tmp
        for log in {loc}/{sample}/tmp/logs/usage_*.tsv; do
            name=$(basename $log)
            IFS=_ read _ stage step <<< "${{name%.tsv}}"
            sed 1d $log | sed "s/^/${{stage}}\t${{step}}\t/" >> {output}.tmp
        done

        awk 'BEGIN {{
            order["stage"]=1; order["assembly"]=2; order["polishing"]=3; order["curation"]=4;
            order["filtering"]=5; order["singletons"]=6; order["coverage"]=7; order["binning"]=8; 
            order["summary"]=9; order["cleanup"]=10;
        }}{{
            key = ($1 in order) ? order[$1] : 9999
            print key "\t" $0
        }}' {output}.tmp | sort -k1,1n | cut -f2- \
        | awk -F'\t' 'BEGIN{{OFS="\t"}} NR==1{{$(NF+1)="wf_name"}} NR>1{{$(NF+1)="{sample}"}} 1' \
        | awk -F'\t' 'BEGIN{{OFS="\t"}} NR==1{{$(NF+1)="wf_read_mode"}} NR>1{{$(NF+1)="{mode}"}} 1' \
        | awk -F'\t' 'BEGIN{{OFS="\t"}} NR==1{{$(NF+1)="wf_bin_mode"}} NR>1{{$(NF+1)="{binmode}"}} 1' \
        | awk -F'\t' 'BEGIN{{OFS="\t"}} NR==1{{$(NF+1)="wf_v"}} NR>1{{$(NF+1)="{wf_v}"}} 1' \
        | awk -F'\t' -v d="$(date +'%Y-%m-%d')" 'BEGIN{{OFS="\t"}} NR==1{{$(NF+1)="wf_date"}} NR>1{{$(NF+1)=d}} 1' > {output}

        rm {output}.tmp
        rsync {output} {loc}/{sample}/results/{sample}_usage.tsv

        if [ "{params.cleanup}" == "TRUE" ]; then
        rsync {output} {loc}/{sample}/tmp/logs/summary_mmlong2-lite_$(date +"%Y%m%d-%Hh%Mm%Ss").tsv 
        rm {loc}/{sample}/tmp/logs/usage_*.tsv; fi
        """

rule Finalise_euk:
    conda: config["env_1"]
    input:
        quast=os.path.join(loc, sample, "tmp/binning/quast.tsv"),
        abund=os.path.join(loc, sample, "tmp/binning/bin_abund.tsv"),
        cov=os.path.join(loc, sample, "tmp/binning/bin_cov.tsv"),
        eukcc=os.path.join(loc, sample, "tmp/binning/round_euk/clust/euks_drep.tsv"),
        usage=os.path.join(loc, sample, "tmp/logs/summary_mmlong2-lite.tsv"),
    output:
        df1=os.path.join(loc, sample, "tmp/binning/bins_euk_mmlong2-lite.tsv"),
        df2=os.path.join(loc, sample, "results", f"{sample}_bins_euk.tsv")
    shell:
        """
        if ! grep -q "r-base" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "r-base " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi
        if [ -f "{loc}/{sample}/results/dependencies.csv" ]; then rm {loc}/{sample}/results/dependencies.csv; fi
        rsync {loc}/{sample}/tmp/dep_mmlong2-lite.csv {loc}/{sample}/results/dependencies.csv
        if [ $(cat {input.eukcc} | wc -l) -lt 1 ]; then touch {output.df1} {output.df2} && exit 0; fi
		
        R --no-echo --silent --args << 'make_df'
        library(dplyr)
        library(stringr)
        options(scipen=999)

        quast <- read.delim("{input.quast}", sep="\t", header=T)
        colnames(quast) <- c("bin","contigs","longest_contig","genome_size","gc","contig_n50","contig_n90","aun","n_per_100kb")
        abund <- read.delim("{input.abund}", sep="\t", header=T)
        colnames(abund) <- c("bin","r_abund")
        cov <- read.delim("{input.cov}", sep="\t", header=T)
        colnames(cov) <- c("bin","cov")

        eukcc <- read.delim("{input.eukcc}", sep="\t", header=T)
        tax_1 <- c("2759"="Eukaryota")
        tax_2 <- c("33154"="Opisthokonta", "33090"="Viridiplantae","2698737"="Sar", "554915"="Amoebozoa", "2611352"="Discoba", "2611341"="Metamonada")
        tax_3 <- c("4751"="Fungi", "3041"="Chlorophyta", "33634"="Stramenopiles", "46625"="Stramenopiles", "543769"="Rhizaria", "2605435"="Evosea", "207245"="Fornicata")

        eukcc$lng <- sapply(strsplit(eukcc$ncbi_lng, "-"), `[`, 5)
        eukcc$name <- tax_3[as.character(eukcc$lng)]
        eukcc$lng <- sapply(strsplit(eukcc$ncbi_lng, "-"), `[`, 4)
        eukcc$name <- ifelse(is.na(eukcc$name),tax_2[as.character(eukcc$lng)],eukcc$name)
        eukcc$lng <- sapply(strsplit(eukcc$ncbi_lng, "-"), `[`, 3)
        eukcc$name <- ifelse(is.na(eukcc$name),tax_1[as.character(eukcc$lng)],eukcc$name)

        eukcc$lng <- NULL
        eukcc$bin_org <- NULL
        colnames(eukcc) <- c("bin","completeness_eukcc","contamination_eukcc","lineage_eukcc","lineage_qc","lineage_name")

        bins <- merge(eukcc,merge(quast,merge(cov,abund,by="bin"),by="bin"), by="bin")
        bins <- bins %>% mutate(id=as.numeric(str_extract(bin, "[^.]+$"))) %>% arrange(id) %>% select(-id)
        bins$wf_name <- "{sample}"
        bins$wf_read_mode <- "{mode}"
        bins$wf_bin_mode <- "{binmode}"
        bins$wf_v <- "{wf_v}"
        bins$wf_date <- Sys.Date()
        write.table(bins, "{output.df1}", quote=F, row.names=FALSE, col.names=TRUE, sep="\t")
        write.table(bins, "{output.df2}", quote=F, row.names=FALSE, col.names=TRUE, sep="\t")
        """

rule Finalise_prok:
    conda: config["env_1"]
    input:
        quast=os.path.join(loc, sample, "tmp/binning/quast.tsv"),
        abund=os.path.join(loc, sample, "tmp/binning/bin_abund.tsv"),
        cov=os.path.join(loc, sample, "tmp/binning/bin_cov.tsv"),
        checkm1=os.path.join(loc, sample, "tmp/binning/checkm1.tsv"),
        checkm2=os.path.join(loc, sample, "tmp/binning/checkm2.tsv"),
        usage=os.path.join(loc, sample, "tmp/logs/summary_mmlong2-lite.tsv"),
        euks=os.path.join(loc, sample, "tmp/binning/bins_euk_mmlong2-lite.tsv") if "TRUE" in config["euk_status"] else [],
    output:
        df1=os.path.join(loc, sample, "tmp/binning/bins_mmlong2-lite.tsv"),
        df2=os.path.join(loc, sample, "results", f"{sample}_bins_prok.tsv")
    shell:
        """
        if ! grep -q "r-base" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "r-base " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi
        if [ -f "{loc}/{sample}/results/dependencies.csv" ]; then rm {loc}/{sample}/results/dependencies.csv; fi
        rsync {loc}/{sample}/tmp/dep_mmlong2-lite.csv {loc}/{sample}/results/dependencies.csv
		
        R --no-echo --silent --args << 'make_df'
        library(dplyr)
        library(stringr)
        options(scipen=999)
        
        quast <- read.delim("{input.quast}", sep="\t", header=T)
        colnames(quast) <- c("bin","contigs","longest_contig","genome_size","gc","contig_n50","contig_n90","aun","n_per_100kb")
        abund <- read.delim("{input.abund}", sep="\t", header=T)
        colnames(abund) <- c("bin","r_abund")
        cov <- read.delim("{input.cov}", sep="\t", header=T)
        colnames(cov) <- c("bin","cov")
        checkm1 <- read.delim("{input.checkm1}", sep="\t", header=T)
        checkm1 <- checkm1[, c("Bin.Id","Completeness","Contamination","Strain.heterogeneity")]
        colnames(checkm1) <- c("bin","completeness_checkm1","contamination_checkm1","strain_heterogeneity_checkm1")
        checkm2 <- read.delim("{input.checkm2}", sep="\t", header=T)
        colnames(checkm2) <- c("bin","completeness_checkm2","contamination_checkm2","model_checkm2","coding_density_checkm2")
        bins <- merge(checkm1,merge(checkm2,merge(quast,merge(cov,abund,by="bin"),by="bin"), by="bin"),by="bin")
        bins$wf_name <- "{sample}"
        bins$wf_read_mode <- "{mode}"
        bins$wf_bin_mode <- "{binmode}"
        bins$wf_v <- "{wf_v}"
        bins$wf_date <- Sys.Date()

        pat <- c("bin.c", "bin.s", "bin.1", "bin.2", "bin.3", "bin.4")
        bins <- bins %>% mutate(id=as.numeric(str_extract(bin, "[^.]+$"))) %>% arrange(id) %>% select(-id)
        bins <- bins %>% arrange(sapply(bin, \(x) which(sapply(pat, grepl, x=x, fixed=TRUE))))

        write.table(bins, "{output.df1}", quote=F,row.names=FALSE, col.names=TRUE, sep="\t")
        write.table(bins, "{output.df2}", quote=F,row.names=FALSE, col.names=TRUE, sep="\t")
        """
