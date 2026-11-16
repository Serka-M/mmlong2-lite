# DESCRIPTION: Snakemake rules for the binning sections (prokaryotic and eukaryotic bin production) of the mmlong2 pipeline
# AUTHOR: Mantas Sereika (mase@bio.aau.dk)
# LICENSE: GNU General Public License

rule Binning_prep_innit:
    conda: config["env_2"]
    input:
        cm2=os.path.join(loc, sample, "tmp/binning/singl/checkm2.tsv"),
        cov1=os.path.join(loc, sample, "tmp/binning/mapping/cov_all.tsv"),
        cov2=os.path.join(loc, sample, "tmp/binning/mapping/cov_all_sub.tsv"),
    output:
        id1=os.path.join(loc, sample, "tmp/binning/round_1/contigs_w_singl.txt"),
        id2=os.path.join(loc, sample, "tmp/binning/round_1/contigs_wo_singl.txt"),
        id3=os.path.join(loc, sample, "tmp/binning/round_1/contigs.txt"),
        contigs=os.path.join(loc, sample, "tmp/binning/round_1/contigs.fasta"),
        cov1=os.path.join(loc, sample, "tmp/binning/round_1/cov.tsv"),
        cov2=os.path.join(loc, sample, "tmp/binning/round_1/cov_sub.tsv"),
    params:
        cov=config["min_prok_cov"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_binning_prep-R1.tsv")
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_1" ]; then rm -r {loc}/{sample}/tmp/binning/round_1; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_2" ]; then rm -r {loc}/{sample}/tmp/binning/round_2; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_3" ]; then rm -r {loc}/{sample}/tmp/binning/round_3; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_4" ]; then rm -r {loc}/{sample}/tmp/binning/round_4; fi
        mkdir {loc}/{sample}/tmp/binning/round_1
        
        sed 1d {loc}/{sample}/tmp/binning/mapping/cov_all_sub.tsv | awk '$2 >= {params.cov}' | cut -f1 | sort > {output.id1}
        comm -23 <(sort {output.id1}) <(sort {loc}/{sample}/tmp/binning/singl/binned.txt) > {output.id2}
        seqkit grep -f {output.id2} {loc}/{sample}/tmp/filtering/asm_filt_prok.fasta | seqkit sort --by-name - > {output.contigs}
        grep ">" {output.contigs} | cut -c2- > {output.id3}
        
        head -n1 {input.cov1} > {output.cov1}
        awk 'BEGIN{{FS=OFS="\t"}} NR==FNR{{list[$1]; next}} $1 in list' {output.id3} {input.cov1} | sort >> {output.cov1}

        head -n1 {input.cov2} > {output.cov2}
        awk 'BEGIN{{FS=OFS="\t"}} NR==FNR{{list[$1]; next}} $1 in list' {output.id3} {input.cov2} | sort >> {output.cov2}
        """

rule Binning_prep_main:
    conda: config["env_2"]
    wildcard_constraints: round = r"\d+"
    input:
        lambda wildcards: ''.join([str(wildcards.loc), "/", str(wildcards.sample), "/tmp/binning/round_",str(int(wildcards.round)-1),"/checkm2.tsv"])
    output:
        id="{loc}/{sample}/tmp/binning/round_{round}/contigs.txt",
        contigs="{loc}/{sample}/tmp/binning/round_{round}/contigs.fasta",
        cov1="{loc}/{sample}/tmp/binning/round_{round}/cov.tsv",
        cov2="{loc}/{sample}/tmp/binning/round_{round}/cov_sub.tsv",
    params:
        round_curr="{round}",
        round_prev=lambda wildcards: str(int(wildcards.round)-1),
    benchmark:
        "{loc}/{sample}/tmp/logs/usage_binning_prep-R{round}.tsv"
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round_curr}" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round_curr}; fi
        mkdir {loc}/{sample}/tmp/binning/round_{params.round_curr}
        
        comm -23 <(sort {loc}/{sample}/tmp/binning/round_{params.round_prev}/contigs.txt) <(sort {loc}/{sample}/tmp/binning/round_{params.round_prev}/binned.txt) > {output.id}
        seqkit grep -f {output.id} {loc}/{sample}/tmp/binning/round_{params.round_prev}/contigs.fasta | seqkit sort --by-name - > {output.contigs}

        head -n1 {loc}/{sample}/tmp/binning/round_{params.round_prev}/cov.tsv > {output.cov1}
        awk 'BEGIN{{FS=OFS="\t"}} NR==FNR{{list[$1]; next}} $1 in list' {output.id} {loc}/{sample}/tmp/binning/round_{params.round_prev}/cov.tsv | sort >> {output.cov1}

        head -n1 {loc}/{sample}/tmp/binning/round_{params.round_prev}/cov_sub.tsv > {output.cov2}
        awk 'BEGIN{{FS=OFS="\t"}} NR==FNR{{list[$1]; next}} $1 in list' {output.id} {loc}/{sample}/tmp/binning/round_{params.round_prev}/cov_sub.tsv | sort >> {output.cov2}
        """

rule Binning_prep_euk:
    conda: config["env_2"]
    input:
        contigs=os.path.join(loc, sample, "tmp/filtering/asm_filt_euk.fasta"),
        id=os.path.join(loc, sample, "tmp/filtering", get_euk(config["tiara_status"],config["eukfilt_status"])),
        cov1=os.path.join(loc, sample, "tmp/binning/mapping/cov_all.tsv"),
        cov2=os.path.join(loc, sample, "tmp/binning/mapping/cov_all_sub.tsv"),
        links=os.path.join(loc, sample, "tmp/binning/contig_bin_prok.tsv") if "TRUE" in config["prok_status"] else [],
    output:
        id1=os.path.join(loc, sample, "tmp/binning/round_euk/contigs_euk_init.txt"),
        id2=os.path.join(loc, sample, "tmp/binning/round_euk/contigs_euk.txt"),
        contigs=os.path.join(loc, sample, "tmp/binning/round_euk/contigs.fasta"),
        cov1=os.path.join(loc, sample, "tmp/binning/round_euk/cov.tsv"),
        cov2=os.path.join(loc, sample, "tmp/binning/round_euk/cov_sub.tsv"),
    params:
        cov=config["min_euk_cov"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_binning_prep-euk.tsv")
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_euk" ]; then rm -r {loc}/{sample}/tmp/binning/round_euk; fi
        mkdir {loc}/{sample}/tmp/binning/round_euk
        
        if [ -f "{input.links}" ]; then
        cut -f1 {input.links} > {loc}/{sample}/tmp/binning/round_euk/contigs_prok.txt
        comm -23 <(sort {input.id}) <(sort {loc}/{sample}/tmp/binning/round_euk/contigs_prok.txt) > {output.id1}
        else cat {input.id} > {output.id1}; fi

        comm -12 <(sort {output.id1}) <(awk '$2 >= {params.cov}' {input.cov2} | cut -f1 | sort) > {output.id2}
        seqkit grep -f {output.id2} {input.contigs} | seqkit sort --by-name - > {output.contigs}

        head -n1 {input.cov1} > {output.cov1}
        awk 'BEGIN{{FS=OFS="\t"}} NR==FNR{{list[$1]; next}} $1 in list' {output.id2} {input.cov1} | sort >> {output.cov1}

        head -n1 {input.cov2} > {output.cov2}
        awk 'BEGIN{{FS=OFS="\t"}} NR==FNR{{list[$1]; next}} $1 in list' {output.id2} {input.cov2} | sort >> {output.cov2}
        """

rule Comebin_prep_innit:
    conda: config["env_2"]
    input:
        "{loc}/{sample}/tmp/binning/round_{round}/contigs.txt"
    output:
        "{loc}/{sample}/tmp/binning/round_{round}/mapping/header.sam"
    params:
        round="{round}",
        type=get_read_type(mode)
    benchmark:
        "{loc}/{sample}/tmp/logs/usage_binning_prep-comebin-innit-R{round}.tsv"
    resources: usage=100
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/mapping" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/mapping; fi
        mkdir {loc}/{sample}/tmp/binning/round_{params.round}/mapping
        awk 'FNR==NR{{contigs[$1]; next}} /^@SQ/ {{contig=$2; sub(/^SN:/, "", contig); if (contig in contigs) print}}' "{input}" <(samtools view --threads {threads} -H "{loc}/{sample}/tmp/binning/mapping/1-{params.type}.bam") > {output}
        """

rule Comebin_prep_main:
    conda: config["env_2"]
    input:
        contigs="{loc}/{sample}/tmp/binning/round_{round}/contigs.txt",
        header="{loc}/{sample}/tmp/binning/round_{round}/mapping/header.sam",
        bam="{loc}/{sample}/tmp/binning/mapping/{reads}.bam"
    output:
        "{loc}/{sample}/tmp/binning/round_{round}/mapping/{reads}.bam"
    params:
        round="{round}",
    benchmark:
        "{loc}/{sample}/tmp/logs/usage_binning_prep-comebin-{reads}-R{round}.tsv"
    resources: usage=config["mapping_usage"]
    threads: proc * (config["mapping_usage"]/100)
    shell:
        """
        grep -wFf "{input.contigs}" <(samtools view --threads {threads} "{input.bam}") | cat {input.header} - | samtools view --threads {threads} -b -o {output} -
        """

rule Binning_comebin:
    conda: config["env_6"]
    input:
        contigs="{loc}/{sample}/tmp/binning/round_{round}/contigs.fasta",
        bams=expand("{{loc}}/{{sample}}/tmp/binning/round_{{round}}/mapping/{reads}.bam", reads=get_reads(reads_diffcov,fastq,mode))
    output:
        "{loc}/{sample}/tmp/binning/round_{round}/comebin/comebin_res/comebin_res.tsv"
    params:
        round="{round}",
        seed=config["seed"],
        len=config["min_contig_len"],
        epochs=config["comebin_epochs"],
    benchmark:
        "{loc}/{sample}/tmp/logs/usage_binning_comebin-R{round}.tsv"
    resources: usage=100
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/comebin" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/comebin; fi
        n_contigs=$(grep ">" {input.contigs} | wc -l)
        if [ $n_contigs -lt 1024 ]; then batch_size=$n_contigs; else batch_size=1024;fi
        python() {{ [[ "$*" == *train* ]] && set -- "$@" --epochs {params.epochs}; command python "$@"; }}
        export -f python
        run_comebin.sh -a {input.contigs} -o {loc}/{sample}/tmp/binning/round_{params.round}/comebin -t {threads} -p {loc}/{sample}/tmp/binning/round_{params.round}/mapping -b $batch_size -d cpu -s {params.seed}
        unset -f python
        if ! grep -q "comebin" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "comebin " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi
        """

rule Binning_metabat2:
    conda: config["env_2"]
    input:
        "{loc}/{sample}/tmp/binning/round_{round}/contigs.fasta"
    output:
        "{loc}/{sample}/tmp/binning/round_{round}/metabat2/bins_metabat2.MemberMatrix.txt"
    params:
        round="{round}",
        seed=config["seed"],
        len=config["min_contig_len"],
        prok_len=config["min_prok_len"],
        euk_len=config["min_euk_len"],
    benchmark:
        "{loc}/{sample}/tmp/logs/usage_binning_metabat2-R{round}.tsv"
    resources: usage=100
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/metabat2" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/metabat2; fi
        mkdir {loc}/{sample}/tmp/binning/round_{params.round}/metabat2
        if [ "{params.round}" == "euk" ]; then mag_len={params.euk_len}; else mag_len={params.prok_len}; fi
        metabat2 -i {input} -a {loc}/{sample}/tmp/binning/round_{params.round}/cov.tsv -o {loc}/{sample}/tmp/binning/round_{params.round}/metabat2/bins_metabat2 -t {threads} -m {params.len} -s $mag_len --saveCls --seed {params.seed}
        if ls {loc}/{sample}/tmp/binning/round_{params.round}/metabat2/*.fa >/dev/null 2>&1; then sed -i '/^>/ s/[[:space:]].*$//' {loc}/{sample}/tmp/binning/round_{params.round}/metabat2/*.fa; fi
        if ! grep -q "metabat2" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "metabat2 " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi	
        """

rule Binning_vamb:
    conda: config["env_7"]
    input:
        "{loc}/{sample}/tmp/binning/round_{round}/contigs.fasta"
    output:
        "{loc}/{sample}/tmp/binning/round_{round}/vamb/vae_clusters_unsplit.tsv"
    params:
        round="{round}",
        seed=config["seed"],
        len=config["min_contig_len"],
        prok_len=config["min_prok_len"],
        euk_len=config["min_euk_len"],
    benchmark:
        "{loc}/{sample}/tmp/logs/usage_binning_vamb-R{round}.tsv"
    resources: usage=100
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/vamb" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/vamb; fi
        if [ "{params.round}" == "euk" ]; then mag_len={params.euk_len}; else mag_len={params.prok_len}; fi
        vamb bin default --fasta {input} --outdir {loc}/{sample}/tmp/binning/round_{params.round}/vamb -m {params.len} --minfasta $mag_len -p {threads} --abundance_tsv {loc}/{sample}/tmp/binning/round_{params.round}/cov_sub.tsv -o "" --seed {params.seed}
        if ! grep -q "vamb" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "vamb " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi	
        """

rule Binning_semibin2:
    conda: config["env_5"]
    input:
        "{loc}/{sample}/tmp/binning/round_{round}/contigs.fasta"
    output:
        "{loc}/{sample}/tmp/binning/round_{round}/semibin/bins_info.tsv"
    params:
        cov=lambda wildcards: get_semibin_cov(wildcards.loc,wildcards.sample,reads_diffcov,wildcards.round,semibin_mod,config["semibin_train_r1"]),
        round="{round}",
        seed=config["seed"],
        len=config["min_contig_len"],
        prok_len=config["min_prok_len"],
        euk_len=config["min_euk_len"],
        type_prok=config["semibin_type_prok"],
        type_euk=config["semibin_type_euk"],
    benchmark:
        "{loc}/{sample}/tmp/logs/usage_binning_semibin2-R{round}.tsv"
    resources: usage=100
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/semibin" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/semibin; fi
        if [ "{params.round}" == "euk" ]; then type="{params.type_euk} --no-recluster"; else type="{params.type_prok}"; fi
        if [ "{params.round}" == "euk" ]; then mag_len={params.euk_len}; else mag_len={params.prok_len}; fi
        SemiBin2 single_easy_bin -i {input} -b {params.cov} -o {loc}/{sample}/tmp/binning/round_{params.round}/semibin -p {threads} -m {params.len} --minfasta-kbs $((${{mag_len}}/1000)) --self-supervised --sequencing-type $type --engine cpu --orf-finder {semibin_prot} --compression none --random-seed {params.seed}
        if [ "{params.round}" == "euk" ]; then touch {output}; fi
        if ! grep -q "semibin" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "semibin " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi  
        """

rule Binning_quickbin:
    conda: config["env_1"]
    input:
        "{loc}/{sample}/tmp/binning/round_{round}/contigs.fasta"
    output:
        "{loc}/{sample}/tmp/binning/round_{round}/quickbin/quickbin.txt"
    params:
        round="{round}",
        seed=config["seed"],
        len=config["min_contig_len"],
        prok_len=config["min_prok_len"],
        euk_len=config["min_euk_len"],
        minseed=config["qb_minseed"],
        minpentamersize=config["qb_minpentamersize"],
        gcwidth=config["qb_gcwidth"],
        depthwidth=config["qb_depthwidth"],
    benchmark:
        "{loc}/{sample}/tmp/logs/usage_binning_quickbin-R{round}.tsv"
    resources: usage=100
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/quickbin" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/quickbin; fi
        if [ "{params.round}" == "euk" ]; then mag_len={params.euk_len}; else mag_len={params.prok_len}; fi
        quickbin.sh -da in={input} mincontig={params.len} mincluster=$mag_len minseed={params.minseed} minpentamersize={params.minpentamersize} gcwidth={params.gcwidth} depthwidth={params.depthwidth} out={loc}/{sample}/tmp/binning/round_{params.round}/quickbin/QuickBin_%.fa {loc}/{sample}/tmp/binning/mapping/*.bam
        if [ ! -d "{loc}/{sample}/tmp/binning/round_{params.round}/quickbin" ]; then mkdir {loc}/{sample}/tmp/binning/round_{params.round}/quickbin; fi
        touch {output}
        if ! grep -q "quickbin" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then quickbin.sh -h | grep "Last modified" | cut -d' ' -f3- | sed 's/,//g' | sed 's/ /_/g' | sed 's/^/quickbin,/' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi	
        """

rule Binning_lorbin:
    conda: config["env_10"]
    input:
        contigs="{loc}/{sample}/tmp/binning/round_{round}/contigs.fasta",
    output:
        "{loc}/{sample}/tmp/binning/round_{round}/lorbin/label.csv"
    params:
        round="{round}",
        len=config["min_contig_len"],
        mag_len=config["min_prok_len"],
    benchmark:
        "{loc}/{sample}/tmp/logs/usage_binning_lorbin-R{round}.tsv"
    resources: usage=100
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/lorbin" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/lorbin; fi
        LorBin bin --fa {input.contigs} -b {loc}/{sample}/tmp/binning/mapping/*.bam -o {loc}/{sample}/tmp/binning/round_{params.round}/lorbin --num_process {threads}
        if ! grep -q "lorbin" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "lorbin " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi  
        """

rule Binning_remag:
    conda: config["env_7"]
    input:
        contigs="{loc}/{sample}/tmp/binning/round_{round}/contigs.fasta",
    output:
        "{loc}/{sample}/tmp/binning/round_{round}/remag/bins.csv"
    params:
        round="{round}",
        seed=config["seed"],
        len=config["min_contig_len"],
        mag_len=config["min_euk_len"],
    benchmark:
        "{loc}/{sample}/tmp/logs/usage_binning_remag-R{round}.tsv"
    resources: usage=100
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/remag" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/remag; fi
        remag -f {input.contigs} -c {loc}/{sample}/tmp/binning/mapping/*.bam -o {loc}/{sample}/tmp/binning/round_{params.round}/remag -t {threads} --min-contig-length {params.len} --min-bin-size {params.mag_len} --random-seed {params.seed} --skip-bacterial-filter --mode metagenomics
        if ! grep -q "remag" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "remag " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi  
        """

rule Binning_binette:
    conda: config["env_8"]
    input:
        asm="{loc}/{sample}/tmp/binning/round_{round}/contigs.fasta",
        vamb="{loc}/{sample}/tmp/binning/round_{round}/vamb/vae_clusters_unsplit.tsv",
        mb2="{loc}/{sample}/tmp/binning/round_{round}/metabat2/bins_metabat2.MemberMatrix.txt",
        sb2="{loc}/{sample}/tmp/binning/round_{round}/semibin/bins_info.tsv"
    output:
        "{loc}/{sample}/tmp/binning/round_{round}/binette/final_bins_quality_reports.tsv",
    params:
        round="{round}",
        db=config["db_checkm2"],
        weight=config["binette_weight"],
        seed=config["seed"],
        compl=lambda wildcards: get_compl(wildcards.round, config["binmode"], config["min_compl_1"], config["min_compl_2"]),
        cont=lambda wildcards: get_cont(wildcards.round, config["binmode"], config["min_cont_1"], config["min_cont_2"]),
        min_mag_len=config["min_prok_len"],
        max_mag_len=config["max_prok_len"],
    benchmark:
        "{loc}/{sample}/tmp/logs/usage_binning_binette-R{round}.tsv"
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/binette" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/binette; fi
        export PYTHONHASHSEED={params.seed}
        dir={loc}/{sample}/tmp/binning/round_{params.round}
        binette -d $dir/vamb/bins $dir/metabat2 $dir/semibin/output_bins -c {input.asm} -o $dir/binette -t {threads} -w {params.weight} --min_completeness {params.compl} --max_contamination {params.cont} --min_length {params.min_mag_len} --max_length {params.max_mag_len} --checkm2_db {params.db}
        if ! grep -q "binette" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "binette " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi  
        """

rule Binning_checkm2:
    conda: config["env_8"]
    input:
        lambda wildcards: get_qc_in(wildcards.loc,wildcards.sample,wildcards.round,config["binmode"],config["lorbin_status"])
    output:
        "{loc}/{sample}/tmp/binning/round_{round}/checkm2.tsv"
    params:
        db=config["db_checkm2"],
        round="{round}",
        dir=lambda wildcards: get_qc_dir(wildcards.loc,wildcards.sample,wildcards.round,config["binmode"],config["lorbin_status"]),
        compl=lambda wildcards: get_compl(wildcards.round, config["binmode"], config["min_compl_1"], config["min_compl_2"]),
        cont=lambda wildcards: get_cont(wildcards.round, config["binmode"], config["min_cont_1"], config["min_cont_2"]),
        mag_len=config["min_prok_len"],
    benchmark:
        "{loc}/{sample}/tmp/logs/usage_binning_qc-R{round}.tsv"
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/checkm2" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/checkm2; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/bins_innit" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/bins_innit; fi
        if [ -d "{loc}/{sample}/tmp/binning/round_{params.round}/bins" ]; then rm -r {loc}/{sample}/tmp/binning/round_{params.round}/bins; fi

        mkdir {loc}/{sample}/tmp/binning/round_{params.round}/bins_innit
        mkdir {loc}/{sample}/tmp/binning/round_{params.round}/bins
        export CHECKM2DB="{params.db}"
        
        if [ $(basename {input}) == "final_bins_quality_reports.tsv" ];
        then
            if ls {params.dir}/*.fa >/dev/null 2>&1; then 
                rsync {params.dir}/*.fa {loc}/{sample}/tmp/binning/round_{params.round}/bins_innit/.
                for bin in {loc}/{sample}/tmp/binning/round_{params.round}/bins_innit/*.fa; do mv "$bin" "$(dirname "$bin")/{sample}.bin.{params.round}.$(basename "$bin" | sed 's/^binette_bin//')"; done; fi
            cut -f1,5,6,8,11 {input} | sed 's/^binette_bin/{sample}.bin.{params.round}./' > {output};
        else
            n=1 && for i in {params.dir}/*.fa; do rsync $i {loc}/{sample}/tmp/binning/round_{params.round}/bins_innit/{sample}.bin.{params.round}.${{n}}.fa && n=$(($n+1)); done
            checkm2 predict -x .fa -i {loc}/{sample}/tmp/binning/round_{params.round}/bins_innit -o {loc}/{sample}/tmp/binning/round_{params.round}/checkm2 -t {threads}
            awk -F'\t' 'NR==1 || $9 >= {params.mag_len}' {loc}/{sample}/tmp/binning/round_{params.round}/checkm2/quality_report.tsv | cut -f1,2,3,4,6 > {output}
        fi
        
        awk -F "\t" '{{ if (($2 >= {params.compl}) && ($3 <= {params.cont})) {{print $1}} }}' {output} > {loc}/{sample}/tmp/binning/round_{params.round}/bins_keep.txt
        if [ $(cat {loc}/{sample}/tmp/binning/round_{params.round}/bins_keep.txt | wc -l) -ge 1 ]; then
        cat {loc}/{sample}/tmp/binning/round_{params.round}/bins_keep.txt | xargs -i --max-procs=1 bash -c 'rsync {loc}/{sample}/tmp/binning/round_{params.round}/bins_innit/{{}}.fa {loc}/{sample}/tmp/binning/round_{params.round}/bins/.'
        grep ">" {loc}/{sample}/tmp/binning/round_{params.round}/bins/*.fa | cut -f2 -d">" | sort > {loc}/{sample}/tmp/binning/round_{params.round}/binned.txt; fi
        """

rule Binning_aggregate_prok:
    conda: config["env_2"]
    input:
        expand("{loc}/{sample}/tmp/binning/round_{round}/checkm2.tsv",sample=sample,loc=loc,round=get_rounds(config["binmode"])),
    output:
        os.path.join(loc, sample, "tmp/binning/checkm2.tsv"),
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_binning_aggregate.tsv")
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/prok_innit" ]; then rm -r {loc}/{sample}/tmp/binning/prok_innit; fi
        mkdir {loc}/{sample}/tmp/binning/prok_innit
        rsync {loc}/{sample}/tmp/binning/*/bins/*.fa {loc}/{sample}/tmp/binning/prok_innit/.
        sed -i '/^>/ s/[[:space:]].*$//' {loc}/{sample}/tmp/binning/prok_innit/*.fa
        
        head -n1 {loc}/{sample}/tmp/binning/singl/checkm2.tsv > {output}
        for stage in "singl" "round_1" "round_2" "round_3" "round_4" 
        do if [ -d "{loc}/{sample}/tmp/binning/$stage" ]; then if [ $(cat {loc}/{sample}/tmp/binning/$stage/bins_keep.txt | wc -l) -ge 1 ]; then grep -w -f {loc}/{sample}/tmp/binning/$stage/bins_keep.txt {loc}/{sample}/tmp/binning/$stage/checkm2.tsv >> {output}; fi; fi; done
        """

rule Binning_checkm1:
    conda: config["env_6"]
    input:
         os.path.join(loc, sample, "tmp/binning/checkm2.tsv"),
    output:
        checkm1=os.path.join(loc, sample, "tmp/binning/checkm1.tsv"),
        bins=os.path.join(loc, sample, "tmp/binning/bins_keep_prok.txt"),
        links=os.path.join(loc, sample, "tmp/binning/contig_bin_prok.tsv"),
    params:
        compl=config["min_compl_checkm1"],
        cont=config["min_cont_checkm1"],
        score=config["min_score_checkm1"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_binning_qc2.tsv")
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/checkm1" ]; then rm -r {loc}/{sample}/tmp/binning/checkm1; fi
        checkm lineage_wf -x fa -t {threads} -f {output.checkm1} --pplacer_threads $(({threads} / 2)) --reduced_tree --tab_table {loc}/{sample}/tmp/binning/prok_innit {loc}/{sample}/tmp/binning/checkm1
        if ! grep -q "checkm-" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "checkm-genome " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi
        
        if [ -d "{loc}/{sample}/tmp/binning/bins_prok" ]; then rm -r {loc}/{sample}/tmp/binning/bins_prok; fi
        mkdir {loc}/{sample}/tmp/binning/bins_prok
        awk -F "\t" '{{ if (($12 >= {params.compl}) && ($13 <= {params.cont}) && ($12 - $13 * 5 > {params.score})) {{print $1}} }}' {output.checkm1} > {output.bins}
        cat {output.bins} | xargs -i --max-procs=1 bash -c 'rsync {loc}/{sample}/tmp/binning/prok_innit/{{}}.fa {loc}/{sample}/tmp/binning/bins_prok/.'

        if ls {loc}/{sample}/tmp/binning/bins_prok/*.fa >/dev/null 2>&1; then
        grep -H ">" {loc}/{sample}/tmp/binning/bins_prok/*.fa | xargs -n1 basename | sed "s/.fa:>/\t/" | awk -F'\t' '{{print $2 "\t" $1}}' > {output.links}
        else touch {output.links}; fi
        
        if [ -d "{loc}/{sample}/results/bins_prok" ]; then rm -r {loc}/{sample}/results/bins_prok; fi
        rsync -r {loc}/{sample}/tmp/binning/bins_prok {loc}/{sample}/results/.
        """

rule Binning_eukcc_prep:
    conda: config["env_3"]
    input:
         os.path.join(loc, sample, "tmp/binning/round_euk/contigs.fasta"),
         os.path.join(loc, sample, "tmp/binning/round_euk/metabat2/bins_metabat2.MemberMatrix.txt") if "FALSE" in config["skip_metabat2_euk"] else [],
         os.path.join(loc, sample, "tmp/binning/round_euk/vamb/vae_clusters_unsplit.tsv") if "FALSE" in config["skip_vamb_euk"] else [],
         os.path.join(loc, sample, "tmp/binning/round_euk/quickbin/quickbin.txt") if "FALSE" in config["skip_quickbin_euk"] else [],
         os.path.join(loc, sample, "tmp/binning/round_euk/semibin/bins_info.tsv") if "FALSE" in config["skip_semibin2_euk"] else [],
         os.path.join(loc, sample, "tmp/binning/round_euk/remag/bins.csv") if "FALSE" in config["skip_remag_euk"] else [],
    output:
        os.path.join(loc, sample, "tmp/binning/round_euk/bins_innit.txt"),
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_binning_eukcc_prep.tsv")
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_euk/bins_innit" ]; then rm -r {loc}/{sample}/tmp/binning/round_euk/bins_innit; fi
        mkdir {loc}/{sample}/tmp/binning/round_euk/bins_innit

        innit="{loc}/{sample}/tmp/binning/round_euk/bins_innit"
        vamb="{loc}/{sample}/tmp/binning/round_euk/vamb/bins"
        metabat="{loc}/{sample}/tmp/binning/round_euk/metabat2"
        quickbin="{loc}/{sample}/tmp/binning/round_euk/quickbin"
        semibin="{loc}/{sample}/tmp/binning/round_euk/semibin/output_bins"
        remag="{loc}/{sample}/tmp/binning/round_euk/remag/bins"

        if ls $vamb/*.fna >/dev/null 2>&1; then rsync $vamb/*.fna $innit/.; fi
        if ls $innit/*.fna >/dev/null 2>&1; then for f in $innit/*.fna; do mv "$f" "$innit/vamb_$(basename "${{f%.fna}}.fa")"; done; fi
        if ls $metabat/*.fa >/dev/null 2>&1; then rsync $metabat/*.fa $innit/.; fi
        if ls $quickbin/*.fa >/dev/null 2>&1; then rsync $quickbin/*.fa $innit/.; fi
        if ls $semibin/*.fa >/dev/null 2>&1; then rsync $semibin/*.fa $innit/.; fi
        if ls $remag/*.fa >/dev/null 2>&1; then rsync $remag/*.fa $innit/.; fi
        if ls $innit/*.fa >/dev/null 2>&1; then sed -i '/^>/ s/[[:space:]].*$//' $innit/*.fa; fi

        find $innit -maxdepth 1 -type f -name '*.fa' -printf '%f\n' > {output}
        """

rule Binning_eukcc:
    conda: config["env_3"]
    input:
        os.path.join(loc, sample, "tmp/binning/round_euk/bins_innit.txt")
    output:
        os.path.join(loc, sample, "tmp/binning/round_euk/eukcc/eukcc.tsv"),
    params:
        db=config["db_eukcc"],
        epa=config["eukcc_epa"],
        merge_max=config["eukcc_merge_max"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_binning_eukcc.tsv")
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_euk/eukcc" ]; then rm -r {loc}/{sample}/tmp/binning/round_euk/eukcc; fi
        bins=$(wc -l < {input})

        if [ "$bins" -gt "{params.merge_max}" ]; then
        merge="--n_combine 0" && echo "Running EukCC with bin merging disabled"
        else
        merge="" && echo "Running EukCC with bin merging enabled"
        fi

        if ls {loc}/{sample}/tmp/binning/round_euk/bins_innit/*.fa >/dev/null 2>&1; then 
            export EUKCC2_DB="{params.db}"
            eukcc folder --out {loc}/{sample}/tmp/binning/round_euk/eukcc --threads {threads} --threads_epa {params.epa} {loc}/{sample}/tmp/binning/round_euk/bins_innit $merge
            awk -F'\t' 'BEGIN{{OFS="\t"}} NR==1{{$(NF+1)="eukcc_qc"}} NR>1{{$(NF+1)="PASS"}} 1' {loc}/{sample}/tmp/binning/round_euk/eukcc/eukcc.csv > {output}
            sed 's/$/\tFAIL/' {loc}/{sample}/tmp/binning/round_euk/eukcc/bad_quality.csv | sed 1d >> {output};
        else mkdir {loc}/{sample}/tmp/binning/round_euk/eukcc && touch {output}; fi
        if ! grep -q "eukcc" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "eukcc " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi
        """

rule Binning_eukclust:
    conda: config["env_1"]
    input:
        os.path.join(loc, sample, "tmp/binning/round_euk/eukcc/eukcc.tsv"),
    output:
        id=os.path.join(loc, sample, "tmp/binning/round_euk/clust/euks.txt"),
        df1=os.path.join(loc, sample, "tmp/binning/round_euk/clust/euks.tsv"),
        df2=os.path.join(loc, sample, "tmp/binning/round_euk/clust/euks.csv"),
        links=os.path.join(loc, sample, "tmp/binning/round_euk/clust/contig_bin.tsv"),
        clust=os.path.join(loc, sample, "tmp/binning/round_euk/clust/clust.tsv"),
    params:
        compl=config["min_compl_euk"],
        cont=config["min_cont_euk"],
        aln_frac=config["aln_frac"],
        ani_pre=config["ani_pre"],
        ani_clust=config["ani_clust"],
        preclust_method=config["preclust_method"],
        clust_method=config["clust_method"],
        qual_formula=config["qual_formula"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_binning_eukclust.tsv")
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/round_euk/clust" ]; then rm -r {loc}/{sample}/tmp/binning/round_euk/clust; fi
        mkdir {loc}/{sample}/tmp/binning/round_euk/clust
        awk -F "\t" '{{ if (($2 >= {params.compl}) && ($3 <= {params.cont})) {{print "{loc}/{sample}/tmp/binning/round_euk/bins_innit/" $1}} }}' {input} > {output.id}
        awk -F "\t" '{{ if (NR==1 || ($2 >= {params.compl}) && ($3 <= {params.cont})) {{print $0}} }}' {input} | sed 's/\\.fa//g' > {output.df1}
        cut -f1,2,3 {output.df1} | tr '\t' ',' | sed '1s/^bin,completeness/genome,completeness/' > {output.df2}

        if ls {loc}/{sample}/tmp/binning/round_euk/eukcc/merged_bins/*.fa >/dev/null 2>&1; then mv {loc}/{sample}/tmp/binning/round_euk/eukcc/merged_bins/*.fa {loc}/{sample}/tmp/binning/round_euk/bins_innit/.; fi
        xargs -a {output.id} -I{{}} awk '/^>/{{bin=FILENAME; sub(/^.*\//,"",bin); sub(/\.fa$/,"",bin); split($1,a," "); print substr(a[1],2) "\t" bin}}' {{}} > {output.links}

        if [ $(cat {output.id} | wc -l) -lt 1 ]; then touch {output.clust}
        else coverm cluster --genome-fasta-list {output.id} --genome-info {output.df2} --output-cluster-definition {output.clust} --precluster-ani {params.ani_pre} --ani {params.ani_clust} --min-aligned-fraction {params.aln_frac} \
            --precluster-method {params.preclust_method} --cluster-method {params.clust_method} -t {threads} --min-completeness {params.compl} --max-contamination {params.cont} --quality-formula {params.qual_formula}; fi
        if ! grep -q "coverm" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "coverm " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi  
        """

rule Binning_euksel:
    conda: config["env_1"]
    input:
        df=os.path.join(loc, sample, "tmp/binning/round_euk/clust/euks.tsv"),
        links=os.path.join(loc, sample, "tmp/binning/round_euk/clust/contig_bin.tsv"),
        clust=os.path.join(loc, sample, "tmp/binning/round_euk/clust/clust.tsv"),
    output:
        df=os.path.join(loc, sample, "tmp/binning/round_euk/clust/euks_drep.tsv"),
        id=os.path.join(loc, sample, "tmp/binning/round_euk/clust/euks_keep.txt"),
    params:
        weigth=config["binette_weight"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_binning_eukselect.tsv")
    shell:
        """
        if ! grep -q "r-base" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "r-base " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi
        if ! grep -q "r-tidyverse" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "r-tidyverse " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi
        if [ $(cat {input.clust} | wc -l) -lt 1 ]; then touch {output.df} {output.id} && exit 0; fi

        R --no-echo --silent --args << 'get_euks' > {loc}/{sample}/tmp/binning/round_euk/clust/eukselect.log 2>&1
        options(scipen=999)
        library(tidyverse)
        library(dplyr)

        qc <- read.delim("{input.df}", sep="\t", header=T)
        contigs <- read.delim("{input.links}", sep="\t", header=F)
        colnames(contigs) <- c("contig","bin")

        clust <- read.delim("{input.clust}", sep="\t", header=F)
        colnames(clust) <- c("cluster","file")
        clust$bin <- sub(".fa$", "", sub(".*/", "", clust$file))

        # First check
        contigs_dup <- contigs %>% group_by(bin,contig) %>% filter(n()>1)
        qc <- qc[! qc$bin %in% contigs_dup$bin,]
        clust <- clust[! clust$bin %in% contigs_dup$bin,]

        clust <- merge(clust,qc[,c("bin","completeness","contamination")],by="bin")
        clust$q_score <- clust$completeness - {params.weigth} * clust$contamination

        repeat {{
        clust_tmp <- clust[order(clust$q_score,decreasing=TRUE),]
        clust_tmp <- clust_tmp[!duplicated(clust_tmp$cluster),]
        contigs_tmp <- contigs[contigs$bin %in% clust_tmp$bin,]

        if (any(duplicated(contigs_tmp$contig))) {{
            contigs_tmp_dup <- contigs_tmp %>% group_by(contig) %>% filter(n()>1)
            clust_tmp_dup <- clust_tmp[clust_tmp$bin %in% contigs_tmp_dup$bin,]
            clust <- clust[! clust$bin %in% tail(clust_tmp_dup$bin, 1),]
            message("Removed ",tail(clust_tmp_dup$bin, 1)," due to overlapping contigs")
        }} else {{ break }} }}

        qc <- qc[qc$bin %in% clust_tmp$bin,]
        qc$bin_org <- qc$bin
        qc$bin <- paste0("{sample}.bin.e.",seq_len(nrow(qc)))

        write.table(qc,"{output.df}",quote=F,row.names=FALSE,col.names=TRUE,sep="\t")
        write.table(qc$bin,"{output.id}",quote=F,row.names=FALSE,col.names=FALSE,sep="\t")
        """

rule Binning_aggregate_euk:
    conda: config["env_1"]
    input:
        os.path.join(loc, sample, "tmp/binning/round_euk/clust/euks_drep.tsv"),
    output:
        os.path.join(loc, sample, "tmp/binning/contig_bin_euk.tsv"),
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_binning_eukagr.tsv")
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/binning/bins_euk" ]; then rm -r {loc}/{sample}/tmp/binning/bins_euk; fi
        mkdir {loc}/{sample}/tmp/binning/bins_euk

        if [ $(cat {input} | wc -l) -ge 2 ]; then
        awk 'NR>1 {{print $1, $6}}' {input} | xargs -n2 bash -c 'rsync {loc}/{sample}/tmp/binning/round_euk/bins_innit/$1.fa {loc}/{sample}/tmp/binning/bins_euk/$0.fa'
        grep -H ">" {loc}/{sample}/tmp/binning/bins_euk/*.fa | xargs -n1 basename | sed "s/.fa:>/\t/" | awk -F'\t' '{{print $2 "\t" $1}}' > {output}
        else touch {output}; fi

        if [ -d "{loc}/{sample}/results/bins_euk" ]; then rm -r {loc}/{sample}/results/bins_euk; fi
        rsync -r {loc}/{sample}/tmp/binning/bins_euk {loc}/{sample}/results/.
        """
