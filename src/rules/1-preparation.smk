# DESCRIPTION: Snakemake rules for the metagenome preparation sections (assembly, polishing, curation, contig filtering and coverage calculation) of the mmlong2 pipeline
# AUTHOR: Mantas Sereika (mase@bio.aau.dk)
# LICENSE: GNU General Public License

rule Assembly_myloasm:
    conda: config["env_4"]
    input:
        fastq
    output:
        asm=os.path.join(loc, sample, "tmp/assembly/assembly.fa"),
        info=os.path.join(loc, sample, get_asm_info("myloasm",config["assembler"],"FALSE")),
    params:
        circ_prob=config["circ_prob"],
        extra=get_asm_params(config["asm_extra"]),
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_assembly_myloasm.tsv")
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/assembly" ]; then rm -r {loc}/{sample}/tmp/assembly; fi
        if [ {mode} == "Nanopore-simplex" ]; then mylo_opt=""; else mylo_opt="--hifi";fi
        myloasm {input} -o {loc}/{sample}/tmp/assembly -t {threads} $mylo_opt {params.extra}
        sed '/^>/ s/_/ /g' {loc}/{sample}/tmp/assembly/assembly_primary.fa > {output.asm}
        grep ">" {loc}/{sample}/tmp/assembly/assembly_primary.fa | cut -c2- > {loc}/{sample}/tmp/assembly/contigs.txt
        sed -e 's/_/\t/g' -e "s/len-//" -e "s/circular-//" -e "s/yes/Y/" -e "s/no/N/" -e "s/possibly/{params.circ_prob}/" -e "s/depth-//" {loc}/{sample}/tmp/assembly/contigs.txt | sed 's/-.*//' | awk -F'\t' '{{t=$3; $3=$4; $4=t; print}}' OFS='\t' > {output.info}
        if ! grep -q "myloasm" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "myloasm " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi		
        """

rule Assembly_metaFlye:
    conda: config["env_2"]
    input:
        fastq
    output:
        asm=os.path.join(loc, sample, "tmp/assembly/assembly.fasta"),
        info=os.path.join(loc, sample, get_asm_info("metaflye",config["assembler"],"FALSE")),
    params:
        extra=get_asm_params(config["asm_extra"]),
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_assembly_metaflye.tsv")
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/assembly" ]; then rm -r {loc}/{sample}/tmp/assembly; fi
        if [ {mode} == "Nanopore-simplex" ]; then flye_opt="--nano-hq"; fi
        if [ {mode} == "PacBio-HiFi" ]; then flye_opt="--read-error 0.01 --pacbio-hifi"; fi
        flye $flye_opt {input} --out-dir {loc}/{sample}/tmp/assembly --threads {threads} --meta {params.extra}
        sed 1d {loc}/{sample}/tmp/assembly/assembly_info.txt | cut -f1,2,3,4 > {output.info}
        if ! grep -q "flye" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "flye " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi		
        """

rule Assembly_metaMDBG:
    conda: config["env_4"]
    input:
        fastq
    output:
        asm=os.path.join(loc, sample, "tmp/assembly/contigs.fasta"),
        info=os.path.join(loc, sample, get_asm_info("metamdbg",config["assembler"],"FALSE")),
    params:
        extra=get_asm_params(config["asm_extra"]),
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_assembly_metamdbg.tsv")
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/assembly" ]; then rm -r {loc}/{sample}/tmp/assembly; fi
        if [ {mode} == "Nanopore-simplex" ]; then mdbg_opt="--in-ont"; else mdbg_opt="--in-hifi";fi
        metaMDBG asm $mdbg_opt {input} --out-dir {loc}/{sample}/tmp/assembly --threads {threads} {params.extra}
        unpigz {output.asm}.gz
        grep ">" {output.asm} | cut -c2- > {loc}/{sample}/tmp/assembly/contigs.txt
        sed -e "s/length=//" -e "s/coverage=//" -e "s/circular=//" -e "s/yes/Y/" -e "s/no/N/" -e 's/ /\t/g' {loc}/{sample}/tmp/assembly/contigs.txt > {output.info}
        if ! grep -q "metamdbg" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "metamdbg " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi		
        """

rule Assembly_custom:
    conda: config["env_2"]
    input:
        config["custom_assembly"],
    output:
        asm=os.path.join(loc, sample, "tmp/assembly/assembly_custom.fa"),
        info=os.path.join(loc, sample, get_asm_info("custom",config["assembler"],"FALSE")),
    params:
        info=config["assembly_info"],
        circ_prob=config["circ_prob"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_assembly_custom.tsv")
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/assembly" ]; then rm -r {loc}/{sample}/tmp/assembly; fi
        mkdir {loc}/{sample}/tmp/assembly

        if head -c2 {input} | grep -q $'\\x1f\\x8b'
        then unpigz -c {input} > {output.asm}
        else rsync {input} {output.asm}; fi
        grep ">" {output.asm} | cut -c2- > {loc}/{sample}/tmp/assembly/contigs.txt

        if grep -q "circular=" {output.asm}
        then sed  -e "s/length=//" -e "s/coverage=//" -e "s/circular=//" -e "s/yes/Y/" -e "s/no/N/" -e 's/ /\t/g' {loc}/{sample}/tmp/assembly/contigs.txt > {output.info}; fi

        if grep -q "_circular-" {output.asm}; then
        sed -e 's/_/\t/g' -e "s/len-//" -e "s/circular-//" -e "s/yes/Y/" -e "s/no/N/" -e "s/possibly/{params.circ_prob}/" -e "s/depth-//" {loc}/{sample}/tmp/assembly/contigs.txt | sed 's/-.*//' | awk -F'\t' '{{t=$3; $3=$4; $4=t; print}}' OFS='\t' > {output.info}
        sed -i '/^>/ s/_/ /g' {output.asm}; fi

        if [ ! "{params.info}" == "FALSE" ]; then 
        if grep -q '^#' {params.info}; then sed 1d {params.info} | cut -f1-4 > {output.info}; else rsync {params.info} {output.info}; fi; fi

        if [ ! -f "{output.info}" ]; then sed -i 's/^\\(>[^ ]*\\).*/\\1/' {output.asm} && seqkit fx2tab -l -n {output.asm} | awk '{{print $1"\t"$2"\t0\tN"}}' > {output.info}; fi
        """

rule Polishing_prep:
    conda: config["env_9"]
    input:
        os.path.join(loc, sample, get_assembler(config["assembler"]))
    output:
        expand("{loc}/{sample}/tmp/polishing/lin_contigs_id_{split}.txt",sample=sample,loc=loc,split=[f"{i:02}" for i in range(1, get_splits(proc,config["medaka_split"]) + 1)]),
        os.path.join(loc, sample, "tmp/polishing/calls_to_draft.bam")
    params:
        splits=get_splits(proc,config["medaka_split"]),
        ram_ref=config["minimap_ram_ref"],
        asm=get_assembler(config["assembler"]),
        seed=config["seed"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_polishing_prep.tsv")
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/polishing" ]; then rm -r {loc}/{sample}/tmp/polishing; fi
        if [ -f "{loc}/{sample}/{params.asm}.map-ont.mmi" ]; then rm {loc}/{sample}/{params.asm}.map-ont.mmi; fi
        if [ -f "{loc}/{sample}/{params.asm}.fai" ]; then rm {loc}/{sample}/{params.asm}.fai; fi
		mkdir {loc}/{sample}/tmp/polishing
        
        grep ">" {input} | cut -c 2- | awk '{{print $1}}' | shuf --random-source=<(yes {params.seed}) > {loc}/{sample}/tmp/polishing/ids.txt
        split -n l/{params.splits} --numeric-suffixes=1 --additional-suffix=.txt -d {loc}/{sample}/tmp/polishing/ids.txt {loc}/{sample}/tmp/polishing/contigs_id_
        for file in {loc}/{sample}/tmp/polishing/contigs_id_*.txt; do filename=lin_$(basename $file)
        awk 'BEGIN {{ ORS = " " }} {{ print }}' $file > {loc}/{sample}/tmp/polishing/$filename; done
        mini_align -I {params.ram_ref}G -i {fastq} -r {input} -m -p {loc}/{sample}/tmp/polishing/calls_to_draft -t {threads}
        """

rule Polishing_consensus:
    conda: config["env_9"]
    input:
        "{contigs}.txt"
    output:
        "{contigs}.hdf"
    params:
        medaka_model=config["medak_mod_pol"],
        medaka_batch=config["medaka_batch"],
    benchmark:
        "{contigs}.tsv"
    resources: usage=config["medaka_usage"]
    threads: proc_sub
    shell:
        """
        if [ -f "{output}" ]; then rm {output}; fi
        medaka inference {loc}/{sample}/tmp/polishing/calls_to_draft.bam {output} --cpu --model {params.medaka_model} --batch {params.medaka_batch} --threads {threads} --region $(cat {input}) 
        """	

rule Polishing_stitch:
    conda: config["env_9"]
    input:
        hdf=expand("{contigs}.hdf", contigs=get_contigs()),
        asm=os.path.join(loc, sample, get_assembler(config["assembler"])),
    output:
        os.path.join(loc, sample, "tmp/polishing/asm_pol.fasta")
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_polishing_stitch.tsv")
    threads: proc
    shell:
        """
        medaka sequence {input.hdf} {input.asm} {output} --threads {threads}
        if [ -f "{loc}/{sample}/tmp/polishing/lin_contigs_id_01.tsv" ]; then
        mv {loc}/{sample}/tmp/polishing/*.tsv {loc}/{sample}/tmp/logs
        for f in {loc}/{sample}/tmp/logs/lin_contigs_id_*.tsv; do mv "$f" "$(dirname "$f")/usage_polishing_consensus-${{f##*/lin_contigs_id_}}"; done; fi
        if ! grep -q "medaka" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then pip list | grep -w "medaka" | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi        
        """	

rule Curation_map:
    conda: config["env_2"]
    input:
        reads=fastq,
        asm=os.path.join(loc, sample, get_assembly(mode,config["assembler"],config["medaka_status"],"FALSE")),
    output:
        os.path.join(loc, sample, "tmp/curation/mapping.bam")
    params:
        map=lambda wildcards: get_map_mode(get_read_type(mode),config["minimap_np"],config["minimap_pb"],config["minimap_il"]),
        ram_ref=config["minimap_ram_ref"],
        ram_reads=config["minimap_ram_reads"],
        ratio=config["minimap_ratio"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_curation_map.tsv")
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/curation" ]; then rm -r {loc}/{sample}/tmp/curation; fi
        mkdir {loc}/{sample}/tmp/curation
        minimap2 -I {params.ram_ref}G -K {params.ram_reads}G -t {threads} --secondary-seq -p {params.ratio} -ax {params.map} {input.asm} {input.reads} | samtools view --threads $(({threads} / 2)) -Sb -F 4 - | samtools sort --threads $(({threads} / 2)) --write-index -o {output}##idx##{output}.bai -
        if ! grep -q "minimap2" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "minimap2 " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi
        if ! grep -q "samtools" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "samtools " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi        
        """	

rule Curation_screening:
    conda: config["env_2"]
    input:
        os.path.join(loc, sample, "tmp/curation/mapping.bam"),
    output:
        cov=os.path.join(loc, sample, "tmp/curation/errors-zero_cov.txt"),
        cov_filt=os.path.join(loc, sample, "tmp/curation/errors-zero_cov.tsv"),
        clip=os.path.join(loc, sample, "tmp/curation/errors-clipping.txt"),
        clip_filt=os.path.join(loc, sample, "tmp/curation/errors-clipping.tsv"),
    params:
        clip_ratio=config["clip_ratio"],
        min_dist=config["min_dist"],
        min_len_zerocov=config["min_len_zerocov"],
        min_cov_clip=config["min_cov_clip"],
        min_cov_clip_all=config["min_cov_clip_all"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_curation_screening.tsv")
    shell:
        """
        org_dir=$(pwd) && cd {loc}/{sample}/tmp/curation
        bam_error_detector --clipping-ratio {params.clip_ratio} --min-dist-to-end {params.min_dist} {input} errors
        cd $org_dir

        awk -F'\t' '$4 > {params.min_len_zerocov}' {output.cov} > {output.cov_filt}
        awk -F'\t' '$5 >= {params.min_cov_clip_all} && $6 >= {params.min_cov_clip}' {output.clip} > {output.clip_filt}
        if ! grep -q "rust" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "rust " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi      
        if ! grep -q "bam_error_detector" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then bam_error_detector --version | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi    
        """	

rule Curation_aggregate:
    conda: config["env_1"]
    input:
        asm=os.path.join(loc, sample, "tmp/assembly/assembly_info.tsv"),
        cov=os.path.join(loc, sample, "tmp/curation/errors-zero_cov.tsv"),
        clip=os.path.join(loc, sample, "tmp/curation/errors-clipping.tsv"),
    output:
        df=os.path.join(loc, sample, "tmp/curation/contigs.tsv"),
        con=os.path.join(loc, sample, "tmp/curation/contigs.txt"),
        asm=os.path.join(loc, sample, "tmp/curation/assembly_info.tsv"),
    params:
        len=config["min_contig_len"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_curation_aggregate.tsv")
    shell:
        """
        if ! grep -q "r-base" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "r-base " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi
        if ! grep -q "r-tidyverse" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "r-tidyverse " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi  

        R --no-echo --silent --args << 'curation' > {loc}/{sample}/tmp/curation/curation.log 2>&1
        library(tidyverse)
        options(scipen=999)

        cov <- read.delim("{input.cov}", sep="\t", header=T)
        if (nrow(cov) == 0) {{cov <- tibble(contig="contig_placeholder", length=1000, range="900-1000", range_size=100)}}

        cov <- cov[order(cov$range,decreasing=FALSE),]
        cov <- cov[order(cov$contig,decreasing=FALSE),]
        cov_ <- cov %>% separate(range, into = c("start", "end"), sep = "-", convert = TRUE) %>% select(contig, length, start, end)

        cov_ <- cov_ %>% arrange(contig, start, end) %>% group_by(contig, length) %>% mutate(group = cumsum(coalesce(start > lag(end) + 1, TRUE))) %>%
                group_by(contig, length, group) %>% summarise(start = min(start), end = max(end), .groups = "drop")

        cov_ <- cov_ %>% group_by(contig, length) %>% arrange(start) %>% summarise(safe = list({{
                        starts <- c(0, end)
                        ends <- c(start, unique(length))
                        tibble(start = starts, end = ends) %>% filter(start < end) }}), .groups = "drop") %>% unnest(safe) %>% select(contig, start, end)

        clip <- read.delim("{input.clip}", sep="\t", header=T)
        if (nrow(clip) == 0) {{clip <- tibble(contig="contig_placeholder", length=1000, pos=500, relative_pos=0.5, cov=1, clipping=1, clipping_ratio=1)}}

        clip <- clip[order(clip$pos,decreasing=FALSE),]
        clip <- clip[order(clip$contig,decreasing=FALSE),]
        clip_ <- clip %>% group_by(contig) %>% summarise(length = max(length), breakpoints = list(sort(unique(pos))))

        clip_ <- clip_ %>% rowwise() %>% mutate(chunks = list({{
                        bp <- breakpoints
                        starts <- c(0, bp)
                        ends <- c(bp - 1, length)
                        data.frame(start = starts, end = ends) }})) %>% tidyr::unnest(chunks) %>% select(contig, start, end)

        sub <- bind_rows(cov_, clip_) %>% pivot_longer(cols = c(start, end), names_to = "type", values_to = "bp") %>% distinct() %>% arrange(contig, bp) %>%
               group_by(contig) %>% summarise(bps = list(bp), .groups = "drop") %>% mutate(intervals = map(bps, ~{{ tibble(
                        start = head(.x, -1), end = tail(.x, -1)) }})) %>% select(contig, intervals) %>% unnest(intervals)

        regions <- bind_rows(cov_, clip_)
        sub <- sub %>% semi_join(regions, by = "contig") %>% rowwise() %>% filter(any(
                        regions$contig == contig &
                        regions$start <= end &
                        regions$end >= start )) %>% ungroup() %>% filter(abs(end - start) >= {params.len})
        write.table(sub,"{output.df}",quote=F,row.names=FALSE,col.names=FALSE,sep="\t")

        asm <- read.delim("{input.asm}", sep="\t", header=F)
        asm_ <- asm[! asm$V1 %in% regions$contig,]

        sub_ <- sub
        sub_$len <- sub_$end - sub_$start 
        sub_$start <- sub_$start + 1
        sub_$cov <- NA
        sub_$circ <- "N"
        sub_$contig_ <- paste0(sub_$contig,"_",sub_$start,"-",sub_$end)

        colnames(asm_) <- c("contig_","len","cov","circ")
        asm_ <- rbind(asm_,sub_[,c("contig_","len","cov","circ")])
        write.table(asm_,"{output.asm}",quote=F,row.names=FALSE,col.names=FALSE,sep="\t")
        write.table(asm[! asm$V1 %in% regions$contig,]$V1,"{output.con}",quote=F,row.names=FALSE,col.names=FALSE,sep="\t")  
        """	

rule Curation_selection:
    conda: config["env_2"]
    input:
        asm=os.path.join(loc, sample, get_assembly(mode,config["assembler"],config["medaka_status"],"FALSE")),
        df=os.path.join(loc, sample, "tmp/curation/contigs.tsv"),
        con=os.path.join(loc, sample, "tmp/curation/contigs.txt"),
    output:
        os.path.join(loc, sample, "tmp/curation/asm_curated.fasta")
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_curation_selection.tsv")
    shell:
        """
        seqkit grep -f {input.con} {input.asm} > {output}
        seqkit subseq --bed {input.df} {input.asm} | sed '/^>/ s/...$//' >> {output}
        """	

rule Filtering_length:
    conda: config["env_2"]
    input:
        asm=os.path.join(loc, sample, get_assembly(mode,config["assembler"],config["medaka_status"],config["curation_status"])),
        info=os.path.join(loc, sample, get_asm_info("TRUE","TRUE",config["curation_status"])),
    output:
        asm=os.path.join(loc, sample, "tmp/filtering/asm_filt_len.fasta"),
        info=os.path.join(loc, sample, "tmp/filtering/assembly_info.tsv"),
    params:
        len=config["min_contig_len"],
        rename=config["rename_status"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_filtering_length.tsv")
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/filtering" ]; then rm -r {loc}/{sample}/tmp/filtering; fi
        mkdir {loc}/{sample}/tmp/filtering

        if [ {params.rename} == "TRUE" ]; then
        seqkit seq -m {params.len} {input.asm} | seqkit replace -p "\\s.+" | seqkit replace -p ^ -r {sample}_ > {output.asm}
        sed 's/^/{sample}_/' {input.info} > {output.info}
        else
        seqkit seq -m {params.len} {input.asm} | seqkit replace -p "\\s.+" > {output.asm}
        cat {input.info} > {output.info}; fi

        rsync {output.asm} {loc}/{sample}/results/{sample}_assembly.fasta
        if ! grep -q "seqkit" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "seqkit " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi        
        """	

rule Filtering_tiara:
    conda: config["env_3"]
    input:
        os.path.join(loc, sample, "tmp/filtering/asm_filt_len.fasta")
    output:
        prok=os.path.join(loc, sample, "tmp/filtering/contigs_filt_prok.txt"),
        euk=os.path.join(loc, sample, "tmp/filtering/contigs_filt_euk.txt"),
    params:
        len=config["min_contig_len"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_filtering_tiara.tsv")
    threads: proc
    shell:
        """
        tiara -i {input} -t {threads} -m {params.len} -o {loc}/{sample}/tmp/filtering/tiara.tsv
        cut -f1,2 {loc}/{sample}/tmp/filtering/tiara.tsv | (grep -e "prokarya" -e "bacteria" -e "archaea" -e "unknown" - || true) | cut -f1 | sort > {output.prok}
        cut -f1,2 {loc}/{sample}/tmp/filtering/tiara.tsv | (grep -e "eukarya" -e "organelle" -e "unknown" - || true) | cut -f1 | sort > {output.euk}
        if ! grep -q "tiara" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "tiara " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi     
        """

rule Filtering_whokaryote:
    conda: config["env_3"]
    input:
        os.path.join(loc, sample, "tmp/filtering/asm_filt_len.fasta")
    output:
        prok=os.path.join(loc, sample, "tmp/filtering/whokaryote/prokaryote_contig_headers.txt"),
        euk=os.path.join(loc, sample, "tmp/filtering/whokaryote/eukaryote_contig_headers.txt"),
    params:
        len=config["min_contig_len"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_filtering_whokaryote.tsv")
    threads: proc
    shell:
        """
        if [ -d "{loc}/{sample}/tmp/filtering/whokaryote" ]; then rm -r {loc}/{sample}/tmp/filtering/whokaryote; fi
        whokaryote.py --contigs {input} --threads {threads} --minsize {params.len} --outdir {loc}/{sample}/tmp/filtering/whokaryote 
        if ! grep -q "whokaryote" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "whokaryote " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi     
        """

rule Filtering_eukskip:
    conda: config["env_3"]
    input:
        os.path.join(loc, sample, "tmp/filtering/assembly_info.tsv")
    output:
        os.path.join(loc, sample, "tmp/filtering/contigs_filt_all.txt"),
    params:
        len=config["min_contig_len"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_filtering_eukskip.tsv")
    threads: proc
    shell:
        """
        cut -f1 {input} > {output}
        """

rule Filtering_euk:
    conda: config["env_2"]
    input:
        asm=os.path.join(loc, sample, "tmp/filtering/asm_filt_len.fasta"),
        prok=os.path.join(loc, sample, "tmp/filtering", get_prok(config["tiara_status"],config["eukfilt_status"])),
        euk=os.path.join(loc, sample, "tmp/filtering", get_euk(config["tiara_status"],config["eukfilt_status"])),
    output:
        prok=os.path.join(loc, sample, "tmp/filtering/asm_filt_prok.fasta"),
        euk=os.path.join(loc, sample, "tmp/filtering/asm_filt_euk.fasta"),
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_filtering_eukaryotes.tsv")
    shell:
        """
        seqkit grep -f {input.prok} {input.asm} > {output.prok}
        seqkit grep -f {input.euk} {input.asm} > {output.euk} 
        if ! grep -q "seqkit" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "seqkit " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi    
        """

rule Singletons_circ:
    conda: config["env_2"]
    input:
        asm=os.path.join(loc, sample, "tmp/filtering/asm_filt_len.fasta"),
        info=os.path.join(loc, sample, "tmp/filtering/assembly_info.tsv"),
        prok=os.path.join(loc, sample, "tmp/filtering/asm_filt_prok.fasta"),
    output:
        os.path.join(loc, sample, "tmp/binning/singl/contig_c.tsv")
    params:
        mag_len=config["min_prok_len"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_singletons_circ.tsv")
    resources: usage=100
    shell:
        """
        if [ ! -d "{loc}/{sample}/tmp/binning" ]; then mkdir {loc}/{sample}/tmp/binning; fi
        if [ ! -d "{loc}/{sample}/tmp/binning/singl" ]; then mkdir {loc}/{sample}/tmp/binning/singl; fi
        if [ ! -d "{loc}/{sample}/tmp/binning/singl/innit" ]; then mkdir {loc}/{sample}/tmp/binning/singl/innit; fi
        
        awk -F "\t" '{{ if (($2 >= {params.mag_len}) && ($4 == "Y")) {{print $1 "\t" "{sample}.bin.c." ++i; next}} }}' {input.info} > {output}
        find {loc}/{sample}/tmp/binning/singl/innit -name "{sample}.bin.c.*.fa" -type f -delete
        if [ $(cat {output} | wc -l) -ge 1 ]; then cat {output} | xargs -i --max-procs=1 -n 2 bash -c 'samtools faidx {input.asm} $0 > {loc}/{sample}/tmp/binning/singl/innit/$1.fa'; fi;
        
        find {loc}/{sample}/tmp/binning/singl/innit -name "*.fa" -type f -size -2k -delete
        if ! grep -q "samtools" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "samtools " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi  
        """

rule Singletons_lin:
    conda: config["env_2"]
    input:
        asm=os.path.join(loc, sample, "tmp/filtering/asm_filt_len.fasta"),
        info=os.path.join(loc, sample, "tmp/filtering/assembly_info.tsv"),
        circ=os.path.join(loc, sample, "tmp/binning/singl/contig_c.tsv"),
    output:
        os.path.join(loc, sample, "tmp/binning/singl/contig_l.tsv")
    params:
        mag_len=config["min_smag_len"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_singletons_lin.tsv")
    shell:
        """
        awk -F "\t" '{{ if (($2 >= {params.mag_len}) && ($4 == "N")) {{print $1 "\t" "{sample}.bin.s." ++i; next}} }}' {input.info} > {output}
        find {loc}/{sample}/tmp/binning/singl/innit -name "{sample}.bin.s.*.fa" -type f -delete
        if [ $(cat {output} | wc -l) -ge 1 ]; then cat {output} | xargs -i --max-procs=1 -n 2 bash -c 'samtools faidx {input.asm} $0 > {loc}/{sample}/tmp/binning/singl/innit/$1.fa'; 
        elif [ $(cat {input.circ} | wc -l) -lt 1 ]; then seqkit fx2tab -nl {input.asm} | sort -k2,2nr | head -n1 | cut -f1 | seqkit grep -f - {input.asm} > {loc}/{sample}/tmp/binning/singl/innit/{sample}.bin.s.1.fa || true
        grep ">" {loc}/{sample}/tmp/binning/singl/innit/{sample}.bin.s.1.fa | cut -c 2- | awk '{{print $1 "\t" "{sample}.bin.s.1"}}' > {output}; fi
        find {loc}/{sample}/tmp/binning/singl/innit -name "*.fa" -type f -size -2k -delete
        """

rule Singletons_checkm2:
    conda: config["env_8"]
    input:
        os.path.join(loc, sample, "tmp/binning/singl/contig_l.tsv")
    output:
        os.path.join(loc, sample, "tmp/binning/singl/checkm2.tsv")
    params:
        db=config["db_checkm2"],
        compl_c=config["min_compl_circ"],
        compl_l=config["min_compl_lin"],
        cont=config["min_cont_singl"],
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_singletons_qc.tsv")
    threads: proc
    shell:
        """
        export CHECKM2DB="{params.db}"
        
        if [ -d "{loc}/{sample}/tmp/binning/singl/checkm2" ]; then rm -r {loc}/{sample}/tmp/binning/singl/checkm2; fi
        checkm2 predict -x .fa -i {loc}/{sample}/tmp/binning/singl/innit -o {loc}/{sample}/tmp/binning/singl/checkm2 -t {threads}
        cut -f1,2,3,4,6 {loc}/{sample}/tmp/binning/singl/checkm2/quality_report.tsv > {output}
        if ! grep -q "checkm2" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then echo "checkm2,$(checkm2 --version)" >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi 
        
        awk -F "\t" '{{ if ((($1 ~ /bin.c/) && ($2 >= {params.compl_c}) && ($3 <= {params.cont})) || (($1 ~ /bin.s/) && ($2 >= {params.compl_l}) && ($3 <= {params.cont}))) {{print $1}} }}' {output} > {loc}/{sample}/tmp/binning/singl/bins_keep.txt
        if [ -d "{loc}/{sample}/tmp/binning/singl/bins" ]; then rm -r {loc}/{sample}/tmp/binning/singl/bins; fi
        mkdir {loc}/{sample}/tmp/binning/singl/bins
		
        if [ $(cat {loc}/{sample}/tmp/binning/singl/bins_keep.txt | wc -l) -ge 1 ]; then
        cat {loc}/{sample}/tmp/binning/singl/bins_keep.txt | xargs -i --max-procs=1 bash -c 'rsync {loc}/{sample}/tmp/binning/singl/innit/{{}}.fa {loc}/{sample}/tmp/binning/singl/bins/.'
        grep ">" {loc}/{sample}/tmp/binning/singl/bins/*.fa | cut -f2 -d">" | sort > {loc}/{sample}/tmp/binning/singl/binned.txt; fi
        """

rule Coverage_prep:
    conda: config["env_2"]
    input:
        prok=os.path.join(loc, sample, "tmp/filtering/asm_filt_prok.fasta"),
    output:
        reads=os.path.join(loc, sample, "tmp/binning/mapping/reads.csv"),
        links=expand("{loc}/{sample}/tmp/binning/mapping/{reads}.lnk",sample=sample,loc=loc,reads=get_reads(reads_diffcov,fastq,mode)),			
    params:
        reads=fastq,
        type=get_read_type(mode),
        extra=reads_diffcov,
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_coverage_prep.tsv")
    resources: usage=100
    shell:
        """
        if [ ! -d "{loc}/{sample}/tmp/binning" ]; then mkdir {loc}/{sample}/tmp/binning; fi
        if [ -d "{loc}/{sample}/tmp/binning/mapping" ]; then rm -r {loc}/{sample}/tmp/binning/mapping; fi
        mkdir {loc}/{sample}/tmp/binning/mapping
        
        printf "{params.type},{params.reads}\n" > {output.reads}
        if [ -f {params.extra} ]; then cat {params.extra} >> {output.reads}; fi
        
        count=0
        cat {output.reads} | while read line || [ -n "$line" ]; do
        count=$(($count + 1))
        type="$(echo $line | cut -f1 -d",")"
        reads="$(echo $line | cut -f2 -d",")"
        ln -s $reads {loc}/{sample}/tmp/binning/mapping/${{count}}-${{type}}.lnk; done
        """

rule Coverage_map:
    conda: config["env_2"]
    input:
        "{loc}/{sample}/tmp/binning/mapping/{reads}.lnk"
    output:
        bam="{loc}/{sample}/tmp/binning/mapping/{reads}.bam",
        cov="{loc}/{sample}/tmp/binning/mapping/{reads}.tsv",
        tr1="{loc}/{sample}/tmp/binning/mapping/{reads}_tr1.tsv",
        tr2="{loc}/{sample}/tmp/binning/mapping/{reads}_tr2.tsv",
    params:
        map=lambda wildcards: get_map_mode(wildcards.reads,config["minimap_np"],config["minimap_pb"],config["minimap_il"]),
        ram_ref=config["minimap_ram_ref"],
        ram_reads=config["minimap_ram_reads"],
        ident=lambda wildcards: get_map_ident(wildcards.reads,config["np_map_ident"],config["pb_map_ident"],config["il_map_ident"]),
        tr1=lambda wildcards: get_map_cut1(wildcards.reads),
        tr2=lambda wildcards: get_map_cut2(wildcards.reads),
    benchmark:
        "{loc}/{sample}/tmp/logs/usage_coverage_map-{reads}.tsv"
    resources: usage=config["mapping_usage"]
    threads: proc * (config["mapping_usage"]/100)
    shell:
        """
        if [ -f "{output.bam}" ]; then rm {output.bam}; fi
        minimap2 -I {params.ram_ref}G -K {params.ram_reads}G -t {threads} -ax {params.map} {loc}/{sample}/tmp/filtering/asm_filt_len.fasta {input} | samtools view --threads $(({threads} / 2)) -Sb -F 2308 - | samtools sort --threads $(({threads} / 2)) --write-index -o {output.bam}##idx##{output.bam}.bai -
        jgi_summarize_bam_contig_depths {output.bam} --percentIdentity {params.ident} --outputDepth {output.cov}
        cut -f{params.tr1} {output.cov} > {output.tr1}
        cut -f{params.tr2} {output.cov} > {output.tr2}
        """

rule Coverage_aggregate:
    conda: config["env_2"]
    input:
        cov1=expand("{loc}/{sample}/tmp/binning/mapping/{reads}_tr1.tsv",sample=sample,loc=loc,reads=get_reads(reads_diffcov,fastq,mode)),
        cov2=expand("{loc}/{sample}/tmp/binning/mapping/{reads}_tr2.tsv",sample=sample,loc=loc,reads=get_reads(reads_diffcov,fastq,mode)),
    output:
        cov1=os.path.join(loc, sample, "tmp/binning/mapping/cov_all.tsv"),
        cov2=os.path.join(loc, sample, "tmp/binning/mapping/cov_all_sub.tsv"),
    benchmark:
        os.path.join(loc, sample, "tmp/logs/usage_coverage_aggregate.tsv")
    shell:
        """
        paste -d "\t" $(ls -v {input.cov1}) > {output.cov1}
        paste -d "\t" $(ls -v {input.cov2}) | sed "s/contigName/contigname/" > {output.cov2}
        if ! grep -q "minimap2" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "minimap2 " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi
        if ! grep -q "metabat2" {loc}/{sample}/tmp/dep_mmlong2-lite.csv; then conda list | grep -w "metabat2 " | tr -s ' ' | awk '{{print $1","$2}}' >> {loc}/{sample}/tmp/dep_mmlong2-lite.csv; fi		
        """
