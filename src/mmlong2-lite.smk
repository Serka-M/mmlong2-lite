# DESCRIPTION: Snakemake workflow for recovering metagenome assembled genomes (MAGs) with long reads from PacBio or Nanopore sequencing data
# AUTHOR: Mantas Sereika (mase@bio.aau.dk)
# LICENSE: GNU General Public License

import os
import sys
import glob
import re
import pandas as pd
import numpy as np

shell.executable("/bin/bash")
singularity: config["sing"]
wf_v=config["version"]

mode=config["mode"]
modes=['Nanopore-simplex', 'PacBio-HiFi']
assemblers=['myloasm', 'metaflye', 'metamdbg', 'custom']

binmode=config["binmode"]
binmodes=['fast', 'default', 'extended']

proc=config["proc"]
proc_sub=config["proc_sub"]

loc=config["loc"]
sample=config["sample"]

fastq=config["fastq"]
reads_diffcov=config["reads_diffcov"]

semibin_mod=config["semibin_mod"]
semibin_mods=['human_gut', 'dog_gut', 'ocean', 'soil', 'cat_gut', 'human_oral', 'mouse_gut', 'pig_gut', 'built_environment', 'wastewater', 'chicken_caecum', 'global']

semibin_prot=config["semibin_prot"]
semibin_prots=['prodigal', 'fast-naive']

def get_reads(reads_diffcov,fastq,mode):
    reads_diffcov_main = pd.DataFrame([get_read_type(mode),fastq]).T

    if reads_diffcov != "none":
        reads_diffcov_sup = pd.read_csv(reads_diffcov, index_col=None, header=None)
        reads_diffcov_main = pd.concat([reads_diffcov_main,reads_diffcov_sup])

    reads_diffcov_main.columns = ['type', 'reads']
    reads_diffcov_main['pos'] = np.arange(1, len(reads_diffcov_main) + 1).astype(str)
    reads_diffcov_main['id'] = reads_diffcov_main[["pos","type"]].agg("-".join, axis=1)
    reads = reads_diffcov_main['id'].tolist()
    return reads

def get_read_type(mode):
    if mode == "PacBio-HiFi": return "PB"
    else: return "NP"

def get_splits(proc,medaka_split):
    if int(proc/4) < medaka_split: return int(proc/4)
    else: return medaka_split

def get_contigs():	
    contigs_all=expand("{loc}/{sample}/tmp/polishing/lin_contigs_id_{split}",sample=sample,loc=loc,split=[f"{i:02}" for i in range(1, get_splits(proc,config["medaka_split"]) + 1)])
    return contigs_all

def get_assembler(assembler):
    if assembler == "metaflye": return "tmp/assembly/assembly.fasta"
    if assembler == "metamdbg": return "tmp/assembly/contigs.fasta"
    if assembler == "myloasm": return "tmp/assembly/assembly.fa"
    else: return "tmp/assembly/assembly_custom.fa"

def get_assembly(mode,assembler,medaka_status,curation):
    if curation == "TRUE": return "tmp/curation/asm_curated.fasta"
    elif mode == "PacBio-HiFi": return get_assembler(assembler)
    elif medaka_status == "TRUE": return "tmp/polishing/asm_pol.fasta"
    else: return get_assembler(assembler)

def get_asm_info(assembler,assembler_config,curation):
    if curation == "TRUE": return "tmp/curation/assembly_info.tsv"
    elif assembler == assembler_config: return "tmp/assembly/assembly_info.tsv"
    else: return assembler

def get_asm_params(asm_extra):
    if asm_extra == "FALSE": return " "
    else: return re.sub(","," ",asm_extra)

def get_prok(tiara_status,eukfilt_status):
    if eukfilt_status == "FALSE": return "contigs_filt_all.txt"
    if tiara_status == "TRUE": return "contigs_filt_prok.txt"
    else: return "whokaryote/prokaryote_contig_headers.txt"

def get_euk(tiara_status,eukfilt_status):
    if eukfilt_status == "FALSE": return "contigs_filt_all.txt"
    if tiara_status == "TRUE": return "contigs_filt_euk.txt"
    else: return "whokaryote/eukaryote_contig_headers.txt"

def get_map_mode(reads,np_map_mode,pb_map_mode,il_map_mode):
    if "NP" in reads: return np_map_mode
    if "PB" in reads: return pb_map_mode
    if "IL" in reads: return il_map_mode

def get_map_ident(reads,np_map_ident,pb_map_ident,il_map_ident):
    if "NP" in reads: return np_map_ident
    if "PB" in reads: return pb_map_ident
    if "IL" in reads: return il_map_ident

def get_map_cut1(reads):
    if reads.startswith("1-"): return "1,2,3,4,5"
    else: return "4,5"

def get_map_cut2(reads):
    if reads.startswith("1-"): return "1,4"
    else: return "4"

def get_semibin_cov(loc,sample,reads_diffcov,round,semibin_mod,semibin_train_r1):
    if round == "1" and reads_diffcov == "none" and semibin_train_r1 == "FALSE": return ''.join([str(loc), "/", str(sample), "/tmp/binning/mapping/1-*.bam", " --environment ", str(semibin_mod), " "])
    if round in ("3", "euk"): return ''.join([str(loc), "/", str(sample), "/tmp/binning/mapping/1-*.bam", " --environment ", str(semibin_mod), " "])
    else: return ''.join([str(loc), "/", str(sample), "/tmp/binning/mapping/*.bam "])

def get_rounds(binmode):
    if binmode == "extended": return list(range(1, 5))
    else: return 1

def get_qc_in(loc,sample,round, binmode, lorbin):
    if binmode == "fast": return ''.join([str(loc), "/", str(sample), "/tmp/binning/round_1/semibin/bins_info.tsv"])
    if round == "2" and lorbin == "TRUE": return ''.join([str(loc), "/", str(sample), "/tmp/binning/round_", str(round), "/lorbin/label.csv"])
    if round == "2": return ''.join([str(loc), "/", str(sample), "/tmp/binning/round_", str(round), "/comebin/comebin_res/comebin_res.tsv"])
    elif round == "4": return ''.join([str(loc), "/", str(sample), "/tmp/binning/round_", str(round), "/quickbin/quickbin.txt"])
    else: return ''.join([str(loc), "/", str(sample), "/tmp/binning/round_", str(round), "/binette/final_bins_quality_reports.tsv"])

def get_qc_dir(loc,sample,round,binmode, lorbin):
    if binmode == "fast": return ''.join([str(loc), "/", str(sample), "/tmp/binning/round_1/semibin/output_bins"])
    if round == "2" and lorbin == "TRUE": return ''.join([str(loc), "/", str(sample), "/tmp/binning/round_", str(round), "/lorbin/output_bins"])
    if round == "2": return ''.join([str(loc), "/", str(sample), "/tmp/binning/round_", str(round), "/comebin/comebin_res/comebin_res_bins"])
    elif round == "4": return ''.join([str(loc), "/", str(sample), "/tmp/binning/round_", str(round), "/quickbin"])
    else: return ''.join([str(loc), "/", str(sample), "/tmp/binning/round_", str(round), "/binette/final_bins"])

def get_compl(round,binmode,compl1,compl2):
    if round == "1" and binmode == "extended": return compl1
    else: return compl2

def get_cont(round,binmode,cont1,cont2):
    if round == "1" and binmode == "extended": return cont1
    else: return cont2

onstart:
    from snakemake.utils import min_version
    min_version("9.9.0")
    if not os.path.exists(fastq): sys.exit(print("Read input (((",fastq,"))) not found. Aborting..."))
    if not os.path.isdir(loc): sys.exit(print("Provided path for output (((",loc,"))) not found. Aborting..."))
    if reads_diffcov != "none" and not os.path.exists(reads_diffcov): sys.exit(print("Dataframe for differential coverage binning (((",reads_diffcov,"))) not found. Aborting..."))
    if not mode in modes: sys.exit(print("Provided workflow mode (((",mode,"))) not recognised. Aborting..."))
    if not config["assembler"] in assemblers: sys.exit(print("Provided assembler (((",config["assembler"],"))) not recognised. Aborting..."))
    if not binmode in binmodes: sys.exit(print("Provided binning mode (((",binmode,"))) not recognised. Aborting..."))
    if not semibin_mod in semibin_mods: sys.exit(print("Provided model for SemiBin (((",semibin_mod,"))) not recognised. Aborting..."))
    if not semibin_prot in semibin_prots: sys.exit(print("Provided gene predictor for SemiBin (((",semibin_prot,"))) not recognised. Aborting..."))
    if len(os.path.join(loc, sample)) > 85: sys.exit(print("Path for provided output too long: (((",os.path.join(loc, sample),")))\nPlease re-run with different output location or shorter name."))
    if not os.path.exists(os.path.join(loc, sample)): os.makedirs(os.path.join(loc, sample))
    if not os.path.exists(os.path.join(loc, sample, "results")): os.makedirs(os.path.join(loc, sample, "results"))
    if not os.path.exists(os.path.join(loc, sample, "tmp")): os.makedirs(os.path.join(loc, sample, "tmp"))
    if not os.path.exists(os.path.join(loc, sample, "tmp", "logs")): os.makedirs(os.path.join(loc, sample, "tmp", "logs"))
    if not os.path.exists(os.path.join(loc, sample, "tmp", "dep_mmlong2-lite.csv")): 
	    with open(os.path.join(loc, sample, "tmp", "dep_mmlong2-lite.csv"), 'w') as f:
		    f.write("dependency,version\n")

onsuccess:
    from datetime import datetime
    now = datetime.now()
    print("MAG recovery for", sample, "with version", wf_v, "completed at", now.strftime("%Y/%m/%d %H:%M:%S"))
    print("Thank you for using mmlong2-lite")

onerror:
    from datetime import datetime
    now = datetime.now()
    print(f"[{now.strftime('%Y/%m/%d %H:%M:%S')}]: ERROR during MAG recovery for '{sample}'")
    print("Inspect Snakemake log files for troubleshooting")

include: "rules/1-preparation.smk"
include: "rules/2-binning.smk"
include: "rules/3-postprocessing.smk"
