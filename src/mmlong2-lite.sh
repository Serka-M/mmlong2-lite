#!/bin/bash
# DESCRIPTION: Wrapper script for running mmlong2-lite Snakemake workflow
# AUTHOR: Mantas Sereika, mase@bio.aau.dk, 2026
# LICENSE: GNU General Public License

# Pre-set default settings
set -eo pipefail
script="$(dirname $(realpath "$0"))"
config=$script/config/mmlong2-lite.yaml
snakefile=$script/mmlong2-lite.smk
help=$script/help.txt
sing=$script/sing-mmlong2-lite
sing_bind=$(grep sing_bind $config | cut -d" " -f2 | tr -d '\r')
version=$(grep version $config | cut -d" " -f2 | tr -d '\r')
zenodo_id=$(grep zenodo_id $config | cut -d" " -f2 | tr -d '\r')
mode=$(grep mode $config | cut -d" " -f2 | tr -d '\r')
fastq=$(grep fastq $config | cut -d" " -f2 | tr -d '\r')
sample=$(grep sample $config | cut -d" " -f2 | tr -d '\r')
loc=$(grep loc $config | cut -d" " -f2 | tr -d '\r')
proc=$(grep -w proc $config | cut -d" " -f2 | tr -d '\r')
cov=$(grep reads_diffcov $config | cut -d" " -f2 | tr -d '\r')
assembler=$(grep assembler $config | cut -d" " -f2 | tr -d '\r')
custom_assembly=$(grep custom_assembly $config | cut -d" " -f2 | tr -d '\r')
assembly_info=$(grep assembly_info $config | cut -d" " -f2 | tr -d '\r')
asm_extra=$(grep asm_extra $config | cut -d" " -f2 | tr -d '\r')
medaka_status=$(grep medaka_status $config | cut -d" " -f2 | tr -d '\r')
medak_mod_pol=$(grep medak_mod_pol $config | cut -d" " -f2 | tr -d '\r')
curation_status=$(grep curation_status $config | cut -d" " -f2 | tr -d '\r')
rename_status=$(grep rename_status $config | cut -d" " -f2 | tr -d '\r')
tiara_status=$(grep tiara_status $config | cut -d" " -f2 | tr -d '\r')
eukfilt_status=$(grep eukfilt_status $config | cut -d" " -f2 | tr -d '\r')
euk_status=$(grep euk_status $config | cut -d" " -f2 | tr -d '\r')
prok_status=$(grep prok_status $config | cut -d" " -f2 | tr -d '\r')
min_contig_len=$(grep min_contig_len $config | cut -d" " -f2 | tr -d '\r')
min_prok_cov=$(grep min_prok_cov $config | cut -d" " -f2 | tr -d '\r')
min_euk_cov=$(grep min_euk_cov $config | cut -d" " -f2 | tr -d '\r')
min_prok_len=$(grep min_prok_len $config | cut -d" " -f2 | tr -d '\r')
min_euk_len=$(grep min_euk_len $config | cut -d" " -f2 | tr -d '\r')
semibin_mod=$(grep semibin_mod $config | cut -d" " -f2 | tr -d '\r')
min_compl_euk=$(grep min_compl_euk $config | cut -d" " -f2 | tr -d '\r')
skip_metabat2_euk=$(grep skip_metabat2_euk $config | cut -d" " -f2 | tr -d '\r')
skip_quickbin_euk=$(grep skip_quickbin_euk $config | cut -d" " -f2 | tr -d '\r')
skip_vamb_euk=$(grep skip_vamb_euk $config | cut -d" " -f2 | tr -d '\r')
skip_semibin2_euk=$(grep skip_semibin2_euk $config | cut -d" " -f2 | tr -d '\r')
skip_remag_euk=$(grep skip_remag_euk $config | cut -d" " -f2 | tr -d '\r')
lorbin_status=$(grep lorbin_status $config | cut -d" " -f2 | tr -d '\r')
cleanup_status=$(grep cleanup_status $config | cut -d" " -f2 | tr -d '\r')
apptainer_status=$(grep apptainer_status $config | cut -d" " -f2 | tr -d '\r')
env_1=$(grep env_1: $config | cut -d" " -f2 | tr -d '\r')
env_2=$(grep env_2 $config | cut -d" " -f2 | tr -d '\r')
env_3=$(grep env_3 $config | cut -d" " -f2 | tr -d '\r')
env_4=$(grep env_4 $config | cut -d" " -f2 | tr -d '\r')
env_5=$(grep env_5 $config | cut -d" " -f2 | tr -d '\r')
env_6=$(grep env_6 $config | cut -d" " -f2 | tr -d '\r')
env_7=$(grep env_7 $config | cut -d" " -f2 | tr -d '\r')
env_8=$(grep env_8 $config | cut -d" " -f2 | tr -d '\r')
env_9=$(grep env_9 $config | cut -d" " -f2 | tr -d '\r')
env_10=$(grep env_10 $config | cut -d" " -f2 | tr -d '\r')
extra=""
dryrun=""
touch=""
conda_prefix=""
stage="OFF"
stages=("OFF" "assembly" "curation" "filtering" "coverage" "singletons" "prokaryotes" "eukaryotes")
binmode=$(grep binmode $config | cut -d" " -f2 | tr -d '\r')
binmodes=("fast" "default" "extended")
rule="Finalise_prok"
apptainer="apptainer"
tmp=$TMPDIR

# Assign inputs
while test $# -gt 0; do
           case "$1" in
                -h | --help) eval "echo \"$(<"$help")\""; exit 0;;
                -v | --version) printf -- "mmlong2-lite, version $version\n"; exit 0;;
                -n | --dryrun) dryrun="-n"; shift;;	
                -t | --touch) touch="--touch"; shift;;	
                -env | --conda_envs_only) apptainer="" && apptainer_status="FALSE" && conda_prefix="--conda-prefix $script/conda"; shift;;	
                -np | --nanopore_reads) shift; fastq=$1 && mode="Nanopore-simplex"; shift;;
                -pb | --pacbio_reads) shift; fastq=$1 && mode="PacBio-HiFi"; shift;;
                -o | --output_dir) shift; sample=$1; shift;;
                -p | --processes) shift; proc=$1; shift;;
                -med | --use_medaka) medaka_status="TRUE"; shift;;
                -who | --use_whokaryote) tiara_status="FALSE"; shift;;
                -srn | --skip_rename) rename_status="FALSE"; shift;;
                -sef | --skip_eukfilt) eukfilt_status="FALSE"; shift;;
                -euk | --bin_eukaryote) euk_status="TRUE"; shift;;
                -spk | --skip_prokaryote) prok_status="FALSE" && euk_status="TRUE" && rule="Finalise_euk"; shift;;
                -cov | --coverage) shift; cov=$1; shift;;
                -bin | --binning_mode) shift; binmode=$1; shift;;
                -run | --run_until) shift; stage=$1; shift;;
                -tmp | --temporary_dir) shift; tmp=$1; shift;;
                -ca | --custom_assembly) shift; custom_assembly=$1 && assembler="custom"; shift;; 
                -cai | --custom_assembly_info) shift; assembly_info=$1; shift;; 
                -mm | --medaka_model) shift; medak_mod_pol=$1 && medaka_status="TRUE"; shift;;				
                -mcl | --min_contig_len) shift; min_contig_len=$1; shift;;
                -mpc | --min_prok_cov) shift; min_prok_cov=$1; shift;;
                -mec | --min_euk_cov) shift; min_euk_cov=$1; shift;;
                -mlp | --min_len_prok) shift; min_prok_len=$1; shift;;
                -mle | --min_len_euk) shift; min_euk_len=$1; shift;;
                -sem | --semibin_model) shift; semibin_mod=$1; shift;;	
                -fly | --use_metaflye) assembler="metaflye"; shift;;
                -dbg | --use_metamdbg) assembler="metamdbg"; shift;;
                -myl | --use_myloasm) assembler="myloasm"; shift;;
                -lrb | --use_lorbin) lorbin_status="TRUE" && binmode="extended"; shift;;
                -mce | --min_compl_euk) shift; min_compl_euk=$1; shift;;  
                -sme | --skip_metabat2_euk) skip_metabat2_euk="TRUE"; shift;;
                -sqe | --skip_quickbin_euk) skip_quickbin_euk="TRUE"; shift;;
                -sve | --skip_vamb_euk) skip_vamb_euk="TRUE"; shift;;
                -sse | --skip_semibin2_euk) skip_semibin2_euk="TRUE"; shift;;
                -sre | --skip_remag_euk) skip_remag_euk="TRUE"; shift;;
                -scr | --skip_assembly_curation) curation_status="FALSE"; shift;;
                -scl | --skip_cleanup) cleanup_status="FALSE"; shift;;
                -r | --rule) shift; rule=$1; shift;;				
                -x | --extra_inputs) shift; extra=$1; shift;;				
                -xa | --extra_asm) shift; asm_extra=$(printf '%s' "$1" | tr ' ' ','); shift;;
                -xbp | --extra_bind_path) shift; sing_bind=$1; shift;;
                *)
                   printf "ERROR: $1 is not a recognized input! \nType mmlong2-lite --help for more information on workflow usage.\n"
                   exit 1;
                   ;;
          esac
  done  

# Basic checks
if [[ ! " ${stages[@]} " =~ " ${stage} " ]]; then echo "ERROR: provided input for early stoppage not recognized" && exit 1; fi
if [[ ! " ${binmodes[@]} " =~ " ${binmode} " ]]; then echo "ERROR: provided setting for binning mode not recognized" && exit 1; fi
if [ $fastq == "none" ]; then printf -- "$usage_text" && echo "ERROR: Sequencing reads not provided" && exit 1; fi
if [ ! -f $fastq ]; then printf -- "$usage_text" && echo "ERROR: Provided sequencing read file not found" && exit 1; fi
if [ ! $cov == "none" ] && [ ! -f $cov ]; then printf -- "$usage_text" && echo "ERROR: Provided dataframe for differential coverage binning not found" && exit 1; fi
if [ $assembler == "custom" ]; then if [ ! -f $custom_assembly ]; then printf -- "$usage_text" && echo "ERROR: Provided custom assembly not found" && exit 1; fi; fi
if [ ! $assembly_info == "FALSE" ]; then if [ ! -f $assembly_info ]; then printf -- "$usage_text" && echo "ERROR: Provided custom assembly information file not found" && exit 1; fi; fi 
bytes_left="$(df --output=avail -B 1 "$(pwd)" | tail -n 1)" && gigabytes_left=$(($bytes_left/1024**3))
if [ $gigabytes_left -lt 100 ]; then  echo "ERROR: Less than 100 GB of free space left in working directory" && exit 1; else echo "$gigabytes_left GB of free space available in working directory"; fi
if ! command -v conda &> /dev/null; then echo "ERROR: Dependency (((conda))) could not be found" && exit 1; fi
if ! command -v snakemake &> /dev/null; then echo "ERROR: Dependency (((snakemake))) could not be found" && exit 1; fi
if ! command -v singularity &> /dev/null; then echo "ERROR: Dependency (((singularity))) could not be found" && exit 1; fi
if ! command -v dirname &> /dev/null; then echo "ERROR: Dependency (((dirname))) could not be found" && exit 1; fi
if [[ "$euk_status" == "TRUE" ]] && [[ "$skip_metabat2_euk" == "TRUE" ]] && [[ "$skip_vamb_euk" == "TRUE" ]] && [[ "$skip_quickbin_euk" == "TRUE" ]] && [[ "$skip_semibin2_euk" == "TRUE" ]] && [[ "$skip_remag_euk" == "TRUE" ]]; then echo "ERROR: No binners selected for eukaryote genome recovery" && exit 1; fi
if [ $mode == "PacBio-HiFi" ]; then medaka_status="FALSE"; fi

# Setup singularity container
if [[ "$apptainer" == "apptainer" ]] && [[ ! -d  "$script/sing-mmlong2-lite" ]]
then
	echo "Downloading and extracting pre-built Singularity container..."
	sleep 10s
	zenodo_get -r $zenodo_id -o $script -v 3
	pv $script/sing-mmlong2-lite-*.tar.gz | pigz -dc - | tar xf - -C $script/.
	rm $script/sing-mmlong2-lite-*.tar.gz
fi

# Configure early stoppage
if [ $stage == "assembly" ]; then if [ $assembler == "myloasm" ]; then rule="Assembly_myloasm"; fi ; fi
if [ $stage == "assembly" ]; then if [ $assembler == "metaflye" ]; then rule="Assembly_metaFlye"; fi ; fi
if [ $stage == "assembly" ]; then if [ $assembler == "metamdbg" ]; then rule="Assembly_metaMDBG"; fi ; fi
if [ $stage == "assembly" ]; then if [ $assembler == "custom" ]; then rule="Assembly_custom"; fi ; fi
if [ $stage == "curation" ]; then rule="Curation_selection" && curation_status="TRUE"; fi
if [ $stage == "filtering" ]; then rule="Filtering_euk"; fi
if [ $stage == "coverage" ]; then rule="Coverage_aggregate"; fi
if [ $stage == "singletons" ]; then rule="Singletons_checkm2"; fi
if [ $stage == "prokaryotes" ]; then rule="Binning_checkm1" && prok_status="TRUE" && euk_status="FALSE"; fi
if [ $stage == "eukaryotes" ]; then rule="Binning_aggregate_euk" && prok_status="FALSE" && euk_status="TRUE"; fi

# Configure input paths
if [[ "$sample" = *"/"* ]]; then if [[ "$sample" = /* ]]; then loc="$(dirname $sample)" && sample="$(basename $sample)"; else loc="$(pwd)/$(dirname $sample)" && sample="$(basename $sample)"; fi; else loc=$(pwd); fi
if [[ ! "$fastq" = /* ]]; then fastq="$(pwd)/$fastq"; fi
if [[ ! "$tmp" = /* ]]; then tmp="$(pwd)/$tmp"; fi
if [[ ! "$cov" = /* && ! "$cov" = "none" ]]; then cov="$(pwd)/$cov"; fi

# Configure cache, temporary file paths, log verbosity, and locale
if [ ! -d ".pycache" ]; then mkdir .pycache; fi
export PYTHONPYCACHEPREFIX="$(pwd)/.pycache"
export TMPDIR=$tmp
export SNAKEMAKE_SOURCE_CACHE=$tmp
export SNAKEMAKE_OUTPUT_CACHE=$tmp
export XDG_CACHE_HOME=$tmp
export PYTHONWARNINGS="ignore::SyntaxWarning"
export LC_ALL=C.UTF-8
export LANG=C.UTF-8

# Configure environments for conda-only mode
if [[ "$apptainer" == "" ]]
then
export CONDA_CHANNEL_PRIORITY=flexible
env_1=$script/envs/env_1.yaml
env_2=$script/envs/env_2.yaml
env_3=$script/envs/env_3.yaml
env_4=$script/envs/env_4.yaml
env_5=$script/envs/env_5.yaml
env_6=$script/envs/env_6.yaml
env_7=$script/envs/env_7.yaml
env_8=$script/envs/env_8.yaml
env_9=$script/envs/env_9.yaml
env_10=$script/envs/env_10.yaml
fi

# Configure Singularity bindings
singularity_binds="-B $(pwd) -B $TMPDIR -B $(dirname $fastq) -B $loc "
if [ ! $custom_assembly == "FALSE" ]; then singularity_binds=$(echo "$singularity_binds -B $(dirname $custom_assembly)"); fi
if [ ! $assembly_info == "FALSE" ]; then singularity_binds=$(echo "$singularity_binds -B $(dirname $assembly_info)"); fi
if [ ! "$sing_bind" == "none" ]; then singularity_binds=$(echo "$singularity_binds -B $sing_bind"); fi
if [ ! $cov == "none" ]; then singularity_binds=$(echo "$singularity_binds -B $(dirname $cov)"); fi
if [ ! $cov == "none" ]; then while read line || [ -n "$line" ]; do
	type="$(echo $line | cut -f1 -d",")"
	reads="$(echo $line | cut -f2 -d",")"
	if [[ ! -f "$reads" ]]; then printf "ERROR: provided read file for differential coverage binning not found \n$reads\n" && exit 1; fi;
	if [[ ! "$type" = "PB" ]] && [[ ! "$type" = "IL" ]] && [[ ! "$type" = "NP" ]]; then printf "ERROR: unrecognised read type for differential coverage binning\n $type with $reads\n" && exit 1; fi;
	singularity_binds=$(echo "$singularity_binds -B $(dirname $reads)")
done < $cov; fi;

# Run snakemake workflow
echo "Starting Snakemake workflow for $sample with $mode reads and $proc multi-threading instances..."
snakemake $dryrun $touch --software-deployment-method conda $apptainer $conda_prefix --apptainer-args "$singularity_binds" --cores $proc --resources usage=100 --rerun-incomplete --nolock -s $snakefile --configfile $config -R $rule --until $rule --config \
  sing=$sing sample=$sample loc=$loc fastq=$fastq proc=$proc mode=$mode reads_diffcov=$cov assembler=$assembler custom_assembly=$custom_assembly assembly_info=$assembly_info asm_extra=$asm_extra \
  medaka_status=$medaka_status medak_mod_pol=$medak_mod_pol curation_status=$curation_status rename_status=$rename_status eukfilt_status=$eukfilt_status tiara_status=$tiara_status binmode=$binmode min_contig_len=$min_contig_len min_prok_cov=$min_prok_cov min_euk_cov=$min_euk_cov prok_status=$prok_status euk_status=$euk_status min_compl_euk=$min_compl_euk \
	min_prok_len=$min_prok_len min_euk_len=$min_euk_len semibin_mod=$semibin_mod lorbin_status=$lorbin_status skip_metabat2_euk=$skip_metabat2_euk skip_quickbin_euk=$skip_quickbin_euk skip_vamb_euk=$skip_vamb_euk skip_semibin2_euk=$skip_semibin2_euk skip_remag_euk=$skip_remag_euk cleanup_status=$cleanup_status apptainer_status=$apptainer_status \
  env_1=$env_1 env_2=$env_2 env_3=$env_3 env_4=$env_4 env_5=$env_5 env_6=$env_6 env_7=$env_7 env_8=$env_8 env_9=$env_9 env_10=$env_10 $extra 2> >(stdbuf -oL grep -Ev "here-document|While bind mounting|underlay of" >&2)
exit 0
