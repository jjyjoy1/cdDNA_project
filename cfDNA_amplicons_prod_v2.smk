import re
import glob
import shutil
import gzip
import subprocess
import pandas as pd
import configparser
import os

# -------------------------------------
# 1) Read safeseq config file
# -------------------------------------
base_dir = '/home/jjy/Documents/GitHubRepos/cfDNA_cfRNA_project'
config_root = '/home/jjy/Documents/GitHubRepos/cfDNA_cfRNA_project/config'
Assay = 'cfDNA'
NextSeqDir = '/home/jjy/Documents/GitHubRepos/cfDNA_cfRNA_project/data/cfDNA_testset_HNLGFAFX2'     #Suppose NGS run folder with flowcell ID

config_file = config_root + '/' + Assay + '/' + Assay + '_amplicons_config.cfg'

Config = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
Config.read(config_file)
###Config.set("Files", "ConfigPath", os.path.split(config_file)[0])

amplicon_file = Config.get("Files", "amplicon_file")
quantispike_file = Config.get("Files", "quantispike_file")
pseudogene_file = Config.get("Files", "pseudogene_file")

# Hardcoded references/annotation (could be moved to config)
Ref_genome = base_dir + '/Ref/GRCh37/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa'
ref_target_wp = base_dir + '/Ref/' + Assay + '_amplicons_with_primer_length.fasta'

target_bed = base_dir + '/Ref/' + Assay + '_amplicons_with_primer_length.bed'

AmpPriInfo = pd.read_csv(amplicon_file)['#AMPLICONID'].to_list()
COSMIC_DB = base_dir + '/Ref/annotationFiles/CosmicCodingMuts_v97.normalized.uniq.vcf.gz'
ClinVar_DB = base_dir + '/Ref/annotationFiles/clinvar_20211204.vcf.gz'
Common_DB = base_dir + '/Ref/annotationFiles/common_all_20180423.vcf.gz'

# -------------------------------------
# 2) Check individual run fastq name pattern
# -------------------------------------

ProjectID = NextSeqDir.split("/")[-2]
RunID = NextSeqDir.split("/")[-1].split('-')[-1] + '_preprod_v1'
FlowCellID = NextSeqDir.split("/")[-1].split('-')[-1].split('_')[-1]

fastqfiles = glob.glob(NextSeqDir + "/*fastq.gz")
fastqnames = []
for f in fastqfiles:
    fn = f.split('/')[-1].split(".")[0]
    fastqnames.append(fn)

# Exclude Undetermined
fastqnames = [item for item in fastqnames if 'Undetermined' not in item]

# Pattern definitions
pattern = r"\dP(?![C])\d*"  # e.g., matching patient IDs
NTCS_pattern = r"NTCS\d+"
NTCSIDX_pattern = r"NTCSIDX\w"

PatIDs = []         # patient IDs
PatFilenames = []
for l in [elem for elem in fastqnames if re.search(pattern, elem)]:
    i = l.split('_')[2]
    PatIDs.append(i)
    PatFilenames.append(l)

PcIDs = []          # positive control IDs
PcFilenames = []
for l in [elem for elem in fastqnames if "PC" in elem]:
    i = l.split('_')[2]
    PcIDs.append(i)
    PcFilenames.append(l)

NTCSIDs = []        # negative control IDs (?)
NTCSFilenames = []
for l in [elem for elem in fastqnames if re.search(NTCS_pattern, elem)]:
    i = l.split('_')[2]
    NTCSIDs.append(i)
    NTCSFilenames.append(l)

NTCSIDXIDs = []
NTCSIDXFilenames = []
for l in [elem for elem in fastqnames if re.search(NTCSIDX_pattern, elem)]:
    i = l.split('_')[2]
    NTCSIDXIDs.append(i)
    NTCSIDXFilenames.append(l)  # <-- FIXED typo here

# Consolidate sample IDs
SampleID = list(set(PatIDs + PcIDs + NTCSIDs + NTCSIDXIDs))
SampleFilenames = PatFilenames + PcFilenames + NTCSFilenames + NTCSIDXFilenames

# For clarity, define sets:
PP = list(set(PatIDs))  # e.g., just patient IDs
PPFilename = list(set(PatFilenames))

PositiveID = list(set(PcIDs))  # e.g., positive controls
PositiveFilename = list(set(PcFilenames))

PwP = list(set(PcIDs + PatIDs))  # patients + positive
PwPFilename = list(set(PatFilenames + PcFilenames))

#print("Samples detected:", SampleFilenames)

# -------------------------------------
# 3) Rules
# -------------------------------------

rule all:
    input:
        # Step 1 outputs
        expand(RunID + '_step001/{sample}.fastq.gz', sample=SampleFilenames),
        RunID + "_step001/" + FlowCellID + "_fastqStats.tab", 
        # Step 2 (merged alignment) only for patients + positives
        expand(RunID + '_step002/' + FlowCellID + '_{sampleid}_alignments.bam', sampleid=PwP),
        # Step 3 (individual alignment) for all
        expand(RunID + '_step003/{sample}.bam.bai', sample=SampleFilenames),
        expand(RunID + '_step003_1/{sample}.bam.bai', sample=SampleFilenames),
        # Step 4
        expand(RunID + '_step004/{sample}.vcf', sample=SampleFilenames),
        expand(RunID + '_step004/{sample}.fastq.gz', sample=SampleFilenames),
        # Step 5
        RunID + '_step005/consensus_merged_fastq_stats.tab',
        # Step 6
        expand(RunID + '_step006_vardict_merged/{sampleid}_merged_filtered_VarDict.vcf', sampleid=PP),
        expand(RunID + '_step004_bwa_vardict_m/{sampleid}_merged_filtered.tsv', sampleid=PwP),
        # Step 7
        expand(RunID + '_step007_freebayes/{sample}.tsv', sample=SampleFilenames),
        expand(RunID + '_step008_freebayes_merged/{sampleid}_merged_filtered_freebayes.vcf', sampleid=PP),
        # Step 9
        expand(RunID + '_step009_lofreq/{sample}.tsv', sample=SampleFilenames),
        expand(RunID + '_step010_lofreq_merged/{sampleid}_merged_filtered_lofreq.vcf', sampleid=PP),
       
   # run:
   #     print("All pipeline steps are considered complete.")


# ------------------
# STEP 1: Generate consensus fastq
# 
# ------------------
rule generate_consensus_seq:
    input:
        NextSeqDir + '/{sample}.fastq.gz',
    output:
        RunID + '_step001/{sample}.fastq.gz',
    threads: 4
    params:
        inputdir = NextSeqDir,
        outputdir = RunID + '_step001',
        configfile = config_file,
    shell:
        """
        python scripts/cfDNA_generateConsensus_prod_v1.py -i {params.inputdir} -o {params.outputdir} -c {params.configfile}
        """

rule gather_stats:
    input:
        expand(
            RunID + "_step001/{sample}.fastq.gz",
            sample=SampleFilenames
        )
    output:
        RunID + "_step001/" + FlowCellID + "_fastqStats.tab"
    params:
        # Possibly pass in the entire directory
        outdir = RunID + "_step001",
        configfile = config_file
    shell:
        """
        python scripts/generateConsensus_v3.py \
            --stats-only \
            -o {params.outdir} \
            -c {params.configfile}
        """

def get_fastqs(sampleid):
    dir_ = RunID + '_step001'
    fqs = glob.glob(f"{dir_}/*{sampleid}*.fastq.gz")
    return fqs

# ------------------
# STEP 2: Amplicon alignment (individual subsample fastq, for checking purpose)
# ------------------
rule amplicons_alignment:
    input:
        lambda wildcards: get_fastqs(wildcards.sampleid),
    output:
        bam = RunID + '_step002/' + FlowCellID + '_{sampleid}_alignments.bam',
        bai = RunID + '_step002/' + FlowCellID + '_{sampleid}_alignments.bam.bai',
    threads: 2
    params:
        config = config_file,
        assay_prefix = 'BC_RUO1',
        inputdir = RunID + '_step001',
        outputdir = RunID + '_step002',
    shell:
        """
        module purge
        module load Biotools/2024.1.1
        python conseq2bam_v1.py -id {params.inputdir} -od {params.outputdir} -c {params.config} -a {params.assay_prefix}
        module unload Biotools/2024.1.1
        """

# ------------------
# STEP 3: Individual amplicons alignment

# ------------------
rule individual_amplicons_alignment:
    input:
        fqs = RunID + '_step001/{sample}.fastq.gz',
    output:
        bam = RunID + '_step003/{sample}.bam',
        bai = RunID + '_step003/{sample}.bam.bai',
    threads: 2
    params:
        config = config_file,
        assay_prefix = 'BC_RUO1',
        inputdir = RunID + '_step001',
        outputdir = RunID + '_step003',
        sample = lambda wildcards: wildcards.sample,
    shell:
        """
        python SafeSeq_conseq2bam_v3.py -id {params.inputdir} -od {params.outputdir} -c {params.config} -a {params.assay_prefix} -s {params.sample}
        """

rule picard_add_readgroups:
    input:
        bam = RunID + '_step003/{sample}.bam',
        bai = RunID + '_step003/{sample}.bam.bai',
    output:
        bam = RunID + '_step003_1/{sample}.bam',
        bai = RunID + '_step003_1/{sample}.bam.bai',
    params:
        extra = lambda wildcards: f'RGID={wildcards.sample} RGLB=BC_P2v1 RGPL=ILLUMINA RGSM=20 RGPU=unit VALIDATION_STRINGENCY=SILENT',
    shell:
        """
        picard AddOrReplaceReadGroups I={input.bam} O={output.bam} {params.extra}
        samtools index {output.bam}
        """

# ------------------
# STEP 4: VarDict on individual
# ------------------
rule individual_vardict_bwa:
    input:
        bam = RunID + '_step003/{sample}.bam',
    output:
        txt = RunID + '_step004/{sample}.txt',
        vcf = RunID + '_step004/{sample}.vcf',
        vcf2 = RunID + '_step004/{sample}_2.vcf',
        vcf3 = RunID + '_step004/{sample}_3.vcf',
        vcf4 = RunID + '_step004/{sample}_4.vcf',
        tsv = RunID + '_step004/{sample}.tsv',
        fq = RunID + '_step004/{sample}.fastq',
        fqz = RunID + '_step004/{sample}.fastq.gz',
    params:
        name = lambda wildcards: wildcards.sample,
        ref = Ref_genome,
        bed = target_bed,
    threads: 2
    shell:
        """
        module load vardict-java/1.5.1
        vardict-java -th {threads} -G {params.ref} -N {params.name} -f 0.0001 -b {input.bam} -I 10 -v -z -c 1 -S 2 -E 3 -q 20 -p -m 40 {params.bed} \
            | teststrandbias.R > {output.txt}
        module purge
        module load BioTools/2023.01
        if [[ -s {output.txt} ]]; then
            cat {output.txt} | perl ./var2vcf_valid.pl > {output.vcf}
        else
            echo "No variants found for {params.name}."
            touch {output.vcf}
        fi
        vcfallelicprimitives -kg {output.vcf} > {output.vcf2}
        vcfleftalign -r {params.ref} -w 30 {output.vcf} > {output.vcf3}
        vcfallelicprimitives -kg {output.vcf} | vcfleftalign -r {params.ref} -w 30 > {output.vcf4}
        vcf2tsv {output.vcf} > {output.tsv}
        samtools bam2fq {input.bam} > {output.fq}
        gzip -k {output.fq}
        """

# ------------------
# STEP 5: Fastq stat
# ------------------
rule fastq_stat_consensus_merged:
    input:
        fq = expand(RunID + '_step001/{sample}.fastq.gz', sample=SampleFilenames)
    output:
        csv = RunID + '_step005/consensus_merged_fastq_stats.tab',
    params:
        pf = pseudogene_file,
        qf = quantispike_file
    shell:
        """
        python consens_fastq_stat_v5.py -i {input.fq} -o {output.csv} -p {params.pf} -q {params.qf}
        """

# ------------------
# STEP 6: Merge + filter VarDict VCF
# ------------------
rule VarDict_merge_filter_vcf:
    input:
        tab = RunID + '_step005/consensus_merged_fastq_stats.tab',
    output:
        vcf = RunID + '_step006_vardict_merged/{sampleid}_merged_filtered_VarDict.vcf',
    params:
        sample = lambda wildcards: wildcards.sampleid,
        dir = RunID + '_step004',
    shell:
        """
        python merge_filter_calMM_v4.py -s {params.sample} -d {params.dir} -t {input.tab} -o {output}
        """


rule freebayes:
    input:
        bam = RunID + '_step003/{sample}.bam',
    output:
        vcf = RunID + '_step007_freebayes/{sample}.vcf',
        vcf2 = RunID + '_step007_freebayes/{sample}_2.vcf',
        vcf3 = RunID + '_step007_freebayes/{sample}_3.vcf',
        vcf4 = RunID + '_step007_freebayes/{sample}_4.vcf',
        tsv = RunID + '_step007_freebayes/{sample}.tsv',
    params:
        ref = Ref_genome,
        bed = target_bed,
    shell:
        """
        module load freebayes/1.3.6
        freebayes -F 0.0001 -C 2 -q 20 --max-complex-gap 15 -O --strict-vcf \
                  --read-mismatch-limit 15 -f {params.ref} -t {params.bed} {input.bam} > {output.vcf}

        vcfallelicprimitives -kg {output.vcf} > {output.vcf2}
        vcfleftalign -r {params.ref} -w 30 {output.vcf} > {output.vcf3}
        vcfallelicprimitives -kg {output.vcf} | vcfleftalign -r {params.ref} -w 30 > {output.vcf4}
        vcf2tsv {output.vcf} > {output.tsv}
        """

rule merge_filter_freebayes_vcf:
    input:
        tab = RunID + '_step005/consensus_merged_fastq_stats.tab',
    output:
        vcf = RunID + '_step008_freebayes_merged/{sampleid}_merged_filtered_freebayes.vcf',
    params:
        sample = lambda wildcards: wildcards.sampleid,
        dir = RunID + '_step007_freebayes',
    shell:
        """
        python merge_filter_calMM_freebayes_v2.py -s {params.sample} -d {params.dir} -t {input.tab} -o {output}
        """


rule lofreq:
    input:
        bam = RunID + '_step003/{sample}.bam',
    output:
        bam = RunID + '_step009_lofreq/{sample}_lofreq_indelqual.bam',
        vcf = RunID + '_step009_lofreq/{sample}.vcf',
        vcf2 = RunID + '_step009_lofreq/{sample}_2.vcf',
        vcf3 = RunID + '_step009_lofreq/{sample}_3.vcf',
        vcf4 = RunID + '_step009_lofreq/{sample}_4.vcf',
        tsv = RunID + '_step009_lofreq/{sample}.tsv',
    params:
        ref = Ref_genome,
        bed = target_bed,
    shell:
        """
        module load freebayes/1.3.6
        lofreq indelqual --dindel {input.bam} -o {output.bam} -f {params.ref}
        lofreq call -l {params.bed} \
                    --call-indels \
                    --max-depth 400000 \
                    --force-overwrite \
                    -f {params.ref} \
                    --no-default-filter \
                    -o {output.vcf} {output.bam}
        module purge
        module load BioTools/2023.01

        vcfallelicprimitives -kg {output.vcf} > {output.vcf2}
        vcfleftalign -r {params.ref} -w 30 {output.vcf} > {output.vcf3}
        vcfallelicprimitives -kg {output.vcf} | vcfleftalign -r {params.ref} -w 30 > {output.vcf4}
        vcf2tsv {output.vcf} > {output.tsv}
        """

rule merge_filter_lofreq_vcf:
    input:
        tab = RunID + '_step005/consensus_merged_fastq_stats.tab',
    output:
        vcf = RunID + '_step010_lofreq_merged/{sampleid}_merged_filtered_lofreq.vcf',
    params:
        sample = lambda wildcards: wildcards.sampleid,
        dir = RunID + '_step009_lofreq',
    shell:
        """
        python merge_filter_calMM_lofreq_v2.py -s {params.sample} -d {params.dir} -t {input.tab} -o {output}
        """




