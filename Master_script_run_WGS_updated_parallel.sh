#!/bin/bash
#SBATCH --job-name WGSBEM
#SBATCH --cpus-per-task=4
#SBATCH --mem=12gb
#SBATCH --time=6:00:00
#SBATCH --ntasks-per-node=2
#SBATCH --nodes=1



# Before you begin!. 
# Things to understand before running.
# I have created 3 conda environments to deal with dependencies and incompatibilities. 
# BEMWGS, mitobim and samtools1_1_6
# A have saved and exported the environments so you can easily install them. To prepare these for your computer just activate conda and use the following commands 
# 'conda create -n BEMWGS' (create the environment) 'conda activate BEMWGS' (activate the environment) and then  'conda env update --file BEMWGSenv.yaml' and everything will save properly. The only other requirement is that you go through and change the conda environment path to your path
# e.g., conda activate /scratch1/odw014/Conda/install/envs/BEMWGS is mine, so change yours to where conda installs (usually your home directory or project directory), see the
# larger conda instructions i sent to everyone if you are unsure. 
# I have left the old module loads as comments below that you can activate if you don't want to do conda but there are some additional limitations (e.g., you can't use samtools consensus, and the analysis 
# will be slower)


# Change the CPUs per task and tasks per node as you feel resources are avaliable. total CPUs used = cpus-per-task * ntasks-per-node 
# Xmx is the java call for memory used e.g., Xmx24G is allocate 24gb of memory per run (which means if in parallel it is 24* number of tasks). Consider limitations in the memory of your system
# Don't forget to download your reference genome you want (HeHo from DNAzoo is what I used to test this) 
# Don't forget to update the file paths for reads and for the genome. THey are just below the modules and DO_XYZ commands
# Also don't forget to fix the file paths for the mitocondrial analysis (near the bottom of the document)

# Warning, If you don't use conda and use the modules, let me know if the parallel doesn't work. It can be very finicky with issues when embedded as a module with other modules 
# I would need to potentially try and change some code to get it to work (or maybe not, I can't really tell now because I am using conda so just let me know how it goes).


# The whole GATK analysis uses multiple iterations of genotype calling to refine accuracy and base calls and so the pipeline will need to be repeated. 
# There are 3 different steps that help with this
# DO_PRELIM_CALLS is to say whether this is the first run or a repeat run. set to true for first run. 
# DO_BASE_RECALIBRATION and DO_APPLY_BQSR are required for the repeat runs to optimise the calls for the program. set to true when you have already completed at least 1 full run with the 
# prelim calls. 
# This pipe does not keep track of how many times you have run the path and refined iterations you will need to do that! It can only help with the first (preliminary calls) and second iterations
# Double check the vcf file paths and change to your file paths as needed. if you want to keep every iteration for the final vcfs, change the name of the vcfs to save/ copy them over to somewhere else
# after each run as the vcf files are overwritten with the new updated vcf here. 


# The mitochondrial extraction part is only partly finished I didn't get the chance to properly research or build full code for it sorry. Let me know if you want to give it a go yourself or 
# I will just plug away in the background on it for a bit.
# The fastp section works, and the mitobim section seems to work. I tested everything on a tiny subset of reads and the program ran but returned no hits. i think this is because of
# there actually being none in the subset I was working on so it should work for the full data. Still, you may need to play around with the settings as they may be too strict etc. 
# you can always look up the command with 'MITObim.pl -h'




eval "$(conda shell.bash hook)"

conda activate /scratch1/odw014/Conda/install/envs/BEMWGS



# module list WGS

#module load fastqc/0.11.9
#module load bcl2fastq-gcc/2.20.0 
#module load stacks-gcc7/2.41
#module load bwa-gcc/0.7.17 #(v0.7.17-r1188)
#module load samtools-gcc/1.6
#module load parallel/20170722
#module load java/1.8.0_66
#module load gatk/4.0.11.0
#module load qualimap/2.1.3
#module load vcftools-gcc/0.1.13
#module load picard/2.2.2

# module list mitochondrial sequencing

#module load mira/4.0.2
#module load fastp/0.20.1
#module load mitobim/20160526
############### Mitobim is the issue. Everything else installed effectively

#module load python-gcc7/3.6.8

#workpaths and data
READSPATH=/scratch1/odw014/BEMWGS/testreads
WORKPATH=/scratch1/odw014/BEMWGS
#Put genome here# 
GENOMEPATH=/scratch1/odw014/BEMWGS/HeHoGenome

#
FASTQall=(`ls ${READSPATH}/*.fastq.gz`)
FASTQR1=(`ls ${READSPATH}/*R1*.fastq.gz`)
FASTQR2=(`ls ${READSPATH}/*R2*.fastq.gz`)

# SLURM_CPUS_ON_NODE should be kept as a variable in the sbatch and will read it from what you have written up the top.
core=${SLURM_CPUS_PER_TASK} # CPUs put on

# analyses to break down into chunks. Keep false until needed


DO_CONVERT_BCL=false  # Isn't needed				#Step 0
DO_DEMULTIPLEX=false # Don't think are needed			#Step 1
DO_BARCODES=false    # Don't think are needed			#Step 1
DO_FASTQC=false							#Step 2
DO_FASTQ_TO_BAM=false
DO_READGROUP_DATA_GEN=false				#Step 3 A
DO_ADD_READGROUP_DATA=false					#Step 4
DO_MARK_ADAPTORS=false						#Step 5 
DO_BAM_TO_FASTQ=false						#Step 6
DO_GENOME_IDX=false						#Step 7 A-C
DO_BWA=false							#Step 8
DO_MERGE_BAM_ALIGNMENT=false					#Step 9
DO_SORT_BAM=false						#Step 10
DO_QUALIMAP=false						#Step 11
DO_FLAG_PCRDUPS=false						#Step 12
DO_BUILD_BAMINDEX=false						#Step 13
DO_BASE_RECALIBRATION=false					#Step 14
DO_APPLY_BQSR=false						#Step 15
DO_SAMTOOLS_CONSENSUS=false                                      #Step NA. As I mentioned with Kat back in September, this is just a quick and easy way of extracting contiguous genomes and will create a single(in the same segements as your reference) genome
DO_PRELIM_CALLS=false						#Step X. For differentiating between first round recalibartion and final recalibarations. True means running first time. False means running second or later times for haplotype caller
DO_HAPLOTYPE_CALLER=false					#Step 16B
DO_COMBINEGVCFS=false						#Step 17
DO_GENOTYPEGVCFS=false						#Step 18					
DO_VCFTOOLS_SUMMARY=false					#Step 19

# Mitochondrial analysis
# Note the mitochondrial analysis is set up to retrieve the mitochondria of each individual for further analysis
# If you want to sequence the mitochondria and for example publish a complete mitochondria, code to combine inds
# of the same species and create a mitochondrial genome with no gaps would likely be better.

DO_MITOCHONDRIAL_ANALYSIS=true					# Mitochondrial analysis prep Step to create directories
DO_TRIM_MITO=true						#Step M1
DO_MITOBIM_RUN=true						#Step M2
DO_MITO_BIM_COLLATE_RESULT=false				#Step M3


# 0 DO_CONVERT_BC
# Not necessary for novaseq data returned

# 1 Demultiplex, remove barcodes. 
# Not necessary for novaseq data returned


#2 FASTQC

# FASTQC is just a simple program which looks at the raw data from a sequencing run and provides summary statistics on the reads including the ratio of Cs Gs As and Ts
# The presence of adaptor contamination, and the overall quality of the sequence depending on its base pair position (1-150 position as the sequencing done was in 150bp chunks)
# For the command here the FASTQ is a list of all the R1 and R2 files in the study and was called in the starting code above with #FASTQ#=#(`ls $WORKPATH/reads/*/*.fastq.gz`)#
# ls is shorthand for list
# -t is the number of CPU cores to use in fastqc and -o is the output path
# Now updated to parallel

# Relevant variables
#FASTQR1
#FASTQR2


PATH_FASTQC=${WORKPATH}/02_FASTQC

if $DO_FASTQC
then



	if [ ! -d $PATH_FASTQC ]
	then
		echo " $(date) mkdir $PATH_FASTQC "
		mkdir $PATH_FASTQC
	fi


echo " Starting FastQC "


# Two parallel runs to analyse first the forward reads and then the reverse reads (fastQC doesn't take in the pairs well together)
parallel --link -j ${SLURM_NTASKS_PER_NODE} 'r1={1}; r2={2}; namefile=$(basename $r1);'\
		'namefile=("${namefile/R*/}");'\
		'echo "Reads $r1 ";'\
		'echo "Sample name: $namefile ";'\
		'fastqc ${r1} ${r2} -t '$core' -o '$PATH_FASTQC' -f fastq ;'\
	::: ${FASTQR1[@]} \
	::: ${FASTQR2[@]}


echo " Finishing FastQC "

fi
#Inspect results to see how they look. This can tell the general success/failure of a run

# Step 3
# We now want to convert the fastq files into bam files 
# gatk uses java so --java-options -Xmx24G indicates the amount of ram to let java use to analyse the sequences. if the run fails you can increase the amount of ram 
# here -Xmx24G indicates 24GB is devoted but some of the later steps require up to 60GB

# Relevant variables
#FASTQR1
#FASTQR2

PATH_BAM_UNALIGNED=${WORKPATH}/03_BAM_UNALIGNED

if $DO_FASTQ_TO_BAM
then



	if [ ! -d $PATH_BAM_UNALIGNED ]
	then
		echo " $(date) mkdir $PATH_BAM_UNALIGNED "
		mkdir $PATH_BAM_UNALIGNED
	fi

    
echo " Starting fastqtobam "

parallel --link -j ${SLURM_NTASKS_PER_NODE} 'r1={1}; r2={2}; namefile=$(basename $r1);'\
		'namefile=("${namefile/_R*/}");'\
		'echo "Reads $r1 $r2 ";'\
		'echo "Sample name: $namefile ";'\
		'gatk --java-options -Xmx32G FastqToSam' \
		'-F1 ${r1}' \
		'-F2 ${r2}' \
		'-O '$PATH_BAM_UNALIGNED'/${namefile}_unaligned_read_pairs.bam' \
		'-SM ${namefile} --SORT_ORDER queryname' \
		'--MAX_RECORDS_IN_RAM 200000' \
		'2> '$PATH_BAM_UNALIGNED'/${namefile}bamconvert.log;'\
	::: ${FASTQR1[@]} \
	::: ${FASTQR2[@]}

echo " Finished fastqtobam "

fi


#step 3A. generate read group data
# Unneccesary for this data. 

# Step 4. Add new read group data

# In step 4 we are creating more informative information tags to go alongside our sequencing reads. 
# another way we can create informative lists is using the grep function in unix. this works similar to grepl in R
# the | symbol means to pipe into the next command so here we have a list everything command and we pipe it into a grep and export only the bam files into a text
# then we again run the read line for each file. 
# we can extract the name again using sub commands which have the format sub(/thing to remove/, "thing to replace it with") using awk we can save it as a new variable
# for sub command  .* is an expression which means anything in that direction so .*pairs_ means remove anything to the left of pairs as well. 

PATH_BAM_UNALIGNED=${WORKPATH}/03_BAM_UNALIGNED

if $DO_ADD_READGROUP_DATA
then

	if [ ! -d $PATH_BAM_UNALIGNED ]
	then
		echo " $(date) mkdir $PATH_BAM_UNALIGNED "
		mkdir $PATH_BAM_UNALIGNED
	fi


echo " Starting AddOrReplaceReadGroups"

Bams=(`ls ${PATH_BAM_UNALIGNED}/*pairs.bam`)

parallel --link -j ${SLURM_NTASKS_PER_NODE} 'r1={1}; namefile=$(basename $r1);'\
		'namefile=("${namefile/.bam/}");'\
		'echo "Reads $r1 ";'\
		'echo "Sample name: $namefile ";'\
		'gatk --java-options -Xmx24G AddOrReplaceReadGroups' \
		'-I ${r1}' \
		'-O '$PATH_BAM_UNALIGNED'/${namefile}_RGup.bam' \
		'-RGID 4' \
		'-RGLB lib1' \
		'-RGPL ILLUMINA' \
		'-RGPU unit1' \
		'-RGSM ${namefile}' \
		'2> '$PATH_BAM_UNALIGNED'/${namefile}bamRG.log;'\
	::: ${Bams[@]}

echo " Finished AddOrReplaceReadGroups"
fi


# Step 5. Mark adaptors #####################
# This step highlights which sequences are likely to be leftover illumina adaptors. It does not remove them but provides information on how many there are and importantly tags
# them as low quality for downstream removal later on

# So before I've given two examples of how to make loops (one with just text files and sed, and one incorporating ls and sed)
# A better way is one that doesn't produce any extra txt files which I have done in this next step
# the declare command lets you put a list directly into the HPCs working memory but won't save it anywhere once the script finishes. 
# here declare -a means we are declaring the list to be inputted as an array which makes for easy access for what we want. 
# we just need to give the array a name (bamlistRGUP) and tell declare how we are populating the array ( this is the ls function)
# so we just list all bam files and that becomes the array. * just means any number of any character either side
# We can then call the array with a simple for loop
# for bam in ${bamlistRGUP[@]}; do
# with the @ keeping track of what index in the array it's up to in the loop (what bam file its up to)
# All script from here on out are done as simple declare arrays and the previous steps can be written up in arrays to. 
# If you want to try. Have a go at trying to change steps 2-4 to use arrays instead of text files now you know how they work. don't forget to save it as a new document first though 
# in case something else gets accidentally deleted :) 


PATH_BAM_UNALIGNED=${WORKPATH}/03_BAM_UNALIGNED				
PATH_TEMPFILES=$WORKPATH/TEMPFILES

if $DO_MARK_ADAPTORS
then


	if [ ! -d $PATH_BAM_UNALIGNED ]
	then
		echo " $(date) mkdir $PATH_BAM_UNALIGNED"
		mkdir $PATH_BAM_UNALIGNED
	fi
	if [ ! -d $PATH_TEMPFILES ]
	then
		echo " $(date) mkdir $PATH_TEMPFILES "
		mkdir $PATH_TEMPFILES
	fi

echo " Starting MarkIlluminaAdapters"

Bams2=(`ls ${PATH_BAM_UNALIGNED}/*RGup.bam`)

parallel --link -j ${SLURM_NTASKS_PER_NODE} 'r1={1}; namefile=$(basename $r1);'\
		'namefile=("${namefile/_RGup.bam/}");'\
		'echo "Reads $r1 ";'\
		'echo "Sample name: $namefile ";'\
		'gatk --java-options -Xmx28G MarkIlluminaAdapters' \
		'-I ${r1}' \
		'-O '$PATH_BAM_UNALIGNED'/${namefile}_adap.bam' \
		'-M '$PATH_BAM_UNALIGNED'/${namefile}_markilluminaadapters_metrics.txt' \
		'2> '$PATH_BAM_UNALIGNED'/${namefile}bam_adapters_markedRG.log;'\
	::: ${Bams2[@]}

echo " Finished MarkIlluminaAdapters"

fi



#Step 6 bam to fastq
# This converts the bam file to a fastq so it can be aligned to the heho genome
# interleave is because the data is paired end (R1 and R2) but combined into one file which we did earlier. 
# Clipping attributes can be changed but read up about the filters and things in this step if you want to change some things

PATH_FASTQ_UPDATED=$WORKPATH/06_FASTQ_UPDATED	

if $DO_BAM_TO_FASTQ
then

			
	if [ ! -d $PATH_FASTQ_UPDATED ]
	then
		echo " $(date) mkdir $PATH_FASTQ_UPDATED "
		mkdir $PATH_FASTQ_UPDATED
	fi

echo " Starting SamToFastq"

bamlistadap=(`ls ${PATH_BAM_UNALIGNED}/*dap.bam`)

parallel --link -j ${SLURM_NTASKS_PER_NODE} 'r1={1}; namefile=$(basename $r1);'\
		'namefile=("${namefile/_adap.bam/}");'\
		'echo "Reads $r1 ";'\
		'echo "Sample name: $namefile ";'\
		'gatk --java-options -Xmx24G SamToFastq' \
		'-I ${r1}' \
		'--FASTQ '$PATH_FASTQ_UPDATED'/${namefile}_interleaved.fq' \
		'--CLIPPING_ATTRIBUTE XT' \
		'--INTERLEAVE true' \
		'--CLIPPING_ACTION 2' \
		'--TMP_DIR '$PATH_TEMPFILES'' \
		'2> '$PATH_FASTQ_UPDATED'/${namefile}samtofastq.log;'\
	::: ${bamlistadap[@]}

echo " Finishing SamToFastq"

fi


# Step 7 build genome index 

# just need to build up the reference genome into something more easily accessible (download the HiC genome from DNA zoo as a fasta first)

if $DO_GENOME_IDX
then

bwa index ${GENOMEPATH}/HeHo_1.0_HiC.fasta

gatk --java-options -Xmx20G CreateSequenceDictionary -R ${GENOMEPATH}/HeHo_1.0_HiC.fasta

samtools faidx ${GENOMEPATH}/HeHo_1.0_HiC.fasta

fi




#Step 8 bwa alignment
# This is where the sequences are aligned to the reference genome
# there are a bunch of different steps to control how stringent/lax this allignment is. if you don't get as much alignment as you want you can make them more lax
# see the bwa website for more details http://bio-bwa.sourceforge.net/bwa.shtml
# we want bwa mem as it is the most comprehensive for longer read sequence data (<70bp)

PATH_BWA=$WORKPATH/08_BWA_aligned
if $DO_BWA
then		
		
	if [ ! -d $PATH_BWA ]
	then
		echo " $(date) mkdir $PATH_BWA "
		mkdir $PATH_BWA
	fi

echo " Starting bwa mem "


fqlist=(`ls ${PATH_FASTQ_UPDATED}/*interleaved*`)

parallel --link -j ${SLURM_NTASKS_PER_NODE} 'r1={1}; namefile=$(basename $r1);'\
		'namefile=("${namefile/_interleaved.fq/}");'\
		'echo "Reads $r1 ";'\
		'echo "Sample name: $namefile ";'\
		'bwa mem -t '$SLURM_CPUS_PER_TASK' -p -M '$GENOMEPATH'/HeHo_1.0_HiC.fasta ${r1} > '$PATH_BWA'/${namefile}_interleaved.sam;'\
	::: ${fqlist[@]}
	

echo " Finished bwa mem "


fi




#Step 9 Merge raw bam with alignment files.
# After alignment to the genome this step copies over all the information from the aligned sequences and gives it back to the raw bam from step 6
# It then uses that information plus the read pairs to correct the whole set of sequences and chromosome maps to make sure all read pairs are actually next to each other
# It also allows for reads to be removed based on set filters. This preseves all sequence data and leads to overall higher accuracy in downstream analysis

PATH_BAM_MERGED=$WORKPATH/09_BAM_MERGED		
if $DO_MERGE_BAM_ALIGNMENT
then		
PATH_BAM_MERGED=$WORKPATH/09_BAM_MERGED		
	if [ ! -d $PATH_BAM_MERGED ]
	then
		echo " $(date) mkdir $PATH_BAM_MERGED "
		mkdir $PATH_BAM_MERGED
	fi



echo " Starting MergeBamAlignment"


samlist=(`ls ${PATH_BWA}/*interleaved.sam`)
Bamsunaligned=(`ls ${PATH_BAM_UNALIGNED}/*RGup.bam`)

parallel --link -j ${SLURM_NTASKS_PER_NODE} 'r1={1}; r2={2}; namefile=$(basename $r1);'\
		'namefile=("${namefile/_interleaved*/}");'\
		'echo "Sample name: $namefile ";'\
		'gatk --java-options -Xmx20G MergeBamAlignment' \
		'-R '$GENOMEPATH'/HeHo_1.0_HiC.fasta' \
		'--UNMAPPED_BAM ${r2}' \
		'--ALIGNED_BAM ${r1}' \
		'-O '$PATH_BAM_MERGED'/${namefile}_aligned_merged.bam' \
		'--CREATE_INDEX true' \
		'--ADD_MATE_CIGAR true' \
		'--CLIP_ADAPTERS' \
		'--CLIP_OVERLAPPING_READS true' \
		'--INCLUDE_SECONDARY_ALIGNMENTS true' \
		'--MAX_INSERTIONS_OR_DELETIONS -1' \
		'--PRIMARY_ALIGNMENT_STRATEGY MostDistant' \
		'--ATTRIBUTES_TO_RETAIN XS' \
		'--TMP_DIR '$WORKPATH'/TEMPFILES' \
		'2> '$PATH_BAM_MERGED'/${namefile}bamconvert.log;'\
	::: ${samlist[@]} \
	::: ${Bamsunaligned[@]}

echo " Finished MergeBamAlignment"

fi

# next we want to sort the bam based on the coordinates of all sequences (base pair position on the genome)
# If there is trouble in this step upp the max records in Ram as it may require more memory

# Step 10. Sort bam
PATH_BAM_MERGED_SORTED=$WORKPATH/10_BAM_MERGED_SORTED

if $DO_SORT_BAM
then

	if [ ! -d $PATH_BAM_MERGED_SORTED ]
	then
		echo " $(date) mkdir $PATH_BAM_MERGED_SORTED "
		mkdir $PATH_BAM_MERGED_SORTED
	fi

echo " Starting SortSam"

bamlist=(`ls ${PATH_BAM_MERGED}/*merged.bam`)

parallel --link -j ${SLURM_NTASKS_PER_NODE} 'r1={1}; namefile=$(basename $r1);'\
		'namefile=("${namefile/_aligned_merged.bam/}");'\
		'echo "Reads $r1 ";'\
		'echo "Sample name: $namefile ";'\
		'gatk --java-options -Xmx24g SortSam' \
		'-I ${r1} -O '$PATH_BAM_MERGED_SORTED'/${namefile}_srt.bam -SO coordinate --MAX_RECORDS_IN_RAM 200000;'\
	::: ${bamlist[@]}


echo " Finished SortSam"


fi


# Step 11 Qualimap 

# Let's generate some summary statistics of what the sequences look like now
# not paralelised because is fairly fast anyway ( I think. can paralelise this as needed to)

	PATH_QUALIMAP=$WORKPATH/11_QUALIMAP
	
if $DO_QUALIMAP
	then

		if [ ! -d $PATH_QUALIMAP ]
		then
			mkdir $PATH_QUALIMAP
		fi
		echo "$(date)	Start qualimap"

bamsrtlist=(`ls ${WORKPATH}/10_BAM_MERGED_SORTED/*srt*`)

for bam in ${bamsrtlist[@]}; do

name=(`awk '{sub(/_srt.bam/, ""); sub(/.*MERGED_SORTED\//, ""); print }'<<< $bam`)

mkdir $PATH_QUALIMAP/${name}

qualimap --java-mem-size=12000M bamqc -bam $bam -outdir ${PATH_QUALIMAP}/${name} -outformat pdf -outfile ${name} -nt $core

done

echo " Finished Qualimap"

fi

# 

#Step 12. Flag PCR duplicates	

BAM_MERGED_FILTERED=$WORKPATH/12_BAM_MERGED_FILTERED
# This step does as it says, it just flags potential PCR duplicates and marks them similarly to how the illumina adaptors were marked

if $DO_FLAG_PCRDUPS

then



	if [ ! -d $BAM_MERGED_FILTERED ]
	then
		echo " $(date) mkdir $BAM_MERGED_FILTERED "
		mkdir $BAM_MERGED_FILTERED
	fi

echo " Starting MarkDuplicates "
bamsrtlist=(`ls ${WORKPATH}/10_BAM_MERGED_SORTED/*srt*`)

parallel --link -j ${SLURM_NTASKS_PER_NODE} 'r1={1}; namefile=$(basename $r1);'\
		'namefile=("${namefile/_unaligned_read_pairs_srt.bam/}");'\
		'echo "Reads $r1 ";'\
		'echo "Sample name: $namefile ";'\
		'gatk --java-options -Xmx24g MarkDuplicates -I ${r1} -O '$BAM_MERGED_FILTERED'/${namefile}_PCRdup.bam -M '$BAM_MERGED_FILTERED'/${namefile}.txt' \
		'2> '$BAM_MERGED_FILTERED'/${namefile}MarkDuplicates.log;'\
	::: ${bamsrtlist[@]}

echo " Finished MarkDuplicates "
fi



# Step 13. build bam index

# This builds index files for each of the bam sample files so they can be read easier 
# not paralelised because is fairly fast anyway ( I think. can paralelise this as needed to)

if $DO_BUILD_BAMINDEX

then

echo " Starting samtools index "

bamlistsrt=(`ls ${BAM_MERGED_FILTERED}/*PCRdup.bam`)

parallel --link -j ${SLURM_NTASKS_PER_NODE} 'r1={1}; namefile=$(basename $r1);'\
		'namefile=("${namefile/_PCRdup.bam/}");'\
		'echo "Reads $r1 ";'\
		'echo "Sample name: $namefile ";'\
		'samtools index -b ${r1}' \
		'2> '$BAM_MERGED_FILTERED'/${namefile}samtoolsidx.log;'\
	::: ${bamlistsrt[@]}


echo " Finished samtools index "

fi

#Step 14 DO_BASE_RECALIBRATION


	# recalibrate base quality score
	# need vcf to compare to (skip if first round, then repeat second half of analysis)
# This is an important step for increasing accuracy but creates the complicated loop in the workflow figure I sent. What this does is recalibrate the quality scores for sequencing
# based on likely SNPs and based on patterns observed over all the sequences. e.g., if the sequence GTGTCG has a really low quality next base in the sequence or has been flagged in the SNPs 
# identified as less reliable, this program will automatically drop the accuracy of the next base in the raw sequences fed to it which will yield overall more accurate results in the final
# data set. We don't have a list of known SNPs though, so we need to run a prelim analysis, get the vcf results of that, come back to this step, insert this vcf results file, recalibrate, generate
# a new SNP results file (down to Step 18 completed) and rinse and repeat 2 or 3 times.

	if $DO_BASE_RECALIBRATION

	then
# Will need to change this to the correct new folder names but this is just the final vcf called genotypes of the prior run of this (if doing multiple for accuracy)
# As of yet not tested (will need to test after step 19 reached, same applies for step 15)

PATH_VCF_SUMMARY=$WORKPATH/19_VCFTOOLS_SUMMARY
vcfname=(${PATH_VCF_SUMMARY}/cohortminers.recode.vcf)


bamlistsrt=(`ls ${WORKPATH}/12_BAM_MERGED_FILTERED/*PCRdup.bam`)

echo " Starting BaseRecalibrator "

echo " $vcfname "

parallel --link -j ${SLURM_NTASKS_PER_NODE} 'r1={1}; namefile=$(basename $r1);'\
		'namefile=("${namefile/_PCRdup.bam/}");'\
		'echo "Reads $r1 ";'\
		'echo "Sample name: $namefile ";'\
		'gatk --java-options -Xmx24g BaseRecalibrator' \
		'--input ${r1}' \
		'--known-sites '$vcfname''\
		'--output '$BAM_MERGED_FILTERED'/${namefile}_recalctable.txt' \
		'--reference '$GENOMEPATH'/HeHo_1.0_HiC.fasta' \
		'2> '$BAM_MERGED_FILTERED'/${namefile}BaseRecalibratorrd1.log;'\
	::: ${bamlistsrt[@]}

echo " Finished BaseRecalibrator "

fi



#Step 15 DO_APPLY_BQSR
# this applies the newly determined quality scores

PATH_BAM_RECALIBRATED=$WORKPATH/15_BAM_RECALIBRATED

if $DO_APPLY_BQSR

then 

if [ ! -d $PATH_BAM_RECALIBRATED ]
	then
		mkdir $PATH_BAM_RECALIBRATED
	fi

echo " Starting ApplyBQSR "

bamlistsrt=(`ls ${WORKPATH}/12_BAM_MERGED_FILTERED/*PCRdup.bam`)

# Put in vcf save location

parallel --link -j ${SLURM_NTASKS_PER_NODE} 'r1={1}; namefile=$(basename $r1);'\
		'namefile=("${namefile/_PCRdup.bam/}");'\
		'echo "Reads $r1 ";'\
		'echo "Sample name: $namefile ";'\
		'gatk --java-options -Xmx24g ApplyBQSR' \
		'-R '$GENOMEPATH'/HeHo_1.0_HiC.fasta' \
		'-I ${r1}' \
		'--bqsr-recal-file '$BAM_MERGED_FILTERED'/${namefile}_recalctable.txt' \
		'-O '$PATH_BAM_RECALIBRATED'/${namefile}_PCRdup.bam' \
		'2> '$BAM_MERGED_FILTERED'/${namefile}_BaseRecalibratorrd2.log;'\
	::: ${bamlistsrt[@]}


echo " Finished ApplyBQSR "

fi
# Not sure yet, but may need to move the .bai files over to, it may create its own though, we will see. 



# Step X (Will filter based on whether step 14 was completed or skipped temporary)
# this just helps identify where in the step 14-Step 18 loop you are up to
if $DO_PRELIM_CALLS

then

BAMDIRECTORY=${WORKPATH}/12_BAM_MERGED_FILTERED

else 

BAMDIRECTORY=${WORKPATH}/15_BAM_RECALIBRATED

fi

echo " $BAMDIRECTORY "



# This step is the most time consuming. It also doesn't have a good parallel. For some reason additional cores/CPUs can be passed from parallel to any command except GATK program 
# commands. This means using parallel I can only parallelise by splitting across individuals but not multi threading within Haplotype caller (--native-pair-hmm-threads setting)
# Because of the RAM limitations, only 4-6 cores will be used at a time which is quite the waste. 

# One alternative is I can split the analysis up into 900ish runs per sample by splitting at the chromosome level in the HeHo genome running each through with one core and stitching
# them all together at the end. This is only recommended if all chromosomes and scaffolds are expected to be in the correct location. Something to discuss on Thurs. The code is 
# all written up for this in the prior analysis pipeline btw. 

#

# Step 16 B

PATH_GENOTYPING=${WORKPATH}/16_GENOTYPING

if $DO_HAPLOTYPE_CALLER
then
	if [ ! -d $PATH_GENOTYPING ]
	then
		mkdir $PATH_GENOTYPING
	fi
	# Generate GVCF files
	# Call germline SNPs and indels via local re-assembly of haplotypes

		echo " $(date)	Start HaplotypeCaller "



#module load fastqc/0.11.9
#module load bcl2fastq-gcc/2.20.0 
#module load stacks-gcc7/2.41
#module load bwa-gcc/0.7.17 #(v0.7.17-r1188)
#module load samtools-gcc/1.6
#module load parallel/20170722
#module load qualimap/2.1.3
#module load vcftools-gcc/0.1.13
#module load picard/2.2.2


bamlistgenotyping=(`ls ${BAMDIRECTORY}/*PCRdup.bam`)

parallel --link -j ${SLURM_NTASKS_PER_NODE} 'r1={1}; namefile=$(basename $r1);'\
		'namefile=("${namefile/_PCRdup.bam/}");'\
		'echo "Reads $r1 ";'\
		'echo "Sample name: $namefile ";'\
		'gatk --java-options -Xmx12g HaplotypeCaller' \
		'--input ${r1}' \
		'--output '$PATH_GENOTYPING'/${namefile}_genotypes.g.vcf.gz' \
		'--reference '$GENOMEPATH'/HeHo_1.0_HiC.fasta' \
		'--emit-ref-confidence GVCF --heterozygosity-stdev 0.015 --native-pair-hmm-threads '$SLURM_CPUS_PER_TASK' --tmp-dir '$WORKPATH'/TEMPFILES' \
		'2> '$PATH_GENOTYPING'/${namefile}_Haplotypecaller.log;'\
	::: ${bamlistgenotyping[@]}

		echo "$(date) Finished HaplotypeCaller "

fi

# STEP NA. samtools consensus. Note, I have it so it only runs off the already edited bams fro step 15 not from 12. 
# Sort the bam file, then index it and then run consensus. -aA tells to output whole genome length even for NAs so you will get something the correct length even if coverage is low
# in areas. -c 0.65 is for calling consensus with 65% of reads required to call a consensus base pair here but can increase as you want for stringency. -d is the minimum read depth
# which if not hit will return N. --het-fract and --ambig work together. If the second most common SNP is above 25% of the reads report the ambig code for the SNP e.g., IUPAC codes
# for either an A/C is here or T/G etc. 
# Show ins and show del are to just show the presence of insertions/deletions. May not b as useful here given these are two different species. 
# Lastly, samtools consensus strips the file of its names so sed is used to 

PATH_SAM_CONSENSUS=${WORKPATH}/16_01_SAMTOOLS_CONSENSUS

if $DO_SAMTOOLS_CONSENSUS
then



	if [ ! -d PATH_SAM_CONSENSUS ]
	then
		echo " $(date) mkdir $PATH_SAM_CONSENSUS "
		mkdir $PATH_SAM_CONSENSUS
	fi

# Change the call and the environment here#
conda activate /scratch1/odw014/Conda/install/envs/samtools1_1_6

bamsrtlist=(`ls ${WORKPATH}/15_BAM_RECALIBRATED/*PCRdup.bam`)

for bam in ${bamsrtlist[@]}; do

name=(`awk '{sub(/_PCRdup.bam/, ""); sub(/.*RECALIBRATED\//, ""); print }'<<< $bam`)

samtools sort $bam -o ${PATH_SAM_CONSENSUS}/${name}_sort.bam
samtools index ${PATH_SAM_CONSENSUS}/${name}_sort.bam && \
samtools consensus -a -A -c 0.7 -d 5 --ambig --het-fract 0.20 --show-ins --show-del ${PATH_SAM_CONSENSUS}/${name}_sort.bam -o ${PATH_SAM_CONSENSUS}/${name}_consensus.fna

done



conda activate /scratch1/odw014/Conda/install/envs/BEMWGS


fi


# Step 17

# Next we combine the GVCFS of each individual to create a population level GVCF



PATH_GENOTYPEGVCFSCOMBINED=${WORKPATH}/17_GENOTYPEGVCFSCOMBINED
if $DO_COMBINEGVCFS 
then
	if [ ! -d $PATH_GENOTYPEGVCFSCOMBINED ]
	then
		mkdir $PATH_GENOTYPEGVCFSCOMBINED
	fi


echo "$(date) Started CombineGVCFs "

ls ${WORKPATH}/16_GENOTYPING/*genotypes.g.vcf.gz > ${PATH_GENOTYPEGVCFSCOMBINED}/gvcfs.list

cd ${WORKPATH}

 gatk --java-options -Xmx40g CombineGVCFs \
   -R ${GENOMEPATH}/HeHo_1.0_HiC.fasta \
	--variant ${PATH_GENOTYPEGVCFSCOMBINED}/gvcfs.list \
	-O ${WORKPATH}/17_GENOTYPEGVCFSCOMBINED/cohort_miners.g.vcf.gz


echo "$(date) Finished CombineGVCFs "

fi

# Step 18
# Finally we call genotypes for the population GVCF 
PATH_GENOTYPEGVCFS=$WORKPATH/18_GENOTYPEGVCFS_FINAL


if $DO_GENOTYPEGVCFS
then
	if [ ! -d $PATH_GENOTYPEGVCFS ]
	then
		mkdir $PATH_GENOTYPEGVCFS
	fi
	# Generate GVCF files
	# Call germline SNPs and indels via local re-assembly of haplotypes

		echo "$(date)	Start GenotypeGVCFs "

echo "$(date) Started GenotypeGVCFs "


gatk --java-options -Xmx48g GenotypeGVCFs \
	-R ${GENOMEPATH}/HeHo_1.0_HiC.fasta \
	-V ${WORKPATH}/17_GENOTYPEGVCFSCOMBINED/cohort_miners.g.vcf.gz \
	-O $PATH_GENOTYPEGVCFS/cohort_miners_genotypescalled.g.vcf.gz

echo "$(date) Finished GenotypeGVCFs "

fi




# Step 19 VCFtools summary

# Summarise the final gvcf data

PATH_VCF_SUMMARY=$WORKPATH/19_VCFTOOLS_SUMMARY
if $DO_VCFTOOLS_SUMMARY
then
	if [ ! -d $PATH_VCF_SUMMARY ]
	then
		mkdir $PATH_VCF_SUMMARY
	fi

echo "$(date) Starting VCFtools "

vcfname=(${PATH_GENOTYPEGVCFS}/cohort_miners_genotypescalled.g.vcf.gz)


indrawvcfs=(`ls ${WORKPATH}/16_GENOTYPING/*_genotypes.g.vcf.gz`)
# pop summary
# can also remove indels or other components at this step (see below comment). This is important to decide on what you want to do with the data afterwards. You will need to research what is
# optimal design for your project.
#http://vcftools.sourceforge.net/man_latest.html
#--remove-indels



vcftools --gzvcf $vcfname --recode --recode-INFO-all --out ${PATH_VCF_SUMMARY}/cohortminers

 gatk IndexFeatureFile \
     -I ${PATH_VCF_SUMMARY}/cohortminers.recode.vcf


echo "$(date) Finished VCFtools "

# next we want individual specific data summaries from the previous analysis.


parallel --link -j ${SLURM_NTASKS_PER_NODE} 'r1={1}; namefile=$(basename $r1);'\
		'namefile=("${namefile/_genotypes.g.vcf.gz/}");'\
		'echo "Reads $r1 ";'\
		'echo "Sample name: $namefile ";'\
		'vcftools --gzvcf ${r1} --recode --recode-INFO-all --out '$PATH_VCF_SUMMARY'/${namefile}' \
		'2> '$PATH_VCF_SUMMARY'/${namefile}_combinedvcf.log;'\
	::: ${indrawvcfs[@]}


fi

# At this point you can then filter your data based on missing data thresholds for inds/loci, HWE, linkage, or MAF/MAC using simple vcftools and setting the correct variable. Or you 
# can import it into R and do it there using vcfR (See my golden perch bioinformatics of vcfR for examples if you need [Kat has it but it's also on my github] ). I also have code there to
# add the relevant information sections and to convert it to a genlight if that's what you want. 






#Mitochondrial analysis

# This uses the raw reads data. 
# Note, I created code up until the completion of Mitobim but the sequence still needs to be circularised and annotated to be a full mitochondrial genome
# Use the DO_MITOCHONDRIAL_ANALYSIS for these steps always so that the cond environment is activate



PATH_M1_READ_MERGING=$WORKPATH/M1_READ_MERGING
PATH_M2_MITOBIM=$WORKPATH/M2_MITOBIM
PATH_M3_MITO_BIM_RESULTS=$WORKPATH/M3_MITOBIM_RESULTS

if $DO_MITOCHONDRIAL_ANALYSIS
then
	if [ ! -d $PATH_M2_MITOBIM ]
	then
		mkdir $PATH_M1_READ_MERGING
		mkdir ${PATH_M1_READ_MERGING}/failed_reads
		mkdir $PATH_M2_MITOBIM
		mkdir $PATH_M3_MITO_BIM_RESULTS
	fi


conda activate /scratch1/odw014/Conda/install/envs/mitobim

fi

# Step M1. Starts using raw data from the beginning

# unlike the main genome sequencing, we need to actually trim the data, not just mark low quality bits. 
# fastp trims and there are many settings which can be changed here for optimal trimming. right now I set it to quality score of 20 minimum and read pairs are identified 
# based on a minimum of 15bp overlap

seed=(${GENOMEPATH}/Helmeted_honeyeater_Mitochondrial_genome.fasta)

if $DO_TRIM_MITO

then

fastqlistR1=(`ls ${READSPATH}/*R1*.fastq.gz`)
fastqlistR2=(`ls ${READSPATH}/*R2*.fastq.gz`)


#location of fasta file used as seed for MITObim assembly
seed=(${GENOMEPATH}/Helmeted_honeyeater_Mitochondrial_genome.fasta)


# Can change these variables a bit to include more data or less. There is not really
# A big issue though. The mitochondria should have more data than everything else so stringent
# or lax there will be ample quantities to reconstruct the mitochondrial genome


parallel --link -j ${SLURM_NTASKS_PER_NODE} 'r1={1}; r2={2}; namefile=$(basename $r1);'\
		'namefile=("${namefile/R*/}");'\
		'echo "Reads $r1 $r2 ";'\
		'echo "Sample name: $namefile ";'\
		'fastp -i ${r1} -I ${r2}' \
		'--merge --merged_out '$PATH_M1_READ_MERGING'/${namefile}trimmed.fastq' \
		'-j '$PATH_M1_READ_MERGING'/fastp_a_miner${namefile}.json -h '$PATH_M1_READ_MERGING'/fastp_a_miner${namefile}.html' \
		'-q 15' \
		'-w ${SLURM_CPUS_PER_TASK}' \
		'-l 75' \
		'--failed_out '$PATH_M1_READ_MERGING'/failed_reads/${nameonly}_failed_reads.fastq' \
		'--low_complexity_filter' \
		'--complexity_threshold 15' \
		'--overrepresentation_analysis' \
		'-5' \
		'--cut_front_window_size 4' \
		'--cut_front_mean_quality 20' \
		'2> '$PATH_M1_READ_MERGING'/${namefile}trim.log;'\
	::: ${fastqlistR1[@]} \
	::: ${fastqlistR2[@]}

fi




#Step M2 
# Mitobim is set to 100 iterations but it took 3-4 per ind in the trial data so I don't think anywhere near that many is needed (it ends once it converges though so it's not an issue) 
# The program doesn't like writing different names for its iterations or overwritting files. This means that new folders need to be created for each ind in each mitbim run (code does this 
# and labels already)
# Note, the command on the HPC for the module mitobim is MITObim_1.8.pl        After version 1.9 they standardised it so that the command is just MITObim.pl so the below works on the conda
# enviroment but you will need to change it if you use the HPC modules

if $DO_MITOBIM_RUN

then


# The program seems to fail through segmentation for some reason. The only way I could fix this is to purge all the environment and reload only mitobim. It worked though!
# I have no idea why though, there was no parallelisation in the script so i don't understand what was driving the segmentation fault.

#module purge
#module load mitobim/20160526


trimmedlist=(`ls ${PATH_M1_READ_MERGING}/*trimmed.fastq`)
seed=(${GENOMEPATH}/Helmeted_honeyeater_Mitochondrial_genome.fasta)

cd ${PATH_M2_MITOBIM}
for fa in ${trimmedlist[@]}; do

nameonly=(`awk '{ sub(/_trimmed.fastq/, ""); sub(/.*MERGING\//, ""); print }'<<< $fa`)


echo " $nameonly "
echo " $fa " 
echo " $seed "

mkdir ${nameonly}
cd ${nameonly}

MITObim.pl -start 1 -end 100 -sample bird${nameonly} -ref honeyeater_mito_genome -readpool $fa --quick ${seed} --clean --NFS_warn_only > ${nameonly}log

cd ..


done

cd $WORKPATH

fi


if $DO_MITO_BIM_COLLATE_RESULT
then

echo " need to do this part still "

fi




# Let me know if you have any questions or issues

