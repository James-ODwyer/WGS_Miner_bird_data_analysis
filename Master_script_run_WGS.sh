#!/bin/bash
#SBATCH 
#SBATCH --cpus-per-task=16
#SBATCH --mem=MaxMemPerNode
#SBATCH --partition=compute
#SBATCH --time=149:00:00
#SBATCH --mail-user=18088076@students.latrobe.edu.au
#SBATCH --ntasks=1

# Running notes. 
# If you are running the genotyping steps, I would recommend changing to long run. Per sample, I estimate it takes 4 days to run. 
# This means for 10 samples it'll take 40 days!



# module list WGS

module load fastqc/0.11.9
module load bcl2fastq-gcc/2.20.0 
module load stacks-gcc7/2.41
module load bwa-gcc/0.7.17 #(v0.7.17-r1188)
module load samtools-gcc/1.6
module load parallel/20170722
module load java/1.8.0_66
module load gatk/4.0.11.0
module load qualimap/2.1.3
module load vcftools-gcc/0.1.13
module load picard/2.2.2

# module list mitochondrial sequencing

module load mira/4.0.2
module load fastp/0.20.1
module load mitobim/20160526
module load python-gcc7/3.6.8

#workpaths and data
WORKPATH=/data/group/murphylab/project/bem/heho_WGS/Test_2_inds
GENOMEPATH=$WORKPATH/heho_genome

FASTQ=(`ls $WORKPATH/reads/*/*.fastq.gz`)
core=16 # CPUs put on

# analyses to break down into chunks. Keep false until needed


DO_CONVERT_BCL=false  # Isn't needed				#Step 0
DO_DEMULTIPLEX=false # Don't think are needed			#Step 1
DO_BARCODES=false    # Don't think are needed			#Step 1
DO_FASTQC=false							#Step 2
DO_FASTQ_TO_BAM=false
DO_READGROUP_DATA_GEN=false     				#Step 3 A
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
DO_SPLT_CHRS=true						#Step 16A. To allow for parallelisation through chromosomes
DO_PRELIM_CALLS=true						#Step X. For differentiating between first round recalibartion and final recalibarations. True means running first time. False means running second or later times for haplotype caller
DO_HAPLOTYPE_CALLER=true					#Step 16B
DO_COMBINEGVCFS=true						#Step 17
DO_GENOTYPEGVCFS=true						#Step 18					
DO_VCFTOOLS_SUMMARY=true					#Step 19

# Mitochondrial analysis
# Note the mitochondrial analysis is set up to retrieve the mitochondria of each individual for further analysis
# If you want to sequence the mitochondria and for example publish a complete mitochondria, code to combine inds
# of the same species and create a mitochondrial genome with no gaps would likely be better.

DO_MITOCHONDRIAL_ANALYSIS=false					# Mitochondrial analysis prep Step to create directories
DO_TRIM_MITO=false						#Step M1
DO_MITOBIM_RUN=false						#Step M2
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
# Using that and the @ symbol we can call each individual item in the list and make a text file
# we can then use a read loop to analyse the each line in the file we created individually
# the $line represents the row in the FASTQLIST.TXT file you are up to and it goes through every line one after another
# -t is the number of CPU cores to use in fastqc and -o is the output path

if $DO_FASTQC
then

PATH_FASTQC=$WORKPATH/02_FASTQC

	if [ ! -d $PATH_FASTQC ]
	then
		echo " $(date) mkdir $PATH_FASTQC "
		mkdir $PATH_FASTQC
	fi

ls 

ls ${FASTQ[@]} > FASTQLIST.txt


cat FASTQLIST.txt | while read line
do

echo "Starting file $line "

	fastqc  -t $core -o $PATH_FASTQC -f fastq $line
	echo "$(date) Done fastqc"
done

fi



# Step 3
# We now want to convert the fastq files into bam files 
# The tricky part about this is that we have both forward and reverse reads we need to combine
# so we use sed to edit the file names in our text lists so we can call each file correctly. 
# After editing we have removed all the forward and reverse read information so it's just the sample name. 
# the -u just filters to unique 
# then we can define forward and reverse reads using the sample name and the correct add ons to represent the forward/reverse section
# Then we run fastqtosam for each line we have in the txt file and the variables (names of everything) update to the correct ones
# gatk uses java so --java-options -Xmx24G indicates the amount of ram to let java use to analyse the sequences. if the run fails you can increase the amount of ram 
# here -Xmx24G indicates 24GB is devoted but some of the later steps require up to 60GB

if $DO_FASTQ_TO_BAM
then

PATH_BAM_UNALIGNED=$WORKPATH/03_BAM_UNALIGNED

	if [ ! -d $PATH_BAM_UNALIGNED ]
	then
		echo " $(date) mkdir $PATH_BAM_UNALIGNED "
		mkdir $PATH_BAM_UNALIGNED
	fi

	sed 's/_R.*_001.fastq.gz/_R/g' FASTQLIST.txt > FASTQLIST2.txt
	sed 's/.*R..//g' FASTQLIST2.txt > FASTQLIST3.txt
	sort FASTQLIST3.txt -u > FASTQLIST4.txt

cat FASTQLIST4.txt | while read line
do
Forward=${WORKPATH}/reads/R1/${line}1_001.fastq.gz
Reverse=${WORKPATH}/reads/R2/${line}2_001.fastq.gz

echo "${line}"| sed 's/_S.*//g' >sampid
name=(`cat sampid`)

gatk --java-options -Xmx24G FastqToSam \
       -F1 ${Forward} \
       -F2 ${Reverse} \
       -O ${PATH_BAM_UNALIGNED}/unaligned_read_pairs_${line}.bam \
       -SM ${name} --SORT_ORDER queryname \
	--MAX_RECORDS_IN_RAM 200000
    

done

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

if $DO_ADD_READGROUP_DATA
then
PATH_BAM_UNALIGNED=$WORKPATH/03_BAM_UNALIGNED
	if [ ! -d $PATH_BAM_UNALIGNED ]
	then
		echo " $(date) mkdir $PATH_BAM_UNALIGNED "
		mkdir $PATH_BAM_UNALIGNED
	fi



ls $WORKPATH/03_BAM_UNALIGNED | grep '.bam' > bam_list.txt


cat bam_list.txt | while read line

do

name=`awk '{ sub(/.*pairs_/, ""); sub(/_S.*/, ""); print }' <<< $line`
echo " $name "
line2=${line%%_R.bam}
echo " $line2 "

gatk --java-options -Xmx24G AddOrReplaceReadGroups \
       -I ${WORKPATH}/03_BAM_UNALIGNED/${line} \
       -O ${WORKPATH}/03_BAM_UNALIGNED/${line2}_RGup.bam \
       -RGID 4 \
       -RGLB lib1 \
       -RGPL ILLUMINA \
       -RGPU unit1 \
       -RGSM $name

done 
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

if $DO_MARK_ADAPTORS
then

PATH_BAM_UNALIGNED=$WORKPATH/03_BAM_UNALIGNED				
PATH_TEMPFILES=$WORKPATH/TEMPFILES
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

declare -a bamlistRGUP=(`ls ${WORKPATH}/03_BAM_UNALIGNED/*RGup.bam*`)

for bam in ${bamlistRGUP[@]}; do
echo " $bam "
bam2=(`awk '{ sub(/.*unaligned/, "unaligned"); sub(/_RGup.bam/, ""); print }'<<< $bam`)
echo " $bam2 "

gatk --java-options -Xmx28G MarkIlluminaAdapters\
 -I $bam\
 -O ${PATH_BAM_UNALIGNED}/${bam2}_adap.bam\
 -M ${PATH_BAM_UNALIGNED}/${bam2}_markilluminaadapters_metrics.txt

done

fi



#Step 6 bam to fastq
# This converts the bam file to a fastq so it can be aligned to the heho genome
# interleave is because the data is paired end (R1 and R2) but combined into one file which we did earlier. 
# Clipping attributes can be changed but read up about the filters and things in this step if you want to change some things

if $DO_BAM_TO_FASTQ
then

PATH_TEMPFILES=$WORKPATH/TEMPFILES
PATH_FASTQ_UPDATED=$WORKPATH/06_FASTQ_UPDATED				
	if [ ! -d $PATH_FASTQ_UPDATED ]
	then
		echo " $(date) mkdir $PATH_FASTQ_UPDATED "
		mkdir $PATH_FASTQ_UPDATED
	fi

declare -a bamlistadap=(`ls ${WORKPATH}/03_BAM_UNALIGNED/*adap.bam*`)

for bam in ${bamlistadap[@]}; do

bam2=(`awk '{ sub(/.*unaligned_/, "unaligned_"); sub(/_adap.bam/, ""); print }'<<< $bam`)

gatk --java-options -Xmx28G SamToFastq\
 -I $bam\
 --FASTQ ${WORKPATH}/06_FASTQ_UPDATED/${bam2}interleaved.fq\
 --CLIPPING_ATTRIBUTE XT\
 --CLIPPING_ACTION 2\
 --INTERLEAVE true\
 --TMP_DIR ${PATH_TEMPFILES}



done

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


if $DO_BWA
then		
PATH_BWA=$WORKPATH/08_BWA_aligned		
	if [ ! -d $PATH_BWA ]
	then
		echo " $(date) mkdir $PATH_BWA "
		mkdir $PATH_BWA
	fi


declare -a fqlist=(`ls ${WORKPATH}/06_FASTQ_UPDATED/*interleaved*`)

for fq in ${fqlist[@]}; do

fq2=(`awk '{ sub(/.*unaligned_/, "aligned_"); sub(/interleaved.fq/, "interleaved"); print }'<<< $fq`)

echo " $fq2 "

bwa mem -t $core -p -M ${GENOMEPATH}/HeHo_1.0_HiC.fasta $fq > ${PATH_BWA}/${fq2}.sam



done

fi




#Step 9 Merge raw bam with alignment files.
# After alignment to the genome this step copies over all the information from the aligned sequences and gives it back to the raw bam from step 6
# It then uses that information plus the read pairs to correct the whole set of sequences and chromosome maps to make sure all read pairs are actually next to each other
# It also allows for reads to be removed based on set filters. This preseves all sequence data and leads to overall higher accuracy in downstream analysis


if $DO_MERGE_BAM_ALIGNMENT
then		
PATH_BAM_MERGED=$WORKPATH/09_BAM_MERGED		
	if [ ! -d $PATH_BAM_MERGED ]
	then
		echo " $(date) mkdir $PATH_BAM_MERGED "
		mkdir $PATH_BAM_MERGED
	fi


declare -a samlist=(`ls ${WORKPATH}/08_BWA_aligned/*interleaved.sam*`)

for sam in ${samlist[@]}; do

name=(`awk '{ sub(/.*aligned/, ""); sub(/interleaved.sam/, ""); print }'<<< $sam`)

echo " $name "

gatk --java-options -Xmx20G MergeBamAlignment \
-R ${GENOMEPATH}/HeHo_1.0_HiC.fasta \
--UNMAPPED_BAM $WORKPATH/03_BAM_UNALIGNED/unaligned${name}_RGup.bam \
--ALIGNED_BAM $sam \
-O ${PATH_BAM_MERGED}/aligned_merged_${name}.bam \
--CREATE_INDEX true \
--ADD_MATE_CIGAR true \
--CLIP_ADAPTERS \
--CLIP_OVERLAPPING_READS true \
--INCLUDE_SECONDARY_ALIGNMENTS true \
--MAX_INSERTIONS_OR_DELETIONS -1 \
--PRIMARY_ALIGNMENT_STRATEGY MostDistant \
--ATTRIBUTES_TO_RETAIN XS \
--TMP_DIR $WORKPATH/TEMPFILES


done

fi

# next we want to sort the bam based on the coordinates of all sequences (base pair position on the genome)
# If there is trouble in this step upp the max records in Ram as it may require more memory

# Step 10. Sort bam
if $DO_SORT_BAM
then
PATH_BAM_MERGED_SORTED=$WORKPATH/10_BAM_MERGED_SORTED

	if [ ! -d $PATH_BAM_MERGED_SORTED ]
	then
		echo " $(date) mkdir $PATH_BAM_MERGED_SORTED "
		mkdir $PATH_BAM_MERGED_SORTED
	fi

declare -a samlist=(`ls ${WORKPATH}/09_BAM_MERGED/*bam*`)


for sam in ${samlist[@]}; do

name=(`awk '{ sub(/.*aligned_merged__read_pairs/, ""); sub(/bam/, ""); print }'<<< $sam`)

echo " $name "


gatk --java-options -Xmx24g SortSam \
	-I $sam -O ${PATH_BAM_MERGED_SORTED}/aligned_merged${name}srt.bam -SO coordinate --MAX_RECORDS_IN_RAM 200000
			
done

fi


# Step 11 Qualimap 

# Let's generate some summary statistics of what the sequences look like now
# Note, this requiers an older version of java to work so we 
	
if $DO_QUALIMAP
	then

PATH_QUALIMAP=$WORKPATH/11_QUALIMAP


		if [ ! -d $PATH_QUALIMAP ]
		then
			mkdir $PATH_QUALIMAP
		fi
		echo "$(date)	Start qualimap"

declare -a bamsrtlist=(`ls ${WORKPATH}/10_BAM_MERGED_SORTED/*srt*`)

		module unload java java/1.8.0_66
		module load java/1.7.0_51

for bam in ${bamsrtlist[@]}; do

name=(`awk '{ sub(/.*aligned_merged_/, ""); sub(/.srt.bam/, ""); print }'<<< $bam`)

	qualimap --java-mem-size=16000M
	qualimap --java-mem-size=16000M bamqc -bam $bam -outdir ${PATH_QUALIMAP} -outformat pdf -outfile ${name} -nt $core 	

done
	fi



module unload java/1.7.0_51
module load java/1.8.0_66


# 

#Step 12. Flag PCR duplicates	

# This step does as it says, it just flags potential PCR duplicates and marks them similarly to how the illumina adaptors were marked

if $DO_FLAG_PCRDUPS

then
BAM_MERGED_FILTERED=$WORKPATH/12_BAM_MERGED_FILTERED
	if [ ! -d $BAM_MERGED_FILTERED ]
	then
		echo " $(date) mkdir $BAM_MERGED_FILTERED "
		mkdir $BAM_MERGED_FILTERED
	fi

declare -a bamsrtlist=(`ls ${WORKPATH}/10_BAM_MERGED_SORTED/*srt*`)

for bam in ${bamsrtlist[@]}; do

name=(`awk '{ sub(/.*aligned_merged_/, ""); sub(/.srt.bam/, ""); print }'<<< $bam`)
	echo $bam	
	echo $name

gatk --java-options -Xmx24g MarkDuplicates -I ${bam} -O ${BAM_MERGED_FILTERED}/aligned_${name}_PCRdup.bam -M ${BAM_MERGED_FILTERED}/${name}.txt

done

fi




# Step 13. build bam index

# This builds index files for each of the bam sample files so they can be read easier 

if $DO_BUILD_BAMINDEX

then

declare -a bamlistsrt=(`ls ${WORKPATH}/12_BAM_MERGED_FILTERED/*PCR*`)

for bam in ${bamlistsrt[@]}; do

echo " ${bam}"

samtools index -b $bam

done

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

vcfname=(/data/group/murphylab/project/bem/heho_WGS/Test_2_inds/18_GENOTYPEGVCFS_FINAL/cohort_miners_genotypescalled.g.vcf.gz)

declare -a bamlistsrt=(`ls ${WORKPATH}/12_BAM_MERGED_FILTERED/*PCRdup.bam`)

for bam in ${bamlistsrt[@]}; do

name=(`awk '{ sub(/.*aligned_/, ""); sub(/_PCRdup.bam/, ""); print }'<<< $bam`)

gatk --java-options -Xmx40g BaseRecalibrator \
			--input $bam \
			--known-sites $vcfname \
			--output recalctable_${name}.txt \
			--reference ${GENOMEPATH}/HeHo_1.0_HiC.fasta



done 

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



declare -a bamlistsrt=(`ls ${WORKPATH}/12_BAM_MERGED_FILTERED/*PCRdup.bam`)

for bam in ${bamlistsrt[@]}; do

name=(`awk '{ sub(/.*aligned_/, ""); sub(/_PCRdup.bam/, ""); print }'<<< $bam`)
name2=(`awk '{ sub(/.*aligned/, "aligned"); sub(/NAN/, ""); print }'<<< $bam`)


gatk --java-options -Xmx40g ApplyBQSR \
   -R ${GENOMEPATH}/HeHo_1.0_HiC.fasta \
   -I $bam \
   --bqsr-recal-file recalctable_${name}.txt \
   -O ${PATH_BAM_RECALIBRATED}/${name2}PCRdup.bam


# Not sure yet, but may need to move the .bai files over to, it may create its own though, we will see. 

done

fi

# Step 16 A

# Step 16 took forever per individual trying to run with just the internal setting CPUs to run everything
# It took up to 4 days per sample which for 10 samples is 40 days and will be run up to 3 times for 4 months of analysis
# but it can be sped up a lot using another form of paralelisation through analysing each chromosome per ind by itself and combining these. ( ~ 16-20 hours per ind instead, so a 4-6 fold speed increase)
# This splits the genome into each chromosome scaffolds (906 chromosome scaffolds, of which 24 are actually chromosomes, the rest are loose bits of DNA)

# The awk is just collecting the second column of data ($2) (minus the first row which is header NR!=1) from the file HeHo_1.0_HiC.dict which we created when indexing the genome.
# This file is a short hand list of all chromosomes and their sizes so it lets us quickly define every chromosome to split the data into 

if $DO_SPLT_CHRS

then

cd ${GENOMEPATH}

awk '{if(NR!=1)print $2}' HeHo_1.0_HiC.dict > chrlst.txt

awk '{sub(/SN:/, ""); print $1}' chrlst.txt > tmp && mv tmp chrlst.txt


mv chrlst.txt ../

cd ${WORKPATH}

fi

# Step X (Will filter based on whether step 14 was completed or skipped temporary)
# this just helps identify where in the step 14-Step 18 loop you are up to
if $DO_PRELIM_CALLS	

then

BAMDIRECTORY=$WORKPATH/12_BAM_MERGED_FILTERED

else 

BAMDIRECTORY=$WORKPATH/15_BAM_RECALIBRATED

fi

echo " $BAMDIRECTORY "



PATH_GENOTYPING=$WORKPATH/16_GENOTYPING

# Step 16 B

# here we use echo and the chromosome list to create 906 gatk haplotypecaller commands, one for each chromosome
# we then run all of them through parallel using a pipe " cat haplotypecaller_job.list | parallel -j $SLURM_CPUS_ON_NODE " 
# and set the number of CPUs (here 16). This sorts through each chromosome separately
# Then those haplotypes per chromosome are stitched back together using GatherVcfs and saving so much time
# To gather them, they need to be in order of CH1 to CH906 which is not how they are read in a standard ls command. So we need to
# create our own list based on the fist part of the name and then add the correct number from 1-906 but make it generalisable to fit in any genome.
# Note, the HPC keeps calling a segmentation fault through parallel on this part. I'm not sure why it's struggling so much but a solution is to purge the modules and build them again
# This blanks out everything and lets us remove potential factors causing the segmentation fault

if $DO_HAPLOTYPE_CALLER
then
	if [ ! -d $PATH_GENOTYPING ]
	then
		mkdir $PATH_GENOTYPING
	fi
	# Generate GVCF files
	# Call germline SNPs and indels via local re-assembly of haplotypes

		echo "$(date)	Start HaplotypeCaller"



module load fastqc/0.11.9
module load bcl2fastq-gcc/2.20.0 
module load stacks-gcc7/2.41
module load bwa-gcc/0.7.17 #(v0.7.17-r1188)
module load samtools-gcc/1.6
module load parallel/20170722
module load qualimap/2.1.3
module load vcftools-gcc/0.1.13
module load picard/2.2.2




echo " $BAMDIRECTORY "

declare -a bamlistgenotyping=(`ls ${BAMDIRECTORY}/*PCRdup.bam`)



for bam in ${bamlistgenotyping[@]}; do


name=(`awk '{ sub(/.*aligned_/, ""); sub(/_PCRdup.bam/, ""); print }'<<< $bam`)

cat chrlst.txt | while read line
do

echo "gatk --java-options -Xmx16g HaplotypeCaller --input ${bam} --output ${PATH_GENOTYPING}/${name}_${line}.inter.g.vcf --reference ${GENOMEPATH}/HeHo_1.0_HiC.fasta --emit-ref-confidence GVCF --heterozygosity-stdev 0.015 --native-pair-hmm-threads 1 --tmp-dir ${WORKPATH}/TEMPFILES --L ${line}" 


done > haplotypecaller_job.list


module purge


module load bcl2fastq-gcc/2.20.0 
module load stacks-gcc7/2.41
module load bwa-gcc/0.7.17 #(v0.7.17-r1188)
module load samtools-gcc/1.6
module load parallel/20170722
module load java/1.7.0_51
module load gatk/4.0.11.0
module load qualimap/2.1.3
module load vcftools-gcc/0.1.13
module load picard/2.2.2



echo "Starting: $(date)"

cat haplotypecaller_job.list | parallel -j $SLURM_CPUS_ON_NODE

echo "Finished: $(date)"


	
		# recombine idividual chrom files into single gvcf file. # Est time per recombine 2-2.5 hours per ind
			cd $PATH_GENOTYPING
			ls ${WORKPATH}/16_GENOTYPING/*inter.g.vcf >gvcfssingleind.list
			a=(`head -1 gvcfssingleind.list`)
			b=(`wc -l gvcfssingleind.list`)


			
			
			echo " $a "
			echo " $b "

			name1=(`awk '{ sub(/_HiC.*/, ""); sub(/.*GENOTYPING\//, ""); print }'<<< $a`)
			name2=(`awk '{ sub(/scaffold.*/, "scaffold_"); sub(/NAN\//, ""); print }'<<< $a`)

		for (( c=1; c<=$b; c++ ))
		do

			echo "${name2}${c}.inter.g.vcf"

		done > gvcfssingleindordered.list


			echo " $name "

				gatk --java-options -Xmx40g GatherVcfs \
				-R ${GENOMEPATH}/HeHo_1.0_HiC.fasta \
				--INPUT gvcfssingleindordered.list \
				-O ${PATH_GENOTYPING}/${name1}.raw.g.vcf
			rm *inter*
			cd ..


done

		echo "$(date)	Done HaplotypeCaller"


module load mira/4.0.2
module load fastp/0.20.1
module load mitobim/20160526
module load python-gcc7/3.6.8



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


ls ${WORKPATH}/16_GENOTYPING/*raw.g.vcf >gvcfs.list
cd ${WORKPATH}
 gatk --java-options -Xmx40g CombineGVCFs \
   -R ${GENOMEPATH}/HeHo_1.0_HiC.fasta \
	--variant gvcfs.list \
	-O ${WORKPATH}/17_GENOTYPEGVCFSCOMBINED/cohort_miners.g.vcf.gz


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


gatk --java-options -Xmx48g GenotypeGVCFs \
	-R ${GENOMEPATH}/HeHo_1.0_HiC.fasta \
	-V ${WORKPATH}/17_GENOTYPEGVCFSCOMBINED/cohort_miners.g.vcf.gz \
	-O $PATH_GENOTYPEGVCFS/cohort_miners_genotypescalled.g.vcf.gz

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


vcfname=(/data/group/murphylab/project/bem/heho_WGS/Test_2_inds/18_GENOTYPEGVCFS_FINAL/cohort_miners_genotypescalled.g.vcf.gz)


declare -a indrawvcfs=(`ls ${WORKPATH}/16_GENOTYPING/*raw.g.vcf`)
# pop summary
# can also remove indels or other components at this step (see below comment). This is important to decide on what you want to do with the data afterwards. You will need to research what is
# optimal design for your project.
#http://vcftools.sourceforge.net/man_latest.html
#--remove-indels

vcftools --gzvcf $vcfname --recode --recode-INFO-all --out ${PATH_VCF_SUMMARY}/cohortminers_filtered.vcf


# next we want individual specific data summaries from the previous analysis.


for ind in ${indrawvcfs[@]}; do


indname=(`awk '{ sub(/.raw.*/, ""); sub(/.*GENOTYPING\//, ""); print }'<<< $ind `)


vcftools --vcf $ind --recode --recode-INFO-all --out ${PATH_VCF_SUMMARY}/${indname}_filtered.vcf


done

fi

# At this point you can then filter your data based on missing data thresholds for inds/loci, HWE, linkage, or MAF/MAC using simple vcftools and setting the correct variable. Or you 
# can import it into R and do it there using vcfR (See my golden perch bioinformatics of vcfR for examples if you need [Kat has it but it's also on my github] ). I also have code there to
# add the relevant information sections and to convert it to a genlight if that's what you want. 






#Mitochondrial analysis

# This uses the raw reads data. 
# Note, I created code up until the completion of Mitobim but the sequence still needs to be circularised and annotated to be a full mitochondrial genome


PATH_M1_READ_MERGING=$WORKPATH/M1_READ_MERGING
PATH_M2_MITOBIM=$WORKPATH/M2_MITOBIM
PATH_M3_MITO_BIM_RESULTS=$WORKPATH/M3_MITOBIM_RESULTS

if $DO_MITOCHONDRIAL_ANALYSIS
then
	if [ ! -d $PATH_M2_MITOBIM ]
	then
		mkdir $PATH_M1_READ_MERGING
		mkdir $PATH_M2_MITOBIM
		mkdir $PATH_M3_MITO_BIM_RESULTS
	fi



fi

# Step M1. Starts using raw data from the beginning

# unlike the main genome sequencing, we need to actually trim the data, not just mark low quality bits. 
# fastp trims and there are many settings which can be changed here for optimal trimming. right now I set it to quality score of 20 minimum and read pairs are identified 
# based on a minimum of 15bp overlap

if $DO_TRIM_MITO

then

declare -a fastqlistR1=(`ls ${WORKPATH}/reads/R1/*fastq*`)
declare -a fastqlistR2=(`ls ${WORKPATH}/reads/R2/*fastq*`)

for fa in ${fastqlistR1[@]}; do

nameR1=(`awk '{ sub(/.fastq.gz/, ""); sub(/.*reads\/R1\//, ""); print }'<<< $fa`)

nameR2=(`awk '{ sub(/R1/, "R2"); sub(/.*reads\/R2\//, ""); print }'<<< $nameR1`)

nameonly=(`awk '{ sub(/R1_/, ""); sub(/\./, ""); print }'<<< $nameR1`)

# Define variables, note, not all these need to be in the loop. I am aware it will save 
# like 2 seconds to move some of them up.

#sample name
sample=($nameonly)

#location of raw forward Illumina reads
forw=(${WORKPATH}/reads/R1/${nameR1}.fastq.gz)
echo " $forw "#location of raw reverse Illumina reads
reve=(${WORKPATH}/reads/R2/${nameR2}.fastq.gz)
#location of fasta file used as seed for MITObim assembly
seed=(${GENOMEPATH}/Helmeted_honeyeater_Mitochondrial_genome.fasta)

#minimum length for reads after trimming
min_length=100
threads=$core

# Can change these variables a bit to include more data or less. There is not really
# A big issue though. The mitochondria should have more data than everything else so stringent
# or lax there will be ample quantities to reconstruct the mitochondrial genome

fastp -i $forw -I $reve \
	--merge \
	--merged_out $PATH_M1_READ_MERGING/${nameonly}_trimmed.fastq \
	-j fastp_a_miner${nameonly}.json -h fastp_a_miner${nameonly}.html \
	--average_qual 0 \
	--overlap_len_require 15 \
	-q 20
	-w $core


done

fi





#Step M2 
# Mitobim is set to 100 iterations but it took 3-4 per ind in the trial data so I don't think anywhere near that many is needed (it ends once it converges though so it's not an issue) 
# The program doesn't like writing different names for its iterations or overwritting files. This means that new folders need to be created for each ind in each mitbim run (code does this 
# and labels already)


if $DO_MITOBIM_RUN

then


# The program seems to fail through segmentation for some reason. The only way I could fix this is to purge all the environment and reload only mitobim. It worked though!
# I have no idea why though, there was no parallelisation in the script so i don't understand what was driving the segmentation fault.

module purge
module load mitobim/20160526


declare -a trimmedlist=(`ls ${PATH_M1_READ_MERGING}/*trimmed*`)
seed=(${GENOMEPATH}/Helmeted_honeyeater_Mitochondrial_genome.fasta)

cd $PATH_M2_MITOBIM
for fa in ${trimmedlist[@]}; do

nameonly=(`awk '{ sub(/_trimmed.fastq/, ""); sub(/.*MERGING\//, ""); print }'<<< $fa`)


echo " $nameonly "
echo " $fa " 
echo " $seed "

mkdir ${nameonly}
cd ${nameonly}

MITObim_1.8.pl -start 1 -end 100 -sample bird${nameonly} -ref honeyeater_mito_genome -readpool $fa --quick ${seed} --clean --NFS_warn_only > ${nameonly}log

cd ..


done

fi


if $DO_MITO_BIM_COLLATE_RESULT
then

echo " need to do this part still "

fi




# Let me know if you have any questions or issues

