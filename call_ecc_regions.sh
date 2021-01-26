#!/bin/bash
#MIT License
#
#Copyright (c) 2020 Pierre Michel Joubert
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
while getopts m:s:t:b: option
do
case "${option}"
in
m) MAPFILE=${OPTARG};;
s) SAMPLE=${OPTARG};;
t) THREADS=${OPTARG};;
b) FILTERED_BAMFILE=${OPTARG};;
esac
done

## USAGE ##
# call eccDNA forming regions using split reads and opposite facing read pairs
# options:
# -m mapfile with names of contigs of interest, as written in the fasta file used to make bwa genome database
# -s sample name/output prefix
# -t threads
# -b bamfile of mapped reads, filtered to contigs of interest


# divide up reads from bam file based on their orientation to ensure that both sides of split read are mapped in the same direction
# filter to split reads
# filter split reads so that either side of the junction is at least 20 bp
# make sure split reads appear only twice, clearly representing an eccDNA junction
# make sure split reads map to the same chromosome (split reads mapping to different chromosomes would need to be analyzed using a different pipeline because opposite facing read pairs wouldn't make sense)
# make sure split read halves are properly oriented to that they represent eccDNA junctions and not potential introns ( ---> gap --- is an eccDNA junction vs --- gap ---> is an intron)
samtools view -f 81 -F 4 ${FILTERED_BAMFILE} > ${SAMPLE}.reverseread1.sam.tmp
splitread_file="${SAMPLE}.reverseread1.sam"
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0, 2}' ${splitread_file}.tmp > ${splitread_file}.tmp.qualityfiltered
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' ${splitread_file}.tmp.qualityfiltered ${splitread_file}.tmp.qualityfiltered | sort -k1,1 -k18,18n > ${splitread_file}.tmp.samechromosome.exactlytwice.qualityfiltered
awk -v OFS='\t' '{
    prev=$0; f4=$4; f1=$1
    getline
    if ($1 == f1 && f4 > $4) {
        print prev
        print $0
    }
}' ${splitread_file}.tmp.samechromosome.exactlytwice.qualityfiltered | awk -v OFS='\t' '{$NF=""; print $0}' | sed 's/[ \t]\+$//' > ${splitread_file}.tmp.oriented.samechromosome.exactlytwice.qualityfiltered

samtools view -f 145 -F 4 ${FILTERED_BAMFILE} > ${SAMPLE}.reverseread2.sam.tmp
splitread_file="${SAMPLE}.reverseread2.sam"
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0, 2}' ${splitread_file}.tmp > ${splitread_file}.tmp.qualityfiltered
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' ${splitread_file}.tmp.qualityfiltered ${splitread_file}.tmp.qualityfiltered | sort -k1,1 -k18,18n > ${splitread_file}.tmp.samechromosome.exactlytwice.qualityfiltered
awk -v OFS='\t' '{
    prev=$0; f4=$4; f1=$1
    getline
    if ($1 == f1 && f4 > $4) {
        print prev
        print $0
    }
}' ${splitread_file}.tmp.samechromosome.exactlytwice.qualityfiltered | awk -v OFS='\t' '{$NF=""; print $0}' | sed 's/[ \t]\+$//' > ${splitread_file}.tmp.oriented.samechromosome.exactlytwice.qualityfiltered

samtools view -f 65 -F 20 ${FILTERED_BAMFILE} > ${SAMPLE}.forwardread1.sam.tmp
splitread_file="${SAMPLE}.forwardread1.sam"
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0, 2}' ${splitread_file}.tmp > ${splitread_file}.tmp.qualityfiltered
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' ${splitread_file}.tmp.qualityfiltered ${splitread_file}.tmp.qualityfiltered | sort -k1,1 -k18,18n > ${splitread_file}.tmp.samechromosome.exactlytwice.qualityfiltered
awk -v OFS='\t' '{
    prev=$0; f4=$4; f1=$1
    getline
    if ($1 == f1 && f4 > $4) {
        print prev
        print $0
    }
}' ${splitread_file}.tmp.samechromosome.exactlytwice.qualityfiltered | awk -v OFS='\t' '{$NF=""; print $0}' | sed 's/[ \t]\+$//' > ${splitread_file}.tmp.oriented.samechromosome.exactlytwice.qualityfiltered

samtools view -f 129 -F 20 ${FILTERED_BAMFILE} > ${SAMPLE}.forwardread2.sam.tmp
splitread_file="${SAMPLE}.forwardread2.sam"
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0, 2}' ${splitread_file}.tmp > ${splitread_file}.tmp.qualityfiltered
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' ${splitread_file}.tmp.qualityfiltered ${splitread_file}.tmp.qualityfiltered | sort -k1,1 -k18,18n > ${splitread_file}.tmp.samechromosome.exactlytwice.qualityfiltered
awk -v OFS='\t' '{
    prev=$0; f4=$4; f1=$1
    getline
    if ($1 == f1 && f4 > $4) {
        print prev
        print $0
    }
}' ${splitread_file}.tmp.samechromosome.exactlytwice.qualityfiltered | awk -v OFS='\t' '{$NF=""; print $0}' | sed 's/[ \t]\+$//' > ${splitread_file}.tmp.oriented.samechromosome.exactlytwice.qualityfiltered

# I did some merging of reads in the past but this ended up being detrimental
# these next two chunks should be removed eventually
# currently both of these should contain 0 reads
samtools view -f 16 -F 5 ${FILTERED_BAMFILE} > ${SAMPLE}.reversemerged.sam.tmp
splitread_file="${SAMPLE}.reversemerged.sam"
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0, 2}' ${splitread_file}.tmp > ${splitread_file}.tmp.qualityfiltered
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' ${splitread_file}.tmp.qualityfiltered ${splitread_file}.tmp.qualityfiltered | sort -k1,1 -k18,18n > ${splitread_file}.tmp.samechromosome.exactlytwice.qualityfiltered
awk -v OFS='\t' '{
    prev=$0; f4=$4; f1=$1
    getline
    if ($1 == f1 && f4 > $4) {
        print prev
        print $0
    }
}' ${splitread_file}.tmp.samechromosome.exactlytwice.qualityfiltered | awk -v OFS='\t' '{$NF=""; print $0}' | sed 's/[ \t]\+$//' > ${splitread_file}.tmp.oriented.samechromosome.exactlytwice.qualityfiltered

samtools view -F 21 ${FILTERED_BAMFILE} > ${SAMPLE}.forwardmerged.sam.tmp
splitread_file="${SAMPLE}.forwardmerged.sam"
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0, 2}' ${splitread_file}.tmp > ${splitread_file}.tmp.qualityfiltered
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' ${splitread_file}.tmp.qualityfiltered ${splitread_file}.tmp.qualityfiltered | sort -k1,1 -k18,18n > ${splitread_file}.tmp.samechromosome.exactlytwice.qualityfiltered
awk -v OFS='\t' '{
    prev=$0; f4=$4; f1=$1
    getline
    if ($1 == f1 && f4 > $4) {
        print prev
        print $0
    }
}' ${splitread_file}.tmp.samechromosome.exactlytwice.qualityfiltered | awk -v OFS='\t' '{$NF=""; print $0}' | sed 's/[ \t]\+$//' > ${splitread_file}.tmp.oriented.samechromosome.exactlytwice.qualityfiltered

# putting them all back together
cat ${SAMPLE}.reverseread1.sam.tmp.oriented.samechromosome.exactlytwice.qualityfiltered \
    ${SAMPLE}.reverseread2.sam.tmp.oriented.samechromosome.exactlytwice.qualityfiltered \
    ${SAMPLE}.forwardread1.sam.tmp.oriented.samechromosome.exactlytwice.qualityfiltered \
    ${SAMPLE}.forwardread2.sam.tmp.oriented.samechromosome.exactlytwice.qualityfiltered \
    ${SAMPLE}.reversemerged.sam.tmp.oriented.samechromosome.exactlytwice.qualityfiltered \
    ${SAMPLE}.forwardmerged.sam.tmp.oriented.samechromosome.exactlytwice.qualityfiltered > ${SAMPLE}.tmp.oriented.samechromosome.exactlytwice.qualityfiltered.all.sam

# converting to bed file
samtools view -b -h <(cat <(samtools view -H ${FILTERED_BAMFILE}) ${SAMPLE}.tmp.oriented.samechromosome.exactlytwice.qualityfiltered.all.sam) > ${SAMPLE}.tmp.oriented.samechromosome.exactlytwice.qualityfiltered.all.bam
bedtools bamtobed -i ${SAMPLE}.tmp.oriented.samechromosome.exactlytwice.qualityfiltered.all.bam | sort -k4,4 -k2,2n > ${SAMPLE}.splitreads.bed

# merging split read halves into single, putative eccDNA forming regions to be confirmed or rejected
awk -v OFS='\t' '{
    prev=$0; f2=$2; f4=$4
    getline
    if ($4 == f4 && f2 < $2) {
        print $1, f2, $3, $4
    }
}' ${SAMPLE}.splitreads.bed > ${SAMPLE}.merged.splitreads.bed

# length filter because we don't expect eccDNAs to be that big
# could be tweaked potentially but this gets rid of very few split reads
awk -v OFS='\t' '$3-$2<1000000' ${SAMPLE}.merged.splitreads.bed > ${SAMPLE}.lengthfiltered.merged.splitreads.bed

# get outward facing read pairs using sam flags
# convert to bed file
# fix names for filtering
# filter to appearing only exactly twice, meaning that only complete read pairs are present
samtools view ${FILTERED_BAMFILE} | awk '{ if (($2 == 81 || $2 == 83 || $2 == 145 || $2 == 147 ) && $9 > 0) print $0 ; else if (($2 == 97 || $2 == 99 || $2 == 161 || $2 == 163) && $9 <0) print $0}' | cat <(samtools view -H ${FILTERED_BAMFILE}) - | samtools view -b -h - > ${SAMPLE}.tmp.outwardfacing.bam
bedtools bamtobed -i ${SAMPLE}.tmp.outwardfacing.bam | sort -k 4,4 > ${SAMPLE}.tmp.outwardfacing.bed
mv ${SAMPLE}.tmp.outwardfacing.bed ${SAMPLE}.tmp.outwardfacing.bed.old
awk 'BEGIN {OFS="\t"}; {print $1,$2,$3,substr($4, 1, length($4)-2),$5,$6}' ${SAMPLE}.tmp.outwardfacing.bed.old > ${SAMPLE}.tmp.outwardfacing.bed.old.trimmed
awk 'NR==FNR{a[$4]++; next} a[$4]==2' ${SAMPLE}.tmp.outwardfacing.bed.old.trimmed ${SAMPLE}.tmp.outwardfacing.bed.old.trimmed > ${SAMPLE}.outwardfacing.bed

# change names of scaffolds using mapfiles for compatability with any genome
chrom_count=$(wc -l ${MAPFILE} | awk '{print $1}')
for (( i = 1 ; i < ${chrom_count}+1; i++)); do echo $i ; done > ${SAMPLE}.tmp.chrom_count
paste ${SAMPLE}.tmp.chrom_count ${MAPFILE} > ${SAMPLE}.tmp.chrom_count_and_names
awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' ${SAMPLE}.tmp.chrom_count_and_names ${SAMPLE}.outwardfacing.bed > ${SAMPLE}.outwardfacing.renamed.bed
awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' ${SAMPLE}.tmp.chrom_count_and_names ${SAMPLE}.lengthfiltered.merged.splitreads.bed > ${SAMPLE}.lengthfiltered.merged.splitreads.renamed.bed

# merge outward facing read pairs into single lines for confirming using python script
sort -k4,4 -k2,2n ${SAMPLE}.outwardfacing.renamed.bed > ${SAMPLE}.sorted.outwardfacing.renamed.bed
awk -v OFS='\t' '{
    prev=$0; f2=$2; f3=$3; f4=$4
    getline
    if ($4 == f4 && f2 < $2 && f3 <$3) {
        print $1, f2, $3, f3, $2
    }
}' ${SAMPLE}.sorted.outwardfacing.renamed.bed > ${SAMPLE}.sorted.grouped.outwardfacing.renamed.bed

# use GNU parallel to speed things up
# split into chunks first then each thread works on a chunk
# parallel.confirmed are split reads confirmed by opposite facing read pairs
split --number=l/${THREADS} --numeric-suffixes=1 --additional-suffix=.bed ${SAMPLE}.lengthfiltered.merged.splitreads.renamed.bed ${SAMPLE}.lengthfiltered.merged.splitreads.renamed.
parallel -j ${THREADS} --link python ${ECC_CALLER_PYTHON_SCRIPTS}/ecc_caller_anygenome_confirmsrs_numpy_gnuparallel.py ${SAMPLE}.lengthfiltered.merged.splitreads.renamed.{}.bed ${SAMPLE}.sorted.grouped.outwardfacing.renamed.bed ${SAMPLE} ${chrom_count} {} ::: $(seq -w 1 ${THREADS})
cat $(find . -maxdepth 1 -name "${SAMPLE}.parallel.confirmed.*" | xargs -r ls -1 | tr "\n" " ") > ${SAMPLE}.parallel.confirmed

# convert scaffolds to 1 index from 0 index
# rename scaffolds in parallel.confirmed
paste ${MAPFILE} ${SAMPLE}.tmp.chrom_count > ${SAMPLE}.tmp.chrom_names_and_count
awk -v OFS='\t' '{print $1+1, $2, $3}' ${SAMPLE}.parallel.confirmed > ${SAMPLE}.parallel.plusone.confirmed
awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' ${SAMPLE}.tmp.chrom_names_and_count ${SAMPLE}.parallel.plusone.confirmed > ${SAMPLE}.confirmedsplitreads.bed
# clean up
# parallel.confirmed MUST be removed here otherwise it causes major problems with the output
rm ${SAMPLE}.parallel.confirmed*

# rm ${SAMPLE}.tmp

rm ${SAMPLE}.lengthfiltered.merged.splitreads.renamed.*.bed
