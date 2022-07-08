#!/bin/bash
#
# The script requires python3, sqlite3, bedtools
#
##################################################################################################################
#
# Dependencies that have to be edited 
#
##################################################################################################################
#
PYTHON_LOC="/home/matthias/miniconda3/bin/python"
SQLITE3_LOC="/home/matthias/miniconda3/bin/sqlite3"
BEDTOOLS_LOC="/home/matthias/miniconda3/bin/bedtools"
#
##################################################################################################################
#
# Files that have to be available 
#
##################################################################################################################
#
SNPFISHER_FILE="./SNPfisher_FLI+TL+WIK.vcf"   # obtained from https://snpfisher.nichd.nih.gov/snpfisher/tracks.html
GCF_FILE="./GCF_000002035.4_Zv9_genomic.bed"  # obtained from zebrafish v9 annotation project, notice GCRv9!!!!!!!
CORRESPONDENCE_FILE="./correspondence.tsv"    # created based on the GCF file for chromosome name matching
#
###################################################################################################################
#
# 3 replicas each of 
#
# WT      (WT_1, WT_2, WT_3) and 
# ADAR KO (KO_1, KO_2, KO_3) 
#
# fish were analyzed 
# 
# starting point from the pipeline output from REDiscover
# these files are called WT_1_editing, WT_2_editing....
#
##################################################################################################################
#
# Create suitable input SNP input file for sqlite
#
##################################################################################################################
#
grep -v "#" ${SNPFISHER_FILE} | awk '{print $1, $2, $3, $4, $5}' | tr -t " " "\t" > SNP_simple.tsv
#
#*****************************************************************************************************************
#
# Iterate over the samples
#
#*****************************************************************************************************************
#
for i in WT_1 WT_2 WT_3 KO_1 KO_2 KO_3
do
##################################################################################################################
#
# Create suitable input sample input file for sqlite
#
##################################################################################################################
#
grep -v "^#" ${i}_editing > ${i}_editing_no_comments
#
##################################################################################################################
#
# Exclude editing sites that coincide with SNPs, or that have an editing penetrance of 1
#
# output in WT_1.sqlite_result, WT_2.sqlite_result, etc.
#
##################################################################################################################
#
cat <<eof | sed "s/WT_1/${i}/g" | $SQLITE3_LOC
.separator "\t"
CREATE TABLE SNP (chrom text, pos int, dummy text, canonical text, actual text);
.import SNP_simple.tsv SNP
CREATE TABLE WT_1 (chrom text, pos int, what text, number int, pen float);
.import WT_1_editing_no_comments WT_1
CREATE TABLE correspondence (edit_name text, SNP_name text);
.import ${CORRESPONDENCE_FILE} correspondence 
CREATE TABLE WT_1_nice as select correspondence.edit_name, correspondence.SNP_name, WT_1.pos, WT_1.what, WT_1.number, WT_1.pen from correspondence inner join WT_1 on correspondence.edit_name=WT_1.chrom;
CREATE TABLE filtered as select * from WT_1_nice left join SNP on WT_1_nice.SNP_name=SNP.chrom and WT_1_nice.pos=SNP.pos;
CREATE TABLE filtered2 as select * from filtered where canonical is NULL;
CREATE TABLE filtered3 as select * from filtered2 where pen!=1.0;
.mode tabs
.output WT_1.sqlite_result
SELECT * from filtered3;
.quit
eof
#
##################################################################################################################
#
# Apply the clustering criterion: 
#
# use a symmetric window of 9 around the site of interest
# keep the site if at least four in the 9 around in the window (5 including the site of interest) are the same type
#
# this choice seems to be a good compromise between specificity and sensitivity for A--
#
# output in WT_1.filtered_4_9, WT_2.filtered_4_9, etc
#
##################################################################################################################
#
$PYTHON_LOC << eof > ${i}.filtered_4_9
with open("${i}.sqlite_result","r") as fp:
    min9=["x","x","x","x","x","x","x"]
    min8=["x","x","x","x","x","x","x"]
    min7=["x","x","x","x","x","x","x"]
    min6=["x","x","x","x","x","x","x"]
    min5=["x","x","x","x","x","x","x"]
    min4=["x","x","x","x","x","x","x"]
    min3=["x","x","x","x","x","x","x"]
    min2=["x","x","x","x","x","x","x"]
    min1=["x","x","x","x","x","x","x"]
    line="start"
    while(line):
        line=fp.readline()
        if (line==""):
            break
        zero=line.split("\t")
        if (min4[3]!="x"):
            outstring=min8[3]+" "+min7[3]+" "+min6[3]+" "+min5[3]+" "+min4[3]+" "+min3[3]+" "+min2[3]+" "+min1[3]+" "+zero[3]
            counter=outstring.count(min4[3])
            #
            # require a minimum neighbor score of 5
            #
            if (counter>=5):
                print(min4[0],min4[1],min4[2],min4[3],min4[4],min4[5])
        min9=min8
        min8=min7
        min7=min6
        min6=min5
        min5=min4
        min4=min3
        min3=min2
        min2=min1
        min1=zero
eof

##################################################################################################################
#
# Converted filtered editing sites to .bed sites 
# 
# output in WT_1_plus_minus.bed, WT_2_plus_minus.bed, etc.
#
##################################################################################################################

$PYTHON_LOC << eof > ${i}_plus_minus.bed
with open("${i}.filtered_4_9","r") as fp:
    while(True):
        line=fp.readline()
        if (line==""):
            break
        else:
            linevec=(line.strip('\n')).split(" ")
            outstring=linevec[0]+"\t"+linevec[2]+"\t"+linevec[2]+"\t"+linevec[3][0:4]+"\t"+linevec[5]
            print(outstring)
eof


##################################################################################################################
#
# assign the editing sites to genes
# 
# output in WT_1_intersect.bed, WT_2_intersect.bed, etc 
#
##################################################################################################################

bedtools intersect -a ${i}_plus_minus.bed -b ${GCF_FILE} -wa -wb >  ${i}_intersect.bed

##################################################################################################################
#
# fix annotations according to strand 
#
# output in WT_1_intersect_annotate.bed, WT_2_intersect_annotate.bed
#
##################################################################################################################

$PYTHON_LOC << eof > ${i}_intersect_annotate.bed
def complement(x):
	if (x=="A"):
		return "T"
	elif (x=="C"):
		return "G"
	elif (x=="G"):
		return "C"
	elif (x=="T"):
		return "A"
	elif (x=="+"):
		return "-"
	elif (x=="-"):
		return "+"
	else:
		return "N"
		
def reverse(symbol):
	a0=complement(symbol[0])
	a1=symbol[1]
	a2=complement(symbol[2])
	a3=complement(symbol[3])
	return a0+a1+a2+a3
	

with open("${i}_intersect.bed","r") as fp:
	line="start"
	while(line):
		line=fp.readline()
		if (line==""):
			break
		else:
			linevec=line.strip('\n').split('\t')
			if (linevec[10]=="+"):
				print(line.strip('\n')+"\t"+linevec[3])
			else:
				print(line.strip('\n')+"\t"+reverse(linevec[3]))
eof

##################################################################################################################
#
# deal with the problem of multiple transcripts for the same gene and alternative exons/alternative exon boundaries 
#
# output in WT_1_intersect_annotate_uniq.bed, WT_2_intersect_annotate_uniq.bed,etc
#
##################################################################################################################

uniq ${i}_intersect_annotate.bed > ${i}_intersect_annotate_uniq.bed

##################################################################################################################
#
# create a short excerpt used for the gene hitlist  
#
# output in WT_1.short.bed, WT_2.short.bed,etc
#
##################################################################################################################

awk -v OFS='\t' '{if ($12=="A>G-") print $9, $5}' ${i}_intersect_annotate_uniq.bed > ${i}.short.bed

##################################################################################################################
#
# create the list with gene based editing scores (sum of site penetrances)
#
# output in WT_1.hitlist, WT_2.hitlist, etc
#
##################################################################################################################

cat <<eof | sed "s/WT_1/${i}/g" | sqlite3
create table edits (gene int, pen float);
.separator "\t"
.import WT_1.short.bed edits
.mode tabs
.output WT_1.hitlist
select gene, sum(pen) from edits group by gene;
.quit
eof

##################################################################################################################
#
# create a sorted version of the hitlist with most edited genes on top
#
# output in WT_1.hitlist.sorted, WT_2.hitlist.sorted, etc
#
##################################################################################################################

sort -k2,2nr ${i}.hitlist > ${i}.hitlist.sorted
#
##################################################################################################################
#
# do the metaanalysis for the types of changes
#
# output in WT_1_meta.log, WT_2_meta.log, etc
#
##################################################################################################################

input=${i}_intersect_annotate_uniq.bed
output=${i}_meta.log

echo -n "A>C "  > $output; awk '{print $12}' $input | grep "A>C-" | wc -l >> $output
echo -n "A>G "  >> $output; awk '{print $12}' $input | grep "A>G-" | wc -l >> $output
echo -n "A>T "  >> $output; awk '{print $12}' $input | grep "A>T-" | wc -l >> $output
echo ""         >> $output
echo -n "C>A "  >> $output; awk '{print $12}' $input | grep "C>A-" | wc -l >> $output
echo -n "C>G "  >> $output; awk '{print $12}' $input | grep "C>G-" | wc -l >> $output
echo -n "C>T "  >> $output; awk '{print $12}' $input | grep "C>T-" | wc -l >> $output
echo ""         >> $output
echo -n "G>A "  >> $output; awk '{print $12}' $input | grep "G>A-" | wc -l >> $output
echo -n "G>C "  >> $output; awk '{print $12}' $input | grep "G>C-" | wc -l >> $output
echo -n "G>T "  >> $output; awk '{print $12}' $input | grep "G>T-" | wc -l >> $output
echo ""         >> $output
echo -n "T>A "  >> $output; awk '{print $12}' $input | grep "T>A-" | wc -l >> $output
echo -n "T>C "  >> $output; awk '{print $12}' $input | grep "T>C-" | wc -l >> $output
echo -n "T>G "  >> $output; awk '{print $12}' $input | grep "T>G-" | wc -l >> $output
#
#*****************************************************************************************************************
#
# Complete iteration over the samples
#
#*****************************************************************************************************************
#
done






















