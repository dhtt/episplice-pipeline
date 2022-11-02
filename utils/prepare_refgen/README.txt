1. Define flank regions, discard exon bodies using GRanges 
script: /home/dhthutrang/ENCODE/refgen/get_flank_gtf.R 
input: /home/dhthutrang/ENCODE/refgen/GCF_000001405.39_GRCh38.p13_genomic.gtf 
output: /home/dhthutrang/ENCODE/refgen/reference_genome.fl200.2022.gtf

2. Modify NCBI chromosome nomenclature, filter out unconventional chromosomes 
script: /home/dhthutrang/ENCODE/refgen/prepare_refseq_gtf.sh
input: /home/dhthutrang/ENCODE/refgen/reference_genome.fl200.2022.gtf 
output: /home/dhthutrang/ENCODE/refgen/temp5.gtf

3. Export into gtf format:
script:  /home/dhthutrang/ENCODE/refgen/get_gtf_refgen.R
input: /home/dhthutrang/ENCODE/refgen/temp5.gtf
output: /home/dhthutrang/ENCODE/refgen/hg19.ncbiRefSeq.2022.gtf

3. Flatten gtf:
script: dexseq_prepare_annotation.py 
input: hg19.ncbiRefSeq.2022.gtf 
output: hg19.ncbiRefSeq.2022.flattened.gtf 

4. For annotation of flank_exon:
i. Get only TRUE exonic part (remove aggregate_gene)
grep exonic_part /home/dhthutrang/ENCODE/refgen/reference_genome.2021_.corrected.gtf > reference_genome.2021_.corrected.exonic_part.gtf
ii. Get only exonic part (remove aggregate_gene but still overlap with flank)
grep exonic_part refgen_collapsed.gtf > refgen_collapsed.exonic_part.gtf
iii. Overlap flank /home/dhthutrang/ENCODE/refgen/reference_genome.fl200.2021_.gtf and exon+flank in refgen_collapsed.exonic_part.gtf
bedtools intersect -a refgen_collapsed.exonic_part.gtf -b /home/dhthutrang/ENCODE/refgen/reference_genome.fl200.2021_.gtf  -loj  > temp.gtf
bedtools groupby -i temp.gtf -g 1-9 -c 12 -o collapse > temp1.gtf
iv. Get only flanks and exon+flank
grep flank temp1.gtf > temp2.gtf
mv temp2.gtf refgen_collapsed.flank.gtf && rm temp*.gtf