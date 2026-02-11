#!/usr/bin/zsh

## run_paftools.sh by Rohan Maddamsetti.
## align complete whole genome assemblies and call variants using minimap2 and paftools.js.
## See README:
## https://github.com/lh3/minimap2/blob/master/misc/README.md#calling-variants-from-haploid-assemblies

## calling variants from asm-to-ref alignment

## D70 WT
## keeping the PAF file is recommended; --cs required!
minimap2 -cx asm5 -t8 --cs ../data/GK2MVG_results/GK2MVG_3_WT_D0/GK2MVG_3_WT_D0_1_reference.fasta ../data/GK2MVG_results/GK2MVG_2_WT_D70/GK2MVG_2_WT_D70_1_reference.fasta > ../results/GK2MVG_2_WT_D70.paf

## sort by reference start coordinate
sort -k6,6 -k8,8n ../results/GK2MVG_2_WT_D70.paf > ../results/GK2MVG_2_WT_D70.srt.paf

paftools.js call -f ../data/GK2MVG_results/GK2MVG_3_WT_D0/GK2MVG_3_WT_D0_1_reference.fasta -L1000 ../results/GK2MVG_2_WT_D70.srt.paf > ../results/GK2MVG_2_WT_D70.vcf

## convert to gdtools format.
gdtools VCF2GD -o ../results/GK2MVG_2_WT_D70.gd ../results/GK2MVG_2_WT_D70.vcf

## and now annotate as HTML
gdtools ANNOTATE -o ../results/GK2MVG_2_WT_D70.html -r ../data/GK2MVG_results/GK2MVG_3_WT_D0/GK2MVG_3_WT_D0_1_reference.gbk ../results/GK2MVG_2_WT_D70.gd

## annotate as annotated GD
gdtools ANNOTATE -f GD -o ../results/annotated_GK2MVG_2_WT_D70.gd -r ../data/GK2MVG_results/GK2MVG_3_WT_D0/GK2MVG_3_WT_D0_1_reference.gbk ../results/GK2MVG_2_WT_D70.gd

## D70 Cured
minimap2 -cx asm5 -t8 --cs ../data/GK2MVG_results/GK2MVG_3_WT_D0/GK2MVG_3_WT_D0_1_reference.fasta ../data/GK2MVG_results/GK2MVG_1_Cured_D70/GK2MVG_1_Cured_D70_1_reference.fasta > ../results/GK2MVG_1_Cured_D70.paf

## sort by reference start coordinate
sort -k6,6 -k8,8n ../results/GK2MVG_1_Cured_D70.paf > ../results/GK2MVG_1_Cured_D70.srt.paf

paftools.js call -f ../data/GK2MVG_results/GK2MVG_3_WT_D0/GK2MVG_3_WT_D0_1_reference.fasta -L1000 ../results/GK2MVG_1_Cured_D70.srt.paf > ../results/GK2MVG_1_Cured_D70.vcf

## convert to gdtools format.
gdtools VCF2GD -o ../results/GK2MVG_1_Cured_D70.gd ../results/GK2MVG_1_Cured_D70.vcf

## and now annotate as HTML
gdtools ANNOTATE -o ../results/GK2MVG_1_Cured_D70.html -r ../data/GK2MVG_results/GK2MVG_3_WT_D0/GK2MVG_3_WT_D0_1_reference.gbk ../results/GK2MVG_1_Cured_D70.gd

## annotate as annotated GD
gdtools ANNOTATE -f GD -o ../results/annotated_GK2MVG_1_Cured_D70.gd -r ../data/GK2MVG_results/GK2MVG_3_WT_D0/GK2MVG_3_WT_D0_1_reference.gbk ../results/GK2MVG_1_Cured_D70.gd


