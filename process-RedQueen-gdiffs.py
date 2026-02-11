#!/usr/bin/env python

"""
process-RedQueen-gdiffs.py by Rohan Maddamsetti.

I take the evolved genomes, and write out evolved-mutations.csv for
downstream analysis in R.

Usage: python process-RedQueen-gdiffs.py

"""

from os.path import join, dirname, realpath
from os import listdir
import pandas as pd

## add genomediff parser to path.
import sys
sys.path.append("genomediff-python")
import genomediff


def get_mutation_data(pop_labels, day, input_dir, FALSE_POSITIVE_BLACKLIST=[]):

        mutation_data = []
        ## skip hidden files like .DS_Store.
        breseq_results = [x for x in listdir(input_dir) if not x.startswith(".")]
        for sample in breseq_results:
                full_f = join(input_dir, sample, "output", "evidence", "annotated.gd")
                infh = open(full_f, 'r', encoding='utf-8')

                my_row = pop_labels[pop_labels.Sample == sample].iloc[0]
                is_cured = my_row["is_cured"]
                Day = my_row["Day"]
                aTc = str(float(my_row["aTc"]))
                spacer_type = my_row["Spacer_type"]
                colony = str(int(my_row["Colony"]))

                gd = genomediff.GenomeDiff.read(infh)
                muts = gd.mutations
                muttype = ""
                allele = ""
                for rec in muts:
                        pos = str(rec.attributes['position'])
                        mut_category = rec.attributes['mutation_category']
                        gene = rec.attributes['gene_name']
                        ## skip over false positive blacklisted genes.
                        if gene in FALSE_POSITIVE_BLACKLIST: continue
                        if 'new_seq' in rec.attributes:
                                allele = rec.attributes['ref_seq'] + "->" + rec.attributes['new_seq']
                        else:
                                allele = rec.type
                        if rec.type == 'SNP':
                                muttype = rec.attributes['snp_type']
                        else:
                                muttype = rec.type
                        mutation_row_data = ','.join([sample, is_cured, Day, aTc, spacer_type, colony, muttype, mut_category, gene, pos, allele])
                        mutation_data.append(mutation_row_data)
        return mutation_data


def filter_mutation_data(evolved_mutation_data, ancestral_mutation_data):
        """ This solution is not efficient, but gets the job done."""
        filtered_mutation_data = []
        for evolved_row in evolved_mutation_data:
                EV_sample, EV_is_cured, EV_Day, EV_aTc, EV_spacer_type, EV_colony, EV_muttype, EV_mut_category, EV_gene, EV_pos, EV_allele = evolved_row.split(',')
                found_in_ancestor = False
                for ancestral_row in ancestral_mutation_data:
                        ANC_sample, ANC_is_cured, ANC_Day, ANC_aTc, ANC_spacer_type, ANC_colony, ANC_muttype, ANC_mut_category, ANC_gene, ANC_pos, ANC_allele = ancestral_row.split(',')
                        if (EV_spacer_type == ANC_spacer_type) and (EV_colony == ANC_colony) and (EV_muttype == ANC_muttype) and (EV_mut_category == ANC_mut_category) and (EV_gene == ANC_gene) and (EV_pos == ANC_pos) and (EV_allele == ANC_allele):
                                found_in_ancestor = True
                if not found_in_ancestor:
                        filtered_mutation_data.append(evolved_row)
        return filtered_mutation_data


def write_evolved_mutations(pop_labels, outf, FALSE_POSITIVE_BLACKLIST):
    outfh = open(outf,"w")
    outfh.write("Sample,is_cured,Day,aTc,Spacer_type,Colony,Mutation,Mutation_Category,Gene,Position,Allele\n")
    Day1_mutation_data = get_mutation_data(pop_labels, 1, "../data/D1_results", FALSE_POSITIVE_BLACKLIST)
    Day100_mutation_data = get_mutation_data(pop_labels, 100, "../data/D100_results", FALSE_POSITIVE_BLACKLIST)

    ## remove Day 1 mutations from Day 100 genomes.
    filtered_Day100_mutation_data = filter_mutation_data(Day100_mutation_data, Day1_mutation_data)
    for row in filtered_Day100_mutation_data:
        outfh.write(row + "\n")
    return


def main():
    srcdir = dirname(realpath(__file__))
    assert srcdir.endswith("src")
    projdir = dirname(srcdir)
    assert projdir.endswith("Hye-in_RedQueen")

    ## these are in the data directory since I got these results from Hye-in and did not run them myself.
    Day1_results_dir = join(projdir,"data","D1_results")
    Day100_results_dir = join(projdir,"data","D100_results")

    pop_label_f = join(projdir,"data","sample-metadata.csv")  
    pop_labels = pd.read_csv(pop_label_f)

    mut_table_outf = "../results/evolved_mutations.csv"
    FALSE_POSITIVE_BLACKLIST = []
    write_evolved_mutations(pop_labels, mut_table_outf, FALSE_POSITIVE_BLACKLIST)

    
main()
    
