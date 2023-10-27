#!/usr/bin/env python3

# tidyGFF2bed.py
# Convert a tidy-*.gff file to a six-column BED file
# with gene name in place of RNA names in intron and UTR.
# Write to stdout by default. 
# Redirect stdout and sort to save the output.
# Usage:
# python tidyGFF2bed.py tidy-*.gff | sort -u | sort -k1,1 -k2,2n > tidy*.bed
#

import sys

def pair_rna_gene_name(gff3_file):
    '''
    Dictionary of Refseq RNA name:gene name from a NCBI GFF3 file.
    '''
    rna_gene_name = {}
    with open(gff3_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.split("\t")
            if line[2] in ["region", "match"]:
                continue
            if line[8].startswith("ID"):  # attribute
                parent = line[8].split("Parent=")
                if len(parent) == 1:
                    continue
                elif parent[1].split(";")[0][0:3] == "rna":
                    parent = parent[1].split(";")[0]
                elif parent[1].split(";")[0][0:4] == "gene":
                    parent = line[8].split("ID=")[1].split(";")[0]
                elif parent[1].split(";")[0][0:2] == "id":
                    parent = parent[1].split(";")[0].replace("id=", "")
                gene = line[8].split("gene=")
                if len(gene) > 1:
                    gene = gene[1].split(";")[0]
                else:
                    gene = gene[0].split("=")[-1].strip()
                rna_gene_name[parent] = gene                
    return rna_gene_name

def tidy_gff_bed(tidy_gff_file):
    '''
    Convert a tidy-*.gff to a BED file with six cloumns:
    chr, start, end, strand, feature, gene name:
    chr1	6338921	6338980	+	CDS     Rb1cc1
    chr1	6338921	6338980	+	exon	Rb1cc1
    chr1	6338980	6340914	+	intron	Rb1cc1
    '''
    # Build the RNA:gene dictionary.
    rna_gene_name_pair = pair_rna_gene_name(tidy_gff_file)
    # Convert selected columns to BED format.
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.split("\t")
            if line[2] in ["region", "match"]:
                continue
            if line[8].startswith("ID"):
                name = line[8].split("gene=")
                if len(name) > 1:
                    name = name[1].split(";")[0].strip()
                else:
                    name = name[0].split("=")[-1].strip()
            elif line[8].startswith("Parent"):  # canon-gff3 ncbi
                name = line[8].split("=")[1].strip().split(",")  # rna-*
                # Replace rna-* with gene name.
                for i in range(len(name)):  
                    if rna_gene_name_pair.get(name[i]) != None:
                        name[i] = rna_gene_name_pair.get(name[i])
                name = ",".join(set(name))
            print(line[0], int(line[3])-1, line[4], line[6], line[2], name, sep="\t")

gff_file = sys.argv[1]
tidy_gff_bed(gff_file)
