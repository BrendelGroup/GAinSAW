from gainsaw import (cfcheck, gsconf, gsparams, gsscrs, BedWrap, LiftOver, PointSet,
                     get_lopset)
from pybedtools import BedTool
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.max_colwidth', None)
pd.set_option('display.width', None)

if cfcheck('test.gainsaw.conf'): gsconf.read('test.gainsaw.conf')


usage1 = '''
***
  First, we load the pickle filtered_mm39_rn7_randomSNPs4ann.lops.pkl generated
  in test07.py. Then we use pybedtools functionality to interset the point
  intervals with genome features read from the requisite genome annotation bed
  files (here we make use of the shortcuts to the annotation files as recorded
  in the configuration file). The output shows the data structures generated to
  classify the SNPs with respect to overlap with genome features, e.g. CDS,
  introns, or RNA genes.

    ps = get_lopset("filtered_mm39_rn7_randomSNPs4ann")
    gbed_file = "mm39"
    qgann = BedWrap(gbed_file)
    gbed_file = "rn7"
    tgann = BedWrap(gbed_file)
    mycrt = gsparams._replace(qlabel='mm39',tlabel='rn7',slop_size=0)
    ps.annotate_lopset(qgann.bed, tgann.bed)

  We then change the filtering criteria to mismatch_rate 8% and re-display the
  qualified points by the stricter criteria:

    mycrt = mycrt._replace(mismatch_rate=8,setpfilter=2)
    ps.filter_pdata(qGenomic_fasta, tGenomic_fasta, mycrt)
    ps.make_beds(mycrt, "q_stringent_mm39randomSNPs_rn7", "mm39",
                        "t_stringent_mm39randomSNPs_rn7", "rn7")
    mycrt = mycrt._replace(slop_size=50)
    ps.align_seqset(qGenomic_fasta, tGenomic_fasta, mycrt)
    ps.write_alignments(gsscrs.lastz)
    ps.annotate_lopset(mycrt, qgann.bed, tgann.bed)

***
'''
print(usage1)


gbed_file = "mm39"
qgann = BedWrap(gbed_file)
gbed_file = "rn7"
tgann = BedWrap(gbed_file)

#We could check on success of the above as follows:
#
#tgann.check_bed()
#tgann.write_bed("tgannBED")


mycrt = gsparams._replace(qlabel='mm39',tlabel='rn7',slop_size=0)
ps = get_lopset("filtered_mm39_rn7_randomSNPs4ann")
print(ps.qbed)

print("\n\nCHECK: Query index:\n")
print(ps.q_index)
print("\n\nCHECK: Target index:\n")
print(ps.t_index)

mycrt = gsparams._replace(qlabel='mm39',tlabel='rn7',slop_size=0)
print("... now extracting genomic features from the annotation files\n")
ps.annotate_lopset(mycrt, qgann.bed, tgann.bed)

print("\n\nCHECK: Query ftype:\n")
print(ps.q_ftype)

print("\n\nCHECK: Unlifted query ftype:\n")
print(ps.u_ftype)

print("\n\nCHECK: Target ftype:\n")
print(ps.t_ftype)

print("\n\nCHECK: Query annotation:\n")
print(ps.q_annotation)

print("\n\nCHECK: Unlifted query annotation:\n")
print(ps.u_annotation)

print("\n\nCHECK: Target annotation:\n")
print(ps.t_annotation)

print("\n\nCHECK: qfeatures:\n")
for point in ps.qfeatures:
    print(f"Query point\t{point['pkey']}\toverlaps feature\t{point['pfeature']}\t{point['pannotation']}")

print("\n\nCHECK: ufeatures:\n")
for point in ps.ufeatures:
    print(f"Unlifted query point\t{point['pkey']}\toverlaps feature\t{point['pfeature']}\t{point['pannotation']}")

print("\n\nCHECK: tfeatures:\n")
for point in ps.tfeatures:
    print(f"Target point\t{point['pkey']}\toverlaps feature\t{point['pfeature']}\t{point['pannotation']}")

print("\n\nAnnotation data frame:\n")
ps.create_annotation_dataframe()
ps.write_annotation_dataframe()


qGenomic_fasta = 'mm39'
tGenomic_fasta = 'rn7'

print("\n\nCHECK: alignments for all non-outofbound points:\n")
ps.align_seqset(qGenomic_fasta, tGenomic_fasta)
ps.write_alignments(gsscrs.lastz)

print("\n\nCHECK: alignments for high-quality points only:\n")
mycrt = mycrt._replace(slop_size=50,usepfilter=2)
ps.align_seqset(qGenomic_fasta, tGenomic_fasta, mycrt)
ps.write_alignments(gsscrs.lastz)

print("\n\nNOW CHANGING TO 8% mismatch_rate:\n")
mycrt = mycrt._replace(mismatch_rate=8,setpfilter=2)
ps.filter_pdata(qGenomic_fasta, tGenomic_fasta, mycrt)
mycrt = mycrt._replace(slop_size=0)
ps.make_beds(mycrt, "q_stringent_mm39randomSNPs_rn7", "mm39",
                    "t_stringent_mm39randomSNPs_rn7", "rn7")
mycrt = mycrt._replace(slop_size=50)
ps.align_seqset(qGenomic_fasta, tGenomic_fasta, mycrt)
ps.write_alignments(gsscrs.lastz)
ps.annotate_lopset(mycrt, qgann.bed, tgann.bed)

print("\n\nCHECK: qfeatures:\n")
for point in ps.qfeatures:
    print(f"Query point\t{point['pkey']}\toverlaps feature\t{point['pfeature']}\t{point['pannotation']}")

print("\n\nCHECK: ufeatures:\n")
for point in ps.ufeatures:
    print(f"Unlifted query point\t{point['pkey']}\toverlaps feature\t{point['pfeature']}\t{point['pannotation']}")

print("\n\nCHECK: tfeatures:\n")
for point in ps.tfeatures:
    print(f"Target point\t{point['pkey']}\toverlaps feature\t{point['pfeature']}\t{point['pannotation']}")

print("\n\nAnnotation data frame for the stringently filtered points:\n")
ps.create_annotation_dataframe()
ps.write_annotation_dataframe()
ps.pickle_annotation_dataframe("stringent_mm39randomSNPs_rn7")
