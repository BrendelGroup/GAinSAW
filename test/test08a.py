from gainsaw import (cfcheck, gsconf, gsparams, BedWrap, LiftOver, PointSet,
                     get_lopset)
from pybedtools import BedTool
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.max_colwidth', None)
pd.set_option('display.width', None)

if cfcheck('test.gainsaw.conf'): gsconf.read('test.gainsaw.conf')


usage1 = '''
***
  First, we load the pickle mm39_rn7_randomSNPs4ann.lops.pkl generated in
  test07.py. Then we use pybedtools functionality to interset the point
  intervals with genome features read from the requisite genome annotation bed
  files. The output shows the data structures generated to classify the SNPs
  with respect to overlap with genome features, e.g. CDS, introns, or RNA genes.

    ps = get_lopset("mm39_rn7_randomSNPs4ann")
    gbed_file = "../data/genomes/mm39/mm39.tidyann.bed"
    qgann = BedWrap(gbed_file)
    gbed_file = "../data/genomes/rn7/rn7.tidyann.bed"
    tgann = BedWrap(gbed_file)
    ps.annotate_lopset(qgann.bed, tgann.bed)

  Finally, we put all the derived data into a data frame and illustrate how
  useful that is for analyses. Note, the case of a query pseudogene point mapped 
  to a target CDS, as well as the non-conservation of lncRNA annotations.

    ps.create_annotation_dataframe()
    ps.write_annotation_dataframe()
    ps.df.qftype = ps.df.qftype.apply(lambda l: '/'.join(i for i in l))
    ps.df.tftype = ps.df.tftype.apply(lambda l: '/'.join(i for i in l))
    pd.crosstab(ps.df.qftype,ps.df.tftype, margins=True, margins_name="Total")

  The following command saves the data frame as a pickle:
    ps.pickle_annotation_dataframe("test08a")
***
'''
print(usage1)


gbed_file = "../data/genomes/mm39/mm39.tidyann.bed"
qgann = BedWrap(gbed_file)
gbed_file = "../data/genomes/rn7/rn7.tidyann.bed"
tgann = BedWrap(gbed_file)

#We could check on success of the above as follows:
#
#tgann.check_bed()
#tgann.write_bed("tgannBED")


mycrt = gsparams._replace(qlabel='mm39',tlabel='rn7',slop_size=0)
ps = get_lopset("mm39_rn7_randomSNPs4ann")
ps.check_pdata(0)

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

print("\n\nEnough checking. What we really want is the annotation data frame:\n")
ps.create_annotation_dataframe()
ps.write_annotation_dataframe()

print("\n\nWith data frame procssing we can derive useful statistics, for example paired annotation counts:\n")
print(pd.crosstab(ps.df.qftype,ps.df.tftype, margins=True, margins_name="Total"))
print()

ps.pickle_annotation_dataframe("test08a")
ps.pickle_lopset("mm39_rn7_randomSNPs_annotated")
