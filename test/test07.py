from gainsaw import (cfcheck, gsconf, gsparams, BedWrap, LiftOver, PointSet,
                     get_lopset)
from pybedtools import BedTool

if cfcheck('test.gainsaw.conf'): gsconf.read('test.gainsaw.conf')


usage1 = '''
***
  Here we take a random sample of 50 SNP positions from mm39, lift them to rn7,
  and create the corresponding bed files q_mm39randomSNPs_rn7.bed and
  t_mm39randomSNPs_rn7.bed.
  We save all the data in the pickle mm39_rn7_randomSNPs4ann.lops.pkl for use in
  test08a.py.

    lo = LiftOver('mm39chr5and7ToRn7')
    ps = PointSet(lo,"mm39randomSNPs.bed")
    mycrt = gsparams._replace(qlabel='mm39',tlabel='rn7',slop_size=0)
    ps.lopset.make_beds(mycrt,
        "q_mm39randomSNPs_rn7", "../data/genomes/mm39/mm39.chrom.sizes",
        "t_mm39randomSNPs_rn7", "../data/genomes/rn7/rn7.chrom.sizes")
    ps.lopset.pickle_lopset("mm39_rn7_randomSNPs4ann")

  We then apply filter_pdata() to select only high-quality points for further
  study. The screen is for points with 50nt extensions  that are (1) within the
  bounds of their respective chain file alignment blocks, and (2) have a
  mismatch rate of no more than 10% (the default value in gsparams). The
  filtered points are saved in the pickle
  filtered_mm39_rn7_randomSNPs4ann.lops.pkl.

  We can retrieve subsets of points as follows (selected by usepfilter):

    ps = get_lopset("filtered_mm39_rn7_randomSNPs4ann")
    mycrt = mycrt._replace(usepfilter=1)
    ps.get_filtered_pdata(mycrt)
    ps.check_pdata(0)

***
'''
print(usage1)

lo = LiftOver('mm39chr5and7ToRn7')
ps = PointSet(lo,"mm39randomSNPs.bed")
ps.check_pdata(0)

print("\nPickeling a lopset with all input points")
print("  [output: mm39_rn7_randomSNPs4ann.lo.pkl]")
mycrt = gsparams._replace(qlabel='mm39',tlabel='rn7',slop_size=0)
ps.lopset.make_beds(mycrt, "q_mm39randomSNPs_rn7", "../data/genomes/mm39/mm39.chrom.sizes",
                    "t_mm39randomSNPs_rn7", "../data/genomes/rn7/rn7.chrom.sizes")
ps.lopset.pickle_lopset("mm39_rn7_randomSNPs4ann")

print("\nNow filtering pdata for high quality context")
print("  [output: filtered_mm39_rn7_randomSNPs4ann.lo.pkl]")

qGenomic_fasta = 'mm39'
tGenomic_fasta = 'rn7'

mycrt = mycrt._replace(setpfilter=2, usepfilter=2, slop_size=50)
ps.lopset.filter_pdata(qGenomic_fasta, tGenomic_fasta, mycrt)
ps.check_pdata(0)
mycrt = mycrt._replace(slop_size=0)
ps.lopset.make_beds(mycrt, "q_filtered_mm39randomSNPs_rn7", "../data/genomes/mm39/mm39.chrom.sizes",
                    "t_filtered_mm39randomSNPs_rn7", "../data/genomes/rn7/rn7.chrom.sizes")
ps.lopset.pickle_lopset("filtered_mm39_rn7_randomSNPs4ann")

print("pdata of accepted points: \n")
ps.lopset.get_filtered_pdata(mycrt)
ps.check_pdata(0)

print("\npdata of points with higher than set mismatch rate: \n")
ps = get_lopset("filtered_mm39_rn7_randomSNPs4ann")
mycrt = mycrt._replace(usepfilter=1)
ps.get_filtered_pdata(mycrt)
ps.check_pdata(0)

print("\npdata of points out of bounds: \n")
ps = get_lopset("filtered_mm39_rn7_randomSNPs4ann")
mycrt = mycrt._replace(usepfilter=0)
ps.get_filtered_pdata(mycrt)
ps.check_pdata(0)
print()
