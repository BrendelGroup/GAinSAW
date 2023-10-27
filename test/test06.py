from gainsaw import (cfcheck, gsconf, gsparams, BedWrap, LiftOver, PointSet, get_lopset)

if cfcheck('test.gainsaw.conf'): gsconf.read('test.gainsaw.conf')

usage1 = '''
***
  With test05.py we explored the context of query/target points in windows of
  length 51 centered on the points. After looking at the results, we may want
  to change the context size or alignment parameters. Using the saved
  LiftOverPointSet pickle, this is easily done:

    psl = get_lopset("test05")
    psl.write_alignments()

    qGenomic_fasta = 'mm39'
    tGenomic_fasta = 'rn7'
    mycrt = gsparams._replace(qlabel='mm39',tlabel='rn7',slop_size=50)
    psl.align_seqset(qGenomic_fasta, tGenomic_fasta, mycrt)
    psl.write_alignments()

  This code loads the saved pickle; writes the previous alignments (with
  slop_size 25); and generates and displays alignments with slop_size 50.
***
'''
print(usage1)

psl = get_lopset("test05")
psl.write_alignments()

qGenomic_fasta = 'mm39'
tGenomic_fasta = 'rn7'

print("\n\nShowing the default parameters: ", gsparams)
print("\nUsing slop_size=50 instead ...")
mycrt = gsparams._replace(qlabel='mm39',tlabel='rn7',slop_size=50)
print("Showing my changed parameters: ", mycrt, "\n")
psl.align_seqset(qGenomic_fasta, tGenomic_fasta, mycrt)
psl.write_alignments()
