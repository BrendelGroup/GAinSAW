from gainsaw import (cfcheck, gsconf, gsparams, gsscrs, LiftOver, PointSet)

if cfcheck('test.gainsaw.conf'): gsconf.read('test.gainsaw.conf')

lo = LiftOver('mm39chr5and7ToRn7','../data/liftovers')
ps = PointSet(lo,"mm39.Ctbp-points.bed")
print()
ps.check_pdata(0)
print()

usage1 = '''
***
  Now we take a closer look in the vicinity of the points of interest,
  independent of the alignment blocks given by the chainfile. The following code
  will take a window extending 25 positions to the left and right of the points
  and show the alignment of the 51-nucleotide segments. As input point set, we
  take the previous mm39.Ctbp.bed, augmented by the first annotated CDS segment
  in mm39 Ctbp1.
    qGenomic_fasta = 'mm39'
    tGenomic_fasta = 'rn7'
    ps = PointSet(lo,"mm39.Ctbp-points.bed")
    ps.check_pdata(0)

    ps.lopset.align_seqset(qGenomic_fasta, tGenomic_fasta)
    ps.lopset.write_alignments(gsscrs.lastz)

  Note: Of course, specifiying 'mm39' and 'rn7' above will only work if you have
  stored the corresponding genome sequence indices and added them to your
  test.gainsaw.conf file after having run test04.py ...

  Note: The line
    if cfcheck('test.gainsaw.conf'): gsconf.read('test.gainsaw.conf')
  checks on the existence of our config file before reading it. Good practice!

  Note: The parameter gssscrs.lastz to write_alignments specifies the scoring
    scheme for the alignment as the one used by lastz. To show gap-free
    alignments, you can explicitly use gsscrs.default or leave out the argument
    (as gsscrs.default is the, well,  default). You can also specify the
    parameter in the following way (to keep a record of all your parameter
    choices in mycrt).

    mycrt = gsparams._replace(scoring=gsscrs.lastz)
    ps.lopset.write_alignments(mycrt.scoring)

***
'''
print(usage1)

qGenomic_fasta = 'mm39'
tGenomic_fasta = 'rn7'

print("\nUsing default slop_size=25")
ps.lopset.align_seqset(qGenomic_fasta, tGenomic_fasta)
ps.lopset.write_alignments(gsscrs.lastz)


usage2 = '''
***
  The output shows the alignments of the sequence segments centered on the
  query and target points. What should we expect?

  Ctbp1 start and end coordinates are well within chain file reported
  alignment blocks (alb). Thus, we expect good matching without gaps throughout.

  chr5:33404990 is the start of the first alb. Thus, we should see the last 26
  nucleotides of the segment match and be preceded by a gap.

  chr5:33405261 is the end of the first alb. Thus, we should see the first 26
  nucleotides of the segment match and be followed by a gap.
 
  Exercise: Convice yourself that all alignments match your expectation.

  Note: But beware - alignments depend on context. Thus, the slop size matters.
        Not that there is a magic number for this, but a larger slop size may
        well change the alignment considerably, as new identity blocks may
        shift the position of gaps. This is precisely why gainsaw facilitates
        honing in on points of interest and exploring the local context. The
        input chain file necessarily represents just one genome alignment
        solution - good enough for large-scale comparisons and detection of
        synteny. But when the focus is on specific points in the genome, the
        alignment should be re-scrutinized for alternatives.

  Let's save everything for our next task:
    ps.lopset.pickle_lopset("./test05.lops.pkl")
***
'''
print(usage2)

ps.lopset.pickle_lopset("test05")
