from gainsaw import (gsconf, gsparams, LiftOver, PointSet)

usage1 = '''
***
  Solution to test03.py exercise:
    mv ../data/liftovers/mm39chr5and7ToRn7.over.chain.lo.pkl ../data/liftovers/mm39chr5and7ToRn7.lo.pkl
    sed -i -e "/\[liftovers\]/a mm39chr5and7ToRn7 = ../data/liftovers/mm39chr5and7ToRn7.lo.pkl" test.gainsaw.conf
***
'''
print(usage1)

gsconf.read('test.gainsaw.conf')
lo = LiftOver('mm39chr5and7ToRn7','../data/liftovers')
ps = PointSet(lo,"mm39.Ctbp.bed")
print()
ps.check_pdata(0)
print()

usage2 = '''
***
  Great. With our conf file updated, the commands
    gsconf.read('test.gainsaw.conf')
    lo = LiftOver('mm39chr5and7ToRn7','../data/liftovers')
    ps = PointSet(lo,"mm39.Ctbp.bed")
    ps.check_pdata(0)
  give us the mm39 Ctbp gene end points mapped to rn7. How well do the
  corresponding genome segments align?
  Our first answer is to re-align the sequences from the alignment blocks in
  the input chainfile that the liftover relied upon. To do this, we load the
  (indexed) genomic sequences and use the following gainsaw methods:
    qGenomic_fasta = 'mm39'
    tGenomic_fasta = rn7'
    mycrt = gsparams._replace(slop_size=0)
    ps.lopset.align_seqset(qGenomic_fasta, tGenomic_fasta, mycrt)
    ps.lopset.write_alignments()
***
'''
print(usage2)

qGenomic_fasta = 'mm39'
tGenomic_fasta = 'rn7'

mycrt = gsparams._replace(slop_size=0)
ps.lopset.align_seqset(qGenomic_fasta, tGenomic_fasta, mycrt)
ps.lopset.write_alignments()

usage3 = '''
***
  The output shows the alignments which should be gap-free (although in some
  cases, there may be alternative alignments with compensating gaps; clearly,
  this depends on the sequences and alignment algorithm).
***
'''
print(usage3)
