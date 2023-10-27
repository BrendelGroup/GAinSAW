from gainsaw import LiftOver
from gainsaw import PointSet
from gainsaw.config import conf

###conf.read('my.gainsaw.conf')

lo = LiftOver('AtoB.over.chain','.')
ps = PointSet(lo,"test.bed")

pdata = ps.lopset.pdata.values()
ps.check_pdata()

qGenomic_fasta = './A.fna'
tGenomic_fasta = './B.fna'

criteria = ()

print("\nUsing slop=0")
slop = 0
ps.lopset.get_seqset(qGenomic_fasta, tGenomic_fasta, slop, criteria)
ps.lopset.write_alignments()

print("\nUsing slop=4")
slop = 4
ps.lopset.get_seqset(qGenomic_fasta, tGenomic_fasta, slop, criteria)
ps.lopset.write_alignments()
