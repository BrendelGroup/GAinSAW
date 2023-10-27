from gainsaw import (gsconf, LiftOver, PointSet)

gsconf.read('test.gainsaw.conf')
lo = LiftOver('mm10ToMm39','../data/liftovers')
ps = PointSet(lo,"mm10.Ctbp.bed")


usage1 = '''
  For a more useful example of using PointSet methods, we'll save the mm39
  coordinates from our first liftover and then lift those over to rat rn7:
    f = open("lifted.mm10.Ctbp.bed","w")
    ps.output_lifted_points(f)
    lo2 = LiftOver('mm39chr5and7ToRn7.over.chain','../data/liftovers')
    ps2 = PointSet(lo2,"lifted.mm10.Ctbp.bed")
    ps2.check_pdata(0)
'''
print(usage1)
f = open("lifted.mm10.Ctbp.bed","w")
ps.output_lifted_points(f)
lo2 = LiftOver('mm39chr5and7ToRn7.over.chain','../data/liftovers')
ps2 = PointSet(lo2,"lifted.mm10.Ctbp.bed")
ps2.check_pdata(0)

usage2 = '''
  Done. We now have the coordinates of the Ctbp genes for both mm39 and rn7,
  starting from the known locations in mm10.

  Exercise: Add a liftovers entry for mm39chr5and7ToRn7 to test.gainsaw.conf
            for future use of the object.
'''
print(usage2)
