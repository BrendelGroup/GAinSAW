from gainsaw import (LiftOver, PointSet)

usage1 = '''
***
  The file mm10.Ctbp.bed shows the locations of the two Ctbp genes in mouse
  mm10. For Ctbp1, start and end points are given on two separate lines. For
  Ctbp2, start and end are on one line.  The commands
    lo = LiftOver('mm10ToMm39.over.chain.gz','../data/liftovers')
    ps = PointSet(lo,"mm10.Ctbp.bed")
  create, first, the LiftOver object lo from the input chain file, and, second,
  the PointSet object ps (with slot ps.lopset), derived from lo and the input
  bed file.
***
'''
print(usage1)
lo = LiftOver('mm10ToMm39.over.chain.gz','../data/liftovers')
ps = PointSet(lo,"mm10.Ctbp.bed")


usage2 = '''
***
  ps.lopset.pdata has the records of the input and lifted coordinates, which
  can be checked by:
    print(ps.lopset.pdata)
***
'''
print(usage2)
print(ps.lopset.pdata)

usage3 = '''
***
  More clearly:
    ps.check_pdata(0)
***
'''
print(usage3)
ps.check_pdata(0)


usage4 = '''
***
  To follow up on the program's suggestion, we can used the saved LiftOver
  object in a subsequent run of the script. To do this, we would need to
  replace the line
    lo = LiftOver('mm10ToMm39.over.chain.gz','../data/liftovers')
  with
    lo = LiftOver('mm10ToMm39.over.chain.gz.lo.pkl','../data/liftovers')
  Try this on your own.

  To simplify further, we can rename the pickle to something shorter and add an
  entry to our test.gainsaw.conf file:
    cp default.gainsaw.conf test.gainsaw.conf
    mv ../data/liftovers/mm10ToMm39.over.chain.gz.lo.pkl ../data/liftovers/mm10ToMm39.lo.pkl
    sed -i -e "/\[liftovers\]/a mm10ToMm39 = ../data/liftovers/mm10ToMm39.lo.pkl" test.gainsaw.conf

  Please do this and then proceed with test02.py.
***
'''
print(usage4)
