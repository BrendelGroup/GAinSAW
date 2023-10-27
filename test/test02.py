from gainsaw import (gsconf, LiftOver, PointSet)

gsconf.read('test.gainsaw.conf')
lo = LiftOver('mm10ToMm39','../data/liftovers')
ps = PointSet(lo,"mm10.Ctbp.bed")

usage1='''
***
  Here we have added to the code
    from gainsaw.config import gsconf
    gsconf.read('test.gainsaw.conf')
  and loaded the previously pickled LiftOver object:
    lo = LiftOver('mm10ToMm39','../data/liftovers')

  The output is the same, but the code wil run faster for big objects.
***
'''
print(usage1)

ps.check_pdata(0)
