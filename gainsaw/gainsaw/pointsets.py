'''
pointsets.py - from the gainsaw package

Scope
    PointSet takes a LiftOver object and a set of points (from a bed file)
    and stores the query and target point pairs as a LiftOverPointSet object.

Usage
    from gainsaw import (LiftOver, PointSet)
    lo = LiftOver('mm10ToMm39.over.chain.gz','../data/chainfiles')
    ps = PointSet(lo,"mm10.Ctbp.bed")
    ps.check_pdata()

Authors
    Wenshu Chen, email: wenschen@indiana.edu
    Volker Brendel, email: vbrendel@indiana.edu

Homepage 
    https://github.com/BrendelGroup/GAainSAW"

Licensed under MIT license.
'''

import gzip
from .process_pointset import LiftOverPointSet

class PointSet:
    def __init__(self, lo, pointsetfile, data_dir='.'):
        '''
        PointSet arguments:
          lo		LiftOver object
          poinsetfile	a bed-formatted set of points from the query genome;
                         the file may be uncompressed or a gzip-compressed file
          data_dir	(optional) the directory where to look for the pointset
                         file
        '''

        if isinstance(pointsetfile, str):
            do_gzip = 1 if pointsetfile.lower().endswith('.gz') else 0
            if do_gzip:
                f = gzip.open(data_dir+'/'+pointsetfile, 'rb')
            else:
                f = open(data_dir+'/'+pointsetfile, 'rb')
        else:
            print ("PROBLEM")
        self.lopset = LiftOverPointSet(lo, f)
        f.close()

    def check_pdata(self, n = 10):
        '''
        Just checking ...
        '''
        if n == 0:
            for p in self.lopset.pdata.values():
               print(p)
        else:
            i = 1
            for p in self.lopset.pdata.values():
                print(p)
                i = i + 1
                if i > n:
                   break

    def output_lifted_points(self,f):
        '''
        Write the lifted coordinates from the input lopset to the specified file.
        '''
        for p in self.lopset.pdata.values():
          for l in p:
            f.write("{0}\t{1:d}\t{2:d}\n".format(l[3],l[2],l[2]+1))
        f.close()

    def output_alb(self,qoutfile,toutfile):
        '''
        TO BE WRITTEN
        '''
        for p in self.lopset.pdata.values():
           print(p)

