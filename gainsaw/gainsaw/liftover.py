'''
liftover.py - from the gainsaw package

This code is adapted from liftover.py of the pyliftover package:
	https://github.com/konstantint/pyliftover
	Copyright 2013, Konstantin Tretyakov.
	http://kt.era.ee/
Thank you Konstantin!

Scope
    The underlying algorithm and data structures (interval tree) are unchanged
    form the pyliftover package. However, the gainsaw implementation record
    additional data from the processed chain file. The basic data unit is an
    arry of tuples like the following:

    [(33247722, 'chr5', 33405066, 'chr5', '+', 30339671, 38300289, 30497015, 38457633, '6', 14025902023, [1])] 
    
    From left to right: query point; query sequence identifier; target point;
    target sequence identifier; target strand; query alignment block start and
    end; target alignment block start and end; chain identifier; chain score;
    and quality flag (1, include in the analyses; 0, exclude from analyses).
    The alignment blocks refer to the ungapped matching segments containing the
    query and target points, respectively. The array will be of length one
    unless the query point maps to multiple targets.

Usage
    from gainsaw import LiftOver
    lo = LiftOver('mm10ToMm39.over.chain.gz', '../data/chainfiles')
    x = lo.convert_coordinate('chr5', 33247722)
    print("mm10 chr5 ", 33247722, " maps to mm39 ", x, "\n")

Authors
    Wenshu Chen, email: wenschen@indiana.edu
    Volker Brendel, email: vbrendel@indiana.edu

Homepage 
    https://github.com/BrendelGroup/GAainSAW"

Licensed under MIT license.
'''

from os import path
import gzip
import pickle

from .config import gsconf
from .process_chainfile import LiftOverChainFile


class LiftOver:
    def __init__(self, chainfile, data_dir='.'):
        '''
        Arguments:
          chainfile	a chain file in either plain or gzip format
          data_dir	(optional) location of where to look for chainfile

        The program saves LiftOver objects as pickles and seeks first to load
        previously saved pickles rather than build the object de novo. Thus,
	  lo = LiftOver('mm39chr5and7ToRn7')
        will look for a 'liftover' entry mm39chr5and7ToRn7 in the gainsaw config
        file and load the corresponding pickle.
          lo = LiftOver('mm39chr5and7ToRn7.lo.pkl','./')
        will load ./mm39chr5and7ToRn7.lo.pkl. And
          lo = LiftOver('mm39chr5and7ToRn7.over.chain','../data/chainfiles')
        will build a LiftOver object from
        ../data/chainfiles/mm39chr5and7ToRn7.over.chain[.gz] for current use, as
        well as save <>.lo.pkl (where <> is the input file name).
        '''

        if not isinstance(chainfile, str):
            print("not a string")
            sys.exit()
        if gsconf.has_option('liftovers', chainfile) and path.exists(f := gsconf['liftovers'][chainfile]):
            print("... loading existing LiftOver object ", f)
            self.lo = pickle.load(open(f, "rb", -1))
        else:
            load_lo = 1 if chainfile.lower().endswith('.lo.pkl') else 0
            if load_lo:
                print("... loading existing LiftOver object ", f := data_dir+'/'+chainfile)
                self.lo = pickle.load(open(f, "rb", -1))
                print("... done")
            else:
                do_gzip = 1 if chainfile.lower().endswith('.gz') else 0
                if do_gzip:
                    f = gzip.open(data_dir+'/'+chainfile, 'rb')
                else:
                    f = open(data_dir+'/'+chainfile, 'rb')
                print("... building LiftOver object for input ", data_dir+'/'+chainfile)
                self.lo = LiftOverChainFile(f)
                f.close()
                print("... done")
                print("... saving LiftOver object ", data_dir+'/'+chainfile+'.lo.pkl')
                with open(data_dir+'/'+chainfile+'.lo.pkl', "wb") as file_:
                    pickle.dump(self.lo, file_, -1)
                print("... done")
                print("Consider adding the pickle to your gainsaw.conf file.\n")


    def convert_coordinate(self, chromosome, position, strand='+'):
        '''
        If chromosome is not found, None is returned.
        Otherwise, a list of mapped coordinates is recorded as documented above.
        '''

        query_results = self.lo.query(chromosome, position)
        if query_results is None:
            return None
        else:
            results = []
            for (query_start, query_end, data) in query_results:
                target_start, target_end, chain = data
                result_position = target_start + (position - query_start)
                result_name = chain.target_name
                result_strand = chain.target_strand
                if result_strand == '-':
                    result_position = chain.target_size - result_position - 1
                    results.append((position, chromosome, result_position, result_name, result_strand,
                                    query_start, query_end,
                                    chain.target_size - target_end, chain.target_size - target_start,
                                    chain.id, chain.score, [1]))
                else:
                    results.append((position, chromosome, result_position, result_name, result_strand,
                                    query_start, query_end, target_start, target_end, chain.id, chain.score, [1]))
            if len(results) > 1:
                results.sort(key=lambda x: x[-2], reverse=True)
            return results
