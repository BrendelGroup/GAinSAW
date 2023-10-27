'''
process_chainfile.py - from the gainsaw package

This code is derived from chainfile.py of the pyliftover package:
	https://github.com/konstantint/pyliftover
	Copyright 2013, Konstantin Tretyakov.
	http://kt.era.ee/
Thank you Konstantin!

Scope
    This code processes an input alignment chain file into an interval tree
    that facilitates the construction of a LiftOver object.
    Note: We use the terms "query" for the, well, query point and genome, and
          "target" for, well, the target genome onto which the query point gets
          mapped.

Usage
    The LiftOverChainFile class is imported by liftover.py.

Authors
    Wenshu Chen, email: wenschen@indiana.edu
    Volker Brendel, email: vbrendel@indiana.edu

Homepage 
    https://github.com/BrendelGroup/GAainSAW"

Licensed under MIT license.
'''

from os import path
import gzip
import urllib
import shutil
import sys

from .intervaltree import IntervalTree


class LiftOverChainFile:
    '''
    The class, which loads and indexes USCS's .over.chain files.
    Specification of the chain format can be found here:
      http://genome.ucsc.edu/goldenPath/help/chain.html
    '''
    
    def __init__(self, f, show_progress=False):
        '''
        Reads chain data from the file and initializes an interval index.
        f must be a file object open for reading.
        If any errors are detected, an Exception is thrown.
        
        If show_progress == True, a progress bar is shown in the console.
        Requires tqdm to be installed.
        '''
        self.chains = self._load_chains(f, show_progress)
        self.chain_index = self._index_chains(self.chains, show_progress)
        
    @staticmethod
    def _load_chains(f, show_progress=False):
        '''
        Loads all LiftOverChain objects from a file into an array. Returns the result.
        '''
        chains = []
        if show_progress:
            from tqdm import tqdm
            pbar = tqdm(total = float('inf'), desc="Reading file", unit=" chains")
        while True:
            line = f.readline()
            if not line:
                break
            if line.startswith(b'#') or line.startswith(b'\n') or line.startswith(b'\r'):
                continue
            if line.startswith(b'chain'):
                # Read chain
                chains.append(LiftOverChain(line, f))
                if show_progress:
                    pbar.update(1)
                continue
        return chains

    @staticmethod
    def _index_chains(chains, show_progress=False):
        '''
        Given a list of LiftOverChain objects, creates a
         dict: query_name --> 
            IntervalTree: <query_from, query_to> -->
                (target_from, target_to, chain)
        Returns the resulting dict.
        Throws an exception on any errors or inconsistencies among chains (e.g. different sizes specified for the same chromosome in various chains).
        '''
        chain_index = {}
        query_size = {}
        target_size = {}
        if show_progress:
            from tqdm import tqdm
            chains = tqdm(chains, desc="Indexing", unit=" chains")
        for c in chains:
            # Verify that sizes of chromosomes are consistent over all chains
            query_size.setdefault(c.query_name, c.query_size)
            if query_size[c.query_name] != c.query_size:
                raise Exception("Chains have inconsistent specification of query chromosome size for %s (%d vs %d)" % (c.query_name, query_size[c.query_name], c.query_size))
            target_size.setdefault(c.target_name, c.target_size)
            if target_size[c.target_name] != c.target_size:
                raise Exception("Chains have inconsistent specification of target chromosome size for %s (%d vs %d)" % (c.target_name, target_size[c.target_name], c.target_size))
            chain_index.setdefault(c.query_name, IntervalTree(0, c.query_size))
            # Register all blocks from the chain in the corresponding interval tree
            tree = chain_index[c.query_name]
            for (qfrom, qto, tfrom, tto) in c.blocks:
                tree.add_interval(qfrom, qto, (tfrom, tto, c))

        # Sort all interval trees
        for k in chain_index:
            chain_index[k].sort()
        return chain_index

    def query(self, chromosome, position):
        '''
        Given a chromosome and position, returns all matching records from the chain index.
        Each record is an interval (query_from, query_to, data)
        where data = (target_from, target_to, chain). Note that depending on chain.target_strand, the target values may need to be reversed (e.g. pos --> chain.target_size - pos).
        
        If chromosome is not found in the index, None is returned.
        '''
        # A somewhat-ugly hack to allow both 'bytes' and 'str' objects to be used as
        # chromosome names in Python 3. As we store chromosome names as strings,
        # we'll transparently translate the query to a string too.
        if type(chromosome).__name__ == 'bytes':
            chromosome = chromosome.decode('utf-8')
        if chromosome not in self.chain_index:
            return None
        else:
            return self.chain_index[chromosome].query(position)


class LiftOverChain:
    '''
    Represents a single chain from an .over.chain file.
    A chain basically maps a set of intervals from "query" coordinates to corresponding coordinates in "target" coordinates.
    '''
    __slots__ = ['score', 'query_name', 'query_size', 'query_start', 'query_end',
                 'target_name', 'target_size', 'target_strand',
                 'target_start', 'target_end', 'id', 'blocks']

    def __init__(self, chainheader, f):
        '''
        Reads the chain from a stream given the first line and a file opened at all remaining lines.
        On error throws an exception.
        '''
        chainheader = chainheader.decode('utf-8')
        fields = chainheader.split()
        if fields[0] != 'chain' and len(fields) not in [12, 13]:
            raise Exception("Invalid chain format. (%s)" % chainheader)
        # chain 4900 chrY 58368225 + 25985403 25985638 chr5 151006098 - 43257292 43257528 1
        self.score = int(fields[1])       # Alignment score
        self.query_name = fields[2]       # E.g. chrY
        self.query_size = int(fields[3])  # Full length of the chromosome
        query_strand = fields[4]          # Must be +
        if query_strand != '+':
            raise Exception("Source strand in an .over.chain file must be +. (%s)" % chainheader)
        self.query_start = int(fields[5]) # Start of query region
        self.query_end = int(fields[6])   # End of query region
        self.target_name = fields[7]      # E.g. chr5
        self.target_size = int(fields[8]) # Full length of the chromosome
        self.target_strand = fields[9]    # + or -
        if self.target_strand not in ['+', '-']:
            raise Exception("Target strand must be - or +. (%s)" % chainheader)
        self.target_start = int(fields[10])
        self.target_end = int(fields[11])
        self.id = None if len(fields) == 12 else fields[12].strip()
        
        # Now read the alignment chain from the file and store it as a list (query_from, query_to) -> (target_from, target_to)
        qfrom, tfrom = self.query_start, self.target_start
        self.blocks = []
        fields = f.readline().decode('utf-8').split()
        while len(fields) == 3:
            size, sgap, tgap = int(fields[0]), int(fields[1]), int(fields[2])
            self.blocks.append((qfrom, qfrom+size, tfrom, tfrom+size))
            qfrom += size + sgap
            tfrom += size + tgap
            fields = f.readline().split()
        if len(fields) != 1:
            raise Exception("Expecting one number on the last line of alignments block. (%s)" % chainheader)
        size = int(fields[0])
        self.blocks.append((qfrom, qfrom+size, tfrom, tfrom+size))
        if (qfrom + size) != self.query_end  or (tfrom + size) != self.target_end:
            raise Exception("Alignment blocks do not match specified block sizes. (%s)" % chainheader)
