'''
process_pointset.py - from the gainsaw package

Scope
    A LiftOverPointSet captures features of a set of query/target points.

Usage
    The LiftOverPointSet class is imported by process_pointset.py.

Authors
    Wenshu Chen, email: wenschen@indiana.edu
    Volker Brendel, email: vbrendel@indiana.edu

Homepage 
    https://github.com/BrendelGroup/GAainSAW"

Licensed under MIT license.
'''

import sys
import datetime
from os import path
import gzip
import pickle
import pandas as pd
from statistics import mean, median
from collections import Counter
from pybedtools import BedTool
from .config import (gsconf, gsparams)
from .align import get_segment_seqset, get_extended_point_seqset, align_segments, string_mismatch

from .align import (get_fasta_index, get_fasta)


class LiftOverPointSet:
    '''
    '''
    
    def __init__(self, lo, f):
        '''
        Arguments:
          lo	LiftOver object
          f	file pointer to a (pointset) file open for reading
        '''
        self.points = self._load_points(f)
        self.pdata  = self._liftover_points(self.points, lo)
        
    @staticmethod
    def _load_points(f):
        '''
        This function loads the points from file f into the array points.
        '''
        points = []
        while True:
            line = f.readline()
            if not line:
                break
            if line.startswith(b'#') or line.startswith(b'\n') or line.startswith(b'\r'):
                continue
            p = line.decode('utf-8').split()
            points.append(tuple([p[0],int(p[1]),int(p[2])-1]))
        return points

    @staticmethod
    def _liftover_points(points,lo):
        '''
        This funtion lifts points according to the LiftOver object lo and
         returns the dictionary point_index.
        '''
        point_index = {}
        for p in points:
            x = lo.convert_coordinate(p[0],p[1])
            point_index[p[0]+'_'+str(p[1])] = x
            ### Note: if lo is None, x will be []
            if p[2] > p[1]:
                x = lo.convert_coordinate(p[0],p[2])
                point_index[p[0]+'_'+str(p[2])] = x

        # Sort UNCLEAR - WHAT ARE WE SORTING ON
        for k in point_index:
            if point_index[k] != None: point_index[k].sort()
        return point_index

    def check_pdata(self, n = 10):
        '''
        Just checking ...
        '''
        if n == 0:
            for p in self.pdata.values():
               print(p)
        else:
            i = 1
            for p in self.pdata.values():
                print(p)
                i = i + 1
                if i > n:
                   break

    def filter_pdata(self, qGenomic_fasta, tGenomic_fasta, params = gsparams):
        '''
        
        '''
        slop_size = params.slop_size
        points_accepted = 0
        points_outofbounds = 0
        points_mismatched = 0
        qfasta_index = get_fasta_index(qGenomic_fasta)
        tfasta_index = get_fasta_index(tGenomic_fasta)
        for i in self.pdata.values():
            if i == None:
                continue
            for j in i:
                if j[0]-slop_size < j[5] or j[0]+slop_size+1 > j[6] or j[2]-slop_size < j[7] or j[2]+slop_size+1 > j[8]:
                    j[11][0] = 0       # out of bounds
                    points_outofbounds += 1
                else:
                    q_seq = get_fasta(qfasta_index, j[1], j[0]-slop_size, j[0]+slop_size+1, "+")
                    t_seq = get_fasta(tfasta_index, j[3], j[2]-slop_size, j[2]+slop_size+1, j[4])
                    mm = string_mismatch(q_seq, t_seq)
                    if mm <= params.mismatch_rate:
                        j[11][0] = params.setpfilter   # few mismatches
                        points_accepted += 1
                    else:
                        j[11][0] = 1
                        points_mismatched += 1
        print(f"\npoints out of bounds: {points_outofbounds}")
        print(f"points mismatched   : {points_mismatched}")
        print(f"points accepted     : {points_accepted}\n")

    def get_filtered_pdata(self, params = gsparams, mode = "equal"):
        '''
        Get the subset of points with a filtered label after running filter function.
        All existing attributes will be retained.
        '''
        filtered_pdata = {}
        label = params.usepfilter
        for k, i in self.pdata.items():
            if i == None or i == []:
                continue
            if mode == 'equal':
                for j in [g for g in i if g[-1][0] == label ]:
                    if k in filtered_pdata:
                        filtered_pdata[k].append(j)
                    else:
                        filtered_pdata[k] = [j]
            elif mode == 'at_least':
                for j in [g for g in i if g[-1][0] >= label ]:
                    if k in filtered_pdata:
                        filtered_pdata[k].append(j)
                    else:
                        filtered_pdata[k] = [j]
            elif mode == 'at_most':
                for j in [g for g in i if g[-1][0] <= label ]:
                    if k in filtered_pdata:
                        filtered_pdata[k].append(j)
                    else:
                        filtered_pdata[k] = [j]
            else:
                continue
        self.pdata = filtered_pdata

    def pdata_stats(self):
        '''
        Return basic stats of a set of pdata values.
        '''
        empty_target = 0
        sl = []
        for i in self.pdata.values():
            if i == None:
                continue
            elif i == []:
                empty_target += 1
            else:
                for j in [g for g in i if g[-1][0] > 0]:
                    sl.append((j[6]-j[5]))                
        print(f"Target returned: {len(self.pdata)-empty_target}")
        print(f"Query without a target: {empty_target}")
        print(f"Median segment length: {int(median(sl))}")
        print(f"Mean segment length: {int(mean(sl))}")
        print(f"Minimum segment length: {min(sl)}")
        print(f"Max segment length: {max(sl)}")  

    def align_seqset(self, qGenomic_fasta, tGenomic_fasta, params = gsparams):
        '''
        '''
        pdata = self.pdata.values()
        slop_size = params.slop_size
        if slop_size == 0:
            self.sequences = get_segment_seqset(pdata, params, qGenomic_fasta, tGenomic_fasta)          
        else:
            self.sequences = get_extended_point_seqset(pdata, params, qGenomic_fasta, tGenomic_fasta, slop_size)

    def write_alignments(self, scoring = gsparams.scoring):
        '''
        '''
        for i in self.sequences:
            seq_names, qs, ts, strand = i
            # align(ts, qs) because the print labels "target", "query".
            print()
            print("*** Alignment(s) ***")
            print(seq_names, "(query|target)", sep="")
            print("Strand = +/", strand, "(query|target)", sep="")
            align_segments(ts, qs, strand, scoring)

    def make_beds(self, params = gsparams,
                  query = 'qpnts', qgsizef = 'unspecified',
                  target = 'tpnts', tgsizef = 'unspecified',
                  data_dir = '.'):
        '''
        Save query and target point contexts as BedTool objects.
          data_dir	(optional) the directory where to look for the
                         genome (size) files.
        '''
        slop_size = params.slop_size
        if gsconf.has_option(qgsizef,'chr') and path.exists(f := gsconf[qgsizef]['chr']):
            print("... using query genome (size) file", f)
            qgsizef = f
        elif path.exists(f := data_dir+'/'+qgsizef):
            print("... using query genome (size) file", f)
            qgsizef = f
        else:
            sys.stderr.write(f"Problem: specified query genome (size) file {f} not found.")
            sys.exit()
        if gsconf.has_option(tgsizef,'chr') and path.exists(f := gsconf[tgsizef]['chr']):
            print("... using target genome (size) file", f)
            tgsizef = f
        elif path.exists(f := data_dir+'/'+tgsizef):
            print("... using target genome (size) file", f)
            tgsizef = f
        else:
            sys.stderr.write(f"Problem: specified target genome (size) file {f} not found.")
            sys.exit()
        q_points, u_points, t_points =  [], [], []
        self.q_index, self.t_index = {}, {}
        pnic = {}   # points on sequences not in the chain file
        for k, i in self.pdata.items():
            if i == None:
                chr, pos = k.rsplit("_",1)
                pnic[chr] = pnic.get(chr, 0) + 1
                continue
            if i == []:
                chr, pos = k.rsplit("_",1)
                pos = int(pos)
                u_points.append((chr,pos,pos+1,'+'))
                continue
            # [(9999999, 'chr1', 9988320, 'chr1', '+', 5396899, 12954106, 5385220, 12942427, '2',21057807908, [1])]
            for j in [g for g in i if g[-1][0] >= params.usepfilter]:
                q_points.append((j[1], j[0], j[0]+1, '+'))
                t_points.append((j[3], j[2], j[2]+1, j[4]))
                qkey = params.qlabel+'_'+j[1]+'_'+str(j[0])+'+'
                tkey = params.tlabel+'_'+j[3]+'_'+str(j[2])+j[4]
                self.q_index.setdefault(qkey, []).extend([tkey])
                self.t_index.setdefault(tkey, []).extend([qkey])
        BedTool(q_points).slop(b=slop_size, g=qgsizef).sort().saveas('accepted_'+query+'.bed')
        self.qbed = BedTool('accepted_'+query+'.bed')
        BedTool(u_points).slop(b=slop_size, g=qgsizef).sort().saveas('unlifted_'+query+'.bed')
        self.ubed = BedTool('unlifted_'+query+'.bed')
        BedTool(t_points).slop(b=slop_size, g=tgsizef).sort().saveas('accepted_'+target+'.bed')
        self.tbed = BedTool('accepted_'+target+'.bed')
        print("...", query+'.bed', " unlifted_"+query+'.bed', " and ", target+'.bed', " written")
        if pnic != {}:
            for chr in pnic:
                sys.stderr.write(f"\nWarning: The input contains {pnic[chr]} points on {chr}.\n")
                sys.stderr.write(f" However, there are no alignment chains for {chr}.")
            sys.stderr.write("\nPlease check whether this is the expected result.\n\n")

    def pickle_lopset(self,lopset_label):
        '''
        Build a lopset pickle.
        '''
        print("... building a pickle for ", lopset_label)
        with open(lopset_label+'.lops.pkl', "wb") as file_:
            pickle.dump(self, file_, -1)
        print("... done")
        print("Consider adding the pickle to your gainsaw.conf file.\n")

    def get_feature_type(self, feature):
        if feature.startswith('pseudogene'):
            return 'pseudogene'
        elif feature.startswith('CDS'):
            return 'CDS'
        elif feature.startswith('five_prime_UTR'):
            return 'five_prime_UTR'
        elif feature.startswith('three_prime_UTR'):
            return 'three_prime_UTR'
        elif feature == 'intron':
            return 'intron'
        elif feature != 'mRNA' and 'RNA' in feature:
            return 'RNA'
        else:
            return 'intergenic'
    
    def annotate_points(self, params, what, annotation):
        # Check for overlap with query intervals
        annotated_points = []
        if what == 'qbed':
            points = self.qbed
            label = params.qlabel
        elif what == 'ubed':
            points = self.ubed
            label = params.qlabel
        else:
            points = self.tbed
            label = params.tlabel

        p = 0
        a = 0
        for point in points:
            p = p+1
            point_details = {
                'pkey': f"{label}_{point.chrom}_{point.start}{point[3]}",
                'pfeature': ['intergenic'],
                'pannotation': []
            }
            for entry in annotation[a:]:
                if point.chrom != entry.chrom:
                    break
                elif point.start < entry.start:
                    break
                elif point.end < entry.end:
                    break
                else:
                    a = a+1
                    if a%500 == 0:
                        now = datetime.datetime.now()
                        print("... ",a," annotations on ",p,what," points added by\t", now.strftime("%Y-%m-%d %H:%M:%S"), file=sys.stderr)
                        sys.stdout.flush()
                    feature_type = self.get_feature_type(entry.fields[8])
                    if point_details['pfeature'] == ['intergenic']:
                        point_details['pfeature'] = [feature_type]
                        point_details['pannotation'] = [entry.fields[4:]]
                    elif not feature_type in point_details['pfeature']:
                        point_details['pfeature'].append(feature_type)
                        point_details['pannotation'].append(entry.fields[4:])
                    if a == len(annotation):
                        break
                    elif annotation[a].start != entry.start or \
                       annotation[a].end   != entry.end   or \
                       annotation[a].chrom != entry.chrom :
                        break	# ... next point
            if len(point_details['pfeature']) > 1:
                point_details['pfeature'], point_details['pannotation'] = \
                  zip(*sorted(zip(point_details['pfeature'], point_details['pannotation'])))
            annotated_points.append(point_details)
        return annotated_points
    
    def annotate_lopset(self, params, qgann, tgann):
        now = datetime.datetime.now()
        print ("Current date and time entering annotate_lopset    :", now.strftime("%Y-%m-%d %H:%M:%S"))
        sys.stdout.flush()

        # Derive the feature intervals as BedTool objects for the types of interest
        feature_types = ['pseudogene', 'CDS', 'five_prime_UTR', 'three_prime_UTR', 'RNA', 'intron']
        feature_intervals = BedTool([feature for feature in qgann if self.get_feature_type(feature.fields[4]) in feature_types]).sort()

        now = datetime.datetime.now()
        print ("Current date and time done with qfeature_intervals:", now.strftime("%Y-%m-%d %H:%M:%S"))
        sys.stdout.flush()

        annotation = self.qbed.intersect(feature_intervals, sorted=True, wa=True, wb=True)

        now = datetime.datetime.now()
        print ("Current date and time done with qbed intersection :", now.strftime("%Y-%m-%d %H:%M:%S"))
        sys.stdout.flush()

        self.qfeatures = self.annotate_points(params, 'qbed', annotation)

        now = datetime.datetime.now()
        print ("Current date and time done with qfeatures         :", now.strftime("%Y-%m-%d %H:%M:%S"))
        sys.stdout.flush()


        annotation = self.ubed.intersect(feature_intervals, sorted=True, wa=True, wb=True)

        now = datetime.datetime.now()
        print ("Current date and time done with ubed intersection :", now.strftime("%Y-%m-%d %H:%M:%S"))
        sys.stdout.flush()

        self.ufeatures = self.annotate_points(params, 'ubed', annotation)

        now = datetime.datetime.now()
        print ("Current date and time done with ufeatures         :", now.strftime("%Y-%m-%d %H:%M:%S"))
        sys.stdout.flush()


        feature_intervals = BedTool([feature for feature in tgann if self.get_feature_type(feature.fields[4]) in feature_types]).sort()

        now = datetime.datetime.now()
        print ("Current date and time done with feature_intervals :", now.strftime("%Y-%m-%d %H:%M:%S"))
        sys.stdout.flush()

        annotation = self.tbed.intersect(feature_intervals, sorted=True, wa=True, wb=True)

        now = datetime.datetime.now()
        print ("Current date and time done with tbed intersection :", now.strftime("%Y-%m-%d %H:%M:%S"))
        sys.stdout.flush()

        self.tfeatures = self.annotate_points(params, 'tbed', annotation)

        now = datetime.datetime.now()
        print ("Current date and time done with tfeatures         :", now.strftime("%Y-%m-%d %H:%M:%S"))
        sys.stdout.flush()

        self.q_ftype, self.u_ftype, self.t_ftype = {}, {}, {}
        self.q_annotation, self.u_annotation, self.t_annotation = {}, {}, {}
        for point_details in self.qfeatures:
            self.q_ftype.setdefault(point_details['pkey'], []).extend(point_details['pfeature'])
            self.q_annotation.setdefault(point_details['pkey'], []).extend(point_details['pannotation'])
        for point_details in self.ufeatures:
            self.u_ftype.setdefault(point_details['pkey'], []).extend(point_details['pfeature'])
            self.u_annotation.setdefault(point_details['pkey'], []).extend(point_details['pannotation'])
        for point_details in self.tfeatures:
            self.t_ftype.setdefault(point_details['pkey'], []).extend(point_details['pfeature'])
            self.t_annotation.setdefault(point_details['pkey'], []).extend(point_details['pannotation'])

    def create_annotation_dataframe(self):
        '''
        Get data from annotation results into a pandas data frame.
        '''
        # Lifted point data frame, self.df
        def _extract_qgann(row):
            if row['qgann'] != []:
                return [','.join(i)  for i in zip(*row['qgann'])]
            else:
                return [None for i in range(6)]

        def _extract_tgann(row):
            if row['tgann'] != []:
                return [','.join(i)  for i in zip(*row['tgann'])]
            else:
                return [None for i in range(6)]

        data_rows = []
        for qkey, tkeys in self.q_index.items():
            for i, tkey in enumerate(tkeys, start=1):
                row_data = {
                    'qkey': f"{qkey}_{i}",
                    'qftype': ','.join(self.q_ftype[qkey]),
                    'tkey': tkey,
                    'tftype': ','.join(self.t_ftype[tkey]),
                    'qgann': self.q_annotation[qkey],
                    'tgann': self.t_annotation[tkey]
                }
                data_rows.append(row_data)
        self.df = pd.DataFrame(data_rows)
        qgann_columns = ['qgann_chr', 'qgann_from', 'qgann_to', 'qgann_strand', 'qgann_ftype', 'qgann_gene']
        self.df[qgann_columns] = self.df.apply(_extract_qgann, axis=1, result_type='expand')
        tgann_columns = ['tgann_chr', 'tgann_from', 'tgann_to', 'tgann_strand', 'tgann_ftype', 'tgann_gene']
        self.df[tgann_columns] = self.df.apply(_extract_tgann, axis=1, result_type='expand')
        self.df.drop(columns=['qgann', 'tgann'], inplace=True)

        # unlifted point data frame, self.udf
        unlifted_annotation = []
        for point in self.ufeatures:
            feature=','.join(point['pfeature'])
            if point['pannotation'] == []:
                annotation = ("None " * 6).split()
            else:
                annotation = point['pannotation'][0]   
            unlifted_annotation.append([point['pkey']] + [feature] + annotation)
        self.udf = pd.DataFrame(unlifted_annotation, columns = ['ukey', 'uftype', 'ugann_chr', 'ugann_from', 'ugann_to', 'ugann_strand', 'ugann_ftype', 'ugann_gene'])

    def write_annotation_dataframe(self):
        '''
        Print out annotation data frames for inspection.
        '''
        print(self.df)
        # print(self.udf)

    def pickle_annotation_dataframe(self, annot_df_label):
        '''
        Build a annotation data frame pickle.
        '''
        print("... building a pickle for ", annot_df_label)
        self.df.to_pickle(annot_df_label+'.anndf.pkl')
        print("... done")
        print("Consider adding the pickle to your gainsaw.conf file.\n")

    def count_feature(self):
        '''
        Generate feature counts for query and target points.
        '''
        qf = [i for f in self.q_ftype.values() for i in f]
        tf = [i for f in self.t_ftype.values() for i in f]
        qf_count = Counter(qf)
        tf_count = Counter(tf)
        print("\nTotal feature counts for query points:")
        for f, c in qf_count.items():
            print(f, c)
        print("\nTotal feature counts for target points:")
        for f, c in tf_count.items():
            print(f, c)

def get_lopset(lopset_label):
    '''
    Retrieve a lopset pickle.
    '''
    if gsconf.has_option('pointsets',lopset_label) and path.exists(f := gsconf['pointsets'][lopset_label]):
        print("... loading existing LiftOverPointSet object ", f)
        lopset = pickle.load(open(f, "rb", -1))
        print("... done")
    else:
        print("... loading existing LiftOverPointSet object ", lopset_label+'.lops.pkl')
        lopset = pickle.load(open(lopset_label+'.lops.pkl', "rb", -1))
        print("... done")
        print("Consider adding the pickle to your gainsaw.conf file.\n")
    return lopset
