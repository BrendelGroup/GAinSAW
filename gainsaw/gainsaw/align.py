from os import path
from Bio import (Align, SeqIO)
import pickle

from .config import (gsconf, gsparams)


def get_fasta_index(fasta_label):
    '''
    Retrieve or build a sequence dictionary with chr:sequence.
    '''
    if gsconf.has_option(fasta_label,'gdx') and path.exists(f := gsconf[fasta_label]['gdx']):
        print("... loading existing fasta index")
        fasta_index = pickle.load(open(f, "rb", -1))
        print("... done")
    elif gsconf.has_option(fasta_label,'fna') and path.exists(f := gsconf[fasta_label]['fna']):
        print("... building and saving fasta index for ", f)
        fasta_index = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
        with open(gsconf[fasta_label]['gdx'], "wb") as file_:
           pickle.dump(fasta_index, file_, -1)
        print("... done")
        print("Consider adding the pickle to your gainsaw.conf file.\n")
    else:
        print("... building and saving fasta index for ", fasta_label)
        fasta_index = SeqIO.to_dict(SeqIO.parse(fasta_label, "fasta"))
        with open(fasta_label+'.pkl', "wb") as file_:
           pickle.dump(fasta_index, file_, -1)
        print("... done")
        print("Consider adding the pickle to your gainsaw.conf file.\n")
    return fasta_index

def get_fasta(fasta_index, chr, start, end, strand):
    '''
    Extract sequence from a sequence dictionary.
    '''
    start = max(0,start)
    seq = fasta_index[chr][start:end].seq
    if strand == '-':
        seq = seq.reverse_complement()
    return seq

def extract_segment_seq(pdata, params, qGenomic_fasta, tGenomic_fasta):
    '''
    Extract both query and target segment sequences in a pdata set.
    Return a dictionary of merged postions as key:(q_seq, t_seq, strand).
    '''
    qfasta_index = get_fasta_index(qGenomic_fasta)
    tfasta_index = get_fasta_index(tGenomic_fasta)
    segment_seqset = {}
    for i in pdata:
        if i == None:
            continue
        # [(9999999, 'chr1', 9988320, 'chr1', '+', 5396899, 12954106, 5385220, 12942427, '2',21057807908, 1)]
        # query point, query chr, target point, target chr, target strand, 
        # query alignment block, target alignment block, chain id, chain score
        for j in [g for g in i if g[-1][0] >= params.usepfilter]:
            q_seq = get_fasta(qfasta_index, j[1], j[5], j[6], '+')
            t_seq = get_fasta(tfasta_index, j[3], j[7], j[8], j[4])
            key = str(j[0])+"_in_"+j[1]+":"+str(j[5])+"-"+str(j[6])+"|"+j[3]+":"+str(j[7])+"-"+str(j[8])             
            segment_seqset[key] = q_seq, t_seq, j[4]
    return segment_seqset
    
def get_segment_seqset(pdata, params, qGenomic_fasta, tGenomic_fasta):
    '''
    Retrieve sequence data for alignment.
    '''
    segment_seqset = []
    segment_seqs = extract_segment_seq(pdata, params, qGenomic_fasta, tGenomic_fasta)
    print("There are ", len(segment_seqs), " set(s) of segment sequences: ")
    for i in segment_seqs.items():
        seq_names, qs, ts, strand = i[0], i[1][0].upper(), i[1][1].upper(), i[1][2]
        segment_seqset.append((seq_names, qs, ts, strand))
    return segment_seqset

def calculate_mismatch(alignment):
    '''
    Parsing pairwise alignment output in PSL format for percentage of mismatch.
    Here we have equal lengthes, no insertions or deletions."
    '''
    psl = alignment.format("psl").split("\t")
    match, mismatch = int(psl[0]), int(psl[1])
    block_count = int(psl[17])
    if block_count > 1:
        print("Number of gaps =", block_count-1)
    pc_mismatch = round(mismatch/(match + mismatch) * 100, 1)
    return str(mismatch)+"/"+str(match+mismatch)+" = "+str(pc_mismatch)+"%"

def align_segments(seq1, seq2, strand, scoring = gsparams.scoring):
    '''
    Align target and query sequences using a specified scoring scheme.
    Scoring should be (open_gap_score, extend_gap_score, substitution_matrix).
    The default is blastn scoring modified for ungapped alignment.
    default  = (-100, -2, "BLASTN")
    blastn = (-7, -2, "BLASTN")
    lastz = (-400, -30, "HOXD70")
    '''
    aligner = Align.PairwiseAligner(open_gap_score = scoring[0], extend_gap_score = scoring[1],
                    substitution_matrix = Align.substitution_matrices.load(scoring[2]))    
    aln = aligner.align(seq1, seq2)
    print("Score = ", aln.score)
    for i in aln:
        pc_mismatch = calculate_mismatch(i)
        print("Fraction of mismatches =", pc_mismatch)
        print(i.format(""))

def extract_extended_point_seq(pdata, params, qGenomic_fasta, tGenomic_fasta, slop_size):
    '''
    Extract both query and target point sequences in a pdata set after extension.
    Return a dictionary of merged postions as key:(q_seq, t_seq, strand).
    '''
    # There is a danger of slopping without chr length control,
    # if the point is close to either end of a chr.
    # A safer way would be using the bedtools slop output.
    # Or, bedtools getfasta, but we have already prepare sequence dbs.
    qfasta_index = get_fasta_index(qGenomic_fasta)
    tfasta_index = get_fasta_index(tGenomic_fasta)
    extended_point_seqset = {}
    for i in pdata:
        if i == None:
            continue
        # [(9999999, 'chr1', 9988320, 'chr1', '+', 5396899, 12954106, 5385220, 12942427, '2',21057807908, 1)]
        for j in [g for g in i if g[-1][0] >= params.usepfilter]:
            q_seq = get_fasta(qfasta_index, j[1], j[0]-slop_size, j[0]+slop_size+1, "+")
            t_seq = get_fasta(tfasta_index, j[3], j[2]-slop_size, j[2]+slop_size+1, j[4])
            key = str(j[0])+"_in_"+j[1]+":"+str(j[0]-slop_size)+"-"+str(j[0]+slop_size+1)+"|"+j[3]+":"+str(j[2]-slop_size)+"-"+str(j[2]+slop_size+1)             
            extended_point_seqset[key] = q_seq, t_seq, j[4]
    return extended_point_seqset
    
def get_extended_point_seqset(pdata, params, qGenomic_fasta, tGenomic_fasta, slop_size):
    '''
    Retrieve sequence data for alignment.
    '''
    extended_point_seqset = []
    extended_point_seqs = extract_extended_point_seq(pdata, params, qGenomic_fasta, tGenomic_fasta, slop_size)
    print()
    print("There are ", len(extended_point_seqs), " set(s) of extended_point sequences: ")
    for i in extended_point_seqs.items():
        ###VERBOSE###print(i)
        seq_names, qs, ts, strand = i[0], i[1][0].upper(), i[1][1].upper(), i[1][2]
        extended_point_seqset.append((seq_names, qs, ts, strand))
    return extended_point_seqset

def string_mismatch(seq1, seq2):
    '''
    Treat a pair of equal length sequences as strings
    for simple mismatch calculation.
    '''
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    m = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            m += 1
    return round(m/len(seq1) * 100, 1)

