#!/usr/bin/env python

#################################################################################$$
# Copyright (c) 2011,2012, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright notice, this 
#   list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice, 
#   this list of conditions and the following disclaimer in the documentation 
#   and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its contributors 
#   may be used to endorse or promote products derived from this software 
#   without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS CONTRIBUTORS 
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS 
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON 
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#################################################################################$$


import sys
import os
import shutil
import glob
import logging
import tempfile
import subprocess
import pkg_resources
import zlib
import subprocess

from pbcore.util.ToolRunner import PBMultiToolRunner
from pbcore.io.FastaIO import SimpleFastaReader

from pbtools.pbdagcon.aligngraph import *
from pbtools.pbdagcon.utils import constructe_aln_graph_from_fasta 
from pbtools.pbdagcon.utils import sorted_nodes
from pbtools.pbdagcon.utils import best_template_by_blasr
from pbtools.pbdagcon.utils import clustering_read
from pbtools.pbdagcon.utils import get_subset_reads
from pbtools.pbdagcon.utils import read_node_vector
from pbtools.pbdagcon.utils import detect_missing


__p4revision__ = "$Revision: #9 $"
__p4change__ = "$Change: 105511 $"
revNum = int(__p4revision__.strip("$").split(" ")[1].strip("#"))
changeNum = int(__p4change__.strip("$").split(":")[-1])
#__version__ = "%s-r%d-c%d" % ( pkg_resources.require("pbtools.aligngraph")[0].version, revNum, changeNum )
__version__ = "rc"


rmap = dict(zip("ACGTN-","TGCAN-"))

def normalize_fasta(fastaFile, refFile, outFile):
    f = SimpleFastaReader(fastaFile)
    recs = []
    with open(outFile, "w") as of:
        for r in f:
            r_id = "%s" %  hex(zlib.adler32(r.name + r.sequence) & 0xffffffff)
            print >>of, ">"+r_id
            seq = r.sequence.upper()
            print >>of, seq 


    output = subprocess.check_output("blasr -bestn 1 -m 1 %s %s" % ( outFile, refFile ), shell=True)
    direction = {}
    output = output.strip().split("\n")
    for l in output:
        l = l.strip().split()
        rId = l[0].split("/")[0]
        if l[2] != l[3]:
            direction[rId] = "-"
        else:
            direction[rId] = "+"

    f = SimpleFastaReader(outFile)
    outData = []
    for r in f:
        r_id = "%s" % r.name
        outData.append(">"+r_id)
        seq = r.sequence.upper()
        if direction != None:
            if direction.get(r_id, "+") != "+":
                seq = "".join([rmap(c) for c in seq[::-1]])
        outData.append(seq)
    with open(outFile,"w") as of:
        print >>of, "\n".join(outData)


def get_consensus(read_fn, init_ref, consensus_fn, consens_seq_name, min_iteration = 4, hp_correction = True):
    g = constructe_aln_graph_from_fasta(read_fn, init_ref, max_num_reads = 150, remove_in_del = False)
    s,c = g.generate_consensus()
    with open(consensus_fn,"w") as f:
        print >>f, ">"+consens_seq_name
        print >>f, s.upper()
    if min_iteration > 1:
        for i in range(min_iteration-2):
            g = constructe_aln_graph_from_fasta(read_fn, consensus_fn, max_num_reads = 150, remove_in_del = False)
            s,c = g.generate_consensus()
            with open(consensus_fn,"w") as f:
                print >>f, ">"+consens_seq_name
                print >>f, s.upper()

        if hp_correction:
            g = constructe_aln_graph_from_fasta(read_fn, consensus_fn, max_num_reads = 150, remove_in_del = False)
            s = detect_missing(g, entropy_th = 0.65)
            with open(consensus_fn,"w") as f:
                print >>f, ">"+consens_seq_name
                print >>f, s.upper()

        g = constructe_aln_graph_from_fasta(read_fn, consensus_fn, max_num_reads = 150, remove_in_del = False)
        s,c = g.generate_consensus()
        with open(consensus_fn,"w") as f:
            print >>f, ">"+consens_seq_name
            print >>f, s.upper()



def generate_haplotype_consensus(inputFastaName, refFastaName, prefix, consensusName, hpFix):


    normalize_fasta(inputFastaName, refFastaName, "%s_input.fa" % prefix)

    get_consensus("%s_input.fa" % prefix, 
                  "%s_ref.fa" % prefix, 
                  "%s.fa" % prefix, 
                  consensusName,
                  hp_correction = False)

    g = constructe_aln_graph_from_fasta("%s_input.fa" % prefix, 
                               "%s.fa" % prefix, 
                               ref_group=consensusName, 
                               max_num_reads = 300, 
                               remove_in_del = False)

    rv, hen = read_node_vector(g, entropy_th = 0.65)
    cluster, cluster_vec = clustering_read(rv, hen, k_cluster = 2, random_seed = 42)


    with open("%s.log" % prefix, "w") as logf:
        print >>logf, len(rv)

        for k in cluster:
            print >>logf, cluster_vec[k], k, len(cluster[k])

        for k in cluster:
            for r in cluster[k]:
                print >>logf, "".join(rv[r]), k, r
    
    if len(cluster[0]) > 0:
        get_subset_reads("%s_input.fa" % prefix, cluster, 0, "%s_h1_input.fa" % prefix)
        
        rid,s = best_template_by_blasr("%s_h1_input.fa" % prefix)
        print rid, len(s)
        with open("%s_h1_ref.fa" % prefix, "w") as f:
            print >>f ,">%s_h1_ref" % prefix
            print >>f, s
        get_consensus("%s_h1_input.fa" % prefix, 
                      "%s_h1_ref.fa" % prefix, 
                      "%s_h1.fa" % prefix, 
                      "%s_h1" % consensusName,
                      hp_correction = hpFix)
    
    if len(cluster[1]) > 0:
        get_subset_reads("%s_input.fa" % prefix, cluster, 1, "%s_h2_input.fa" % prefix)
        
        rid,s = best_template_by_blasr("%s_h2_input.fa" % prefix)
        print rid, len(s)
        with open("%s_h2_ref.fa" % prefix, "w") as f:
            print >>f ,">%s_h2_ref" % prefix
            print >>f, s
        get_consensus("%s_h2_input.fa" % prefix, 
                      "%s_h2_ref.fa" % prefix, 
                      "%s_h2.fa" % prefix, 
                      "%s_h2" % consensusName,
                      hp_correction = hpFix)

class HapConsensus(PBMultiToolRunner):

    def __init__(self):
        desc = ["""Making haplotype consensus sequences from a group reads storing in a fasta file.\
                All of the reads are expected to cover most of the target templates that gets sequenced.\
                If the reads have broad read length distribution and not all of them cover the same region of\
                a template, this code won't generat correct result.\
                This code is designed for getting consensus up to the length of reads (~10k).\
                It is not optimized for getting consensus for larger templates."""]
        super(HapConsensus, self).__init__('\n'.join(desc))

        subparsers = self.getSubParsers()

        desc = ['using self-self alignment for getting the seed sequence for constructing consensus']
        parser_d = subparsers.add_parser('d', help = "generate the seed sequence for consensus using alignments",
                                         description = "\n".join(desc), parents = [self.parser])
        parser_d.add_argument('input', metavar = 'input.fasta',
                              help = 'an input fasta file')
        
        desc = ['using a reference file as seed for consensus']
        parser_r = subparsers.add_parser('r', help = "using a reference fasta as the seed sequence",
                                         description = "\n".join(desc), parents = [self.parser])
        parser_r.add_argument('input', metavar = 'input.fasta',
                              help = 'an input fasta file')
        parser_r.add_argument('ref', metavar = 'ref.fasta',
                              help = 'a reference fasta file')

        for subp in (parser_r, parser_d): 
            subp.add_argument('--output', metavar = 'file-name', dest = 'outFileName', default = "ag_consensus", 
                               help = 'consensus output filename')
            subp.add_argument('--outputDir', metavar = 'directory-name', dest = 'outDirName', default = "./", 
                               help = 'consensus output working directory')
            subp.add_argument('--cname', metavar = 'consensus-seq-name', dest = 'consensusSeqName', default = "consensus", 
                               help = 'consensus sequence name')
            subp.add_argument('--disable_hp_correction', action='store_true', default = False, dest="disable_hp_corr",
                               help = 'disable aggressive homopolymer missing errot detection and correction')
            subp.add_argument('--hp_correction_th', dest = 'entropy_th', default = 0.65,
                               help = 'homopolymer missing correction entropy threshold')
            subp.add_argument('--n_iter', default = 4, dest = 'niter',
                              help = 'number of iteration of consensus correction')
            subp.add_argument('--max_n_reads', default = 150, dest = 'maxNReads',
                              help = 'the maximum number of reads used for consensus')
                    
    def getVersion(self):
        return __version__
    
    def denovoConsensus(self):
        inputFastaName = self.args.input 
        rid,s =best_template_by_blasr(inputFastaName)
        prefix = self.args.outFileName.split(".")
        if len(prefix) > 1:
            prefix = ".".join(prefix[:-1])
        else:
            prefix = ".".join(prefix)
        full_prefix = os.path.join(self.args.outDirName, prefix)
        with open("%s_ref.fa" % full_prefix, "w") as f:
            print >>f ,">%s_ref" % self.args.consensusSeqName
            print >>f, s
        hp_corr = False if self.args.disable_hp_corr else True
        generate_haplotype_consensus(inputFastaName, "%s_ref.fa" % full_prefix, full_prefix, self.args.consensusSeqName, hpFix = hp_corr)

    def refConsensus(self):
        inputFastaName = self.args.input 
        prefix = self.args.outFileName.split(".")
        if len(prefix) > 1:
            prefix = ".".join(prefix[:-1])
        else:
            prefix = ".".join(prefix)
        full_prefix = os.path.join(self.args.outDirName, prefix)
        hp_corr = False if self.args.disable_hp_corr else True
        generate_haplotype_consensus(inputFastaName, self.args.ref, full_prefix, self.args.consensusSeqName, hpFix = hp_corr)

    def run(self):
        if self.args.subName == 'd':
            self.denovoConsensus()
        elif self.args.subName == 'r':
            self.refConsensus()

if __name__ == '__main__':    
    sys.exit(HapConsensus().start())

