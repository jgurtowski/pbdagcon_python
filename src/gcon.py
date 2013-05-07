#!/usr/bin/env python

#################################################################################$$
# Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted (subject to the limitations in the
# disclaimer below) provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright
#  notice, this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above
#  copyright notice, this list of conditions and the following
#  disclaimer in the documentation and/or other materials provided
#  with the distribution.
#
#  * Neither the name of Pacific Biosciences nor the names of its
#  contributors may be used to endorse or promote products derived
#  from this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
# GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
# BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
# USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.
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
from pbcore.io import FastaReader

from pbtools.pbdagcon.aligngraph import *
from pbtools.pbdagcon.utils import construct_aln_graph_from_fasta 
from pbtools.pbdagcon.utils import sorted_nodes
from pbtools.pbdagcon.utils import best_template_by_blasr
from pbtools.pbdagcon.utils import clustering_read
from pbtools.pbdagcon.utils import get_subset_reads
from pbtools.pbdagcon.utils import read_node_vector
from pbtools.pbdagcon.utils import detect_missing
from pbtools.pbdagcon.utils import mark_lower_case_base

try:
    __p4revision__ = "$Revision$"
    __p4change__ = "$Change$"
    revNum = int(__p4revision__.strip("$").split(" ")[1].strip("#"))
    changeNum = int(__p4change__.strip("$").split(":")[-1])
    __version__ = "%s-r%d-c%d" % ( pkg_resources.require("pbtools.pbdagcon")[0].version, revNum, changeNum )
except:
    __version__ = "pbtools.pbdagcon-github"

rmap = dict(zip("ACGTN-","TGCAN-"))

def normalize_fasta(fastaFile, refFile, outFile):
    f = FastaReader(fastaFile)
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

    f = FastaReader(outFile)
    outData = []
    for r in f:
        r_id = "%s" % r.name
        outData.append(">"+r_id)
        seq = r.sequence.upper()
        if direction != None:
            if direction.get(r_id, "+") != "+":
                seq = "".join([rmap[c] for c in seq[::-1]])
        outData.append(seq)
    with open(outFile,"w") as of:
        print >>of, "\n".join(outData)


def get_consensus(read_fn, init_ref, consensus_fn, consens_seq_name, 
                  hp_correction = True, 
                  min_iteration = 4, 
                  max_num_reads = 150,
                  entropy_th = 0.65,
                  min_cov = 8):

    g = construct_aln_graph_from_fasta(read_fn, init_ref, max_num_reads = max_num_reads, remove_in_del = False)
    s,c = g.generate_consensus(min_cov = min_cov)
    with open(consensus_fn,"w") as f:
        print >>f, ">"+consens_seq_name
        print >>f, s.upper()

    if min_iteration > 1:
        for j in range(2):
            for i in range(min_iteration-2):
                g = construct_aln_graph_from_fasta(read_fn, consensus_fn, max_num_reads = max_num_reads, remove_in_del = False)
                s,c = g.generate_consensus(min_cov = min_cov)
                with open(consensus_fn,"w") as f:
                    print >>f, ">"+consens_seq_name
                    print >>f, s.upper()

            if hp_correction:
                g = construct_aln_graph_from_fasta(read_fn, consensus_fn, max_num_reads = max_num_reads, remove_in_del = False)
                s = detect_missing(g, entropy_th = entropy_th)
                with open(consensus_fn,"w") as f:
                    print >>f, ">"+consens_seq_name
                    print >>f, s.upper()

            g = construct_aln_graph_from_fasta(read_fn, consensus_fn, max_num_reads = max_num_reads, remove_in_del = False)
            s,c = g.generate_consensus(min_cov = min_cov)
            s = mark_lower_case_base(g, entropy_th = entropy_th)
            with open(consensus_fn,"w") as f:
                print >>f, ">"+consens_seq_name
                print >>f, s
        return g


def output_dag_info(aln_graph, out_file_name):

    with open(out_file_name,"w") as f:
        read_ids = []
        read_ids_set = set()
        for n in sorted_nodes(aln_graph):
            for r in n.info:
                if r not in read_ids_set:
                    read_ids.append(r)
                    read_ids_set.add(r)


        read_id_to_pos = dict(( (x[1],x[0]) for x in enumerate(list(read_ids))) )
        for read_id in read_ids:
            pileup_pos = read_id_to_pos[read_id]
            print >>f, "\t".join( [ "R", "%d" % pileup_pos, read_id ] )

        backbone_node_to_pos =  aln_graph.backbone_node_to_pos 
        ne, hne = aln_graph.get_high_entropy_nodes(coverage_th=0)
        node_to_entropy = dict( [ (v[1],v[2]) for v in ne ] ) 

        data = []
        consensus_pos = 0
        for n in sorted_nodes(aln_graph):
            s = ["."] * len(read_id_to_pos)
            for r in n.info:
                s[read_id_to_pos[r]] = n.base

            if n.base not in ["B", "E"]:
                bpos = backbone_node_to_pos[n.backbone_node]
            entropy_th = 0
            ent = node_to_entropy[n] if n in node_to_entropy else 0
            if n.base not in ["B","E"] and ent >= entropy_th:
                data.append ( ( n.ID, consensus_pos, backbone_node_to_pos[n.backbone_node],\
                      "+" if n in aln_graph.consensus_path else "-",\
                      "+" if n.is_backbone == True else "-",\
                      n.base, "".join(s), len(n.info),\
                      n.backbone_node.coverage,\
                      node_to_entropy[n] if n in node_to_entropy else "-" ) )
                if n in aln_graph.consensus_path:
                    consensus_pos += 1

        for l in data:
            print >>f, "N"+"\t"+"\t".join([str(c) for c in l])

        for e_id in aln_graph.edges:
            e = aln_graph.edges[e_id]
            print >> f, "\t".join( ["E", "%d" % e_id, "%d" % e.in_node.ID, "%d" % e.out_node.ID, "%d" % e.count ] )

def generate_consensus(input_fasta_name, 
                       ref_fasta_name, 
                       prefix, 
                       consensus_name, 
                       hpFix, 
                       min_iteration,
                       max_num_reads,
                       entropy_th,
                       dump_dag_info,
                       min_cov):

    normalize_fasta(input_fasta_name, ref_fasta_name, "%s_input.fa" % prefix)

    g = get_consensus("%s_input.fa" % prefix, 
                      ref_fasta_name, 
                      "%s.fa" % prefix, 
                      consensus_name,
                      hp_correction = hpFix,
                      min_iteration = min_iteration,
                      max_num_reads = max_num_reads,
                      entropy_th = entropy_th,
                      min_cov = min_cov)
    if dump_dag_info == True:
        output_dag_info(g, "%s_dag.dat" % prefix)

class Consensus(PBMultiToolRunner):

    def __init__(self):
        desc = ["Making consensus sequence from a group reads storing in a fasta file. "
                "All of the reads are expected to cover most of the target templates that gets sequenced. "
                "If the reads have broad read length distribution and not all of them cover the same region of "
                "a template, this code won't generate correct result. "
                "This code is designed for getting consensus up to the length of reads (~10k). "
                "It is not optimized for getting consensus for larger templates."]
        super(Consensus, self).__init__('\n'.join(desc))

        subparsers = self.subParsers

        desc = ['using self-self alignment for getting the seed sequence for constructing consensus']
        parser_d = subparsers.add_parser('d', help = "generate consensus using the alignments between the input sequneces to find seed sequence",
                                         description = "\n".join(desc))
        parser_d.add_argument('input', metavar = 'input.fasta',
                              help = 'an input fasta file')
        
        desc = ['using a reference file as seed for consensus']
        parser_r = subparsers.add_parser('r', help = "using a reference fasta as the seed sequence",
                                         description = "\n".join(desc))
        parser_r.add_argument('input', metavar = 'input.fasta',
                              help = 'an input fasta file')
        parser_r.add_argument('ref', metavar = 'ref.fasta',
                              help = 'a reference fasta file')

        for subp in (parser_r, parser_d): 
            subp.add_argument('-o', '--output', metavar = 'file-name', dest = 'out_file_name', default = "g_consensus", 
                               help = 'consensus output filename')
            subp.add_argument('-d', '--output_dir', metavar = 'directory-name', dest = 'out_dir_name', default = "./", 
                               help = 'consensus output working directory')
            subp.add_argument('--cname', metavar = 'consensus-seq-name', dest = 'consensus_seq_name', default = "consensus", 
                               help = 'consensus sequence name')
            subp.add_argument('--disable_hp_correction', action='store_true', default = False, dest="disable_hp_corr",
                               help = 'disable aggressive homopolymer missing errot detection and correction')
            subp.add_argument('--hp_correction_th', dest = 'entropy_th', default = 0.65,
                               help = 'homopolymer missing correction entropy threshold')
            subp.add_argument('--n_iter', default = 4, dest = 'niter',
                              help = 'number of iteration of consensus correction')
            subp.add_argument('--min_cov', default = 8, dest = 'min_cov',
                              help = 'minimum coverage for generate consensus')
            subp.add_argument('--max_n_reads', default = 150, dest = 'max_num_reads',
                              help = 'the maximum number of reads used for consensus')
            subp.add_argument('--dump_dag_info', action='store_true', default = False, dest="dump_dag_info",
                               help = 'dump the information of the dag, including a pileup view of the alignments')
                    
    def getVersion(self):
        return __version__
    
    def denovoConsensus(self):
        prefix = self.args.out_file_name.split(".")
        input_fasta_name = self.args.input 

        rid,s =best_template_by_blasr(input_fasta_name)
        if len(prefix) > 1:
            prefix = ".".join(prefix[:-1])
        else:
            prefix = ".".join(prefix)
        full_prefix = os.path.join(self.args.out_dir_name, prefix)
        with open("%s_ref.fa" % full_prefix, "w") as f:
            print >>f ,">%s_ref" % self.args.consensus_seq_name
            print >>f, s
        hp_corr = False if self.args.disable_hp_corr else True
        generate_consensus(input_fasta_name, "%s_ref.fa" % full_prefix, full_prefix, self.args.consensus_seq_name, 
                           hp_corr, int(self.args.niter), 
                           int(self.args.max_num_reads), 
                           float(self.args.entropy_th),
                           self.args.dump_dag_info,
                           int(self.args.min_cov))

    def refConsensus(self):
        input_fasta_name = self.args.input 
        prefix = self.args.out_file_name.split(".")
        if len(prefix) > 1:
            prefix = ".".join(prefix[:-1])
        else:
            prefix = ".".join(prefix)
        full_prefix = os.path.join(self.args.out_dir_name, prefix)
        hp_corr = False if self.args.disable_hp_corr else True
        generate_consensus(input_fasta_name, self.args.ref, full_prefix, self.args.consensus_seq_name,
                           hp_corr, int(self.args.niter), 
                           int(self.args.max_num_reads), 
                           float(self.args.entropy_th),
                           self.args.dump_dag_info,
                           int(self.args.min_cov))

    def run(self):
        logging.debug("Arguments" + str(self.args))
        if self.args.subCommand == 'd':
            self.denovoConsensus()
        elif self.args.subCommand == 'r':
            self.refConsensus()

if __name__ == '__main__':    
    sys.exit(Consensus().start())

