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
import logging
import pkg_resources

from pbcore.util.ToolRunner import PBMultiToolRunner
from pbcore.io import FastaReader
from pbtools.pbdagcon.q_sense import *

try:
    __p4revision__ = "$Revision: #15 $"
    __p4change__ = "$Change: 115421 $"
    revNum = int(__p4revision__.strip("$").split(" ")[1].strip("#"))
    changeNum = int(__p4change__.strip("$").split(":")[-1])
    __version__ = "%s-r%d-c%d" % ( pkg_resources.require("pbtools.pbdagcon")[0].version, revNum, changeNum )
except:
    __version__ = "pbtools.pbdagcon-github"


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
            subp.add_argument('--enable_hp_correction', action='store_true', default = False, dest="enable_hp_corr",
                               help = 'enable aggressive homopolymer missing errot detection and correction')
            subp.add_argument('--mark_lower_case', action='store_true', default = False, dest="mark_lower_case",
                               help = 'mark low quality consensus base with lower case letter')
            subp.add_argument('--hp_correction_th', dest = 'entropy_th', default = 0.65,
                               help = 'homopolymer missing correction entropy threshold')
            subp.add_argument('--n_iter', default = 4, dest = 'niter',
                              help = 'number of iteration of consensus correction')
            subp.add_argument('--min_cov', default = 8, dest = 'min_cov',
                              help = 'minimum coverage for generate consensus')
            subp.add_argument('--max_cov', default = 60, dest = 'max_cov',
                              help = 'maximum coverage for generate consensus')
            subp.add_argument('--max_n_reads', default = 150, dest = 'max_num_reads',
                              help = 'the maximum number of reads used for consensus')
            subp.add_argument('--nproc', default = 4, dest = 'nproc',
                              help = 'the number cpu core used by blasr, default to 4 cores')
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
        hp_corr = True if self.args.enable_hp_corr else False
        mark_lower_case = True if self.args.mark_lower_case else False
        generate_consensus(input_fasta_name, "%s_ref.fa" % full_prefix, full_prefix, self.args.consensus_seq_name, 
                           hp_corr, int(self.args.niter), 
                           int(self.args.max_num_reads), 
                           float(self.args.entropy_th),
                           self.args.dump_dag_info,
                           int(self.args.min_cov),
                           int(self.args.max_cov),
                           mark_lower_case,
                           int(self.args.nproc))

    def refConsensus(self):
        input_fasta_name = self.args.input 
        prefix = self.args.out_file_name.split(".")
        if len(prefix) > 1:
            prefix = ".".join(prefix[:-1])
        else:
            prefix = ".".join(prefix)
        full_prefix = os.path.join(self.args.out_dir_name, prefix)
        hp_corr = True if self.args.enable_hp_corr else False
        mark_lower_case = True if self.args.mark_lower_case else False
        generate_consensus(input_fasta_name, self.args.ref, full_prefix, self.args.consensus_seq_name,
                           hp_corr, int(self.args.niter), 
                           int(self.args.max_num_reads), 
                           float(self.args.entropy_th),
                           self.args.dump_dag_info,
                           int(self.args.min_cov),
                           int(self.args.max_cov),
                           mark_lower_case,
                           int(self.args.nproc))

    def run(self):
        logging.debug("Arguments" + str(self.args))
        if self.args.subName == 'd':
            self.denovoConsensus()
        elif self.args.subName == 'r':
            self.refConsensus()

if __name__ == '__main__':    
    sys.exit(Consensus().start())

