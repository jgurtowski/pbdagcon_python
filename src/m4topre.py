#!/usr/bin/env python

# Super-simple converter from blasr m4 alignments to pbdagcon 'pre'
# alignments. For use in the pre-assembler dagcon workflow.
import sys
import csv
import string
import heapq
from collections import namedtuple, defaultdict
from pbcore.io.FastaIO import FastaReader


M4Record = namedtuple('M4Record',
     ('qname tname score pctsimilarity qstrand qstart qend qseqlength '
     'tstrand tstart tend tseqlength mapqv'))

rc = string.maketrans('actgACTG', 'tgacTGAC')


class TopAlignments(object):
    bestn = 10

    def __init__(self):
        self.empty = (0, 0, None)
        self.top = [self.empty for x in xrange(0, TopAlignments.bestn)]

    def __call__(self):
        return  # noop

    def remove(self, i):
        self.top[i] = self.empty


def main():
    myM4 = sys.argv[1]
    allM4 = sys.argv[2]
    reads = sys.argv[3]
    bestn = int(sys.argv[4])

    TopAlignments.bestn = bestn

    myTargets = defaultdict(list)
    myQueries = defaultdict(TopAlignments)

    # load my m4 chunk
    for m in map(M4Record._make, csv.reader(open(myM4), delimiter=' ')):
        myTargets[m.tname].append(m)
        subreadId = m.qname[:m.qname.rfind('/')]
        myQueries[subreadId]

    # if we're chunked locate relevant alignments.  remove alignments from
    # the target that fall outside bestn
    if myM4 != allM4:
        for m in map(M4Record._make, csv.reader(open(allM4), delimiter=' ')):
            subreadId = m.qname[:m.qname.rfind('/')]
            if subreadId in myQueries:
                score = -int(m.score)
                alen = int(m.tend) - int(m.tstart)
                heapq.heappushpop(myQueries[subreadId].top, (score, alen, m))

        for tid, alnList in myTargets.iteritems():
            for m in alnList:
                subreadId = m.qname[:m.qname.rfind('/')]
                score = -int(m.score)
                alen = int(m.tend) - int(m.tstart)
                if (score, alen, m) not in myQueries[subreadId].top:
                    alnList.remove(m)

    # load related sequences in all m4 files
    seqs = {}
    f = FastaReader(reads)
    for e in f:
        if e.name in myTargets or e.name in myQueries:
            seqs[e.name] = e.sequence

    # generate pre-alignments
    for target, alignments in myTargets.iteritems():
        for m in alignments:
            qname = m.qname[:m.qname.rfind('/')]
            qs = int(m.qstart)
            qe = int(m.qend)
            qseq = seqs[qname][qs:qe]
            strand = '-' if m.tstrand == '1' else '+'
            ts = int(m.tstart)
            te = int(m.tend)
            if strand == '+':
                tseq = seqs[m.tname][ts:te]
            else:
                tseq = seqs[m.tname].translate(rc)[::-1][ts:te]

            print ' '.join([qname, m.tname, strand,
                m.tseqlength, m.tstart, m.tend, qseq, tseq])

if __name__ == '__main__':
    sys.exit(main())
