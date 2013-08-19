#!/usr/bin/env python

# Super-simple converter from blasr m4 alignments to pbdagcon 'pre'
# alignments. For use in the pre-assembler dagcon workflow.
import sys
import csv
import string
import heapq
from itertools import ifilter
from collections import namedtuple, defaultdict
from pbcore.io.FastaIO import FastaReader


M4Record = namedtuple('M4Record',
     ('qname tname score pctsimilarity qstrand qstart qend qseqlength '
     'tstrand tstart tend tseqlength mapqv'))

# alias to avoid loop eval
tuplfy = M4Record._make

# dna compliment
rc = string.maketrans('actgACTG', 'tgacTGAC')


def sortTargScore(rec):
    f = rec.split()
    return (f[1], int(f[2]))


def bestnTrue(rec, myq):
    m = tuplfy(rec.split())
    score = -int(m.score)
    alen = int(m.tend) - int(m.tstart)
    return (score, alen, rec) in myq[m.qname].top


class AlnLimiter(object):
    def __init__(self, limit=76):
        self.count = 0
        self.target = ''
        self.limit = limit

    def __call__(self, rec):
        target = rec.split()[1]
        if target != self.target:
            self.count = 0
            self.target = target
        self.count += 1
        return self.count < self.limit


class TopAlignments(object):
    bestn = 10

    def __init__(self):
        self.empty = (0, 0, None)
        self.top = [self.empty for x in xrange(0, TopAlignments.bestn)]

    def __call__(self):
        return  # noop

    def add(self, aln):
        heapq.heappushpop(self.top, aln)


def main():
    myM4 = sys.argv[1]
    allM4 = sys.argv[2]
    reads = sys.argv[3]
    bestn = int(sys.argv[4])

    TopAlignments.bestn = bestn

    # tracks bestn
    myQueries = defaultdict(TopAlignments)

    myM4Recs = set()

    # load my m4 chunk
    myM4Hndl = open(myM4)
    recAdd = myM4Recs.add
    for rec in myM4Hndl:
        recAdd(rec)
        m = tuplfy(rec.split())
        score = -int(m.score)
        alen = int(m.tend) - int(m.tstart)
        myQueries[m.qname].add((score, alen, rec))

    myM4Hndl.close()

    # if we're chunked locate relevant alignments
    if myM4 != allM4:
        # assuming fofn here
        m4files = [x.rstrip() for x in open(allM4) if x.rstrip() != myM4]
        for m4 in m4files:
            m4Hndl = open(m4)
            for rec in m4Hndl:
                m = tuplfy(rec.split())
                if m.qname in myQueries:
                    score = -int(m.score)
                    alen = int(m.tend) - int(m.tstart)
                    myQueries[m.qname].add((score, alen, rec))
            m4Hndl.close()

        # remove alignments that fall outside of bestn
        myM4Recs = [x for x in myM4Recs if bestnTrue(x, myQueries)]

    # sort by target name/score
    myM4Recs.sort(key=sortTargScore)

    # take a max number of alignments for each target
    limiter = AlnLimiter()
    myM4Recs = [x for x in ifilter(limiter, myM4Recs)]

    # load only related sequences
    seqs = {}
    f = FastaReader(reads)
    for e in f:
        if e.name in myQueries:
            seqs[e.name] = e.sequence

    # generate pre-alignments
    for rec in myM4Recs:
        m = tuplfy(rec.split())
        qs = int(m.qstart)
        qe = int(m.qend)
        qseq = seqs[m.qname][qs:qe]
        strand = '-' if m.tstrand == '1' else '+'
        ts = int(m.tstart)
        te = int(m.tend)
        if strand == '+':
            tseq = seqs[m.tname][ts:te]
        else:
            tseq = seqs[m.tname].translate(rc)[::-1][ts:te]

        print ' '.join([m.qname, m.tname, strand,
            m.tseqlength, str(ts), str(te), qseq, tseq])

if __name__ == '__main__':
    sys.exit(main())
