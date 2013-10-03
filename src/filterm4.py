#!/usr/bin/env python

# Filters for unique, highest scoring subread query/target pairs from an m4
# file. Helps get rid of chimeras, at the cost of some yield.

import sys
from collections import namedtuple

M4Record = namedtuple('M4Record', ('qname tname score pctsimilarity qstrand '
                                   'qstart qend qseqlength tstrand tstart '
                                   'tend tseqlength mapqv'))


class Count(object):
    """Tracks record count for original and filtered"""
    def __init__(self):
        self.orig = 0
        self.filt = 0

    def __repr__(self):
        return "Record count: original=%i, filtered=%i\n" % \
            (self.orig, self.filt)


def printUniq(qgroup, count):
    top = dict()
    for q in qgroup:
        m = M4Record._make(q.split())
        k = "%s%s" % (m.qname, m.tname)
        if k in top:
            n = M4Record._make(top[k].split())
            if int(m.score) < int(n.score):
                top[k] = q
        else:
            top[k] = q

    for r in top.values():
        count.filt += 1
        print r,

    qgroup[:] = []


def main():
    m4file = sys.argv[1]
    m4Hndl = open(m4file)
    qgroup = []
    curr = ''
    count = Count()
    for rec in m4Hndl:
        count.orig += 1
        m = M4Record._make(rec.split())
        if curr != m.qname:
            printUniq(qgroup, count)
            qgroup.append(rec)
            curr = m.qname
        else:
            qgroup.append(rec)

    printUniq(qgroup, count)
    sys.stderr.write(str(count))

if __name__ == '__main__':
    sys.exit(main())
