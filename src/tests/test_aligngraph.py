from nose.tools import assert_equal
from nose import SkipTest
import random
from pbtools.pbdagcon.aligngraph import *

def generate_simulated_reads(pi=None, pd=None, n = 4):
    import random
    random.seed(42)
    seq   = "ATATTTGGC"
    seq1  = "ATAGCCGGC"
    seq2  = "ATACCCGGC"
    seq3  = "ATATCCGGC"
    seq4  = "ATATCGGC"
    if pi == None:
        pi = 0.03
    if pd == None:
        pd = 0.03
    out_seq = []
    for i in range(n):
        c = 0
        s = []
        if i % 4 == 0:
            ss = seq1
        elif i % 4 == 1:
            ss = seq2
        elif i % 4 == 2:
            ss = seq3
        else:
            ss = seq4
        while 1:
            if random.uniform(0,1) < pi:
                s.append(random.choice( ("A","G","C","T") ) )
                continue
            if random.uniform(0,1) < pd:
                c += 1
                continue
            if c < len(ss):
                s.append(ss[c])
                c += 1
            else:
                break
        out_seq.append( "".join(s) )
    return seq, out_seq

class TestPhiCoeff:
    def test_phi_coeff(self):
        # assert_equal(expected, phi_coeff(xvec, yvec))
        raise SkipTest # TODO: implement your test here

class TestConvertMismatches:
    def test_convert_mismatches(self):
        assert_equal( ('C-AC', 'CG-C'), convert_mismatches("CAC","CGC") ) 
        assert_equal( ('CAACAT', 'CAA--T'), convert_mismatches("CAACAT","C-A-AT" ) )
        assert_equal( ('CCG--T', 'CCGACT'), convert_mismatches("-C--CGT","CCGAC-T") ) 

class TestAlnEdge:
    def test___init__(self):
        # aln_edge = AlnEdge(in_node, out_node)
        raise SkipTest # TODO: implement your test here

    def test___repr__(self):
        # aln_edge = AlnEdge(in_node, out_node)
        # assert_equal(expected, aln_edge.__repr__())
        raise SkipTest # TODO: implement your test here

    def test_add_to_score(self):
        # aln_edge = AlnEdge(in_node, out_node)
        # assert_equal(expected, aln_edge.add_to_score(s))
        raise SkipTest # TODO: implement your test here

    def test_increase_count(self):
        # aln_edge = AlnEdge(in_node, out_node)
        # assert_equal(expected, aln_edge.increase_count())
        raise SkipTest # TODO: implement your test here

    def test_set_score(self):
        # aln_edge = AlnEdge(in_node, out_node)
        # assert_equal(expected, aln_edge.set_score(s))
        raise SkipTest # TODO: implement your test here

class TestAlnNode:
    def test___init__(self):
        # aln_node = AlnNode(base)
        raise SkipTest # TODO: implement your test here

    def test___repr__(self):
        # aln_node = AlnNode(base)
        # assert_equal(expected, aln_node.__repr__())
        raise SkipTest # TODO: implement your test here

    def test_add_in_edge(self):
        # aln_node = AlnNode(base)
        # assert_equal(expected, aln_node.add_in_edge(in_edge))
        raise SkipTest # TODO: implement your test here

    def test_addout_edge(self):
        # aln_node = AlnNode(base)
        # assert_equal(expected, aln_node.addout_edge(out_edge))
        raise SkipTest # TODO: implement your test here

    def test_increase_weight(self):
        # aln_node = AlnNode(base)
        # assert_equal(expected, aln_node.increase_weight(w))
        raise SkipTest # TODO: implement your test here

class TestAlnGraph:
    def test___init__(self):
        backbone_seq, reads = generate_simulated_reads()
        aln_graph = AlnGraph(backbone_seq)
        assert len(aln_graph.nodes) == len(backbone_seq) + 2

    def test_add_alignment(self):
        aln_graph = AlnGraph("ATATTAGGC")
        alns = [((0,  9, 'A-TAGCCGGC'),   (2, 9, 'ATTA---GGC')), 
                ((0, 10, 'ATA-TACCGAG-'), (0, 9, 'ATATTA--G-GC')), 
                ((0, 10, 'ATCATCC--GGC'), (0, 9, 'AT-AT--TAGGC')), 
                ((0,  9, 'ATA-TACGGC'),   (0, 9, 'ATATTA-GGC'))]
        for aln in alns:
            aln_graph.add_alignment( aln )
        assert len(aln_graph.nodes) != 0
        assert len(aln_graph.edges) != 0


    def test_add_edge(self):
        # aln_graph = AlnGraph(backbone_seq)
        # assert_equal(expected, aln_graph.add_edge(edge))
        raise SkipTest # TODO: implement your test here

    def test_add_node(self):
        # aln_graph = AlnGraph(backbone_seq)
        # assert_equal(expected, aln_graph.add_node(node))
        raise SkipTest # TODO: implement your test here

    def test_delete_edge(self):
        # aln_graph = AlnGraph(backbone_seq)
        # assert_equal(expected, aln_graph.delete_edge(edge))
        raise SkipTest # TODO: implement your test here

    def test_delete_node(self):
        # aln_graph = AlnGraph(backbone_seq)
        # assert_equal(expected, aln_graph.delete_node(node))
        raise SkipTest # TODO: implement your test here

    def test_find_best_path(self):
        # aln_graph = AlnGraph(backbone_seq)
        # assert_equal(expected, aln_graph.find_best_path())
        raise SkipTest # TODO: implement your test here

    def test_generate_consensus(self):
        # aln_graph = AlnGraph(backbone_seq)
        # assert_equal(expected, aln_graph.generate_consensus())
        raise SkipTest # TODO: implement your test here

    def test_merge_in_nodes(self):
        # aln_graph = AlnGraph(backbone_seq)
        # assert_equal(expected, aln_graph.merge_in_nodes(nodes, node))
        raise SkipTest # TODO: implement your test here

    def test_merge_nodes(self):
        aln_graph = AlnGraph("ATAATTGGC")
        alns = [((0, 9, 'ATAG--CTGGC'), (0, 9, 'ATA-AT-TGGC')), 
                ((0, 9, 'ATAG--CTGGC'), (0, 9, 'ATA-AT-TGGC')), 
                ((0, 9, 'ATAG-TTGGC'), (0, 9, 'ATA-ATTGGC')), 
                ((0, 9, 'ATAG-TTGGC'), (0, 9, 'ATA-ATTGGC')), 
                ((0, 9, 'ATAG--CTGGC'), (0, 9, 'ATA-AT-TGGC'))]
        for aln in alns:
            aln_graph.add_alignment( aln )
        aln_graph.merge_nodes()

    def test_merge_out_nodes(self):
        # aln_graph = AlnGraph(backbone_seq)
        # assert_equal(expected, aln_graph.merge_out_nodes(node, nodes))
        raise SkipTest # TODO: implement your test here

    def test_output_consensus_fasta(self):
        # aln_graph = AlnGraph(backbone_seq)
        # assert_equal(expected, aln_graph.output_consensus_fasta(fn, rID))
        raise SkipTest # TODO: implement your test here

    def test_track_path(self):
        # aln_graph = AlnGraph(backbone_seq)
        # assert_equal(expected, aln_graph.track_path(seq, node))
        raise SkipTest # TODO: implement your test here

class TestOutputDot:
    def test_output_dot(self):
        # assert_equal(expected, output_dot(aln_graph, fn, r))
        raise SkipTest # TODO: implement your test here

class TestOutputDot2:
    def test_output_dot_2(self):
        # assert_equal(expected, output_dot_2(aln_graph, fn))
        raise SkipTest # TODO: implement your test here

class TestGenerateSimulatedReads:
    def test_generate_simulated_reads(self):
        # assert_equal(expected, generate_simulated_reads())
        raise SkipTest # TODO: implement your test here

class TestSimpleTest:
    def test_simple_test(self):
        # assert_equal(expected, simple_test())
        raise SkipTest # TODO: implement your test here

