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


"""
Classes for representing alignments as graph
"""
import re

def convert_mismatches(alnQ, alnT):
    qs = []
    ts = []
    assert len(alnQ) == len(alnT)
    alnQ = list(alnQ)
    alnT = list(alnT)
    for i in range(len(alnQ)):
        qb = alnQ[i]
        tb = alnT[i]
        if qb != tb and qb != "-" and tb != "-" :
            qs.append("-")
            ts.append(tb)
            qs.append(qb)
            ts.append("-")
        else:
            qs.append(qb)
            ts.append(tb)
    alnQ = qs
    alnT = ts

    #push gaps to the right
    for i in range(len(alnQ)-1):
        if alnT[i] == "-":
            j = i
            while 1:
                j += 1
                c = alnT[j]
                if c != "-" or j >= len(alnQ)-1:
                    break
            if c == alnQ[i]:
                alnT[i] = c
                alnT[j] = "-"

        if alnQ[i] == "-":
            j = i
            while 1:
                j += 1
                c = alnQ[j]
                if c != "-" or j >= len(alnT)-1:
                    break
            if c == alnT[i]:
                alnQ[i] = c
                alnQ[j] = "-"
    qs2 = []
    ts2 = []
    for i in range(len(alnQ)):
        if alnQ[i] != '-' or alnT[i] != '-':
            qs2.append(alnQ[i])
            ts2.append(alnT[i])

    return "".join(qs2), "".join(ts2)

class AlnEdge(object):

    def __init__(self, in_node, out_node):
        self.ID = id(self) 
        self.in_node = in_node
        self.out_node = out_node
        in_node.addout_edge(self)
        out_node.add_in_edge(self)
        self.visited = False
        self.count = 0
        self.score = 0

    def increase_count(self):
        self.count += 1

    def set_score(self, s):
        self.score = s

    def add_to_score(self, s):
        self.score += s

    def __repr__(self):
        return "(edge_ID:%d, in_node:%s, out_node:%s)" % (self.ID, self.in_node.__repr__(), self.out_node.__repr__() )

class AlnNode(object):
    def __init__(self, base):
        self.ID = id(self) 
        self.base = base
        self._in_edges = []
        self._out_edges = []
        self.is_backbone = False
        self.coverage = 0
        self.weight = 0
        self.backbone_node = None
        self.best_in_edge = None
        self.best_out_edge = None
        self.info = []

    def add_in_edge(self, in_edge):
        self._in_edges.append(in_edge)

    def addout_edge(self, out_edge):
        self._out_edges.append(out_edge)
    
    def increase_weight(self, w = 1):
        self.weight += w

    def __repr__(self):
        return "(node_id:%d, base:%s, b:%s, w:%d, c:%d)" % (self.ID, self.base, self.is_backbone, self.weight, self.coverage )

class AlnGraph(object):

    def __init__(self, backbone_seq):
        """
        AlnGraph is instantiated by giving a "backbone" sequence
        >>> g = AlnGraph("ACGTCTAT")
        """

        self.nodes = {}
        self.edges = {}
        self.backbone_nodes = {} 
        self.backbone_to_reads = {} 
        self.nodes_to_edge = {}
        self._max_node_id = 0
        self._max_edge_id = 0
        self.consensus_path = None
        self.consensus_str = None

        self.begin_node = AlnNode("B")
        self.begin_node.backbone_node = self.begin_node
        self.add_node(self.begin_node)
        self.read_range = {}
        last_node = self.begin_node

        # set up a barebone linear graph for the seed sequence 
        for pos in range(len(backbone_seq)):
            node = AlnNode(backbone_seq[pos])
            self.backbone_nodes[pos] = node
            node.is_backbone = True
            node.backbone_node = node
            node.increase_weight()
            self.add_node(node)

            edge = AlnEdge( last_node, node)
            edge.set_score(0)
            edge.increase_count()
            self.add_edge(edge)
            last_node = node

        self.end_node = AlnNode("E")
        self.end_node.backbone_node = self.end_node
        self.add_node(self.end_node)
        edge = AlnEdge( last_node, self.end_node)
        edge.set_score(0)
        edge.increase_count()
        self.add_edge(edge)

        self.backbone_node_to_pos = dict( zip( self.backbone_nodes.values(), 
                                               self.backbone_nodes.keys() ) ) 

    def add_alignment(self, aln, rId = None):

        read_pos, read_end_pos, read_aln_seq = aln[0]
        backbone_pos, backbone_end_pos, backbone_aln_seq = aln[1]

        if rId != None:
            self.read_range[rId] = (backbone_pos, backbone_end_pos)

        last_node = self.begin_node

        for aln_pos in range( len( read_aln_seq ) ):
            read_base = read_aln_seq[aln_pos]
            backbone_base = backbone_aln_seq[aln_pos]
            if read_base == backbone_base:
                node = self.backbone_nodes[backbone_pos]
                node.increase_weight()
                node.backbone_node.coverage += 1

                if (last_node.ID, node.ID) not in self.nodes_to_edge:
                    edge = AlnEdge( last_node, node)
                    self.add_edge(edge)
                    edge.increase_count()
                else:
                    self.nodes_to_edge[ (last_node.ID, node.ID) ].increase_count()

                read_pos += 1
                backbone_pos += 1
                last_node = node
                if rId != None:
                    node.info.append(rId)

            elif read_base == "-" and backbone_base != "-":
                node = self.backbone_nodes[ backbone_pos ]
                node.backbone_node.coverage += 1
                backbone_pos += 1

            elif read_base != "-" and backbone_base == "-":
                node = AlnNode(read_base)
                node.increase_weight()
                node.backbone_node = self.backbone_nodes[ backbone_pos ]

                self.add_node(node)

                edge = AlnEdge( last_node, node)
                self.add_edge(edge)
                edge.increase_count()

                read_pos += 1
                last_node = node
                if rId != None:
                    node.info.append(rId)

            self.backbone_to_reads.setdefault(backbone_pos,set())

            if rId != None:
                self.backbone_to_reads[backbone_pos].add(rId)

        if (last_node.ID, self.end_node.ID) not in self.nodes_to_edge:
            edge = AlnEdge( last_node, self.end_node)
            edge.set_score(0)
            edge.increase_count()
            self.add_edge(edge)

    def add_edge(self, edge):
        edge.ID = id(edge)
        self.edges[edge.ID] = edge
        self.nodes_to_edge[ (edge.in_node.ID, edge.out_node.ID) ] = edge

    def add_node(self, node):
        self.nodes[node.ID] = node

    def delete_edge(self, edge):
        n1 = edge.in_node
        n2 = edge.out_node
        n1._out_edges.remove(edge)
        n2._in_edges.remove(edge)
        del self.edges[edge.ID]
        del edge

    def delete_node(self, node):
        for in_edge in node._in_edges:
            self.delete_edge(in_edge)
        for out_edge in node._out_edges:
            self.delete_edge(out_edge)
        del self.nodes[node.ID]
        del node

    def merge_in_nodes(self, node):

        node_groups ={}
        for in_edge in node._in_edges:
            in_node = in_edge.in_node
            if len(in_node._out_edges) == 1:
                node_groups.setdefault( in_node.base, [])
                node_groups[in_node.base].append(in_node)

        for b in node_groups:
            if len(node_groups[b]) <= 1:
                continue

            nodes = node_groups[b]

            for n in nodes[1:]:
                nodes[0]._out_edges[0].count += n._out_edges[0].count
                nodes[0].increase_weight(n.weight)
                nodes[0].info.extend(n.info)

            for n in nodes[1:]:
                for in_edge in n._in_edges:
                    n1 = in_edge.in_node
                    if (n1.ID, nodes[0].ID) not in self.nodes_to_edge:
                        e = AlnEdge(n1, nodes[0])
                        self.add_edge(e)
                        e.count = self.nodes_to_edge[ (n1.ID, n.ID) ].count
                        e.visited = self.nodes_to_edge[ (n1.ID, n.ID) ].visited 
                    else:
                        e = self.nodes_to_edge[ (n1.ID, nodes[0].ID) ]
                        e.count += self.nodes_to_edge[ (n1.ID, n.ID) ].count
                        e.add_to_score( self.nodes_to_edge[ (n1.ID, n.ID) ].score )
                self.delete_node(n)

            self.merge_in_nodes(nodes[0])


    def merge_out_nodes(self, node):

        node_groups ={}

        for out_edge in node._out_edges:
            out_node = out_edge.out_node
            if len(out_node._in_edges) == 1:
                node_groups.setdefault( out_node.base, [])
                node_groups[out_node.base].append(out_node)

        for b in node_groups:
            if len(node_groups[b]) <= 1:
                continue

            nodes = node_groups[b]

            for n in nodes[1:]:
                nodes[0]._in_edges[0].count += n._in_edges[0].count
                nodes[0].increase_weight(n.weight)
                nodes[0].info.extend(n.info)

            for n in nodes[1:]:
                for out_edge in n._out_edges:
                    n2 = out_edge.out_node
                    if (nodes[0].ID, n2.ID) not in self.nodes_to_edge:
                        e = AlnEdge(nodes[0], n2)
                        self.add_edge(e)
                        e.count = self.nodes_to_edge[ (n.ID, n2.ID) ].count 
                        e.visited = self.nodes_to_edge[ (n.ID, n2.ID) ].visited 
                    else:
                        e = self.nodes_to_edge[ (nodes[0].ID, n2.ID) ]
                        e.count += self.nodes_to_edge[ (n.ID, n2.ID) ].count
                        e.add_to_score( self.nodes_to_edge[ (n.ID, n2.ID) ].score )

                self.delete_node(n)

    def merge_nodes(self):
        for e in self.edges.values():
            e.visited = False

        seed_nodes = []

        for node_id, node in self.nodes.items():
            if len(node._in_edges) == 0:
                seed_nodes.append(node)
        
        while 1:
            if len(seed_nodes) == 0: 
                break
            node = seed_nodes.pop(0)
            self.merge_in_nodes(node)
            self.merge_out_nodes(node)
           
            for out_edge in node._out_edges:
                out_node = out_edge.out_node
                out_edge.visited = True
                in_edges = [ e for e in out_node._in_edges if e.visited == False ]
                if len(in_edges) == 0:
                    seed_nodes.append(out_node)

    def find_best_path(self):
        seed_nodes = []
        node_score = {}
        best_node_score_edge = {}

        for node_id, node in self.nodes.items():
            if len(node._out_edges) == 0:
                seed_nodes.append(node)
                node_score[node.ID] = 0 
                node.best_out_edge = None
                node.best_in_edge = None

        for edge in self.edges.values():
            edge.visited = False

        while 1:
            if len(seed_nodes) == 0: 
                break
            
            node = seed_nodes.pop(0)
            best_edge = None
            best_score = None
            for out_edge in node._out_edges:
                out_node = out_edge.out_node

                score = node_score[out_node.ID]  

                backbone_node_coverage = out_node.backbone_node.coverage 
                if out_node.is_backbone == True and out_node.weight == 1:
                    new_score = score  - 10
                else:
                    new_score = out_edge.count  - backbone_node_coverage * 0.5 + score
                if new_score > best_score or best_score == None:
                    best_edge = out_edge
                    best_score = new_score
                    
            if best_edge != None:
                best_node_score_edge[node.ID] = best_edge
                node_score[node.ID] = best_score

            for in_edge in node._in_edges:
                in_node = in_edge.in_node
                in_edge.visited = True
                out_edges = [ e for e in in_node._out_edges if e.visited == False ]
                if len(out_edges) == 0:
                    seed_nodes.append(in_node)

        node = self.begin_node
        consensus_path = []

        while 1:
            consensus_path.append(node)
            if node.ID not in best_node_score_edge:
                break
            else:
                best_out_edge = best_node_score_edge[node.ID]
                node.best_out_edge = best_out_edge
                best_node_score_edge[node.ID].out_node.best_in_edge = best_node_score_edge[node.ID]
                node = best_node_score_edge[node.ID].out_node
                node.best_in_edge = best_out_edge

        return consensus_path


    def generate_consensus(self, min_cov = 8):

        self.merge_nodes()

        if self.consensus_path == None:
            self.consensus_path = self.find_best_path()

        s = []
        c = []
        cov_good = []
        for n in self.consensus_path:
            if n not in [self.begin_node, self.end_node]:
                s.append(n.base)
                if len(n.info) >= min_cov:
                    #print len(n.info), n.info
                    cov_good.append("1")
                else:
                    cov_good.append("0")

                rn = n.weight
                if n.is_backbone == True:
                    rn -= 1
                if n.best_out_edge != None:
                    en = n.best_out_edge.count
                else:
                    en = 0
                if n.best_in_edge != None:
                    en2 = n.best_in_edge.count
                else:
                    en2 = 0
                c.append( (rn, en2, en) )

        s = "".join(s)
        cov_good = "".join(cov_good)

        p = re.compile("1+")
        best_range = (0, 0)
        best_l = 0
        for m in p.finditer(cov_good):
            b, e = m.start(), m.end()
            if e-b > best_l:
                best_range = (b,e)
                best_l = e - b
        b,e = best_range
        if e-b != 0:
            self.consensus_str = s[b:e]
            c = c[b:e]
        else:
            self.consensus_str = ""
            c = []
        return self.consensus_str, c


    def get_high_entropy_nodes(self, 
                               ignore_backbone = False, 
                               coverage_th = 20,
                               overlap_th = 10,
                               entropy_th = 0.5):
        from math import log 
        if self.consensus_path == None:
            self.generate_consensus()

        backbone_node_to_pos = self.backbone_node_to_pos  

        node_entropy = []

        for node_id in self.nodes:

            node = self.nodes[node_id]

            if node.backbone_node.coverage < coverage_th: continue
            if ignore_backbone and node.is_backbone: continue
                
            p = 1.0*(len(node.info)+1)/(node.backbone_node.coverage+1)
            if abs(p-1) < 1e-5 or p > 1:
                ent = 0
            else:
                ent = - p*log(p) - (1-p)*log(1-p)

            node_entropy.append( [ node_id, node, ent ] )

        node_entropy.sort( key = lambda x:-x[2] )

        high_entropy_nodes = [ n for n in node_entropy if n[2] > entropy_th]

        return node_entropy, high_entropy_nodes

    def get_sorted_nodes(self):
        seed_nodes = []
        for node_id, node in self.nodes.items():
            if len(node._in_edges) == 0:
                seed_nodes.append(node)

        for edge in self.edges.values():
            edge.visited = False

        sorted_nodes = []

        while 1:
            if len(seed_nodes) == 0:  
                break

            node = seed_nodes.pop(0)

            for in_edge in node._in_edges:
                in_node = in_edge.in_node

            sorted_nodes.append(node)

            for out_edge in node._out_edges:
                out_node = out_edge.out_node
                out_edge.visited = True
                in_edges = [ e for e in out_node._in_edges if e.visited == False ]
                if len(in_edges) == 0:
                    seed_nodes.append(out_node)

        return sorted_nodes

    def get_read_node_vector(self, entropy_th = 0.5):
        
        ne, hne = self.get_high_entropy_nodes(coverage_th = 0, entropy_th = entropy_th)
        node_to_entropy = dict( [ (v[1],v[2]) for v in ne ] )
        read_ids = set()
        
        for n in self.nodes.values():
            for r in n.info:
                read_ids.add(r)
                
        backbone_node_to_pos = self.backbone_node_to_pos  

        sn = self.get_sorted_nodes()
        high_entropy_nodes = [ (n, 
                                backbone_node_to_pos[n.backbone_node], 
                                node_to_entropy[n]) for n in sn if node_to_entropy[n] > entropy_th]
        
        read_to_nodes = {}
        for r_id in read_ids:
            read_to_nodes[r_id] = ["-"] * len(high_entropy_nodes)
            for i,n in enumerate(high_entropy_nodes):
                if r_id in n[0].info:
                    read_to_nodes[r_id][i] = n[0].base
            

        return read_to_nodes, high_entropy_nodes

    def output_consensus_fasta(self,  fn, rID):
        f = open(fn, "a")
        print >>f,">%s" % rID
        for i in range(0, len(self.consensus_str), 60):
            print >>f, self.consensus_str[i:i+60]
        f.close()

