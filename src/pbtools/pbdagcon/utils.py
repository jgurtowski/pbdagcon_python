from math import log, sqrt
import numpy as np
from .aligngraph import convert_mismatches, AlnGraph

class AlignGraphUtilError(Exception):
    pass

def phi_coeff(xvec, yvec):
    nx1 = sum( xvec == 1 )
    nx0 = len(xvec) - nx1
    ny1 = sum( yvec == 1 )
    ny0 = len(yvec) - ny1
    xyvec = (xvec << 1) + yvec
    n11 = sum( xyvec == 3 )
    n10 = sum( xyvec == 2 )
    n01 = sum( xyvec == 1 )
    n00 = sum( xyvec == 0 )
    return (n11 * n00 - n10 * n01) / sqrt( nx1 * nx0 * ny1 * nx0 + 1)

class Simple_Alignment_Hit(object):
    """ A simple class to wrap the output of the blasr "-m 5" option """
    def __init__(self, rm5_line):
        rm5_line = rm5_line.strip().split()
        self.query_id = rm5_line[0]
        self.query_length = int(rm5_line[1])
        self.query_start = int(rm5_line[2])
        self.query_end = int(rm5_line[3])
        self.query_strand = rm5_line[4]
        self.target_id = rm5_line[5]
        self.target_length = int(rm5_line[6])
        self.target_start = int(rm5_line[7])
        self.target_end = int(rm5_line[8])
        self.target_strand = rm5_line[9]
        self.alignedQuery = rm5_line[16]
        self.alignedTarget = rm5_line[18]

def simple_align_hit_iterator(rm1_fn, ref_group = None):
    with open(rm1_fn) as f:
        for l in f:
            ll = l.strip().split()
            if ref_group != None and ll[5] != ref_group:
                continue
            yield Simple_Alignment_Hit(l)

def get_aln_array(aln_iter,
                  max_num_reads = None, 
                  remove_in_del = False, 
                  min_length = None):

    rMap = dict(zip("ACGTacgtNn-", "TGCAtgcaNn-"))

    alns = []
    nread = 0
    i = 0
    bp = {}
    reads = {}

    for h in aln_iter:

        if min_length != None:
            if h.target_end - h.target_start < min_length:
                continue

        nread += 1

        if max_num_reads != None:
            if nread > max_num_reads:
                break

        read_id = h.query_id 

        if read_id in reads: continue #avoid duplicate readId in the data
        
        ts = h.target_start
        te = h.target_end
        qs = h.query_start
        qe = h.query_end

        alnT, alnQ = h.alignedTarget.upper(), h.alignedQuery.upper()
        dqPos = 1
        
        if h.target_strand == '-':
            alnT = "".join([rMap[c] for c in alnT[::-1]])
            alnQ = "".join([rMap[c] for c in alnQ[::-1]])
            dqPos = -1
            qs, qe = qe - 1, qs - 1

        if remove_in_del:
            aln_pair = zip(alnQ, alnT)
            new_aln_pair = []
            for b1, b2 in aln_pair:
                if b1 == "-":
                    b1 = b2
                if b2 != "-":
                    new_aln_pair.append( (b1,b2) )
            alnQ, alnT = zip(*new_aln_pair)
            alnQ = "".join(alnQ)
            alnT = "".join(alnT)

        alnQ, alnT = convert_mismatches(alnQ,alnT)

        if alnQ[0] == "-" or alnT[0] == "-":
            continue
        if alnQ[-1] == "-" or alnT[-1] == "-":
            continue

        #print h.target_strand 
        #print "r %10d %10d" % (qs, qe), alnQ
        #print "t %10d %10d" % (ts, te), alnT
        #print

        bp[read_id] = ( (qs, qe), ts)
        alns.append(  ( ( qs, qe, alnQ ), ( ts, te, alnT ) , read_id ) )
        tPos = ts
        qPos = qs
        assert len(alnQ) == len(alnT)
        reads[read_id] = alnQ.replace("-","")
        i += 1

    return alns

def constructe_aln_graph_from_cmph5(cmpH5_fn, 
                                    backbone_fasta_fn, 
                                    max_num_reads = None, 
                                    remove_in_del = False, 
                                    ref_group = None,
                                    min_length = None):

    from pbcore.io import cmph5
    from pbcore.io.cmph5 import CmpH5

    cmpType = CmpH5.SFCmpH5
    cmph5f = cmph5.factory.create( cmpH5_fn, "r", cmpType=cmpType )

    aln_hit_iterator = None

    if ref_group != None:
        aln_hit_iterator = cmph5f["/%s" % ref_group].alnHitIterator()
    else:
        aln_hit_iterator = cmph5f.alnHitIterator()

    alns = get_aln_array(aln_hit_iterator, 
                         max_num_reads = max_num_reads, 
                         remove_in_del = remove_in_del, 
                         min_length = min_length)
    
    backboneSeq = open(backbone_fasta_fn).read()
    backboneSeq = "".join(backboneSeq.split()[1:])

    g = AlnGraph(backboneSeq)

    i = 0
    for aln in alns:
        rId = aln[2]
        aln = aln[0:2]
        g.add_alignment( aln, "%s" % rId)

    cmph5.factory.close(cmph5f)
    return g

def constructe_aln_graph_from_fasta(read_fasta_fn, 
                                    backbone_fasta_fn, 
                                    max_num_reads = None, 
                                    remove_in_del = False, 
                                    ref_group = None,
                                    min_length = None):

    import os

    os.system("blasr %s %s -m 5 -nproc 16 -out %s" % (read_fasta_fn, backbone_fasta_fn, read_fasta_fn+".aln_unsorted"))
    os.system("cat %s | sort > %s" % (read_fasta_fn+".aln_unsorted", read_fasta_fn+".aln" ))
    
    aln_hit_iterator = simple_align_hit_iterator(read_fasta_fn+".aln", ref_group = ref_group)

    alns = get_aln_array(aln_hit_iterator, 
                         max_num_reads = max_num_reads, 
                         remove_in_del = remove_in_del, 
                         min_length = min_length)
        
    backboneSeq = open(backbone_fasta_fn).read()
    backboneSeq = "".join(backboneSeq.split()[1:])

    g = AlnGraph(backboneSeq)

    i = 0
    for aln in alns:
        rId = aln[2]
        rId = rId.split("/")[0]
        aln = aln[0:2]
        g.add_alignment( aln, "%s" % rId)

    return g

def sorted_nodes(g):
    return g.get_sorted_nodes()

def read_node_vector(g, entropy_th = 2.5):
    read_to_nodes, high_entropy_nodes = g.get_read_node_vector(entropy_th=entropy_th)
    return read_to_nodes, high_entropy_nodes

def clustering_read( read_to_nodes, high_entropy_nodes, k_cluster = 2, random_seed = 42, cleanup_th = 0.5):
    
    import random
    random.seed(random_seed)

    cluster = {}
    read_to_binary_vector = {}
    count = 0
    for r in read_to_nodes:   
        read_to_binary_vector[r] = np.array([ 1 if c != "-" else -1 for c in read_to_nodes[r] ])
        k = random.randint(0, k_cluster-1)
        #k = count % k_cluster
        cluster.setdefault(k,[])
        cluster[k].append(r)
        count += 1
    
    n_iteration = 0
    cluster_vec = {}
    while 1:
        for k in range(k_cluster):
            new_vec = np.zeros(len(high_entropy_nodes),dtype=np.float)
            #print k, len(cluster[k])
            for r in cluster[k]:
                new_vec += read_to_binary_vector[r]
            new_vec /=  (len(cluster[k])+1)
            cluster_vec[k] = np.array([1 if v>=0 else -1 for v in new_vec])

        #avoid degnerated cluster    
        for k in range(k_cluster):
            for j in range(k+1, k_cluster):
                if sum(cluster_vec[j] == cluster_vec[k]) > 0.5*len(high_entropy_nodes):
                    cluster_vec[j] = np.array([random.choice([-1,1]) for v in new_vec])

        cluster = {}
        for r in read_to_nodes:   
            distances = []
            for k in range(k_cluster):
                cluster.setdefault(k,[])
                distances.append( (sum(read_to_binary_vector[r] * cluster_vec[k]), k) )

            distances.sort()

            cluster[distances[-1][1]].append(r)
        
        n_iteration += 1
        if n_iteration  > 10:
            break
            
    cluster = {}
    for r in read_to_nodes:   
        distances = []
        for k in range(k_cluster):
            cluster.setdefault(k,[])
            distances.append( (sum(read_to_binary_vector[r] * cluster_vec[k]), k) )

        distances.sort()
        if distances[-1][0] > cleanup_th*len(high_entropy_nodes):
            cluster[distances[-1][1]].append(r)
    return cluster, cluster_vec

def find_best_aligned_read(cmpH5_fn, ref_group=None):

    cmpType = CmpH5.SFCmpH5
    cmph5f = cmph5.factory.create( cmpH5_fn, "r", cmpType=cmpType )
    rMap = dict(zip("ACGTacgtNn-", "TGCAtgcaNn-"))
    read_data = []
    all_read = []
    aln_hit_iterator = None
    if ref_group != None:
        aln_hit_iterator = cmph5f["/%s" % ref_group].alnHitIterator
    else:
        aln_hit_iterator = cmph5f.alnHitIterator
    for h in aln_hit_iterator():
        ref_len = h.target_end - h.target_start
        q_len = h.query_end - h.query_start
        rId = "/".join(h.query_id.split("/")[:2])
        read = h.alignedQuery.replace("-","").upper()
        #if h.target_strand == '-':
        #    read = "".join([rMap[c] for c in read[::-1]])
        
        
        read_data.append( ( ref_len/10, 100*h.nMatch/len(h.alignedTarget), q_len/10, rId, read) )
        all_read.append( (rId, read ) )
    read_data.sort()
 
    cmph5.factory.close(cmph5f)
    return read_data[-1], all_read
        
def best_template_by_blasr(fasta_fn, len_threshold = 300, min_number_reads = 1):
    import os
    import numpy as np
    from pbcore.io import FastaIO
    f = FastaIO.SimpleFastaReader(fasta_fn)
    read_dict = {}
    for r in f:
        r_id = r.name.split("/")[0]
        read_dict[r_id] = r.sequence

    
    rtn = os.system("blasr -nproc 8 -bestn 10 -nCandidates 50 -m 1 %s %s -out %s" % (fasta_fn, fasta_fn, fasta_fn+".saln"))
    if rtn != 0:
        return None
    
    scores = {}
    with open(fasta_fn + ".saln") as f:
        for l in f:
            l = l.strip().split()
            l = l[:9]
            r1, r2 = l[:2]
            r1 = r1.split("/")[0]
            r2 = r2.split("/")[0]
            if r1 == r2:
                continue
            if int(l[7]) - int(l[6]) < len_threshold : continue
            scores.setdefault(r1,[])
            #scores[r1].append(int(l[4]))
            scores[r1].append(-float(l[5]))
    score_array = []
    for r, s in scores.items():
        score_array.append( (np.mean(s), r) )
    score_array.sort()

    if min_number_reads != None:
        if len(score_array) < min_number_reads:
            raise AlignGraphUtilError("not enough number of reads in best_template_by_blasr")
        else:
            return score_array[0][1], read_dict[score_array[0][1]]
    else:
        return r_id, read_dict[r_id]

def sorted_node_data(aln_graph, entropy_th = 0.5, interval = None):
    ne, hne = aln_graph.get_high_entropy_nodes(coverage_th=0)
    node_to_entropy = dict( [ (v[1],v[2]) for v in ne ] ) 
    read_ids = set()
    for n in aln_graph.nodes.values():
        for r in n.info:
            read_ids.add(r)
    read_id_to_pos = dict(( (x[1],x[0]) for x in enumerate(list(read_ids))) )

    backbone_node_to_pos =  aln_graph.backbone_node_to_pos 

    data = []
    for n in sorted_nodes(aln_graph):
        s = [" "] * len(read_id_to_pos)
        for r in n.info:
            s[read_id_to_pos[r]] = n.base

        if n.base not in ["B", "E"]:
            bpos = backbone_node_to_pos[n.backbone_node]
            if interval != None and (bpos < interval[0] or bpos > interval[1]):
                continue
        ent = node_to_entropy[n] if n in node_to_entropy else 0
        if n.base not in ["B","E"] and ent >= entropy_th:
            data.append ( (  backbone_node_to_pos[n.backbone_node],\
                  "+" if n in aln_graph.consensus_path else "-",\
                  "+" if n.is_backbone == True else "-",\
                  n.base, "".join(s), len(n.info),\
                  n.backbone_node.coverage,\
                  node_to_entropy[n] if n in node_to_entropy else "-" ) )
    return data

def detect_missing(aln_graph, entropy_th = 0.66, interval = None):
    data = sorted_node_data(aln_graph, entropy_th = 0, interval = interval)
    s = []
    for d in data:
        if d[7] > entropy_th and d[2] == "-":
            s.append(d[3].lower())
            continue
        if d[2] =="+" and d[7] > entropy_th:
            #pass
            s.append(d[3].lower())
        elif d[2] == "+":
            s.append(d[3].upper())
    return "".join(s)

def get_subset_reads(fasta_fn ,cluster_dict, cluster_index, out_file_name):
    from pbcore.io import FastaIO
    f = FastaIO.SimpleFastaReader(fasta_fn)
    
    with open(out_file_name,"w") as out_f:
        for r in f:
            read_id = r.name
            read_seq = r.sequence.upper()

            if read_id in cluster_dict[cluster_index]:
                print >>out_f, ">"+r.name
                print >>out_f, r.sequence
