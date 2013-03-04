from utils import sorted_node_data

def output_dot(aln_graph, fn, r = None):
    f = open(fn, "w")
    try:
        c_path, c_str = aln_graph.consensus_path, aln_graph.consensus_str
    except:
        c_path = []

    backbone_node_to_pos = aln_graph.backbone_node_to_pos 

    print >>f, "digraph X { rankdir=LR;"
    for node_id, node in aln_graph.nodes.items():
        x = backbone_node_to_pos.get(node.backbone_node, -1)
        #if node.weight < 10: continue
        if r != None:
            if x < r[0] or x > r[1]:
                continue

        if node in c_path:
            if node in aln_graph.backbone_nodes.values():
                fcol = "yellow"
                col = "blue"
            else:
                fcol = "yellow"
                col = "green"
        else:
            if node in aln_graph.backbone_nodes.values():
                fcol = "white"
                col = "blue"
            else:
                fcol = "white"
                col = "green"
        if node in aln_graph.backbone_nodes.values():
            print >>f,'%d [label="%d%s,%d" color=%s fillcolor=%s style=filled]' % (node.ID, 
                                        node.weight, node.base, x, col, fcol) ,";"
        else:
            if node.base not in ["B", "E"]:
                x = backbone_node_to_pos[node.backbone_node]
                print >>f,'%d [label="%d%s,%d" color=%s fillcolor=%s style=filled]' % (node.ID, node.weight, node.base, x, col, fcol ) ,";"
            else:
                print >>f,'%d [label="%d%s" color=%s fillcolor=%s, style=filled]' % (node.ID, node.weight, node.base, col, fcol ) ,";"
        

    for edge_ID, edge in aln_graph.edges.items():
        in_node  = edge.in_node
        out_node = edge.out_node
        x = backbone_node_to_pos.get( in_node.backbone_node, -1)
        if r != None:
            if x < r[0] or x > r[1]:
                continue
        x = backbone_node_to_pos.get(out_node.backbone_node, -1)
        if r != None:
            if x < r[0] or x > r[1]:
                continue

        count = edge.count
        score = edge.count * 0.2
        if score < 1: score = 1
        if in_node.best_out_edge != None and out_node == in_node.best_out_edge.out_node:
            print >>f, '%d -> %d [label="%d", weight=%f, penwidth=%d, color=black];' % (in_node.ID, out_node.ID, edge.count, edge.score, score)
        else:
            print >>f, '%d -> %d [label="%d", weight=%f, penwidth=%d, color=grey];' % (in_node.ID, out_node.ID, edge.count, edge.score, score)
    print >>f,"}"
    f.close()

def output_dot_2(aln_graph, fn):
    f = open(fn, "w")
    c_path, c_str = aln_graph.consensus_path, aln_graph.consensus_str
    backbone_node_to_pos = aln_graph
    nodeEntropy = []
    for nodeId in aln_graph.nodes:
        node = aln_graph.nodes[nodeId]
        p = 1.0*(len(node.info)+1)/(node.backbone_node.coverage+1)
        if abs(p-1) < 1e-5 or p > 1:
            ent = 0 
        else:
            ent = - p*log(p) - (1-p)*log(1-p)
        nodeEntropy.append( [ nodeId, node, ent ] )
    nodeEntropy.sort( key = lambda x:-x[2] )

    high_entropy_nodes = [ n for n in nodeEntropy if n[2] > 0.5]
    print high_entropy_nodes
    print >>f, "graph X {"
    for node_id, node, ent in high_entropy_nodes:
        x = backbone_node_to_pos[node.backbone_node]
        if node.is_backbone:
            print >>f,'%d [label="%d,%d%s,%d,%f", color="red" ];' % (node.ID, node.ID, node.weight, node.base, x, ent) ,";"
        else:
            print >>f,'%d [label="%d,%d%s,%d,%f" ];' % (node.ID, node.ID, node.weight, node.base, x, ent) ,";"
    
    for node_id_1, node_1, ent_1 in high_entropy_nodes:
        for node_id_2, node_2, ent_2 in high_entropy_nodes:
            if node_id_1 == node_id_2: continue
            if node_id_1 > node_id_2: continue
            s1 = set(node_1.info)
            s2 = set(node_2.info)
            s3 = aln_graph.backbone_to_reads[backbone_node_to_pos[node_1.backbone_node]]
            s4 = aln_graph.backbone_to_reads[backbone_node_to_pos[node_2.backbone_node]]
            s34 = list( s3 & s4 )
            xvec = np.array( [ 1 if r in s1 else 0 for r in s34], dtype=np.uint32 )
            yvec = np.array( [ 1 if r in s2 else 0 for r in s34], dtype=np.uint32 )
            distance = phi_coeff(xvec, yvec) 

            if distance > 0.5 :
                print >>f,"%d -- %d [color=black, len=%f, penwidth=%d];" % (node_1.ID, node_2.ID, distance, 4*distance) 
        
    print >>f,"}"

def output_read_network(aln_graph, 
                        dot_fn = None, 
                        read_graph_fn = None, 
                        consensus_fn = None, 
                        ignore_backbone = False, 
                        coverage_th = 20,
                        overlap_th = 10,
                        entropy_th = 0.5):

    node_entropy, high_entropy_nodes = \
    aln_graph.get_high_entropy_nodes(ignore_backbone = ignore_backbone, 
                                     coverage_th = coverage_th,
                                     overlap_th = overlap_th,
                                     entropy_th = entropy_th) 

    s = aln_graph.consensus_str

    if consensus_fn != None:
        with open(consensus_fn, "w") as con_f:
            print >>con_f, ">consensus1"
            print >>con_f, s
            
    if dot_fn != None:
        f = open(dot_fn, "w")
        print >>f, "graph X {rankdir=LR;"
        r2n = {}
        for node_id, node, ent in high_entropy_nodes:
            for r in node.info:
                r2n.setdefault(r, set())
                r2n[r].add(node_id)
                col = "blue"
                r_short = r.split("/")[1]
                print >>f,'%s [label="%s" color=%s];' % (r_short, r_short, col)

        for r1 in r2n:
            for r2 in r2n:
                if r1 <= r2: continue
                s1 = r2n[r1]
                s2 = r2n[r2]
                overlap = len(s1 & s2)
                r1_short = r1.split("/")[1]
                r2_short = r2.split("/")[1]
                if overlap > overlap_th:
                    print >>f,"%s -- %s [color=black];" % (r1_short, r2_short)
        print >>f,"}"
        f.close()

    if read_graph_fn != None:
        f2 = open(read_graph_fn, "w")
        r2n = {}
        for node_id, node, ent in high_entropy_nodes:
            for r in node.info:
                r2n.setdefault(r, set())
                r2n[r].add(node_id)

        for r1 in r2n:
            for r2 in r2n:
                if r1 <= r2: continue
                s1 = r2n[r1]
                s2 = r2n[r2]
                overlap = len(s1 & s2)
                print >>f2, r1, r2, overlap
        f2.close()

def dump_sorted_node_data(g, entropy_th = 0.5, interval = None):
    data = sorted_node_data(g, entropy_th = entropy_th, interval = interval)
    for d in data:
        print "%06d" % d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]

