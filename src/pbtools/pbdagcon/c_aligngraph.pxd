cdef class AlnNode

cdef class AlnEdge:
    cdef int count
    cdef int score
    cdef unsigned long int ID
    cdef AlnNode in_node
    cdef AlnNode out_node
    cdef bint visited

cdef class AlnNode:
    cdef int weight
    cdef base
    cdef unsigned long int ID
    cdef list _in_edges
    cdef list _out_edges
    cdef bint is_backbone
    cdef int coverage
    cdef AlnNode backbone_node
    cdef AlnEdge best_in_edge
    cdef AlnEdge best_out_edge
    cdef list info


cdef class AlnGraph:
    cdef dict nodes
    cdef dict edges
    cdef dict backbone_nodes
    cdef dict backbone_to_reads
    cdef dict nodes_to_edge
    cdef dict backbone_node_to_pos
    cdef AlnNode begin_node
    cdef AlnNode end_node
    cdef list consensus_path
    cdef int _max_node_id
    cdef int _max_edge_id
    cdef consensus_str
    cdef read_range
