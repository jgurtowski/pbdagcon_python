
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
