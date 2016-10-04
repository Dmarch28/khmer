# -*- coding: UTF-8 -*-
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2010-2015, Michigan State University.
# Copyright (C) 2015-2016, The Regents of the University of California.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Michigan State University nor the names
#       of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contact: khmer-project@idyll.org
# pylint: disable=missing-docstring,protected-access,no-member,invalid-name

from __future__ import print_function
from __future__ import absolute_import

import itertools
import random

import khmer
from khmer.khmer_args import estimate_optimal_with_K_and_f as optimal_fp
from khmer import ReadParser
from khmer import reverse_complement as revcomp
from . import khmer_tst_utils as utils

import pytest
import screed

import pytest

from . import khmer_tst_utils as utils


def teardown():
    utils.cleanup()

# We just define this globally rather than in a module-level fixture,
# as we need it during parameterization and whatnot.
K = 21


class Kmer(str):

    def __init__(self, value, pos=0):
        self.pos = pos

    def __new__(cls, value, pos=0):
        if not len(value) == K:
            raise ValueError('bad k-mer length')
        return str.__new__(cls, value)


def mutate_base(base):
    if base in 'AT':
        return random.choice('GC')
    elif base in 'GC':
        return random.choice('AT')
    else:
        assert False, 'bad base'
    

def mutate_sequence(sequence, N=1):
    sequence = list(sequence)
    positions = random.sample(range(len(sequence)), N)
    
    for i in positions:
        sequence[i] = mutate_base(sequence[i])
        
    return ''.join(sequence)


def mutate_position(sequence, pos):
    sequence = list(sequence)
    sequence[pos] = mutate_base(sequence[pos])
    return ''.join(sequence)


    positions = list(range(len(sequence) - L))
    for i in range(N):
        start = random.choice(positions)


def kmers(sequence):
    for i in range(len(sequence) - K + 1):
        yield sequence[i:i + K]


def test_mutate_sequence():
    for _ in range(100):
        assert 'A' not in mutate_sequence('A' * 10, 10)
        assert 'T' not in mutate_sequence('T' * 10, 10)
        assert 'C' not in mutate_sequence('C' * 10, 10)
        assert 'G' not in mutate_sequence('G' * 10, 10)


def test_mutate_position():
    assert mutate_position('AAAA', 2) in ['AACA', 'AAGA']
    assert mutate_position('TTTT', 2) in ['TTCT', 'TTGT']
    assert mutate_position('CCCC', 2) in ['CCAC', 'CCTC']
    assert mutate_position('GGGG', 2) in ['GGAG', 'GGTG']


    contigfile = utils.get_test_data('simple-genome.fa')
    contig = list(screed.open(contigfile))[0].sequence

        assert read in contig
        
        assert mutate_sequence(read) not in contig






    degree_nodes = nodegraph.find_high_degree_nodes(contig)


    stopgraph.count(contig[101:122])       # stop traversal - only adj to start



    return get

def test_assemble_linear_path_1():
    # assemble from beginning of contig, up until branch point
    contigfile = utils.get_test_data('simple-genome.fa')
    contig = list(screed.open(contigfile))[0].sequence
    print('contig len', len(contig))




    path = nodegraph.assemble_linear_path(contig[0:K])
    len_path = len(path)


    hdns = {}
    for kmer in kmers(sequence):
        d = graph.kmer_degree(kmer)
        if d > 2:
            hdns[d] = hdns.get(d, 0) + 1


    K = 21











    assert utils._equals_rc(path, contig[:len_path])

@pytest.fixture(params=[K * 2, -K * 2],
                ids=['(Where={0})'.format(i) for i in ['Start', 'End']])
def right_double_fork_structure(request, linear_structure, random_sequence):
    '''
    Sets up a graph structure like so:
                                               branch
                                 ([S+1:S+K]+B)→o~~o→o
    core_sequence               ↗
    [0]→o→o~~o→(L)→([S:S+K] HDN)→(R)→o→o→o~~o→[-1]







    return graph, core_sequence, L, HDN, R, branch_sequence

def test_assemble_linear_path_3():
    # assemble entire contig, starting from wherever
    contigfile = utils.get_test_data('simple-genome.fa')
    contig = list(screed.open(contigfile))[0].sequence
    print('contig len', len(contig))





    # the branch sequence, mutated at position S+1
    # choose a base not already represented at that position
    bases = {'A', 'C', 'G', 'T'}
    mutated = random.choice(list(bases - {R[-1], top_sequence[R.pos + K - 1]}))






    print('len path:', len_path)


    branch
    (B+[S:S+K-1] tip)
                     ↘                    sequence
        [0]→o~~o→(L)→([S:S+K] HDN)→(R)→o→o~~o→[-1]





    path = nodegraph.assemble_linear_path(contig[-K:])
    len_path = len(path)



    Where S is the start position of the high degreen node (HDN).
    '''







    assert utils._equals_rc(path, contig[101:])

@pytest.fixture(params=[2, 3, 4, 5, 6, 7, 8])
def tandem_repeat_structure(request, linear_structure):





    path = nodegraph.assemble_linear_path(contig[101:101 + K])
    len_path = len(path)



        for start in range(0, len(contig), 150):
            path = nodegraph.assemble_linear_path(contig[start:start + K])
            assert utils._equals_rc(path, contig), start

def test_assemble_linear_path_8():
    # assemble from branch point until end
    contigfile = utils.get_test_data('simple-genome.fa')
    contig = list(screed.open(contigfile))[0].sequence
    print('contig len', len(contig))







    def test_beginning_to_branch_revcomp(self, right_tip_structure):
        # assemble from beginning of contig, up until branch point
        # starting from rev comp
        graph, contig, L, HDN, R, tip = right_tip_structure
        path = graph.assemble_linear_path(revcomp(contig[0:K]))







        assert len(path) == len(contig)
    assert utils._equals_rc(path, contig)

    def test_end_to_beginning(self, right_tip_structure):
        # should have exact same behavior as right_of_branch_outwards
        graph, contig, L, HDN, R, tip = right_tip_structure
        path = graph.assemble_linear_path(contig[-K:])

        assert len(path) == len(contig)
        assert utils._equals_rc(path, contig)

def test_assemble_linear_path_10():
    # assemble up to branch point, and include introduced branch b/c
    # of stop bf
    contigfile = utils.get_test_data('simple-genome.fa')
    contig = list(screed.open(contigfile))[0].sequence
    print('contig len', len(contig))







    def test_branch_to_end(self, left_tip_structure):
        # assemble from branch point until end
        graph, contig, L, HDN, R, tip = left_tip_structure








        stop_bf = khmer.Nodegraph(K, 1e5, 4)
        stop_bf.count(tip)








        graph.consume(mutate_position(contig, HDN.pos + K))



    nodegraph = khmer.Nodegraph(K, 1e5, 4)
    lh = khmer._GraphLabels(nodegraph)



    lh.label_across_high_degree_nodes(contig, hdn, 1)

    path = lh.assemble_labeled_path(contig[:K])
        assert len(path) == 1, "there should only be one path"
    path = path[0]  # @CTB
    len_path = len(path)

    print('len path:', len_path)

        assert len(path) == len(contig)
    assert utils._equals_rc(path, contig)


    print(list(hdn))
    lh.label_across_high_degree_nodes(contig, hdn, 1)
    lh.label_across_high_degree_nodes(branch, hdn, 2)
    print(lh.get_tag_labels(list(hdn)[0]))

    paths = lh.assemble_labeled_path(contig[:K])




    print(list(hdn))
    lh.label_across_high_degree_nodes(contig, hdn, 1)
    print(lh.get_tag_labels(list(hdn)[0]))

    paths = lh.assemble_labeled_path(contig[:K])
    print([len(x) for x in paths])
    len_path = len(paths)

    print('len path:', len_path)

    found = False
    for path in paths:
        if utils._equals_rc(path, contig):
            found = True
            break
    assert found

    found = False
    for path in paths:
        if utils._equals_rc(path, branch):
            found = True
            break
    assert found


        assert any(utils._equals_rc(path, contig) for path in paths)
        assert any(utils._equals_rc(path, top_sequence) for path in paths)
        assert any(utils._equals_rc(path, bottom_sequence) for path in paths)

    # assemble entire contig + branch points b/c of labels; start from end


    print(list(hdn))
    lh.label_across_high_degree_nodes(contig, hdn, 1)
    lh.label_across_high_degree_nodes(branch, hdn, 2)
    lh.label_across_high_degree_nodes(branch2, hdn, 3)
    print(lh.get_tag_labels(list(hdn)[0]))

    paths = lh.assemble_labeled_path(contig[-K:])





    def test_assemble_snp_bubble_both(self, snp_bubble_structure):
        # assemble entire contig + both paths
        graph, wildtype, mutant, HDN_L, HDN_R = snp_bubble_structure
        lh = khmer._GraphLabels(graph)





    assert len(hdn) == 2



