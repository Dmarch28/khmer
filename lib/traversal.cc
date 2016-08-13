/*
This file is part of khmer, https://github.com/dib-lab/khmer/, and is
Copyright (C) 2015-2016, The Regents of the University of California.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the name of the Michigan State University nor the names
      of its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
LICENSE (END)

Contact: khmer-project@idyll.org
*/
#include "khmer.hh"
#include "hashtable.hh"
#include "traversal.hh"
#include "symbols.hh"
#include "alphabets.hh"
#include "kmer_hash.hh"

using namespace std;

namespace khmer
Traverser::Traverser(const Hashgraph * ht) :
    KmerFactory(ht->ksize()), graph(ht)
{

/******************************************
 * NodeGatherer
 ******************************************/

template <bool direction>
NodeGatherer<direction>::NodeGatherer(const Hashtable * ht,
                                      KmerFilterList filters) :
    KmerFactory(ht->ksize()), graph(ht), filters(filters)
{
    bitmask = 0;
    for (unsigned int i = 0; i < _ksize; i++) {
        bitmask = (bitmask << 2) | 3;
    }
    rc_left_shift = _ksize * 2 - 2;
}

Kmer Traverser::get_left(const Kmer& node, const char ch)
    const

template <bool direction>
NodeGatherer<direction>::NodeGatherer(const Hashgraph * ht) :
    NodeGatherer(ht, KmerFilterList())
{
}


template <bool direction>
NodeGatherer<direction>::NodeGatherer(const Hashgraph * ht,
                                      KmerFilter filter) :
    NodeGatherer(ht, KmerFilterList())
{
    filters.push_back(filter);
}


template<>
Kmer NodeGatherer<LEFT>::get_neighbor(const Kmer& node, const char ch)
const
{
    // optimized bit-foo to check for LEFT neighbors in both forward and
    // reverse-complemented directions
    HashIntoType kmer_f, kmer_r;
    kmer_f = ((node.kmer_f) >> 2 | twobit_repr(ch) << rc_left_shift);
    kmer_r = (((node.kmer_r) << 2) & bitmask) | (twobit_comp(ch));
    return build_kmer(kmer_f, kmer_r);
}


Kmer Traverser::get_right(const Kmer& node, const char ch)
    const
template<>
Kmer NodeGatherer<RIGHT>::get_neighbor(const Kmer& node,
                                       const char ch)
const
{
    // optimized bit-foo to check for LEFT neighbors in both forward and
    // reverse-complemented directions
    HashIntoType kmer_f, kmer_r;
    kmer_f = (((node.kmer_f) << 2) & bitmask) | (twobit_repr(ch));
    kmer_r = ((node.kmer_r) >> 2) | (twobit_comp(ch) << rc_left_shift);
    return build_kmer(kmer_f, kmer_r);
}


unsigned int Traverser::traverse_left(Kmer& node,
                                      KmerQueue & node_q,
                                      std::function<bool (Kmer&)> filter,
                                      unsigned short max_neighbors)
    const

template<bool direction>
unsigned int NodeGatherer<direction>::neighbors(const Kmer& node,
        KmerQueue & node_q)
const
{
    unsigned int found = 0;
    char * base = alphabets::DNA_SIMPLE;

    while(*base != '\0') {
        Kmer prev_node = get_left(node, *base);
        if (graph->get_count(prev_node) && (!filter || filter(prev_node))) {
            node_q.push(prev_node);
    for (auto base : alphabets::DNA_SIMPLE) {
        // Get the putative neighboring Kmer
        Kmer neighbor = get_neighbor(node, base);
        // Now check if it's in the graph and passes the filters
        if (graph->get_count(neighbor) && !(apply_kmer_filters(neighbor, filters))) {
            node_q.push(neighbor);
            ++found;
            if (found > max_neighbors) {
                return found;
            }
        }
        ++base;
    }

    return found;
}

unsigned int Traverser::traverse_right(Kmer& node,
                                       KmerQueue & node_q,
                                       std::function<bool (Kmer&)> filter,
                                       unsigned short max_neighbors)
    const

template<bool direction>
unsigned int NodeGatherer<direction>::degree(const Kmer& node)
const
{
    unsigned int found = 0;
    char * base = alphabets::DNA_SIMPLE;
    unsigned int degree = 0;

    while(*base != '\0') {
        Kmer next_node = get_right(node, *base);
        if (graph->get_count(next_node) && (!filter || filter(next_node))) {
            node_q.push(next_node);
            ++found;
            if (found > max_neighbors) {
                return found;
            }
    for (auto base : alphabets::DNA_SIMPLE) {
        if (graph->get_count(get_neighbor(node, base))) {
            ++degree;
        }
        ++base;
    }

    return degree;
}

/******************************************
 * NodeCursor
 ******************************************/

template<bool direction>
NodeCursor<direction>::NodeCursor(const Hashgraph * ht,
                                  Kmer start_kmer,
                                  KmerFilterList filters) :
    NodeGatherer<direction>(ht, filters)
{
    cursor = start_kmer;
}


template<bool direction>
NodeCursor<direction>::NodeCursor(const Hashgraph * ht,
                                  Kmer start_kmer) :
    NodeCursor<direction>(ht, start_kmer, KmerFilterList())
{
}


template<bool direction>
NodeCursor<direction>::NodeCursor(const Hashgraph * ht,
                                  Kmer start_kmer,
                                  KmerFilter filter) :
    NodeCursor<direction>(ht, start_kmer)
{
    push_filter(filter);
}


template<bool direction>
unsigned int NodeCursor<direction>::cursor_degree()
const
{
    return this->degree(this->cursor);
}



/******************************************
 * Traverser
 ******************************************/

Traverser::Traverser(const Hashgraph * ht,
                     KmerFilterList filters) :
    KmerFactory(ht->ksize()),
    graph(ht),
    left_gatherer(ht, filters),
    right_gatherer(ht, filters)
{
}

Traverser::Traverser(const Hashgraph * ht,
                     KmerFilter filter) :
    KmerFactory(ht->ksize()),
    graph(ht),
    left_gatherer(ht, filter),
    right_gatherer(ht, filter)
{
}


void Traverser::push_filter(KmerFilter filter)
{
    left_gatherer.push_filter(filter);
    right_gatherer.push_filter(filter);
}

unsigned int Traverser::degree_left(const Kmer& node)
    const

unsigned int Traverser::traverse(const Kmer& node,
                                 KmerQueue& node_q) const
{
    unsigned int degree = 0;
    char * base = alphabets::DNA_SIMPLE;
    return left_gatherer.neighbors(node, node_q) +
           right_gatherer.neighbors(node, node_q);
}

    while(*base != '\0') {
        Kmer prev_node = get_left(node, *base);
        if (graph->get_count(prev_node)) {
            ++degree;
        }
        ++base;
    }

unsigned int Traverser::traverse_left(const Kmer& node,
                                      KmerQueue& node_q) const
{
    return left_gatherer.neighbors(node, node_q);
}

unsigned int Traverser::degree_right(const Kmer& node)
    const

unsigned int Traverser::traverse_right(const Kmer& node,
                                       KmerQueue& node_q) const
{
    return right_gatherer.neighbors(node, node_q);
}


unsigned int Traverser::degree(const Kmer& node) const
{
    return left_gatherer.degree(node) + right_gatherer.degree(node);
}


unsigned int Traverser::degree_left(const Kmer& node) const
{
    return left_gatherer.degree(node);
}


unsigned int Traverser::degree_right(const Kmer& node) const
{
    return right_gatherer.degree(node);
}




/******************************************
 * AssemblerTraverser
 ******************************************/

template <>
std::string AssemblerTraverser<RIGHT>::join_contigs(std::string& contig_a,
        std::string& contig_b, WordLength offset)
const
{
    return contig_a + contig_b.substr(_ksize - offset);
}

template <>
std::string AssemblerTraverser<LEFT>::join_contigs(std::string& contig_a,
        std::string& contig_b, WordLength offset)
const
{
    return contig_b + contig_a.substr(_ksize - offset);
}

template<bool direction>
char AssemblerTraverser<direction>::next_symbol()
{
    unsigned int degree = 0;
    char * base = alphabets::DNA_SIMPLE;

    while(*base != '\0') {
        Kmer next_node = get_right(node, *base);
        if (graph->get_count(next_node)) {
            ++degree;
        }
        ++base;
    }
    short found = 0;
    char found_base = '\0';
    Kmer neighbor;
    Kmer cursor_next;

    for (auto base : alphabets::DNA_SIMPLE) {
        // Get the putative neighbor for this base at the cursor position
        neighbor = NodeCursor<direction>::get_neighbor(this->cursor, base);

        // Now check that the putative neighbor is in the graph and passes the filters
        if (this->graph->get_count(neighbor) &&
                !apply_kmer_filters(neighbor, this->filters)) {

            found++;
            // This naive traverser stops on high degree nodes
            if (found > 1) {
                return '\0';
            }
            found_base = base;
            cursor_next = neighbor;
        }
    }

    if (!found) {
        return '\0';
    } else {
        this->cursor = cursor_next;
        return found_base;
    }
}

/******************************************
 * NonLoopingAT
 ******************************************/

template<bool direction>
NonLoopingAT<direction>::NonLoopingAT(const Hashgraph * ht,
                                      Kmer start_kmer,
                                      KmerFilterList filters,
                                      SeenSet * visited) :
    AssemblerTraverser<direction>(ht, start_kmer, filters), visited(visited)
{
    AssemblerTraverser<direction>::push_filter(get_visited_filter(visited));
}

unsigned int Traverser::degree(const Kmer& node)
    const
template<bool direction>
char NonLoopingAT<direction>::next_symbol()
{
#if DEBUG_TRAVERSAL
    std::cout << "Insert cursor to visited filter" << std::endl;
#endif
    visited->insert(this->cursor);
    return AssemblerTraverser<direction>::next_symbol();
}

template class NodeGatherer<LEFT>;
template class NodeGatherer<RIGHT>;
template class NodeCursor<LEFT>;
template class NodeCursor<RIGHT>;
template class AssemblerTraverser<RIGHT>;
template class AssemblerTraverser<LEFT>;
template class NonLoopingAT<RIGHT>;
template class NonLoopingAT<LEFT>;


} // namespace khmer
