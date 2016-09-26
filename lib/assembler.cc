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

//#include "khmer.hh"
//#include "hashtable.hh"
//#include "traversal.hh"
#include "hashtable.hh"
#include "traversal.hh"

#include "assembler.hh"

#include <algorithm>
#include <iostream>

#define DEBUG 0
#define DEBUG_AT 0

using namespace khmer;
#include <algorithm>
#include <iostream>

using namespace std;

namespace khmer
{

template<bool direction>
AssemblerTraverser<direction>::AssemblerTraverser(const Hashtable * ht,
                                 Kmer start_kmer,
                                 KmerFilterList filters) :
    Traverser(ht), filters(filters)
                                 KmerFilterList filters,
                                 bool direction) :
    Traverser(ht), filters(filters), direction(direction)
                                 bool traverse_right) :
    Traverser(ht), filters(filters), traverse_right(traverse_right)
Assembler::Assembler(const Hashtable * ht) :
    Traverser(ht)
namespace khmer
{
    cursor = start_kmer;
}

template<>
Kmer AssemblerTraverser<LEFT>::get_neighbor(Kmer& node,
                                            const char symbol) {
    return get_left(node, symbol);
}

template<>
Kmer AssemblerTraverser<RIGHT>::get_neighbor(Kmer& node,
                                             const char symbol) {
    return get_right(node, symbol);
}

template<>
unsigned int AssemblerTraverser<LEFT>::cursor_degree()
    const
{
    return degree_left(cursor);
}

template<>
unsigned int AssemblerTraverser<RIGHT>::cursor_degree()
    const
{
    return degree_right(cursor);
}

template <>
std::string AssemblerTraverser<RIGHT>::join_contigs(std::string& contig_a,
                                                           std::string& contig_b)
    const
{
    return contig_a + contig_b.substr(_ksize);
}

template <>
std::string AssemblerTraverser<LEFT>::join_contigs(std::string& contig_a,
                                                         std::string& contig_b)
    const
{
    return contig_b + contig_a.substr(_ksize);
}

template<bool direction>
char AssemblerTraverser<direction>::next_symbol()
{
    char * symbol_ptr = alphabets::DNA_SIMPLE;
    char base;
    short found = 0;
    Kmer neighbor;
    Kmer cursor_next;

    while(*symbol_ptr != '\0') {
        neighbor = get_neighbor(cursor, *symbol_ptr);

        if (graph->get_count(neighbor) &&
            !apply_kmer_filters(neighbor, filters)) {

            found++;
            if (found > 1) {
                return '\0';
            }
            base = *symbol_ptr;
            cursor_next = neighbor;
        }
        symbol_ptr++;
    }

    if (!found) {
        return '\0';
    } else {
        cursor = cursor_next;
        return base;
    }
}
            if(traverse_right) { // NOTE: hoping this gets optimized out because const
                neighbor = get_right(cursor, *base);
            } else {
                neighbor = get_left(cursor, *base);
            }
            std::cout << "Try: " << (char)*base << " (" << neighbor << ")" << std::endl;
            if (graph->get_count(neighbor) &&
                !apply_kmer_filters(neighbor, filters)) {
                std::cout << "Found " << (char)*base << std::endl;
                found++;
                if (found > 1) {
                    break;
                }
            }
            base++;
        }
        std::cout << "Found: " << found << std::endl;
        if (found != 1) {
            break;
        } else {
            contig_kmers.push_back(neighbor);
            cursor = neighbor;
        }
    }
}
/********************************
 * Simple Linear Assembly
 ********************************/


template<bool direction>
bool AssemblerTraverser<direction>::set_cursor(Kmer& node)
LinearAssembler::LinearAssembler(const Hashgraph * ht) :
    graph(ht), _ksize(ht->ksize())
{

<template bool direction>
}

unsigned int AssemblerTraverser::get_path_length()
    const
{
    return _ksize + (contig_kmers.size() - 1);
}

bool AssemblerTraverser::set_cursor(Kmer& node)
{
    if(!apply_kmer_filters(node, filters)) {
        cursor = node;
        return true;
    }
    return false;
}

template<bool direction>
void AssemblerTraverser<direction>::push_filter(KmerFilter filter)
{
    filters.push_back(filter);
}

template<bool direction>
KmerFilter AssemblerTraverser<direction>::pop_filter()
{
    KmerFilter back = filters.back();
    filters.pop_back();
    return back;
}

template<bool direction>
NonLoopingAT<direction>::NonLoopingAT(const Hashtable * ht,
                                      Kmer start_kmer,
                                      KmerFilterList filters,
                                      const SeenSet * visited) :
    AssemblerTraverser<direction>(ht, start_kmer, filters), visited(visited)
{
    AssemblerTraverser<direction>::push_filter(get_visited_filter(visited));
}

template<bool direction>
char NonLoopingAT<direction>::next_symbol()
{
    visited->insert(AssemblerTraverser<direction>::cursor);
    #if DEBUG
    std::cout << "next_symbol; visited " << visited->size() << std::endl;
    #endif
    return AssemblerTraverser<direction>::next_symbol();
}

/********************************
 * Simple Linear Assembly
 ********************************/

LinearAssembler::LinearAssembler(const Hashtable * ht) :
    graph(ht), _ksize(ht->ksize())
{

}

// Starting from the given seed k-mer, assemble the maximal linear path in
// both directions.
//
// No guarantees on direction, of course - this may return the reverse
// complement of the input sequence.
std::string LinearAssembler::assemble(const Kmer seed_kmer,
                                      const Hashtable * stop_bf)
const
                                const Hashtable * stop_bf)
std::string Assembler::assemble_linear_path(const Kmer seed_kmer,
                                            const Hashtable * stop_bf)
std::string LinearAssembler::assemble(const Kmer seed_kmer,
                                      const Hashgraph * stop_bf)
    const
{
    std::string right_contig = assemble_right(seed_kmer, stop_bf);
    std::string left_contig = assemble_left(seed_kmer, stop_bf);
    std::list<KmerFilter> node_filters;
    if (stop_bf) {
        auto stop_bf_filter = [&] (Kmer& n) {
            return stop_bf->get_count(n);
        };
        node_filters.push_back(stop_bf_filter);
    }

#if DEBUG
    std::string right_contig;
    AssemblerTraverser<RIGHT> right_cursor(graph, seed_kmer, node_filters);
    assemble_right(right_contig, right_cursor);
    AssemblerTraverser<RIGHT> right_cursor(graph, start_kmer, node_filters);
    assemble_right(seed_kmer, right_contig, node_filters);
    std::string right = assemble_right(seed_kmer, node_filters);
    std::string left = assemble_left(seed_kmer, node_filters);
    AssemblerTraverser right_traverser(graph, seed_kmer, node_filters);
    AssemblerTraverser left_traverser(graph, seed_kmer, node_filters, 0);

    //std::string start_kmer = seed_kmer.get_string_rep(_ksize);

    //std::string right = _assemble_right(start_kmer, node_filters);
    //std::string left = _assemble_left(start_kmer, node_filters);
    std::string start_kmer = seed_kmer.get_string_rep(_ksize);
    std::string right = _assemble_right(start_kmer.c_str(), stop_bf);
    if (graph->get_count(seed_kmer) == 0) {
        // If the seed k-mer is not in the de Bruijn graph, stop trying to make
        // something happen. It's not going to happen!
        return "";
    }
    std::string right_contig = assemble_right(seed_kmer, stop_bf);
    std::string left_contig = assemble_left(seed_kmer, stop_bf);

    std::string left_contig;
    AssemblerTraverser<LEFT> left_cursor(graph, seed_kmer, node_filters);
    assemble_left(left_contig, left_cursor);

    #if DEBUG
    std::cout << "Left: " << left_contig << std::endl;
    std::cout << "Right: " << right_contig << std::endl;
#endif
    #endif
    std::string right = right_traverser.assemble();
    std::string left = left_traverser.assemble();
    std::string right = _assemble_right(start_kmer, node_filters);
    std::string left = _assemble_left(start_kmer, node_filters);
    std::string right = _assemble_right(start_kmer, stop_bf);
    std::string left = _assemble_left(start_kmer, stop_bf);
    start_kmer = _revcomp(start_kmer);
    std::string left = _assemble_right(start_kmer.c_str(), stop_bf);
#if DEBUG_ASSEMBLY
    std::cout << "Left: " << left_contig << std::endl;
    std::cout << "Right: " << right_contig << std::endl;
#endif

    right_contig = right_contig.substr(_ksize);
    return left_contig + right_contig;
}


std::string LinearAssembler::assemble_right(const Kmer seed_kmer,
        const Hashtable * stop_bf)
const
{
    std::list<KmerFilter> node_filters;
    if (stop_bf) {
        node_filters.push_back(get_stop_bf_filter(stop_bf));
    }

    AssemblerTraverser<RIGHT> cursor(graph, seed_kmer, node_filters);
    return _assemble_directed<RIGHT>(cursor);
}


std::string LinearAssembler::assemble_left(const Kmer seed_kmer,
        const Hashtable * stop_bf)
const
{
    std::list<KmerFilter> node_filters;
    if (stop_bf) {
        node_filters.push_back(get_stop_bf_filter(stop_bf));
    }

    AssemblerTraverser<LEFT> cursor(graph, seed_kmer, node_filters);
    return _assemble_directed<LEFT>(cursor);
}

template <>
std::string LinearAssembler::_assemble_directed<LEFT>(AssemblerTraverser<LEFT>&
        cursor)
const
{
    std::string contig = cursor.cursor.get_string_rep(_ksize);
    if (!cursor.cursor.is_forward()) {
        contig = _revcomp(contig);
    }

#if DEBUG
    std::cout << "## assemble_left\nStart Contig: " << contig << std::endl;
#endif

    reverse(contig.begin(), contig.end());
    char next_base;

    while ((next_base = cursor.next_symbol()) != '\0') {
        contig += next_base;
    }

    reverse(contig.begin(), contig.end());

    return contig;
    return cursor.get_cursor();
    return contig;
    std::string contig = _assemble_right(_revcomp(start_kmer), node_filters);
    return _revcomp(contig);
    left = left.substr(_ksize);
    return _revcomp(left) + right;
    right_contig = right_contig.substr(_ksize);
    return left_contig + right_contig;
}

template<>
std::string LinearAssembler::_assemble_directed<RIGHT>
(AssemblerTraverser<RIGHT>& cursor)
const
std::string LinearAssembler::_assemble_directed<RIGHT>(AssemblerTraverser<RIGHT>& cursor)

std::string LinearAssembler::_assemble_directed(AssemblerTraverser<RIGHT>& cursor)
Kmer LinearAssembler::assemble_right(std::string& contig,
                                     AssemblerTraverser<RIGHT>& cursor)
                                     AssemblerTraverser<LEFT>& cursor)
std::string LinearAssembler::assemble_right(const Kmer start_kmer,
                                       std::list<KmerFilter>& node_filters)
                                       const Hashtable * stop_bf)
std::string LinearAssembler::assemble_right(const Kmer seed_kmer,
        const Hashgraph * stop_bf)
    const
{
    std::string contig = cursor.cursor.get_string_rep(_ksize);
    char next_base;
    std::string kmer = start_kmer;
    std::string contig = kmer;
    std::list<KmerFilter> node_filters;
    if (stop_bf) {
        node_filters.push_back(get_stop_bf_filter(stop_bf));
    }

    while (1) {
        char * base = alphabets::DNA_SIMPLE;
        bool found = false;
        char found_base;
        bool found2 = false;
    AssemblerTraverser<RIGHT> cursor(graph, seed_kmer, node_filters);
    return _assemble_directed<RIGHT>(cursor);
}

        while(*base != 0) {
            std::string try_kmer = kmer.substr(1) + (char) *base;
            Kmer try_hashed = build_kmer(try_kmer);

            // a hit!
            if (graph->get_count(try_hashed) &&
                !apply_kmer_filters(try_hashed, node_filters)) {

#if DEBUG
    std::cout << "## assemble_right\nContig: " << contig << std::endl;
#endif

    while ((next_base = cursor.next_symbol()) != '\0') {
        contig += next_base;
    }

    return contig;
}
                if (found) {
                    found2 = true;
                    break;
std::string LinearAssembler::assemble_left(const Kmer seed_kmer,
        const Hashgraph * stop_bf)
const
{
    std::list<KmerFilter> node_filters;
    if (stop_bf) {
        node_filters.push_back(get_stop_bf_filter(stop_bf));
                }
                found_base = (char) *base;
                found = true;


/********************************
 * Labeled Assembly
 ********************************/

LabeledLinearAssembler::LabeledLinearAssembler(const LabelHash * lh) :
    graph(lh->graph), lh(lh), _ksize(lh->graph->ksize())
{
    linear_asm = new LinearAssembler(graph);
}


// Starting from the given seed k-mer, assemble all maximal linear paths in
// both directions, using labels to skip over tricky bits.
StringVector LabeledLinearAssembler::assemble(const Kmer seed_kmer,
        const Hashtable * stop_bf)
const
{
#if DEBUG
    std::cout << "Assemble Labeled: " << seed_kmer.repr(_ksize) << std::endl;
#endif

    KmerFilterList node_filters;
    if (stop_bf) {
        node_filters.push_back(get_stop_bf_filter(stop_bf));
    }

    SeenSet * visited = new SeenSet();

#if DEBUG
    std::cout << "Assemble Labeled RIGHT: " << seed_kmer.repr(_ksize) << std::endl;
#endif
    StringVector right_paths;
    NonLoopingAT<RIGHT> rcursor(graph, seed_kmer, node_filters, visited);
    _assemble_directed<RIGHT>(rcursor, right_paths);

#if DEBUG
    std::cout << "Assemble Labeled LEFT: " << seed_kmer.repr(_ksize) << std::endl;
#endif
    StringVector left_paths;
    NonLoopingAT<LEFT> lcursor(graph, seed_kmer, node_filters, visited);
    _assemble_directed<LEFT>(lcursor, left_paths);

    StringVector paths;
    for (unsigned int i = 0; i < left_paths.size(); i++) {
        for (unsigned int j = 0; j < right_paths.size(); j++) {
            std::string right = right_paths[j];
            right = right.substr(_ksize);
            std::string contig = left_paths[i] + right;
            paths.push_back(contig);
        }
    }

    visited->clear();
    return paths;
}

template <bool direction>
void LabeledLinearAssembler::_assemble_directed(NonLoopingAT<direction>&
        start_cursor,
        StringVector& paths)
const
{

    std::string root_contig = linear_asm->_assemble_directed<direction>
                              (start_cursor);
    Kmer end_kmer = start_cursor.cursor;

    if (start_cursor.cursor_degree() > 1) {               // hit a HDN
#if DEBUG
        std::cout << "Contig thus far: " << root_contig << std::endl;
        std::cout << "HDN: " << end_kmer.repr(_ksize) << "\n";
#endif // DEBUG

        LabelSet labels;
        lh->get_tag_labels(end_kmer, labels);

        /* For each label, we try to find spanning paths. We create
         * a new cursor starting at the end k-mer, with the existing node
         * filters; then we give that to the label spanning function.
         *
         * NOTE: This implies that there may be some non-deterministic
         * behavior: the ordering of the labels could vary, and decides
         * the recursion termination.
         */
        if(labels.size() == 0) {
            // if no labels are found there's nothing to be done, return
#if DEBUG
            std::cout << "no labels" << std::endl;
#endif
            paths.push_back(root_contig);
            return;
        } else {
#if DEBUG
            std::cout << "num labels: " << labels.size() << std::endl;
#endif
            for (Label label : labels) {
                /* Copy the current cursor at end_cursor for the spanning function.
                 * We'll now assemble, following the given label, as far as we can.
                 * We add an extra filter to the list: now, if we find no labels, we
                 * continue assembling; if we find labels and ours is included, we
                 * continue; and if we find labels and ours is not included, we stop.
                 * This cursor should also have the filters for visited k-mers and
                 * the stop bloom filter already.
                 */
#if DEBUG
                std::cout << "label: " << label << std::endl;
#endif
                NonLoopingAT<direction> span_cursor(start_cursor);
                span_cursor.push_filter(get_label_filter(label, lh));
                std::string spanning_contig = linear_asm->_assemble_directed<direction>
                                              (span_cursor);

                if(spanning_contig.length() == _ksize) {
#if DEBUG
                    std::cout << "zero length spanning contig" << spanning_contig << std::endl;
#endif
                    paths.push_back(root_contig);
                    continue;
                }

                // Remove the label filter
                span_cursor.pop_filter();

                // Recurse and gather paths
                StringVector continue_contigs;
                _assemble_directed<direction>(span_cursor, continue_contigs);

                if (continue_contigs.size() == 0) {
                    paths.push_back(span_cursor.join_contigs(root_contig,
                                    spanning_contig));
                } else {
                    for (auto continue_contig : continue_contigs) {
                        std::string full_contig = span_cursor.join_contigs(root_contig,
                                                  spanning_contig);
                        paths.push_back(span_cursor.join_contigs(full_contig, continue_contig));
                    }
                }
            } //end for
        }
    } else {
        paths.push_back(root_contig);
    }
    } else {
        paths.push_back(contig);
    }
}

//TODO make member function to create this filter
std::string LabeledLinearAssembler::_assemble_across_labels(AssemblerTraverser& start_cursor,
                                                            const Label label)
    const
{
    /* We'll now assemble, following the given label, as far as we can.
     * We add an extra filter to the list: now, if we find no labels, we
     * continue assembling; if we find labels and ours is included, we
     * continue; and if we find labels and ours is not included, we stop.
     * This cursor should also have the filters for visited k-mers and
     * the stop bloom filter already.
     */
    KmerFilter label_filter = [&] (Kmer& node) {
        LabelSet ls;
        lh->get_tag_labels(node, ls);
        if (ls.size() == 0) {
            return false;
        } else {
            return !set_contains(ls, label);
        }
    };
    start_cursor.add_filter(label_filter);

    return linear_asm._assemble_directed(start_cursor);
    AssemblerTraverser<LEFT> cursor(graph, seed_kmer, node_filters);
    return _assemble_directed<LEFT>(cursor);
}

template <>
std::string LinearAssembler::_assemble_directed<LEFT>(AssemblerTraverser<LEFT>&
        cursor)
const
{
    std::string contig = cursor.cursor.get_string_rep(_ksize);
    if (!cursor.cursor.is_forward()) {
        contig = _revcomp(contig);
    }

#if DEBUG_ASSEMBLY
    std::cout << "## assemble_linear_left[start] at " << contig << std::endl;
#endif

    reverse(contig.begin(), contig.end());
    char next_base;
    unsigned int found = 0;

    while ((next_base = cursor.next_symbol()) != '\0') {
        contig += next_base;
        found++;
            }
            base++;

    reverse(contig.begin(), contig.end());
#if DEBUG_ASSEMBLY
    std::cout << "## assemble_linear_left[end] found " << found << std::endl;
#endif

    return contig;
        }
        if (!found || found2) {
            break;
        } else {
            contig += found_base;
            kmer = kmer.substr(1) + found_base;
            found = true;

template<>
std::string LinearAssembler::_assemble_directed<RIGHT>
(AssemblerTraverser<RIGHT>& cursor)
const
{
    std::string contig = cursor.cursor.get_string_rep(_ksize);
    if (!cursor.cursor.is_forward()) {
        contig = _revcomp(contig);
        }
    char next_base;
    unsigned int found = 0;

#if DEBUG_ASSEMBLY
    std::cout << "## assemble_linear_right[start] at " << contig << std::endl;
#endif

    while ((next_base = cursor.next_symbol()) != '\0') {
        contig += next_base;
        found++;
    }

    return cursor.get_cursor();
}
#if DEBUG_ASSEMBLY
    std::cout << "## assemble_linear_right[end] found " << found << std::endl;
#endif
    return contig;
}

/*
std::string Assembler::_assemble_directed(const char * start_kmer,
                                       const Hashtable * stop_bf,
                                       const bool assemble_left)

/********************************
 * Labeled Assembly
 ********************************/

SimpleLabeledAssembler::SimpleLabeledAssembler(const LabelHash * lh) :
    graph(lh->graph), lh(lh), _ksize(lh->graph->ksize())
{
    linear_asm = new LinearAssembler(graph);
}


// Starting from the given seed k-mer, assemble all maximal linear paths in
// both directions, using labels to skip over tricky bits.
StringVector SimpleLabeledAssembler::assemble(const Kmer seed_kmer,
        const Hashgraph * stop_bf)
const
{
#if DEBUG_ASSEMBLY
    std::cout << "Assemble Labeled: " << seed_kmer.repr(_ksize) << std::endl;
#endif

    KmerFilterList node_filters;
    if (stop_bf) {
        node_filters.push_back(get_stop_bf_filter(stop_bf));
    }

    SeenSet * visited = new SeenSet();

#if DEBUG_ASSEMBLY
    std::cout << "Assemble Labeled RIGHT: " << seed_kmer.repr(_ksize) << std::endl;
#endif
    StringVector right_paths;
    NonLoopingAT<RIGHT> rcursor(graph, seed_kmer, node_filters, visited);
    _assemble_directed<RIGHT>(rcursor, right_paths);

#if DEBUG_ASSEMBLY
    std::cout << "Assemble Labeled LEFT: " << seed_kmer.repr(_ksize) << std::endl;
#endif
    StringVector left_paths;
    NonLoopingAT<LEFT> lcursor(graph, seed_kmer, node_filters, visited);
    _assemble_directed<LEFT>(lcursor, left_paths);

    StringVector paths;
    for (unsigned int i = 0; i < left_paths.size(); i++) {
        for (unsigned int j = 0; j < right_paths.size(); j++) {
            std::string right = right_paths[j];
            right = right.substr(_ksize);
            std::string contig = left_paths[i] + right;
            paths.push_back(contig);
        }
    }

    visited->clear();
    return paths;
}

template <bool direction>
void SimpleLabeledAssembler::_assemble_directed(NonLoopingAT<direction>&
        start_cursor,
        StringVector& paths)
    const
{
    Kmer kmer = this->build_kmer(start_kmer);
    std::cout << "starting on kmer " << kmer.get_string_rep(_ksize) << std::endl;
    std::string contig = start_kmer;
    Traverser traverser(this);
    KmerQueue neighbors;
    unsigned short found;
#if DEBUG_ASSEMBLY
    std::cout << "## assemble_labeled_directed_" << direction << " [start] at " <<
              start_cursor.cursor.repr(_ksize) << std::endl;
#endif

    auto keep_func = [&] (Kmer& node) {
        return !stop_bf || !stop_bf->get_count(node);
    };
    // prime the traversal with the first linear segment
    std::string root_contig = linear_asm->_assemble_directed<direction>
                              (start_cursor);
#if DEBUG_ASSEMBLY
    std::cout << "Primed: " << root_contig << std::endl;
    std::cout << "Cursor: " << start_cursor.cursor.repr(_ksize) << std::endl;
#endif
    StringVector segments;
    std::vector< NonLoopingAT<direction> > cursors;

    std::cout << "start loop" << std::endl;
    while (1) {
        if (assemble_left) {
            found = traverser.traverse_left(kmer, neighbors, keep_func, 1);
    segments.push_back(root_contig);
    cursors.push_back(start_cursor);

    while(segments.size() != 0) {

        std::string segment = segments.back();
        NonLoopingAT<direction> cursor = cursors.back();
#if DEBUG_ASSEMBLY
        std::cout << "Pop: " << segments.size() << " segments on stack." << std::endl;
        std::cout << "Segment: " << segment << std::endl;
        std::cout << "Cursor: " << cursor.cursor.repr(_ksize) << std::endl;
        std::cout << "n_filters: " << cursor.n_filters() << std::endl;
#endif
        segments.pop_back();
        cursors.pop_back();

        // check if the cursor has hit a HDN or reached a dead end
        if (cursor.cursor_degree() > 1) {

            LabelSet labels;
            lh->get_tag_labels(cursor.cursor, labels);

            if(labels.size() == 0) {
                // if no labels are found we can do nothing; gather this contig
#if DEBUG_ASSEMBLY
                std::cout << "no-label dead-end" << std::endl;
#endif
                paths.push_back(segment);
                continue;
        } else {
            found = traverser.traverse_right(kmer, neighbors, keep_func, 1);
                // if there are labels, try to hop the HDN.
                // first, get a label filter
                cursor.push_filter(get_simple_label_intersect_filter(labels, lh));
                KmerQueue branch_starts;
                // now get neighbors that pass the filter
                cursor.neighbors(branch_starts);
                // remove the filter
                cursor.pop_filter();

                // no neighbors found; done with this path
                if (branch_starts.empty()) {
#if DEBUG_ASSEMBLY
                    std::cout << "no-neighbors dead-end" << std::endl;
#endif
                    paths.push_back(segment);
                    continue;
                }

                // found some neighbors; extend them
                while(!branch_starts.empty()) {
                    // spin off a cursor for the new branch
                    NonLoopingAT<direction> branch_cursor(cursor);
                    branch_cursor.cursor = branch_starts.front();
                    branch_starts.pop();

#if DEBUG_ASSEMBLY
                    std::cout << "Branch cursor: " << branch_cursor.cursor.repr(
                                  _ksize) << std::endl;
#endif

                    // assemble linearly as far as possible
                    std::string branch = linear_asm->_assemble_directed<direction>(branch_cursor);
                    // create a new segment with the branch
                    std::string new_segment = branch_cursor.join_contigs(segment, branch, 1);
#if DEBUG_ASSEMBLY
                    std::cout << "Push segment: " << new_segment << std::endl;
                    std::cout << "Push cursor: " << branch_cursor.cursor.repr(_ksize) << std::endl;
#endif
                    segments.push_back(new_segment);
                    cursors.push_back(branch_cursor);
                }
        }
        //std::cout << "check traverser result" << std::endl;
        if (found == 0) {
            std::cout << "no neighbors, break" << std::endl;
            break;
        } else if (found > 1) {
            KmerQueue().swap(neighbors); // clear queue
            std::cout << "break" << std::endl;
            break;
        } else {
            //std::cout << "put base on contig" << std::endl;
            //contig += revtwobit_repr(kmer & 3);
            //std::cout << "get new kmer" << std::endl;
            kmer = neighbors.front();
            if (assemble_left) {
                contig += (kmer.get_string_rep(this->_ksize)[0]);
            // this segment is a dead-end; keep the contig
#if DEBUG_ASSEMBLY
            std::cout << "degree-1 dead-end" << std::endl;
#endif
            paths.push_back(segment);
            continue;
        }
    }
}



/*
template <bool direction>
void SimpleLabeledAssembler::_assemble_directed(NonLoopingAT<direction>&
        start_cursor,
        StringVector& paths)
const
{
#if DEBUG_ASSEMBLY
    std::cout << "## assemble_labeled_directed_" << direction << " [start] at " <<
        start_cursor.cursor.repr(_ksize) << std::endl;
#endif
    std::string root_contig = linear_asm->_assemble_directed<direction>
                              (start_cursor);

    if (start_cursor.cursor_degree() > 1) {               // hit a HDN
#if DEBUG_ASSEMBLY
        std::cout << "Root contig: " << root_contig << std::endl;
        std::cout << "HDN: " << start_cursor.cursor.repr(_ksize) << "\n";
        std::cout << start_cursor.cursor_degree() << std::endl;
#endif // DEBUG_ASSEMBLY

        LabelSet labels;
        lh->get_tag_labels(start_cursor.cursor, labels);


        if(labels.size() == 0) {
            // if no labels are found there's nothing to be done, return

            paths.push_back(root_contig);
            return;
            } else {
                contig += (kmer.get_string_rep(this->_ksize)[this->_ksize-1]);
#if DEBUG_ASSEMBLY
            std::cout << "Found " << labels.size() << " labels" << std::endl;
#endif



            start_cursor.push_filter(get_simple_label_intersect_filter(labels, lh));
            KmerQueue branch_starts;
            start_cursor.neighbors(branch_starts);
            start_cursor.pop_filter();
#if DEBUG_ASSEMBLY
                std::cout << branch_starts.size() << " neighbors found" << std::endl;
#endif
            if (branch_starts.size() == 0) {
                paths.push_back(root_contig);
                return;
            }

            StringVector branch_contigs;
            while(!branch_starts.empty()) { // TODO: change from queue
                NonLoopingAT<direction> branch_cursor(start_cursor);
                branch_cursor.cursor = branch_starts.front();
                branch_starts.pop();

                _assemble_directed<direction>(branch_cursor, branch_contigs);
            }

            //std::cout << "pop!" << std::endl;
            neighbors.pop();
            for (auto branch_contig : branch_contigs) {
                std::string full_contig = start_cursor.join_contigs(root_contig,
                                                                   branch_contig,
                                                                   1);
                paths.push_back(full_contig);
            }
        }
    } else {
        paths.push_back(root_contig);
    }
    return contig;
}

*/

}

}
