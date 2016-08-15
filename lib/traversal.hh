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
#ifndef TRAVERSAL_HH
#define TRAVERSAL_HH

#include <queue>
#include <functional>

#include "khmer.hh"
#include "hashtable.hh"
#include "kmer_hash.hh"
#include "kmer_filters.hh"

namespace khmer
{

/*
#ifndef LEFT
#define LEFT 0
#endif
#ifndef RIGHT
#define RIGHT 1
#endif
*/

#ifndef LEFT
#define LEFT 0
#endif
#ifndef RIGHT
#define RIGHT 1
#endif

class Hashtable;
class LabelHash;
class Hashgraph;
class LabelHash;

/**
 * @brief Gather neighbors from a given node.
 *
 * The most basic traversal utility. Stores a list of KmerFilter functions, and given
 * a Kmer, finds all its neighbors that pass the filter function.s
 *
 * @tparam direction The direction in the graph to gather nodes from.
 */
template<bool direction>
class NodeGatherer: public KmerFactory
{
    friend class Hashgraph;

protected:

    KmerFilterList filters;
    HashIntoType bitmask;
    unsigned int rc_left_shift;
    const Hashtable * graph;

public:

    explicit NodeGatherer(const Hashtable * ht,
                          KmerFilterList filters);

    explicit NodeGatherer(const Hashtable * ht);

    explicit NodeGatherer(const Hashtable * ht, KmerFilter filter);

    /**
     * @brief Push a new filter on to the filter stack.
     */
    void push_filter(KmerFilter filter)
    {
        filters.push_back(filter);
    }

    /**
     * @brief Pop a filter off the stack.
     *
     * @return The filter.
     */
    KmerFilter pop_filter()
    {
        KmerFilter back = this->filters.back();
        this->filters.pop_back();
        return back;
    }

    unsigned int n_filters()
    {
        return filters.size();
    }

    /**
     * @brief Build the Kmer for the potential neighbor node of the given Kmer.
     *
     * When templated for RIGHT, will return the Kmer built from the length K-1 suffix of the
     * input Kmer with the new base appended; when LEFT, the length K-1 prefix of the input Kmer
     * with the new base prepended.
     *
     * @param node The starting node.
     * @param ch The new base to build from.
     *
     * @return The new Kmer.
     */
    Kmer get_neighbor(const Kmer& node, const char ch) const;

    /**
     * @brief Get all neighbors which are present in the graph and pass the filters.
     *
     * @param node The Kmer to start at.
     * @param node_q To collect the results.
     *
     * @return Number of neighbors found.
     */
    unsigned int neighbors(const Kmer& node,
                           KmerQueue &node_q) const;

    /**
     * @brief Get the degree of the given Kmer in the templated direction.
     *
     * @param node The Kmer to check.
     *
     * @return The degree.
     */
    unsigned int degree(const Kmer& node) const;
};


/**
 * @brief A stateful NodeGatherer. Stores its current position.
 *
 * @tparam direction The direction to gather nodes from.
 */
template <bool direction>
class NodeCursor: public NodeGatherer<direction>
{

public:

    // The current position.
    Kmer cursor;
    using NodeGatherer<direction>::push_filter;

    explicit NodeCursor(const Hashtable * ht,
                        Kmer start_kmer,
                        KmerFilterList filters);

    explicit NodeCursor(const Hashtable * ht,
                        Kmer start_kmer);

    explicit NodeCursor(const Hashtable * ht,
                        Kmer start_kmer,
                        KmerFilter filter);

    /**
     * @brief Get the neighbors for the current position.

     *
     * @param node_q To collection the results.
     *
     * @return Number of neighbors found.
     */
    unsigned int neighbors(KmerQueue& node_q) const
    {
        return NodeGatherer<direction>::neighbors(cursor, node_q);
    }

    /**
     * @return Degree of the current cursor position and direction.
     */
    unsigned int cursor_degree() const;

};


/**
 * @brief Wraps a LEFT and RIGHT NodeGatherer.
 */
class Traverser: public KmerFactory
{

protected:

    const Hashtable * graph;
    const Hashgraph * graph;
    NodeGatherer<LEFT> left_gatherer;
    NodeGatherer<RIGHT> right_gatherer;

public:
    explicit Traverser(const Hashgraph * ht);

    explicit Traverser(const Hashtable * ht,
                       KmerFilterList filters);

    explicit Traverser(const Hashtable * ht) : Traverser(ht, KmerFilterList()) {}

    explicit Traverser(const Hashtable * ht,
                       KmerFilter filter);

    void push_filter(KmerFilter filter);

    unsigned int traverse(const Kmer& node,
                          KmerQueue& node_q) const;

    unsigned int traverse_left(const Kmer& node,
                               KmerQueue& node_q) const;

    unsigned int traverse_right(const Kmer& node,
                                KmerQueue& node_q) const;

    unsigned int degree(const Kmer& node) const;
    unsigned int degree_left(const Kmer& node) const;
    unsigned int degree_right(const Kmer& node) const;

};


/**
 * @brief A NodeCursor specialized for assembling contigs.
 *
 * @tparam direction The direction to assemble.
 */
template <bool direction>
class AssemblerTraverser: public NodeCursor<direction>
{

    Kmer get_left(const Kmer& node, const char ch) const;
    Kmer get_right(const Kmer& node, const char ch) const;
public:
    using NodeCursor<direction>::NodeCursor;

    /**
     * @brief Get the next symbol.
     *
     * Finds the next symbol which passes the filters, so long as there is only
     * one branch. Does not return a new symbol if there are multiple potential neighbors.
     *
     * @return A member of alphabets::DNA_SIMPLE if a neighbor is found; '\0' otherwise.
     */
    virtual char next_symbol();

    unsigned int traverse_left(Kmer& node,
                               KmerQueue &node_q,
                               KmerFilter filter=0,
                               unsigned short max_neighbors=4) const;
    unsigned int traverse_right(Kmer& node,
                                KmerQueue &node_q,
                                KmerFilter filter=0,
                                unsigned short max_neighbors=4) const;
    unsigned int traverse(Kmer& node,
                          KmerQueue &node_q,
                          KmerFilter filter=0) const {
        unsigned int found;
        found = traverse_left(node, node_q, filter);
        found += traverse_right(node, node_q, filter);
        return found;
    };
    /**
     * @brief Utility function to join two overlapping contigs with proper directionality.
     *
     *  By default, assumes the two contigs overlap by length K. This can be reduced via the
     *  offset parameter.
     *
     * @param contig_a
     * @param contig_b
     * @param offset Number of bases to subtract from K when joining.
     *
     * @return The joined contig.
     */
    std::string join_contigs(std::string& contig_a,
                             std::string& contig_b,
                             WordLength offset = 0) const;
};

    unsigned int degree_left(const Kmer& node) const;
    unsigned int degree_right(const Kmer& node) const;
    unsigned int degree(const Kmer& node) const;

/**
 * @brief An AssemblerTraverser which does not traverse to Kmers it has already encountered.
 *
 * Simply adds a new filter to check if the Kmer has been seen, and adds the Kmer to the set
 * of seen Kmers after calling ::next_symbol.
 *
 * @tparam direction The direction to assemble.
 */
template<bool direction>
class NonLoopingAT: public AssemblerTraverser<direction>
{
protected:

    SeenSet * visited;

public:

    explicit NonLoopingAT(const Hashtable * ht,
                          Kmer start_kmer,
                          KmerFilterList filters,
                          SeenSet * visited);
    virtual char next_symbol();
};


template<bool direction>
class AssemblerTraverser: public Traverser
{

protected:

    KmerFilterList filters;

private:

    Kmer get_neighbor(Kmer& node, const char symbol);

public:

    Kmer cursor;

    explicit AssemblerTraverser(const Hashtable * ht,
                             Kmer start_kmer,
                             KmerFilterList filters);

    char next_symbol();
    bool set_cursor(Kmer& node);
    void push_filter(KmerFilter filter);
    KmerFilter pop_filter();
    unsigned int cursor_degree() const;

    std::string join_contigs(std::string& contig_a, std::string& contig_b) const;
};


template<bool direction>
class NonLoopingAT: public AssemblerTraverser<direction>
{
protected:

    const SeenSet * visited;

public:

    explicit NonLoopingAT(const Hashtable * ht,
                          Kmer start_kmer,
                          KmerFilterList filters,
                          const SeenSet * visited);
    char next_symbol();
};

}
#endif
