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
#ifndef FILTER_HH
#define FILTER_HH

#include <functional>

#include "khmer.hh"

#include "kmer_hash.hh"
#include "hashtable.hh"
#include "labelhash.hh"

#define DEBUG_FILTERS 0

namespace khmer
{

class Hashtable;
class LabelHash;

// A function which takes a Kmer and returns true if it
// is to be filtered / ignored
typedef std::function<bool (Kmer&)> KmerFilter;
typedef std::list<KmerFilter> KmerFilterList;

bool apply_kmer_filters(const Kmer& node, const KmerFilterList& filters);

inline bool apply_kmer_filters(Kmer& node, std::list<KmerFilter>& filters)
{
    if (!filters.size()) {
        return false;
    }
KmerFilter get_label_filter(const Label label, const LabelHash * lh);

    for(auto filter : filters) {
        if (filter(node)) {
            return true;
        }
    }
KmerFilter get_simple_label_intersect_filter(const LabelSet& src_labels,
        const LabelHash * lh,
        const unsigned int min_cov = 5);

    return false;
}
KmerFilter get_stop_bf_filter(const Hashtable * stop_bf);


inline KmerFilter get_label_filter(const Label label, const LabelHash * lh)
{
    KmerFilter filter = [=] (Kmer& node) {
        LabelSet ls;
        lh->get_tag_labels(node, ls);
        #if DEBUG_FILTERS
        if (ls.size() == 0) {
            std::cout << "no labels to jump to!" << std::endl;
        }
        #endif

        return !set_contains(ls, label);

    };
    return filter;
}


inline KmerFilter get_stop_bf_filter(const Hashtable * stop_bf)
{
    KmerFilter filter = [=] (Kmer& n) {
        return stop_bf->get_count(n);
    };
    return filter;
}


inline KmerFilter get_visited_filter(const SeenSet * visited)
{
    KmerFilter filter = [=] (Kmer& node) {
        #if DEBUG_FILTERS
        if(set_contains(*visited, node)) {
            std::cout << "loop!" << std::endl;
        }
        #endif
        return set_contains(*visited, node);
    };
    return filter;
}
KmerFilter get_visited_filter(const SeenSet * visited);


}
#endif
