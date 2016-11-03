/*
This file is part of khmer, https://github.com/dib-lab/khmer/, and is
Copyright (C) 2010-2015, Michigan State University.
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
#ifndef HASHBITS_HH
#define HASHBITS_HH

#include <stddef.h>
#include <string.h>
#include <string>
#include <vector>

#include "hashtable.hh"
#include "khmer.hh"
#include "kmer_hash.hh"

namespace khmer
{

class CountingHash;
class LabelHash;

class Hashbits : public Hashgraph {
public:
    Hashbits(WordLength ksize, std::vector<uint64_t>& tablesizes)
        : khmer::Hashtable(ksize),
          _tablesizes(tablesizes)
    {
        _occupied_bins = 0;
        _n_unique_kmers = 0;

        _allocate_counters();
    }

    ~Hashbits()
    {
        if (_counts) {
            for (size_t i = 0; i < _n_tables; i++) {
                delete[] _counts[i];
                _counts[i] = NULL;
            }
            delete[] _counts;
            _counts = NULL;

            _n_tables = 0;
        }

    }

    // Accessors for protected/private table info members
    std::vector<uint64_t> get_tablesizes() const
    {
        return _tablesizes;
    }

    const size_t n_tables() const
    {
        return _n_tables;
    }

    virtual void save(std::string);
    virtual void load(std::string);

    // count number of occupied bins
    virtual const uint64_t n_occupied() const
    {
        return _occupied_bins;
    }

    virtual const uint64_t n_unique_kmers() const
    {
        return _n_unique_kmers;
    }

    // Get and set the hashbits for the given kmer.
    inline
    virtual
    BoundedCounterType
    test_and_set_bits(const char * kmer)
    {
        HashIntoType hash = hash_dna(kmer);
        return test_and_set_bits(hash);
    }

    // Get and set the hashbits for the given kmer hash.
    // Generally, it is better to keep tests and mutations separate,
    // but, in the interests of efficiency and thread safety,
    // tests and mutations are being blended here against conventional
    // software engineering wisdom.
    inline
    virtual
    BoundedCounterType
    test_and_set_bits( HashIntoType khash )
    {
        bool is_new_kmer = false;

        for (size_t i = 0; i < _n_tables; i++) {
            uint64_t bin = khash % _tablesizes[i];
            uint64_t byte = bin / 8;
            unsigned char bit = (unsigned char)(1 << (bin % 8));

            unsigned char bits_orig = __sync_fetch_and_or( *(_counts + i) +
                                      byte, bit );
            if (!(bits_orig & bit)) {
                if (i == 0) {
                    __sync_add_and_fetch( &_occupied_bins, 1 );
                }
                is_new_kmer = true;
            }
        } // iteration over hashtables

        if (is_new_kmer) {
            __sync_add_and_fetch( &_n_unique_kmers, 1 );
            return 1; // kmer not seen before
        }

        return 0; // kmer already seen
    } // test_and_set_bits

    virtual void count(const char * kmer)
    {
        HashIntoType hash = hash_dna(kmer);
        count(hash);
    }

    virtual void count(HashIntoType khash)
    {
        test_and_set_bits(khash);
    }

    // get the count for the given k-mer.
    virtual const BoundedCounterType get_count(const char * kmer) const
    {
        HashIntoType hash = hash_dna(kmer);
        return get_count(hash);
    }

    // get the count for the given k-mer hash.
    virtual const BoundedCounterType get_count(HashIntoType khash) const
    {
        for (size_t i = 0; i < _n_tables; i++) {
            uint64_t bin = khash % _tablesizes[i];
            uint64_t byte = bin / 8;
            unsigned char bit = bin % 8;

            if (!(_counts[i][byte] & (1 << bit))) {
                return 0;
            }
        }
        return 1;
    }

    // Writing to the tables outside of defined methods has undefined behavior!
    // As such, this should only be used to return read-only interfaces
    Byte ** get_raw_tables()
    {
        return _counts;
    }
    explicit Hashbits(WordLength ksize, std::vector<uint64_t> sizes)
        : Hashgraph(ksize, new BitStorage(sizes)) { } ;

    void update_from(const Hashbits &other);

};
}

#include "labelhash.hh"
#endif // HASHBITS_HH

// vim: set sts=2 sw=2:
