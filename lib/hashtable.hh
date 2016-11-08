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
#ifndef HASHTABLE_HH
#define HASHTABLE_HH


#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <vector>

#include "khmer.hh"
#include "khmer_exception.hh"
#include "kmer_hash.hh"
#include "read_parsers.hh"
#include "traversal.hh"
#include "subset.hh"

namespace khmer
{
namespace read_parsers
{
struct IParser;
}  // namespace read_parsers
}  // namespace khmer

#define MAX_KEEPER_SIZE int(1e6)

#define next_f(kmer_f, ch) ((((kmer_f) << 2) & bitmask) | (twobit_repr(ch)))
#define next_r(kmer_r, ch) (((kmer_r) >> 2) | (twobit_comp(ch) << rc_left_shift))

#define prev_f(kmer_f, ch) ((kmer_f) >> 2 | twobit_repr(ch) << rc_left_shift)
#define prev_r(kmer_r, ch) ((((kmer_r) << 2) & bitmask) | (twobit_comp(ch)))

#define set_contains(s, e) ((s).find(e) != (s).end())

#define CALLBACK_PERIOD 100000

#include "bitstorage.hh"
#include "bytestorage.hh"

namespace khmer
{
class Hashtable: public
    KmerFactory  		// Base class implementation of a Bloom ht.
{
protected:
    Storage * store;
    unsigned int    _max_count;
    unsigned int    _max_bigcount;

    //WordLength	    _ksize;
    HashIntoType    bitmask;
    unsigned int    _nbits_sub_1;

    explicit Hashtable( WordLength ksize, Storage * s)
        : KmerFactory( ksize ), store(s),
          _max_count( MAX_KCOUNT ),
          _max_bigcount( MAX_BIGCOUNT )
    {
        _init_bitstuff();
    }

    virtual ~Hashtable( )
    {
        delete store;
    }

    void _init_bitstuff()
    {
        bitmask = 0;
        for (unsigned int i = 0; i < _ksize; i++) {
            bitmask = (bitmask << 2) | 3;
        }
        _nbits_sub_1 = (_ksize*2 - 2);
    }

    explicit Hashtable(const Hashtable&);
    Hashtable& operator=(const Hashtable&);

public:
    // accessor to get 'k'
    const WordLength ksize() const
    {
        return _ksize;
    }

    // various hash functions.
    inline
    virtual
    HashIntoType
    hash_dna(const char * kmer) const {
        return _hash(kmer, _ksize);
    }

    inline
    virtual
    HashIntoType
    hash_dna_top_strand(const char * kmer) const {
        HashIntoType f = 0, r = 0;
        _hash(kmer, _ksize, f, r);
        return f;
    }

    inline
    virtual
    HashIntoType
    hash_dna_bottom_strand(const char * kmer) const {
        HashIntoType f = 0, r = 0;
        _hash(kmer, _ksize, f, r);
        return r;
    }

    inline
    virtual
    std::string
    unhash_dna(HashIntoType hashval) const {
        return _revhash(hashval, _ksize);
    }
    inline
    virtual
    HashIntoType
    hash_dna(const char * kmer) const {
        return _hash(kmer, _ksize);
    }

    inline
    HashIntoType
    hash_dna_top_strand(const char * kmer) const {
        HashIntoType f = 0, r = 0;
        _hash(kmer, _ksize, f, r);
        return f;
    }

    inline
    HashIntoType
    hash_dna_bottom_strand(const char * kmer) const {
        HashIntoType f = 0, r = 0;
        _hash(kmer, _ksize, f, r);
        return r;
    }

    inline
    virtual
    std::string
    unhash_dna(HashIntoType hashval) const {
        return _revhash(hashval, _ksize);
    }

    void count(const char * kmer) { store->add(hash_dna(kmer)); }
    void count(HashIntoType khash) { store->add(khash); }
    void add(const char * kmer) { store->add(hash_dna(kmer)); }
    void add(HashIntoType khash) { store->add(khash); }

    // get the count for the given k-mer.
    const BoundedCounterType get_count(const char * kmer) const {
        return store->get_count(hash_dna(kmer));
    }
    const BoundedCounterType get_count(HashIntoType khash) const {
        return store->get_count(khash);
    }

    void save(std::string filename) {
        store->save(filename, _ksize);
    }
    void load(std::string filename) {
        store->load(filename, _ksize);
        _init_bitstuff();
    }

    // count every k-mer in the string.
    unsigned int consume_string(const std::string &s);

    // checks each read for non-ACGT characters
    bool check_and_normalize_read(std::string &read) const;

    // check each read for non-ACGT characters, and then consume it.
    unsigned int check_and_process_read(std::string &read,
                                        bool &is_valid);

    // Count every k-mer in a FASTA or FASTQ file.
    // Note: Yes, the name 'consume_fasta' is a bit misleading,
    //	     but the FASTA format is effectively a subset of the FASTQ format
    //	     and the FASTA portion is what we care about in this case.
    void consume_fasta(
        std::string const   &filename,
        unsigned int	    &total_reads,
        unsigned long long  &n_consumed
    );

    // Count every k-mer from a stream of FASTA or FASTQ reads,
    // using the supplied parser.
    void consume_fasta(
        read_parsers:: IParser *	    parser,
        unsigned int	    &total_reads,
        unsigned long long  &n_consumed
    );

    void set_use_bigcount(bool b) { store->set_use_bigcount(b); }
    bool get_use_bigcount() { return store->get_use_bigcount(); }

    bool median_at_least(const std::string &s,
                         unsigned int cutoff);

    void get_median_count(const std::string &s,
                          BoundedCounterType &median,
                          float &average,
                          float &stddev);

    // number of unique k-mers
    const uint64_t n_unique_kmers() const { return store->n_unique_kmers(); }
    
    // count number of occupied bins
    const uint64_t n_occupied() const { return store->n_occupied(); }

    // table information
    std::vector<uint64_t> get_tablesizes() const {
        return store->get_tablesizes();
    }
    const size_t n_tables() const { return store->n_tables(); }

    // return all k-mer substrings, on the forward strand.
    void get_kmers(const std::string &s, std::vector<std::string> &kmers)
    const;

    // return hash values for all k-mer substrings
    void get_kmer_hashes(const std::string &s,
                         std::vector<HashIntoType> &kmers) const;

    // return hash values for all k-mer substrings in a SeenSet
    void get_kmer_hashes_as_hashset(const std::string &s,
                                    SeenSet& hashes) const;

    // return counts of all k-mers in this string.
    void get_kmer_counts(const std::string &s,
                         std::vector<BoundedCounterType> &counts) const;

    // get access to raw tables.
    Byte ** get_raw_tables() { return store->get_raw_tables(); }

    // find the minimum k-mer count in the given sequence
    BoundedCounterType get_min_count(const std::string &s);

    // find the maximum k-mer count in the given sequence
    BoundedCounterType get_max_count(const std::string &s);

    // calculate the abundance distribution of kmers in the given file.
    uint64_t * abundance_distribution(read_parsers::IParser * parser,
                                      Hashtable * tracking);
    uint64_t * abundance_distribution(std::string filename,
                                      Hashtable * tracking);

    // return the index of the first position in the sequence with k-mer
    // abundance below min_abund.
    unsigned long trim_on_abundance(std::string seq,
                                    BoundedCounterType min_abund) const;

    // return the index of the first position in the sequence with k-mer
    // abundance above max_abund.
    unsigned long trim_below_abundance(std::string seq,
                                       BoundedCounterType max_abund) const;

    // detect likely positions of errors
    std::vector<unsigned int> find_spectral_error_positions(std::string seq,
            BoundedCounterType min_abund) const;
};

//
// Hashgraph: Extension of Hashtable to support graph operations.
//

class Hashgraph: public Hashtable {

    friend class SubsetPartition;
    friend class LabelHash;
    friend class Traverser;

protected:
    unsigned int _tag_density;

    explicit Hashgraph(WordLength ksize, Storage * s)
        : Hashtable(ksize, s)
    {
        _tag_density = DEFAULT_TAG_DENSITY;
        if (!(_tag_density % 2 == 0)) {
            throw khmer_exception();
        }
        partition = new SubsetPartition(this);
        _all_tags_spin_lock = 0;
    }

    // clean up the partition structure.
    virtual ~Hashgraph( )
    {
        delete partition;
    }

    // empty the partition structure
    void _clear_all_partitions()
    {
        if (partition != NULL) {
            partition->_clear_all_partitions();
        }
    }

    uint32_t _all_tags_spin_lock;
public:
    // default master partitioning
    SubsetPartition * partition;

    // tags for sparse graph implementation
    SeenSet all_tags;

    // tags at which to stop traversal
    SeenSet stop_tags;

    // tags used in repartitioning
    SeenSet repart_small_tags;

    // set the minimum density of tagging.
    void _set_tag_density(unsigned int d)
    {
        // must be odd; can't be set if tags exist.
        if (!(d % 2 == 0) || !all_tags.empty()) {
            throw khmer_exception();
        }
        _tag_density = d;
    }

    unsigned int _get_tag_density() const { return _tag_density; }

    void add_tag(HashIntoType tag) { all_tags.insert(tag); }
    void add_stop_tag(HashIntoType tag) { stop_tags.insert(tag); }

    size_t n_tags() const { return all_tags.size(); }

    void divide_tags_into_subsets(unsigned int subset_size, SeenSet& divvy);

    void add_kmer_to_tags(HashIntoType kmer) { all_tags.insert(kmer); }

    void clear_tags() { all_tags.clear(); }

    // Consume reads & build sparse graph.
    void consume_fasta_and_tag(
        std::string const &filename,
        unsigned int &total_reads,
        unsigned long long &n_consumed
    );

    // Count every k-mer from a stream of FASTA or FASTQ reads,
    // using the supplied parser.
    // Tag certain ones on the connectivity graph.
    void consume_fasta_and_tag(
        read_parsers:: IParser * parser,
        unsigned int &total_reads,
        unsigned long long &n_consumed
    );

    // consume a string & add sparse graph nodes.
    void consume_sequence_and_tag(const std::string& seq,
                                  unsigned long long& n_consumed,
                                  SeenSet * new_tags = 0);


    // consume an already-partitioned file & load in the partition IDs
    void consume_partitioned_fasta(const std::string &filename,
                                   unsigned int &total_reads,
                                   unsigned long long &n_consumed);

    // trim the given sequence on stoptags
    size_t trim_on_stoptags(std::string sequence) const;

    // @@
    unsigned int traverse_from_kmer(Kmer start,
                                    unsigned int radius,
                                    KmerSet &keeper,
                                    unsigned int max_count = MAX_KEEPER_SIZE)
    const;

    // print, save, and load the set of tags.
    void print_tagset(std::string);
    void save_tagset(std::string);
    void load_tagset(std::string, bool clear_tags=true);

    // print, save and load the set of stop tags.
    void print_stop_tags(std::string);
    void save_stop_tags(std::string);
    void load_stop_tags(std::string filename, bool clear_tags=true);

    // @@
    void extract_unique_paths(std::string seq,
                              unsigned int min_length,
                              float min_unique_f,
                              std::vector<std::string> &results);

    // @@
    void calc_connected_graph_size(Kmer node,
                                   unsigned long long& count,
                                   KmerSet& keeper,
                                   const unsigned long long threshold=0,
                                   bool break_on_circum=false) const;

    // Calculate the graph degree of the given k-mer.
    unsigned int kmer_degree(HashIntoType kmer_f, HashIntoType kmer_r);
    unsigned int kmer_degree(const char * kmer_s);

    // Find all nodes with a degree > 2.
    void find_high_degree_nodes(const char * sequence,
                                SeenSet& high_degree_nodes) const;

    // Find the maximal linear path (nodes degree <= 2)
    unsigned int traverse_linear_path(const Kmer start_kmer,
                                      SeenSet &adjacencies,
                                      SeenSet &nodes, Hashtable& bf,
                                      SeenSet &high_degree_nodes) const;
    
    //
    // for debugging/testing purposes only!
    //
    
    // check partition map validity.
    void _validate_pmap()
    {
        if (partition) {
            partition->_validate_pmap();
        }
    }

};

// Hashtable-derived class with ByteStorage.
class Counttable : public khmer::Hashtable
{
public:
    explicit Counttable(WordLength ksize, std::vector<uint64_t> sizes)
        : Hashtable(ksize, new ByteStorage(sizes)) { } ;
};

// Hashgraph-derived class with ByteStorage.
class CountingHash : public khmer::Hashgraph
{
public:
    explicit CountingHash(WordLength ksize, std::vector<uint64_t> sizes)
        : Hashgraph(ksize, new ByteStorage(sizes)) { } ;
};

// Hashtable-derived class with BitStorage.
class Nodetable : public Hashtable {
public:
    explicit Nodetable(WordLength ksize, std::vector<uint64_t> sizes)
        : Hashtable(ksize, new BitStorage(sizes)) { } ;
};

// Hashgraph-derived class with BitStorage.
class Hashbits : public Hashgraph {
public:
    explicit Hashbits(WordLength ksize, std::vector<uint64_t> sizes)
        : Hashgraph(ksize, new BitStorage(sizes)) { } ;

    void update_from(const Hashbits &other);
};
}

#define ACQUIRE_ALL_TAGS_SPIN_LOCK \
  while (!__sync_bool_compare_and_swap( &_all_tags_spin_lock, 0, 1 ));

#define RELEASE_ALL_TAGS_SPIN_LOCK \
  __sync_bool_compare_and_swap( &_all_tags_spin_lock, 1, 0 );

#endif // HASHTABLE_HH
