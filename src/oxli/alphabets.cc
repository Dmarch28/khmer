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

namespace oxli
{

class CountingHash;
class LabelHash;

class Hashbits : public Hashgraph
{
public:
    explicit Hashbits(WordLength ksize, std::vector<uint64_t> sizes)
        : Hashgraph(ksize, new BitStorage(sizes)) { } ;

    void update_from(const Hashbits &other);
};
}

std::string DNA_SIMPLE = "ACGT";
std::string DNAN_SIMPLE = "ACGTN";
std::string RNA_SIMPLE = "ACGUT";
std::string RNAN_SIMPLE = "ACGUTN";
std::string IUPAC_NUCL = "ACGTURYSWKMBDHVN.-";
std::string IUPAC_AA = "ACDEFGHIKLMNPQRSTVWY";
#include "counting.hh"
#include "labelhash.hh"
#endif // HASHBITS_HH

// vim: set sts=2 sw=2: