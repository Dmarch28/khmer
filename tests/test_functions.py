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
# pylint: disable=missing-docstring,invalid-name,no-member
from __future__ import print_function
from __future__ import absolute_import

import screed
import khmer
import os
import sys
import collections
import pytest
from khmer.utils import check_is_pair, broken_paired_reader
from . import khmer_tst_utils as utils
from khmer.utils import (check_is_pair, broken_paired_reader, check_is_left,
                         check_is_right)
from khmer.kfile import check_input_files, get_file_writer
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


def test_forward_hash():
    assert khmer.forward_hash('AAAA', 4) == 0
    assert khmer.forward_hash('TTTT', 4) == 0
    assert khmer.forward_hash('CCCC', 4) == 170
    assert khmer.forward_hash('GGGG', 4) == 170

    h = 13607885392109549066
    assert khmer.forward_hash('GGTTGACGGGGCTCAGGGGGCGGCTGACTCCG', 32) == h


def test_get_file_writer_fail():
    somefile = utils.get_temp_filename("potato")
    somefile = open(somefile, "w")
    stopped = True
    try:
        get_file_writer(somefile, True, True)
        stopped = False
    except ValueError as err:
        assert "Cannot specify both bzip and gzip" in str(err), str(err)

    assert stopped, "Expected exception"


def test_forward_hash_no_rc():
    h = khmer.forward_hash_no_rc('AAAA', 4)
    assert h == 0, h

    h = khmer.forward_hash_no_rc('TTTT', 4)
    assert h == 85, h

    h = khmer.forward_hash_no_rc('CCCC', 4)
    assert h == 170, h

    h = khmer.forward_hash_no_rc('GGGG', 4)
    assert h == 255, h


def test_reverse_hash():
    s = khmer.reverse_hash(0, 4)
    assert s == "AAAA"

    s = khmer.reverse_hash(85, 4)
    assert s == "TTTT"

    s = khmer.reverse_hash(170, 4)
    assert s == "CCCC"

    s = khmer.reverse_hash(255, 4)
    assert s == "GGGG"


def test_reverse_hash_longs():
    # test explicitly with long integers, only needed for python2
    # the builtin `long` exists in the global scope only
    global long # pylint: disable=global-variable-undefined
    if sys.version_info > (3,):
        long = int

    s = khmer.reverse_hash(long(0), 4)
    assert s == "AAAA"

    s = khmer.reverse_hash(long(85), 4)
    assert s == "TTTT"

    s = khmer.reverse_hash(long(170), 4)
    assert s == "CCCC"

    s = khmer.reverse_hash(long(255), 4)
    assert s == "GGGG"


def test_reverse_hash_raises():
    with pytest.raises(TypeError) as excinfo:
        khmer.reverse_hash('2345', 4)

    assert 'must be int' in str(excinfo.value)


def test_hash_murmur3():
    assert khmer.hash_murmur3('AAAA') == 526240128537019279
    assert khmer.hash_murmur3('TTTT') == 526240128537019279
    assert khmer.hash_murmur3('CCCC') == 14391997331386449225
    assert khmer.hash_murmur3('GGGG') == 14391997331386449225


def test_hash_no_rc_murmur3():
    h = khmer.hash_no_rc_murmur3('AAAA')
    assert h == 5231866503566620412, h

    h = khmer.hash_no_rc_murmur3('TTTT')
    assert h == 5753003579327329651, h

    h = khmer.hash_no_rc_murmur3('CCCC')
    assert h == 3789793362494378039, h

    h = khmer.hash_no_rc_murmur3('GGGG')
    assert h == 17519752047064575358, h


def test_get_primes():
    primes = khmer.get_n_primes_near_x(7, 20)

    assert primes == [19, 17, 13, 11, 7, 5, 3]


def test_get_primes_fal():
    try:
        khmer.get_n_primes_near_x(5, 5)
        assert 0, "previous statement should fail"
    except RuntimeError as err:
        assert "unable to find 5 prime numbers < 5" in str(err)


def test_extract_countgraph_info_badfile():
    try:
        khmer.extract_countgraph_info(
            utils.get_test_data('test-abund-read-2.fa'))
        assert 0, 'this should fail'
    except ValueError:
        pass


def test_extract_countgraph_info():
    fn = utils.get_temp_filename('test_extract_counting.ct')
    for size in [1e6, 2e6, 5e6, 1e7]:
        ht = khmer.Countgraph(25, size, 4)
        ht.save(fn)

        try:
            info = khmer.extract_countgraph_info(fn)
        except ValueError as err:
            assert 0, 'Should not throw a ValueErorr: ' + str(err)
        ksize, table_size, n_tables, _, _, _, _ = info
        print(ksize, table_size, n_tables)

        assert(ksize) == 25
        assert table_size == size
        assert n_tables == 4

        try:
            os.remove(fn)
        except OSError as err:
            assert 0, '...failed to remove ' + fn + str(err)


def test_extract_nodegraph_info_badfile():
    try:
        khmer.extract_nodegraph_info(
            utils.get_test_data('test-abund-read-2.fa'))
        assert 0, 'this should fail'
    except ValueError:
        pass


def test_extract_nodegraph_info():
    fn = utils.get_temp_filename('test_extract_nodegraph.pt')
    for size in [1e6, 2e6, 5e6, 1e7]:
        ht = khmer.Nodegraph(25, size, 4)
        ht.save(fn)

        info = khmer.extract_nodegraph_info(fn)
        ksize, table_size, n_tables, _, _, _ = info
        print(ksize, table_size, n_tables)

        assert(ksize) == 25
        assert table_size == size, table_size
        assert n_tables == 4

        try:
            os.remove(fn)
        except OSError as err:
            print('...failed to remove {fn}'.format(fn) + str(err),
                  file=sys.stderr)


def test_check_file_status_kfile():
    fn = utils.get_temp_filename('thisfiledoesnotexist')

    old_stderr = sys.stderr
    sys.stderr = capture = StringIO()

    try:
        check_input_files(fn, False)
    except SystemExit:
        assert "does not exist" in capture.getvalue(), capture.getvalue()
    finally:
        sys.stderr = old_stderr


def test_check_file_status_kfile_force():
    fn = utils.get_temp_filename('thisfiledoesnotexist')

    old_stderr = sys.stderr
    sys.stderr = capture = StringIO()

    try:
        check_input_files(fn, True)
    except OSError:
        assert False
    finally:
        sys.stderr = old_stderr

    assert "does not exist" in capture.getvalue(), capture.getvalue()


def test_check_is_pair_1():
    read1 = screed.Record(name='seq', quality='###', sequence='AAA')
    read2 = screed.Record(name='seq2', quality='###', sequence='AAA')

    assert not check_is_pair(read1, read2)


def test_check_is_pair_2():
    read1 = screed.Record(name='seq/1', quality='###', sequence='AAA')
    read2 = screed.Record(name='seq/2', quality='###', sequence='AAA')

    assert check_is_pair(read1, read2)


def test_check_is_pair_3_fq():
    read1 = screed.Record(name='seq 1::', quality='###', sequence='AAA')
    read2 = screed.Record(name='seq 2::', quality='###', sequence='AAA')

    assert check_is_pair(read1, read2)


def test_check_is_pair_3_broken_fq_1():
    read1 = screed.Record(name='seq', quality='###', sequence='AAA')
    read2 = screed.Record(name='seq 2::', quality='###', sequence='AAA')

    assert not check_is_pair(read1, read2)


def test_check_is_pair_3_broken_fq_2():
    read1 = screed.Record(name='seq 1::', quality='###', sequence='AAA')
    read2 = screed.Record(name='seq', quality='###', sequence='AAA')

    assert not check_is_pair(read1, read2)


def test_check_is_pair_3_fa():
    read1 = screed.Record(name='seq 1::', sequence='AAA')
    read2 = screed.Record(name='seq 2::', sequence='AAA')

    assert check_is_pair(read1, read2)


def test_check_is_pair_4():
    read1 = screed.Record(name='seq/1', quality='###', sequence='AAA')
    read2 = screed.Record(name='seq/2', sequence='AAA')

    try:
        check_is_pair(read1, read2)
        assert False                    # check_is_pair should fail here.
    except ValueError:
        pass


def test_check_is_pair_4b():
    read1 = screed.Record(name='seq/1', sequence='AAA')
    read2 = screed.Record(name='seq/2', quality='###', sequence='AAA')

    try:
        check_is_pair(read1, read2)
        assert False                    # check_is_pair should fail here.
    except ValueError:
        pass


def test_check_is_pair_5():
    read1 = screed.Record(name='seq/1', sequence='AAA')
    read2 = screed.Record(name='seq/2', sequence='AAA')

    assert check_is_pair(read1, read2)


def test_check_is_pair_6():
    read1 = screed.Record(name='seq1', sequence='AAA')
    read2 = screed.Record(name='seq2', sequence='AAA')

    assert not check_is_pair(read1, read2)


def test_check_is_pair_7():
    read1 = screed.Record(name='seq/2', sequence='AAA')
    read2 = screed.Record(name='seq/1', sequence='AAA')

    assert not check_is_pair(read1, read2)


class Test_BrokenPairedReader(object):
    stream = [FakeFastaRead(name='seq1/1', sequence='A'*5),
              FakeFastaRead(name='seq1/2', sequence='A'*4),
              FakeFastaRead(name='seq2/1', sequence='A'*5),
              FakeFastaRead(name='seq3/1', sequence='A'*3),
              FakeFastaRead(name='seq3/2', sequence='A'*5)]

    def gather(self, **kw):
        iter = broken_paired_reader(self.stream, **kw)

        x = []
        for n, is_pair, read1, read2 in iter:
            if is_pair:
                x.append((read1.name, read2.name))
            else:
                x.append((read1.name, None))

        return x

    def testDefault(self):
        x = self.gather(min_length=1)

        expected = [('seq1/1', 'seq1/2'),
                    ('seq2/1', None),
                    ('seq3/1', 'seq3/2')]
        assert x == expected, x

    def testMinLength(self):
        x = self.gather(min_length=3)

        expected = [('seq1/1', 'seq1/2'),
                    ('seq2/1', None),
                    ('seq3/1', 'seq3/2')]
        assert x == expected, x

    def testMinLength_2(self):
        x = self.gather(min_length=4)

        expected = [('seq1/1', 'seq1/2'),
                    ('seq2/1', None),
                    ('seq3/2', None)]
        assert x == expected, x

    def testForceSingle(self):
        x = self.gather(force_single=True)

        expected = [('seq1/1', None),
                    ('seq1/2', None),
                    ('seq2/1', None),
                    ('seq3/1', None),
                    ('seq3/2', None)]
        assert x == expected, x

    def testForceSingleAndMinLength(self):
        x = self.gather(min_length=5, force_single=True)

        expected = [('seq1/1', None),
                    ('seq2/1', None),
                    ('seq3/2', None)]
        assert x == expected, x


def test_check_is_right():
    assert not check_is_right('seq1/1')
    assert not check_is_right('seq1 1::N')
    assert check_is_right('seq1/2')
    assert check_is_right('seq1 2::N')

    assert not check_is_right('seq')
    assert not check_is_right('seq 2')


def test_check_is_left():
    assert check_is_left('seq1/1')
    assert check_is_left('seq1 1::N')
    assert not check_is_left('seq1/2')
    assert not check_is_left('seq1 2::N')

    assert not check_is_left('seq')
    assert not check_is_left('seq 1')

    assert check_is_left(
        '@HWI-ST412:261:d15khacxx:8:1101:3149:2157 1:N:0:ATCACG')


def gather(stream, **kw):
    itr = broken_paired_reader(stream, **kw)

    x = []
    m = 0
    num = 0
    for num, is_pair, read1, read2 in itr:
        if is_pair:
            x.append((read1.name, read2.name))
        else:
            x.append((read1.name, None))
        m += 1

    return x, num, m


class Test_BrokenPairedReader(object):
    stream = [screed.Record(name='seq1/1', sequence='A' * 5),
              screed.Record(name='seq1/2', sequence='A' * 4),
              screed.Record(name='seq2/1', sequence='A' * 5),
              screed.Record(name='seq3/1', sequence='A' * 3),
              screed.Record(name='seq3/2', sequence='A' * 5)]

    def testDefault(self):
        x, n, m = gather(self.stream, min_length=1)

        expected = [('seq1/1', 'seq1/2'),
                    ('seq2/1', None),
                    ('seq3/1', 'seq3/2')]
        assert x == expected, x
        assert m == 3
        assert n == 3, n

    def testMinLength(self):
        x, n, m = gather(self.stream, min_length=3)

        expected = [('seq1/1', 'seq1/2'),
                    ('seq2/1', None),
                    ('seq3/1', 'seq3/2')]
        assert x == expected, x
        assert m == 3
        assert n == 3, n

    def testMinLength_2(self):
        x, n, m = gather(self.stream, min_length=4)

        expected = [('seq1/1', 'seq1/2'),
                    ('seq2/1', None),
                    ('seq3/2', None)]
        assert x == expected, x
        assert m == 3
        assert n == 3, n

    def testForceSingle(self):
        x, n, m = gather(self.stream, force_single=True)

        expected = [('seq1/1', None),
                    ('seq1/2', None),
                    ('seq2/1', None),
                    ('seq3/1', None),
                    ('seq3/2', None)]
        assert x == expected, x
        assert m == 5
        assert n == 4, n

    def testForceSingleAndMinLength(self):
        x, n, m = gather(self.stream, min_length=5, force_single=True)

        expected = [('seq1/1', None),
                    ('seq2/1', None),
                    ('seq3/2', None)]
        assert x == expected, x
        assert m == 3, m
        assert n == 2, n


def test_BrokenPairedReader_OnPairs():
    stream = [FakeFastaRead(name='seq1/1', sequence='A' * 5),
              FakeFastaRead(name='seq1/2', sequence='A' * 4),
              FakeFastaRead(name='seq3/1', sequence='A' * 3),
              FakeFastaRead(name='seq3/2', sequence='A' * 5)]

    x, n, m = gather(stream, min_length=4, require_paired=True)

    expected = [('seq1/1', 'seq1/2')]
    assert x == expected, x
    assert m == 1
    assert n == 0, n


def test_BrokenPairedReader_OnPairs_2():
    stream = [FakeFastaRead(name='seq1/1', sequence='A' * 5),
              FakeFastaRead(name='seq1/2', sequence='A' * 4),
              FakeFastaRead(name='seq3/1', sequence='A' * 5),   # switched
              FakeFastaRead(name='seq3/2', sequence='A' * 3)]   # wrt previous

    x, n, m = gather(stream, min_length=4, require_paired=True)

    expected = [('seq1/1', 'seq1/2')]
    assert x == expected, x
    assert m == 1
    assert n == 0, n


def test_BrokenPairedReader_OnPairs_3():
    stream = [FakeFastaRead(name='seq1/1', sequence='A' * 5),
              FakeFastaRead(name='seq1/2', sequence='A' * 4),
              FakeFastaRead(name='seq3/1', sequence='A' * 3),   # both short
              FakeFastaRead(name='seq3/2', sequence='A' * 3)]

    x, n, m = gather(stream, min_length=4, require_paired=True)

    expected = [('seq1/1', 'seq1/2')]
    assert x == expected, x
    assert m == 1
    assert n == 0, n


def test_BrokenPairedReader_OnPairs_4():
    stream = [FakeFastaRead(name='seq1/1', sequence='A' * 3),  # too short
              FakeFastaRead(name='seq1/2', sequence='A' * 4),
              FakeFastaRead(name='seq3/1', sequence='A' * 4),
              FakeFastaRead(name='seq3/2', sequence='A' * 5)]

    x, n, m = gather(stream, min_length=4, require_paired=True)

    expected = [('seq3/1', 'seq3/2')]
    assert x == expected, x
    assert m == 1
    assert n == 0, n
