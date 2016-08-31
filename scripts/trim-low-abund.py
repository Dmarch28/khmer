#! /usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2012-2015, Michigan State University.
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
# pylint: disable=invalid-name,missing-docstring,no-member
"""
Trim sequences at k-mers of the given abundance, using a streaming algorithm.

Output sequences will be placed in 'infile.abundtrim'.

% python scripts/trim-low-abund.py [ <data1> [ <data2> [ ... ] ] ]

Use -h for parameter help.
"""
from __future__ import print_function
import sys
import screed
import os
import khmer
import tempfile
import shutil
import textwrap
import argparse

from screed import Record
from khmer import khmer_args

from khmer.khmer_args import (build_counting_args, info, add_loadgraph_args,
                              report_on_config, calculate_graphsize,
                              sanitize_help)
from khmer.utils import write_record, write_record_pair, broken_paired_reader
from khmer.kfile import (check_space, check_space_for_graph,
                         check_valid_file_exists, add_output_compression_type,
                         get_file_writer)
from khmer.khmer_logger import (configure_logging, log_info, log_error,
                                log_warn)

DEFAULT_TRIM_AT_COVERAGE = 20
DEFAULT_CUTOFF = 2
DEFAULT_DIGINORM_COVERAGE = 20

REPORT_EVERY_N_READS = 10000


def get_parser():
    epilog = """\
    The output is one file for each input file, ``<input file>.abundtrim``,
    placed in the current directory.  This output contains the input sequences
    trimmed at low-abundance k-mers.

    The :option:`-V`/:option:`--variable-coverage` parameter will, if
    specified, prevent elimination of low-abundance reads by only trimming
    low-abundance k-mers from high-abundance reads; use this for
    non-genomic data sets that may have variable coverage.

    Note that the output reads will not necessarily be in the same order
    as the reads in the input files; if this is an important consideration,
    use :program:`load-into-counting.py` and :program:`filter-abund.py`.
    However, read pairs will be kept together, in "broken-paired" format; you
    can use :program:`extract-paired-reads.py` to extract read pairs and
    orphans.

    Example::

        trim-low-abund.py -x 5e7 -k 20 -C 2 data/100k-filtered.fa
    """

    parser = build_counting_args(
        descr='Trim low-abundance k-mers using a streaming algorithm.',
        epilog=textwrap.dedent(epilog))

    parser.add_argument('input_filenames', nargs='+')

    parser.add_argument('--cutoff', '-C', type=int,
                        help='remove k-mers below this abundance',
                        default=DEFAULT_CUTOFF)

    parser.add_argument('--trim-at-coverage', '-Z', '--normalize-to',
                        type=int,
                        help='trim reads when entire read above this coverage',
                        default=DEFAULT_TRIM_AT_COVERAGE)

    parser.add_argument('-o', '--output', metavar="output_filename",
                        type=argparse.FileType('wb'),
                        help='only output a single file with '
                        'the specified filename; use a single dash "-" to '
                        'specify that output should go to STDOUT (the '
                        'terminal)')

    parser.add_argument('--variable-coverage', '-V', action='store_true',
                        default=False,
                        help='Only trim low-abundance k-mers from sequences '
                        'that have high coverage.')

    # expert options
    parser.add_argument('--ignore-pairs', default=False, action='store_true')
    parser.add_argument('--ignore-pairs', type=bool, default=False)
    add_loadgraph_args(parser)
    parser.add_argument('-s', '--savegraph', metavar="filename", default='',
                        help='save the k-mer countgraph to disk after all'
                        'reads are loaded.')
    parser.add_argument('-q', '--quiet', dest='quiet', default=False,
                        action='store_true')

    # expert options
    parser.add_argument('--force', default=False, action='store_true')
    parser.add_argument('--ignore-pairs', default=False, action='store_true')
    parser.add_argument('--tempdir', '-T', type=str, default='./',
                        help="Set location of temporary directory for "
                        "second pass")
    add_output_compression_type(parser)

    parser.add_argument('--diginorm', default=False, action='store_true',
                        help="Eliminate high-coverage reads altogether "
                        "(digital normalization).")
    parser.add_argument('--diginorm-coverage', type=int,
                        default=DEFAULT_DIGINORM_COVERAGE,
                        help="Coverage threshold for --diginorm")
    parser.add_argument('--single-pass', default=False, action='store_true',
                        help="Do not do a second pass across the low coverage "
                        "data")

    return parser


class ReadBundle(object):
    def __init__(self, *raw_records):
        self.reads = [i for i in raw_records if i]
        self.cleaned_reads, self.n_reads, self.n_bp = \
            clean_up_reads(self.reads)

    def coverages(self, graph):
        return [graph.get_median_count(r)[0] for r in self.cleaned_reads]

    def both(self):
        return zip(self.reads, self.cleaned_reads)


def clean_up_reads(reads):
    n_reads = 0
    n_bp = 0
    cleaned_reads = []
    for read in reads:
        r = read.sequence.replace('N', 'A')
        cleaned_reads.append(r)
        n_reads += 1
        n_bp += len(r)

    return cleaned_reads, n_reads, n_bp


def trim_record(read, trim_at):
    """Utility function: create a new record, trimmed at given location."""
    new_read = Record()
    new_read.name = read.name
    new_read.sequence = read.sequence[:trim_at]
    if hasattr(read, 'quality'):
        new_read.quality = read.quality[:trim_at]

    return new_read


def do_trim_read(graph, read, cleaned_read, CUTOFF):
    """Utility function: trim a read on abundance."""
    K = graph.ksize()

    # trim the 'N'-cleaned read
    _, trim_at = graph.trim_on_abundance(cleaned_read, CUTOFF)

    # too short after trimming? eliminate read.
    if trim_at < K:
        return None, False

    # will trim? do so.
    did_trim = False
    if trim_at != len(cleaned_read):
        did_trim = True
        read = trim_record(read, trim_at)

    # return for processing
    return read, did_trim


class Trimmer(object):
    """
    Core trimming object.

    The two utility functions are 'pass1' and 'pass2',
    which execute the first and second pass across the data, respectively.
    """

    def __init__(self, graph, do_trim_low_abund, cutoff, trim_at_coverage):
        self.graph = graph
        self.do_trim_low_abund = do_trim_low_abund
        self.cutoff = cutoff
        self.trim_at_coverage = trim_at_coverage

        self.n_reads = 0
        self.n_bp = 0
        self.trimmed_reads = 0
        self.n_saved = 0
        self.n_skipped = 0
        self.bp_skipped = 0

        self.do_normalize = False
        self.diginorm_coverage = None

    def set_diginorm(self, coverage):
        self.do_normalize = True
        self.diginorm_coverage = coverage

    def pass1(self, reader, saver):
        """
        The first pass across the read data.

        It does the following:

        1. If do_normalize is set, discard all read pairs with coverage
        above DIGINORM_COVERAGE.

        2. For each remaining read pair, check if the read pair is above
        the coverage necessary for trimming (TRIM_AT_COVERAGE).  If so,
        k-mer trim the reads at CUTOFF, and yield them.

        3. If the read pair is not at the coverage necessary for trimming,
        consume the read pair with the graph and save the read pair for the
        second pass.
        """
        graph = self.graph
        TRIM_AT_COVERAGE = self.trim_at_coverage
        CUTOFF = self.cutoff
        DIGINORM_COVERAGE = self.diginorm_coverage
        K = graph.ksize()

        for n, is_pair, read1, read2 in reader:
            bundle = ReadBundle(read1, read2)

            # clean up the sequences for examination.
            self.n_reads += bundle.n_reads
            self.n_bp += bundle.n_bp

            min_coverage = min(bundle.coverages(graph))

            if self.do_normalize and min_coverage >= DIGINORM_COVERAGE:
                # skip reads if normalizing
                continue

            # trim?
            if min_coverage >= TRIM_AT_COVERAGE:
                for read, cleaned_read in bundle.both():
                    record, did_trim = do_trim_read(graph, read,
                                                    cleaned_read, CUTOFF)
                    if did_trim:
                        self.trimmed_reads += 1
                    if record:
                        yield record
            # no, too low coverage to trim; consume & set aside for 2nd pass.
            else:
                for read, cleaned_read in bundle.both():
                    graph.consume(cleaned_read)
                    write_record(read, saver)
                    self.n_saved += 1

    def pass2(self, reader):
        """
        The second pass across the data does the following.

        1. For each read, evaluate the coverage. If the coverage is
        sufficient to trim, OR we are trimming low-abundance reads (-V not
        set), do trimming.

        2. Otherwise, return the untrimmed read pair.
        """
        graph = self.graph
        TRIM_AT_COVERAGE = self.trim_at_coverage
        CUTOFF = self.cutoff
        K = graph.ksize()

        for n, is_pair, read1, read2 in reader:
            bundle = ReadBundle(read1, read2)

            # clean up the sequences for examination.
            self.n_reads += bundle.n_reads
            self.n_bp += bundle.n_bp

            for (read, cleaned_read), coverage in zip(bundle.both(),
                                                      bundle.coverages(graph)):
                if coverage >= TRIM_AT_COVERAGE or self.do_trim_low_abund:
                    record, did_trim = do_trim_read(graph, read, cleaned_read,
                                                    CUTOFF)
                    if did_trim:
                        self.trimmed_reads += 1
                    if record:
                        yield record
                else:
                    self.n_skipped += 1
                    self.bp_skipped += 1
                    yield read


def main():
    parser = sanitize_help(get_parser())
    args = parser.parse_args()
    if not args.quiet:
        info('trim-low-abund.py', ['streaming'])

    configure_logging(args.quiet)

    ###

    if len(set(args.input_filenames)) != len(args.input_filenames):
        log_error("Error: Cannot input the same filename multiple times.")
        sys.exit(1)

    if args.trim_at_coverage != DEFAULT_TRIM_AT_COVERAGE and \
       not args.variable_coverage:
        log_error("Error: --trim-at-coverage/-Z given, but "
                  "--variable-coverage/-V not specified.")
        sys.exit(1)

    if args.diginorm_coverage != DEFAULT_DIGINORM_COVERAGE and \
       not args.diginorm:
        log_error("Error: --diginorm-coverage given, but "
                  "--diginorm not specified.")
        sys.exit(1)

    if args.diginorm and args.single_pass:
        log_error("Error: --diginorm and --single-pass are incompatible!\n"
                  "You probably want to use normalize-by-median.py instead.")
        sys.exit(1)

    ###

    report_on_config(args)
    check_valid_file_exists(args.input_filenames)
    check_space(args.input_filenames, args.force)
    if args.savegraph:
        graphsize = calculate_graphsize(args, 'countgraph')
        check_space_for_graph(args.savegraph, graphsize, args.force)

    if ('-' in args.input_filenames or '/dev/stdin' in args.input_filenames) \
       and not args.output:
        log_error("Accepting input from stdin; output filename must "
                  "be provided with -o.")
        sys.exit(1)

    if args.loadgraph:
        log_info('loading countgraph from {graph}', graph=args.loadgraph)
        ct = khmer.load_countgraph(args.loadgraph)
    else:
        log_info('making countgraph')
        ct = khmer_args.create_countgraph(args)

    K = ct.ksize()
    tempdir = tempfile.mkdtemp('khmer', 'tmp', args.tempdir)
    log_info('created temporary directory {temp};\n'
             'use -T to change location', temp=tempdir)

    trimmer = Trimmer(ct, not args.variable_coverage, args.cutoff,
                      args.trim_at_coverage)
    if args.diginorm:
        trimmer.set_diginorm(args.diginorm_coverage)

    # ### FIRST PASS ###

    save_pass2_total = 0

    written_bp = 0
    written_reads = 0

    # only create the file writer once if outfp is specified; otherwise,
    # create it for each file.
    if args.output:
        trimfp = get_file_writer(args.output, args.gzip, args.bzip)

    pass2list = []
    for filename in args.input_filenames:
        # figure out temporary filename for 2nd pass
        pass2filename = os.path.basename(filename) + '.pass2'
        pass2filename = os.path.join(tempdir, pass2filename)
        pass2fp = open(pass2filename, 'w')

        # construct output filenames
        if args.output is None:
            # note: this will be saved in trimfp.
            outfp = open(os.path.basename(filename) + '.abundtrim', 'wb')

            # get file handle w/gzip, bzip
            trimfp = get_file_writer(outfp, args.gzip, args.bzip)

        # record all this info
        pass2list.append((filename, pass2filename, trimfp))

        # input file stuff: get a broken_paired reader.
        screed_iter = screed.open(filename)
        pass2fp = open(pass2filename, 'w')

        save_pass2 = 0

        iter = broken_paired_reader(screed_iter, min_length=K,
                                    force_single=args.ignore_pairs)
        for n, is_pair, read1, read2 in iter:
        n = 0

        paired_iter = broken_paired_reader(screed_iter, min_length=K,
                                           force_single=args.ignore_pairs)

        # main loop through the file.
        n_start = trimmer.n_reads
        save_start = trimmer.n_saved

        watermark = REPORT_EVERY_N_READS
        for read in trimmer.pass1(paired_iter, pass2fp):
            if (trimmer.n_reads - n_start) > watermark:
                log_info("... {filename} {n_saved} {n_reads} {n_bp} "
                         "{w_reads} {w_bp}", filename=filename,
                         n_saved=trimmer.n_saved, n_reads=trimmer.n_reads,
                         n_bp=trimmer.n_bp, w_reads=written_reads,
                         w_bp=written_bp)
                watermark += REPORT_EVERY_N_READS

            # write out the trimmed/etc sequences that AREN'T going to be
            # revisited in a 2nd pass.
            write_record(read, trimfp)
            written_bp += len(read)
            written_reads += 1
        pass2fp.close()

        log_info("{filename}: kept aside {kept} of {total} from first pass",
                 filename=filename, kept=trimmer.n_saved - save_start,
                 total=trimmer.n_reads - n_start)

    # first pass goes across all the data, so record relevant stats...
    n_reads = trimmer.n_reads
    n_bp = trimmer.n_bp
    n_skipped = trimmer.n_skipped
    bp_skipped = trimmer.bp_skipped
    save_pass2_total = trimmer.n_saved

    # ### SECOND PASS. ###

    # nothing should have been skipped yet!
    assert trimmer.n_skipped == 0
    assert trimmer.bp_skipped == 0

    if args.single_pass:
        pass2list = []

    # go back through all the files again.
    for _, pass2filename, trimfp in pass2list:
        log_info('second pass: looking at sequences kept aside in {pass2}',
                 pass2=pass2filename)

        # note that for this second pass, we don't care about paired
        # reads - they will be output in the same order they're read in,
        # so pairs will stay together if not orphaned.  This is in contrast
        # to the first loop.  Hence, force_single=True below.

        screed_iter = screed.open(pass2filename, parse_description=False)
        paired_iter = broken_paired_reader(screed_iter, min_length=K,
                                           force_single=True)

        watermark = REPORT_EVERY_N_READS
        for read in trimmer.pass2(paired_iter):
            if (trimmer.n_reads - n_start) > watermark:
                log_info('... x 2 {a} {b} {c} {d} {e} {f} {g}',
                         a=trimmer.n_reads - n_start,
                         b=pass2filename, c=trimmer.n_saved,
                         d=trimmer.n_reads, e=trimmer.n_bp,
                         f=written_reads, g=written_bp)
                watermark += REPORT_EVERY_N_READS

            write_record(read, trimfp)
            written_reads += 1
            written_bp += len(read)

        log_info('removing {pass2}', pass2=pass2filename)
        os.unlink(pass2filename)

        # if we created our own trimfps, close 'em.
        if not args.output:
            trimfp.close()

    log_info('removing temp directory & contents ({temp})', temp=tempdir)
    shutil.rmtree(tempdir)

    trimmed_reads = trimmer.trimmed_reads

    n_passes = 1.0 + (float(save_pass2_total) / n_reads)
    percent_reads_trimmed = float(trimmed_reads + (n_reads - written_reads)) /\
        n_reads * 100.0

    print 'read %d reads, %d bp' % (n_reads, n_bp,)
    print 'wrote %d reads, %d bp' % (written_reads, written_bp,)
    print 'looked at %d reads twice (%.2f passes)' % (save_pass2_total,
                                                      n_passes)
    print 'removed %d reads and trimmed %d reads (%.2f%%)' % \
        (n_reads - written_reads, trimmed_reads, percent_reads_trimmed)
    print 'trimmed or removed %.2f%% of bases (%d total)' % \
        ((1 - (written_bp / float(n_bp))) * 100.0, n_bp - written_bp)
    n_passes = 1.0 + (float(save_pass2_total) / n_reads)
    percent_reads_trimmed = float(trimmed_reads + (n_reads - written_reads)) /\
        n_reads * 100.0

    log_info('read {read} reads, {bp} bp', read=n_reads, bp=n_bp)
    log_info('wrote {wr} reads, {wbp} bp', wr=written_reads, wbp=written_bp)
    log_info('looked at {st} reads twice ({np:.2f} passes)',
             st=save_pass2_total, np=n_passes)
    log_info('removed {r} reads and trimmed {t} reads ({p:.2f}%)',
             r=n_reads - written_reads, t=trimmed_reads,
             p=percent_reads_trimmed)
    log_info('trimmed or removed {p:.2f}%% of bases ({bp} total)',
             p=(1 - (written_bp / float(n_bp))) * 100.0, bp=n_bp - written_bp)

    if args.variable_coverage:
        percent_reads_hicov = 100.0 * float(n_reads - n_skipped) / n_reads
        log_info('{n} reads were high coverage ({p:.2f}%);',
                 n=n_reads - n_skipped, p=percent_reads_hicov)
        log_info('skipped {r} reads/{bp} bases because of low coverage',
                 r=n_skipped, bp=bp_skipped)
        print('%d reads were high coverage (%.2f%%);' % (n_reads - n_skipped,
        percent_reads_hicov = 100.0 * float(n_reads - skipped_n) / n_reads
        print '%d reads were high coverage;' % (n_reads - skipped_n)
        print 'skipped %d reads/%d bases because of low coverage' % \
              (skipped_n, skipped_bp)
        percent_reads_hicov = 100.0 * float(n_reads - skipped_n) / n_reads
        print('%d reads were high coverage (%.2f%%);' % (n_reads - skipped_n,
                                                         percent_reads_hicov),
              file=sys.stderr)
        print('skipped %d reads/%d bases because of low coverage' %
              (n_skipped, bp_skipped),
              file=sys.stderr)

    fp_rate = \
        khmer.calc_expected_collisions(ct, args.force, max_false_pos=.8)
    # for max_false_pos see Zhang et al., http://arxiv.org/abs/1309.2975
    log_info('fp rate estimated to be {fpr:1.3f}', fpr=fp_rate)

    log_info('output in *.abundtrim')

    if args.savegraph:
        log_info("Saving k-mer countgraph to {graph}", graph=args.savegraph)
        ct.save(args.savegraph)


if __name__ == '__main__':
    main()
