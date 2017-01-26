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
#ifndef _CPY_HASHGRAPH_HH
#define _CPY_HASHGRAPH_HH

#include <Python.h>
#include "_cpy_utils.hh"
#include "_cpy_hashtable.hh"
#include "hashgraph.hh"

namespace khmer {


typedef struct {
    khmer_KHashtable_Object khashtable;
    Hashgraph * hashgraph;
} khmer_KHashgraph_Object;


extern PyTypeObject khmer_KHashgraph_Type
CPYCHECKER_TYPE_OBJECT_FOR_TYPEDEF("khmer_KHashgraph_Object");


extern PyMethodDef khmer_hashgraph_methods[];


//
// Method definitions
//


PyObject *
hashgraph_find_high_degree_nodes(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_neighbors(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_traverse_linear_path(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_assemble_linear_path(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_n_tags(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_print_stop_tags(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_print_tagset(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_load_stop_tags(khmer_KHashgraph_Object * me, PyObject * args);



PyObject *
hashgraph_save_stop_tags(khmer_KHashgraph_Object * me, PyObject * args);



PyObject *
hashgraph_repartition_largest_partition(khmer_KHashgraph_Object * me,
                                        PyObject * args);


PyObject *
hashgraph_calc_connected_graph_size(khmer_KHashgraph_Object * me,
                                    PyObject * args);


PyObject *
hashgraph_kmer_degree(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_trim_on_stoptags(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_do_subset_partition(khmer_KHashgraph_Object * me, PyObject * args);

    if (!PyArg_ParseTuple(args, "|OOOO", &start_kmer_obj, &end_kmer_obj,
                          &break_on_stop_tags_o,
                          &stop_big_traversals_o)) {
        return NULL;
    }
    if (!ht_convert_PyObject_to_HashIntoType(start_kmer_obj, start_kmer,
            hashgraph)) {
        return NULL;
    }
    if (!ht_convert_PyObject_to_HashIntoType(end_kmer_obj, end_kmer,
            hashgraph)) {
        return NULL;
    }


PyObject *
hashgraph_merge_subset(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_merge_from_disk(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_consume_fasta_and_tag_with_reads_parser(khmer_KHashgraph_Object * me,
        PyObject * args);
hashgraph_consume_seqfile_and_tag_with_reads_parser(khmer_KHashgraph_Object * me,
        PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    python::khmer_ReadParser_Object * rparser_obj = NULL;

    if (!PyArg_ParseTuple( args, "O!", &python::khmer_ReadParser_Type,
                           &rparser_obj)) {
        return NULL;
    }

    FastxParserPtr& rparser = rparser_obj->parser;

    // call the C++ function, and trap signals => Python
    const char         *value_exception = NULL;
    const char         *file_exception  = NULL;
    unsigned long long  n_consumed      = 0;
    unsigned int        total_reads     = 0;
    std::string exc_string;

    Py_BEGIN_ALLOW_THREADS
    try {
        hashgraph->consume_seqfile_and_tag<FastxReader>(rparser, total_reads, n_consumed);
    } catch (khmer_file_exception &exc) {
        exc_string = exc.what();
        file_exception = exc_string.c_str();
    } catch (khmer_value_exception &exc) {
        exc_string = exc.what();
        value_exception = exc_string.c_str();
    }
    Py_END_ALLOW_THREADS

PyObject *
hashgraph_consume_partitioned_fasta(khmer_KHashgraph_Object * me,
                                    PyObject * args);


    try {
        hashgraph->consume_partitioned_fasta<FastxReader>(filename, total_reads, n_consumed);
    } catch (khmer_file_exception &exc) {
        PyErr_SetString(PyExc_OSError, exc.what());
        return NULL;
    } catch (khmer_value_exception &exc) {
        PyErr_SetString(PyExc_ValueError, exc.what());
        return NULL;
    }

    return Py_BuildValue("IK", total_reads, n_consumed);
}

static
PyObject *
hashgraph_find_all_tags(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_assign_partition_id(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_add_tag(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_add_stop_tag(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_get_stop_tags(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_get_tagset(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_output_partitions(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_save_partitionmap(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_load_partitionmap(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph__validate_partitionmap(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_count_partitions(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_subset_count_partitions(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_subset_partition_size_distribution(khmer_KHashgraph_Object * me,
        PyObject * args);


PyObject *
hashgraph_load_tagset(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_save_tagset(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_save_subset_partitionmap(khmer_KHashgraph_Object * me,
                                   PyObject * args);


PyObject *
hashgraph_load_subset_partitionmap(khmer_KHashgraph_Object * me,
                                   PyObject * args);


PyObject *
hashgraph__set_tag_density(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph__get_tag_density(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph__validate_subset_partitionmap(khmer_KHashgraph_Object * me,
                                        PyObject * args);


PyObject *
hashgraph_set_partition_id(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_join_partitions(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_get_partition_id(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_divide_tags_into_subsets(khmer_KHashgraph_Object * me,
                                   PyObject * args);


PyObject *
hashgraph_count_kmers_within_radius(khmer_KHashgraph_Object * me,
                                    PyObject * args);


PyObject *
hashgraph_extract_unique_paths(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_consume_and_tag(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_get_tags_and_positions(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_find_all_tags_list(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_consume_fasta_and_tag(khmer_KHashgraph_Object * me, PyObject * args);
hashgraph_consume_seqfile_and_tag(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    const char * filename;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    // call the C++ function, and trap signals => Python

    unsigned long long n_consumed;
    unsigned int total_reads;

    try {
        hashgraph->consume_seqfile_and_tag<FastxReader>(filename, total_reads, n_consumed);
    } catch (khmer_file_exception &exc) {
        PyErr_SetString(PyExc_OSError, exc.what());
        return NULL;
    } catch (khmer_value_exception &exc) {
        PyErr_SetString(PyExc_ValueError, exc.what());
        return NULL;
    }

    return Py_BuildValue("IK", total_reads, n_consumed);
}

static PyMethodDef khmer_hashgraph_methods[] = {
    //
    // graph/traversal functionality
    //

    {
        "neighbors",
        (PyCFunction)hashgraph_neighbors, METH_VARARGS,
        "Get a list of neighbor nodes for this k-mer.",
    },
    {
        "calc_connected_graph_size",
        (PyCFunction)hashgraph_calc_connected_graph_size, METH_VARARGS, ""
    },
    {
        "kmer_degree",
        (PyCFunction)hashgraph_kmer_degree, METH_VARARGS,
        "Calculate the number of immediate neighbors this k-mer has in "
        "the graph."
    },
    {
        "count_kmers_within_radius",
        (PyCFunction)hashgraph_count_kmers_within_radius, METH_VARARGS,
        "Calculate the number of neighbors with given radius in the graph."
    },
    {
        "find_high_degree_nodes",
        (PyCFunction)hashgraph_find_high_degree_nodes, METH_VARARGS,
        "Examine the given sequence for degree > 2 nodes and add to  "
        "list; used in graph contraction.",
    },
    {
        "traverse_linear_path",
        (PyCFunction)hashgraph_traverse_linear_path, METH_VARARGS,
        "Traverse the path through the graph starting with the given "
        "k-mer and avoiding high-degree nodes, finding (and returning) "
        "traversed k-mers and any encountered high-degree nodes.",
    },
    {
        "assemble_linear_path",
        (PyCFunction)hashgraph_assemble_linear_path, METH_VARARGS,
        "Assemble a purely linear path starting with the given "
        "k-mer, returning traversed k-mers and any encountered high-degree "
        "nodes.",
    },

    //
    // tagging / sparse graph functionality
    //

    {
        "consume_and_tag",
        (PyCFunction)hashgraph_consume_and_tag, METH_VARARGS,
        "Consume a sequence and tag it."
    },
    {
        "get_tags_and_positions",
        (PyCFunction)hashgraph_get_tags_and_positions, METH_VARARGS,
        "Retrieve tags and their positions in a sequence."
    },
    {
        "find_all_tags_list",
        (PyCFunction)hashgraph_find_all_tags_list, METH_VARARGS,
        "Find all tags within range of the given k-mer, return as list"
    },
    {
        "consume_seqfile_and_tag",
        (PyCFunction)hashgraph_consume_seqfile_and_tag, METH_VARARGS,
        "Consume all sequences in a FASTA/FASTQ file and tag the resulting "
        "graph."
    },
    {
        "extract_unique_paths",
        (PyCFunction)hashgraph_extract_unique_paths, METH_VARARGS,
        "@CTB remove."
    },
    {
        "print_tagset",
        (PyCFunction)hashgraph_print_tagset, METH_VARARGS,
        "Print out all of the tags."
    },
    {
        "add_tag",
        (PyCFunction)hashgraph_add_tag, METH_VARARGS,
        "Add a k-mer to the tagset."
    },
    {
        "get_tagset",
        (PyCFunction)hashgraph_get_tagset, METH_VARARGS,
        "Get all tagged k-mers as DNA strings."
    },
    {
        "load_tagset",
        (PyCFunction)hashgraph_load_tagset, METH_VARARGS,
        "Load tags from a file."
    },
    {
        "save_tagset",
        (PyCFunction)hashgraph_save_tagset, METH_VARARGS,
        "Save tags to a file."
    },
    {
        "n_tags",
        (PyCFunction)hashgraph_n_tags, METH_VARARGS,
        "Return the count of all tags."
    },
    {
        "divide_tags_into_subsets",
        (PyCFunction)hashgraph_divide_tags_into_subsets, METH_VARARGS,
        "Divide tags equally up into subsets of given size."
    },
    {
        "_get_tag_density",
        (PyCFunction)hashgraph__get_tag_density, METH_VARARGS,
        "Get the tagging density."
    },
    {
        "_set_tag_density",
        (PyCFunction)hashgraph__set_tag_density, METH_VARARGS,
        "Set the tagging density."
    },

    //
    // partitioning
    //
    {
        "do_subset_partition",
        (PyCFunction)hashgraph_do_subset_partition, METH_VARARGS,
        "Partition the graph starting from a given subset of tags."
    },
    {
        "find_all_tags",
        (PyCFunction)hashgraph_find_all_tags, METH_VARARGS,
        "Starting from the given k-mer, find all closely connected tags."
    },
    {
        "assign_partition_id",
        (PyCFunction)hashgraph_assign_partition_id, METH_VARARGS,
        "Assign a partition ID to a given tag."
    },
    {
        "output_partitions",
        (PyCFunction)hashgraph_output_partitions, METH_VARARGS,
        "Write out sequences in given filename to another file, annotating "
        "with partition IDs."
    },
    {
        "load_partitionmap",
        (PyCFunction)hashgraph_load_partitionmap, METH_VARARGS,
        "Load a partitionmap for a given subset."
    },
    {
        "save_partitionmap",
        (PyCFunction)hashgraph_save_partitionmap, METH_VARARGS,
        "Save a partitionmap for the given subset."
    },
    {
        "_validate_partitionmap",
        (PyCFunction)hashgraph__validate_partitionmap, METH_VARARGS,
        "Run internal validation checks."
    },
    {
        "consume_seqfile_and_tag_with_reads_parser",
        (PyCFunction)hashgraph_consume_seqfile_and_tag_with_reads_parser,
        METH_VARARGS,
        "Count all k-mers using the given reads parser"
    },
    {
        "consume_partitioned_fasta",
        (PyCFunction)hashgraph_consume_partitioned_fasta, METH_VARARGS,
        "Count all k-mers in a given file"
    },
    {
        "merge_subset",
        (PyCFunction)hashgraph_merge_subset, METH_VARARGS,
        "Merge the given subset into this one."
    },
    {
        "merge_subset_from_disk",
        (PyCFunction)hashgraph_merge_from_disk, METH_VARARGS,
        "Merge the given subset (filename) into this one."
    },
    {
        "count_partitions",
        (PyCFunction)hashgraph_count_partitions, METH_VARARGS,
        "Count the number of partitions in the master partitionmap."
    },
    {
        "subset_count_partitions",
        (PyCFunction)hashgraph_subset_count_partitions, METH_VARARGS,
        "Count the number of partitions in this subset partitionmap."
    },
    {
        "subset_partition_size_distribution",
        (PyCFunction)hashgraph_subset_partition_size_distribution,
        METH_VARARGS,
        "Get the size distribution of partitions in this subset."
    },
    {
        "save_subset_partitionmap",
        (PyCFunction)hashgraph_save_subset_partitionmap, METH_VARARGS,
        "Save the partition map for this subset."
    },
    {
        "load_subset_partitionmap",
        (PyCFunction)hashgraph_load_subset_partitionmap, METH_VARARGS,
        "Save the partition map for this subset."
    },
    {
        "_validate_subset_partitionmap",
        (PyCFunction)hashgraph__validate_subset_partitionmap, METH_VARARGS,
        "Run internal validation checks on this subset."
    },
    {
        "set_partition_id",
        (PyCFunction)hashgraph_set_partition_id, METH_VARARGS,
        "Set the partition ID for this tag."
    },
    {
        "join_partitions",
        (PyCFunction)hashgraph_join_partitions, METH_VARARGS,
        "Join the partitions of these two tags."
    },
    {
        "get_partition_id",
        (PyCFunction)hashgraph_get_partition_id, METH_VARARGS,
        "Get the partition ID of this tag."
    },
    {
        "repartition_largest_partition",
        (PyCFunction)hashgraph_repartition_largest_partition, METH_VARARGS,
        "Repartition the largest partition (in the face of stop tags)."
    },

    // stop tags
    {
        "load_stop_tags",
        (PyCFunction)hashgraph_load_stop_tags, METH_VARARGS,
        "Load the set of stop tags."
    },
    {
        "save_stop_tags",
        (PyCFunction)hashgraph_save_stop_tags, METH_VARARGS,
        "Save the set of stop tags."
    },
    {
        "print_stop_tags",
        (PyCFunction)hashgraph_print_stop_tags, METH_VARARGS,
        "Print out the set of stop tags."
    },
    {
        "trim_on_stoptags",
        (PyCFunction)hashgraph_trim_on_stoptags, METH_VARARGS,
        "Trim the reads on the given stop tags."
    },
    {
        "add_stop_tag",
        (PyCFunction)hashgraph_add_stop_tag, METH_VARARGS,
        "Add this k-mer as a stop tag."
    },
    {
        "get_stop_tags",
        (PyCFunction)hashgraph_get_stop_tags, METH_VARARGS,
        "Return a DNA list of all of the stop tags."
    },
    {NULL, NULL, 0, NULL}           /* sentinel */
};

}

#endif
