#! /usr/bin/env python3
# $ -S $HOME/.pyenv/shims/python3
# $ -l s_vmem=32G -l mem_req=32G
# $ -cwd
# $ -o $HOME/ugelogs/
# $ -e $HOME/ugelogs/

"""
Convert positions from genome coordinates to transcript coordinates

Usage:
  genome2tx [options] <gtf> <bam>

Options:
  --strandness TYPE : none/rf/fr [default: none]
  <gtf>             : GTF formatted gene annotation file
  <bam>             : BAM (sorted) formatted alignment file

"""

import os
import sys
import itertools
import copy
from typing import Iterable
import time

from docopt import docopt
import pysam
import yaml
import pandas as pd
from gtfparse import read_gtf
from sqlite3 import Row
from sqlalchemy import create_engine
from numba import jit


def to_ranges(positions: Iterable[int]):
    try:
        if len(positions) != len(set(positions)):
            sys.stderr.write(str(positions))
            raise Exception('Duplicated position!')
    except Exception as e:
        sys.stderr.write(e)
        sys.exit(1)

    if type(positions) is list:
        positions.sort()

    for _, g in itertools.groupby(enumerate(positions), lambda x: x[1] - x[0]):
        groups = list(g)
        range_ = range(groups[0][1], -~groups[-1][1])
        yield range_


@jit
def strand_str(alignment, strandness):
    # TODO: Stranded sequence
    if strandness == 'rf':
        if alignment.is_read1:
            return '+' if alignment.is_reverse else '-'
        else:
            return '-' if alignment.is_reverse else '+'
    elif strandness == 'fr':
        if alignment.is_read1:
            return '-' if alignment.is_reverse else '+'
        else:
            return '+' if alignment.is_reverse else '-'
    else:
        return ('+', '-')


def pickle_annotations(gtf_path, columns: list, features=['exon']):
    annotations = read_gtf(
        gtf_path,
        features=features).filter(
            columns).sort_values(by=['seqname', 'transcript_id', 'exon_number'])

    annotations.to_pickle("annotations.pkl")
    sys.exit()


def load_annotations(gtf_path, columns: list, features=['exon']):
    if not __debug__:
        # load from pickled object
        annotations = pd.read_pickle(gtf_path)
        return annotations

    # NOTE: GENCODE data format;
    # https://www.gencodegenes.org/pages/data_format.html
    annotations = read_gtf(gtf_path
                           ).query("feature == {}".format(features)).filter(
        columns).sort_values(by=['seqname', 'transcript_id', 'exon_number'])

    return annotations


def load_alignments(bam_path):
    # NOTE: pysam.AlignedSegment;
    # https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment
    THREADS = 1

    try:
        alignments = pysam.AlignmentFile(bam_path, mode='rb', threads=THREADS)
        alignments.check_index()
    except (ValueError, AttributeError) as e:
        sys.stderr.write("Index is not found: {}\n".format(e))
        sys.stderr.write('Create index...')
        pysam.index(bam_path)
        try:
            alignments = pysam.AlignmentFile(bam_path, mode='rb', threads=THREADS)
            alignments.check_index()
        except Exception as e:
            sys.stderr.write("Alignments file open was failed: {}\n".format(e))
            sys.exit()

    return alignments


@jit
def offsets(annotations):
    offsets = [None] * len(annotations)
    transcript_id_previous = end_previous = offset = None

    transcript_ids = list(annotations['transcript_id'])
    starts = list(annotations['start'])
    ends = list(annotations['end'])

    for i, (t, s, e) in enumerate(zip(transcript_ids, starts, ends)):
        if transcript_id_previous != t:
            offset = ~-s
            transcript_id_previous = t
        else:
            offset = offset + ~-s - end_previous

        end_previous = e
        offsets[i] = offset

    return offsets


def create_db_engine(
        annotations,
        chromosomes,
        table='exons',
        columns_index=['strand', 'transcript_id', 'start', 'end', ('strand', 'transcript_id'),  ('strand', 'start'), ('strand', 'end'), ('strand', 'start', 'end')]):
    db_engine = create_engine('sqlite://', echo=False)
    db_engine.row_factory = Row
    annotations.query(
        "seqname == {}".format(chromosomes)
    ).to_sql(table, con=db_engine, if_exists='replace')

    for c in columns_index:
        if type(c) is str:
            c = (c, )
        query = "CREATE INDEX ix_{0}_{1} on {0}({2})".format(
            table,
            '_'.join(c),
            ', '.join(c))
        db_engine.execute(query)

    return db_engine


# @jit
# def _build_query(strand, positions, table, columns):
#     columns_str = ', '.join(columns)
#     query = ''

#     # FIXME:
#     if strand == ('+', '-'):
#         where_clauses = [
#             "start >= {1} AND start <= {2}",
#             "end >= {1} AND end <= {2}",
#             "start <= {1} AND end >= {2}"
#         ]
#     else:
#         where_clauses = [
#             "strand = '{0}' AND start >= {1} AND start <= {2}",
#             "strand = '{0}' AND end >= {1} AND end <= {2}",
#             "strand = '{0}' AND start <= {1} AND end >= {2}"
#         ]

#     for r in to_ranges(positions):
#         for w in where_clauses:
#             if query:
#                 query += ' UNION ALL '
#             query += "SELECT {} FROM {} WHERE {}".format(
#                 columns_str,
#                 table,
#                 w.format(strand, r.start + 1, r.stop - 2))

#     print(query)
#     return query


def build_query(strand, transcript_id, table, columns):
    columns_str = ', '.join(columns)

    if strand == ('+', '-'):
        where_clause = "transcript_id == '{}'".format(transcript_id)
    else:
        where_clause = "strand = '{}' AND transcript_id == '{}'".format(strand, transcript_id)

    query = "SELECT {} FROM {} WHERE {}".format(
                columns_str,
                table,
                where_clause)

    return query


@jit
def overlaps(
        strand,
        transcript_id,
        annotations,
        table='exons',
        columns=['start', 'end']):
    results = annotations.execute(build_query(
        strand, transcript_id, table, columns)).fetchall()
    return results


# @jit
# def overlaps(
#         strand,
#         positions,
#         annotations,
#         table='exons',
#         columns=['transcript_id']):
#     results = annotations.execute(build_query(
#         strand, positions, table, columns)).fetchall()
#     return results


def transcript_id_derived_from(query_name):
    # NOTE: polyester generated reads are named according to the format:
    # 'readX/transcript_id'
    # e.g. read1/ENST00000501122.2

    try:
        transcript_id = query_name.split('/')[1].split(';')[0]
        return transcript_id
    except Exception:
        sys.stderr.write("Query name ({}) is invalid.\n".format(query_name))
        return query_name


def init_stats(keys_, default=None):
    stats = {}

    for key in keys_:
        if type(key) is dict:
            for k in key.keys():
                stats[k] = init_stats(key[k])
        else:
            stats[key] = [] if not default else default

    return stats


def counted_up_dict(dict_):
    for k, v in dict_.items():
        if type(v) is dict:
            dict_[k] = counted_up_dict(v)
        else:
            dict_[k] = len(v)

    return dict_


def print_results(stats: dict):
    print('-' * 80)
    print(yaml.dump(stats))
    print(counted_up_dict(copy.deepcopy(stats)))


# TODO: Exact match
# def list_relpositions(positions: Iterable[int], exons):
#     list_ = []

#     transcript_ids = list(exons['transcript_id'])
#     starts = list(exons['start'])
#     ends = list(exons['end'])
#     offsets = list(exons['offset'])

#     for i, (t, s, e, o) in enumerate(zip(transcript_ids, starts, ends, offsets)):
#         exon_positions = range(s, -~e)
#         positions_overlapped = sorted(set(positions).intersection(set(exon_positions)))
#         relpositions = [p - o for p in positions_overlapped]
#         list_[i] = (transcript_ids[i], relpositions)

#     return list_


# TODO: Exact match
# def transcript_relregions(list_relpositions: List[List[int]]):
#     list_relpositions.sort(key=lambda x: x[0])
#     for transcript_id, g in itertools.groupby(list_relpositions, lambda x: x[0]):
#         relregions = list(to_ranges([p for positions in g for p in positions[1]]))

#         yield transcript_id, relregions


def main():
    start = time.time()

    options = docopt(__doc__)
    strandness = options['--strandness']
    gtf_path = options['<gtf>']
    bam_path = options['<bam>']

    for path in [gtf_path, bam_path]:
        try:
            if not os.path.exists(path):
                raise Exception("File not exists!: {}.".format(path))
        except Exception as e:
            sys.stderr.write(e.args)
            sys.exit(1)

    columns_select = ['seqname', 'strand', 'start', 'end',
                      'transcript_id', 'exon_number']
    annotations = load_annotations(gtf_path, columns_select)

    print('gtfparse:', time.time() - start)

    annotations['offset'] = offsets(annotations)
    columns_select.append('offset')

    print('offsets:', time.time() - start)

    alignments = load_alignments(bam_path)

    print('load_alignments:', time.time() - start)

    keys = [{'mapped': ['true', 'false']}, 'unmapped']
    stats = init_stats(keys)
    chromosome = None
    for i, alignment in enumerate(alignments.fetch()):

        if i % 100000 == 0:
            print('alignment-{}:'.format(i), time.time() - start)

        if alignment.is_unmapped:
            mapped_or_unmapped = 'unmapped'
            stats[mapped_or_unmapped].append(i)
            continue

        mapped_or_unmapped = 'mapped'

        # OPTIMIZE: Total or partial?
        reference_name = alignment.reference_name
        if chromosome != reference_name:
            chromosome = reference_name
            annotations_chrXX = create_db_engine(annotations, [chromosome])

        strand = strand_str(alignment, strandness)
        # NOTE: pysam return 0 based postions
        positions = alignment.get_reference_positions()

        transcript_id = transcript_id_derived_from(alignment.query_name)
        exons_overlapped = overlaps(strand, transcript_id, annotations_chrXX)
        # NOTE: ajust to 0 based positions
        positions_exons = set().union(*[set(range(e[0] - 1, e[1])) for e in exons_overlapped])

        # CASE: TRUE POSITIVE
        if len(positions_exons.intersection(positions)) > 1:
            true_or_false = 'true'
        else:
            true_or_false = 'false'

        stats[mapped_or_unmapped][true_or_false].append(i)

        # print("{} query_name: {} {} {}".format(i, alignment.query_name, str(alignment.is_reverse), true_or_false))
        # print("{}\n{}".format(list(to_ranges(positions)), list(to_ranges(positions_exons))))

    annotations_chrXX = None
    alignments.close()

    print('close:', time.time() - start)

    print_results(stats)


if __name__ == '__main__':
    start = time.time()
    main()
    print('main:', time.time() - start)
