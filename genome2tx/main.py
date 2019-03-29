#! /usr/bin/env python3
#$ -S $HOME/.pyenv/shims/python3
#$ -l s_vmem=32G -l mem_req=32G
#$ -cwd
#$ -o $HOME/ugelogs/
#$ -e $HOME/ugelogs/

"""
Convert positions from genome coordinates to transcript coordinates

Usage:
  genome2tx <gtf> <bam>

"""

import os
import sys
import itertools
import copy
from typing import Iterable, List
import time

from docopt import docopt
import pysam
import yaml
from gtfparse import read_gtf
import sqlite3
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


@jit
def create_db_engine(annotations,
                     chromosomes,
                     table='exons',
                     columns_index=['strand', 'start', 'end']):
    db_engine = create_engine('sqlite://', echo=False)
    db_engine.row_factory = sqlite3.Row
    annotations.query(
        "seqname == {}".format(chromosomes)
    ).to_sql(table, con=db_engine, if_exists='replace')

    # OPTIMIZE: Reconsider target columns
    for c in columns_index:
        q = "CREATE INDEX ix_{0}_{1} on {0}({1});".format(table, c)
        _ = db_engine.execute(q)

    return db_engine


@jit
# @profile
def exons_overlapped(annotations, strand: str,
                     positions: List[int],
                     columns='transcript_id'):
    columns_str = ', '.join(columns)
    query = ''

    # HACK: Consider query building
    for r in to_ranges(positions):
        where_ = "(strand = '{0}' AND ("\
                 "start >= {1} AND start <= {2} OR "\
                 "end >= {1} AND end <= {2} OR "\
                 "start <= {1} AND end >= {2})"\
                 ")".format(
                     strand, r.start, ~-r.stop)

        if query:
            query += ' UNION '
        query += "SELECT {} FROM exons WHERE {}".format(columns_str, where_)

    exons = annotations.execute(query).fetchall()
    return exons


def transcript_id_derived_from(query_name):
    try:
        transcript_id = query_name.split('/')[1]
        return transcript_id
    except Exception:
        sys.stderr.write("Query name ({}) is invalid.".format(query_name))
        return query_name


def init_stats(keys_=[{'mapped': [{'true': ['paired', 'unpaired'], 'false': ['paired', 'unpaired']}]}, 'unmapped'],
               default=None):
    stats = {}

    for key in keys_:
        if type(key) is dict:
            for k in key.keys():
                stats[k] = init_stats(key[k])
        else:
            stats[key] = [] if not default else default

    return stats


def count_up_dict(dict_):
    for k, v in dict_.items():
        if type(v) is dict:
            dict_[k] = count_up_dict(v)
        else:
            dict_[k] = len(v)

    return dict_


def print_results(stats: dict):
    print('-' * 80)
    print(yaml.dump(stats))
    print(count_up_dict(copy.deepcopy(stats)))


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


# @profile
def main():
    start = time.time()

    options = docopt(__doc__)
    gtf_path = options['<gtf>']
    bam_path = options['<bam>']

    for path in [gtf_path, bam_path]:
        try:
            if not os.path.exists(path):
                raise Exception("File not exists!: {}.".format(path))
        except Exception as e:
            sys.stderr.write(e)
            sys.exit(1)

    # NOTE: GENCODE data format;
    # https://www.gencodegenes.org/pages/data_format.html
    columns_select = ['seqname', 'strand', 'start', 'end',
                      'transcript_id', 'exon_number']
    annotations = read_gtf(gtf_path,
                           features=['exon']).filter(
                               columns_select).sort_values(
                                   by=['seqname', 'transcript_id', 'exon_number'])

    print('gtfparse:', time.time() - start)

    annotations['offset'] = offsets(annotations)
    columns_select.append('offset')

    print('offsets:', time.time() - start)

    # NOTE: pysam.AlignedSegment;
    # https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment
    THREADS = 2
    try:
        alignments = pysam.AlignmentFile(bam_path, mode='rb', threads=THREADS)
        alignments.check_index()
    except (ValueError, AttributeError) as e:
        sys.stderr.write("Index is not found: {}".format(e))
        sys.stderr.write('Create index...')
        pysam.index(bam_path)
        alignments.check_index()
    except Exception as e:
        sys.stderr.write("Alignments file open was failed: {}".format(e))
        sys.exit()

    print('load_alignments:', time.time() - start)
    stats = init_stats()
    chromosome = None
    for i, alignment in enumerate(alignments.fetch(contig='chr1')):

        if i % 1000 == 0:
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

        positions = alignment.get_reference_positions()
        strand = '-' if alignment.is_reverse else '+'

        exons = exons_overlapped(annotations_chrXX, strand, positions, columns_select)
        transcript_ids = {e['transcript_id'] for e in exons}
        paired_or_unpaired = 'paired' if alignment.is_paired else 'unpaired'

        if len(transcript_ids) < 1:
            true_or_false = 'false'
            stats[mapped_or_unmapped][true_or_false][paired_or_unpaired].append(i)
            continue

        true_or_false = 'true' if transcript_id_derived_from(alignment.query_name) in transcript_ids else 'false'
        stats[mapped_or_unmapped][true_or_false][paired_or_unpaired].append(i)

        # TODO: Exact match
        # for transcript_id, regions in transcript_relregions(
        #         list_relpositions(positions, exons)):
        #             pass

    annotations_chrXX = None
    alignments.close()

    print('close:', time.time() - start)

    print_results(stats)


if __name__ == '__main__':
    start = time.time()
    main()
    print('main:', time.time() - start)
