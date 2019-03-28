#! /usr/bin/env python3

"""
Convert positions from genome coordinates to transcript coordinates

Usage:
  genome2tx <gtf> <bam>

"""

import os
import sys
import itertools
from typing import Iterable, List

from docopt import docopt
import pysam
import yaml
from gtfparse import read_gtf


def to_ranges(positions: Iterable[int]):
    try:
        if len(positions) != len(set(positions)):
            raise Exception('Duplicated position!')

        if type(positions) is list:
            positions.sort()

        enumerate(positions)
    except Exception:
        sys.stderr.write(str(positions))
        sys.exit(1)

    for _, g in itertools.groupby(enumerate(positions), lambda x: x[1] - x[0]):
        groups = list(g)
        yield range(groups[0][1], -~groups[-1][1])


# @profile
def offsets(annotations):
    list_ = [None] * len(annotations)
    transcript_id_previous = end_previous = None

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
        list_[i] = offset

    return list_


# @profile
def exons_overlapped(annotations, strand: str, positions: List[int]):
    query = '|'.join(["strand == '{0}' & "
                      "(start >= {1} & start <= {2} | end >= {1} & end <= {2} | start <= {1} & end >= {2})".format(
                          strand, r.start, ~-r.stop) for r in to_ranges(positions)])

    exons = annotations.query(query)
    return exons


# @profile
def list_relpositions(positions: List[int], exons):
    list_ = []

    transcript_ids = list(exons['transcript_id'])
    starts = list(exons['start'])
    ends = list(exons['end'])
    offsets = list(exons['offset'])

    for i, (t, s, e, o) in enumerate(zip(transcript_ids, starts, ends, offsets)):
        exon_positions = range(s, -~e)
        positions_overlapped = sorted(set(positions).intersection(set(exon_positions)))
        relpositions = [p - o for p in positions_overlapped]
        list_[i] = (transcript_ids[i], relpositions)

    return list_


# @profile
def transcript_relregions(list_relpositions: List[List[int]]):
    list_relpositions.sort(key=lambda x: x[0])
    for transcript_id, g in itertools.groupby(list_relpositions, lambda x: x[0]):
        relregions = list(to_ranges([p for positions in g for p in positions[1]]))

        yield transcript_id, relregions


# @profile
def derived_from(query_name):
    # Query name format must be below:
    # readXX/transcript_id
    try:
        transcript_id = query_name.split('/')[1]
        return transcript_id
    except Exception:
        sys.stderr.write("Query name ({}) is invalid.".format(query_name))
        return query_name


# @profile
def main():
    options = docopt(__doc__)
    gtf_path = options['<gtf>']
    bam_path = options['<bam>']

    for path in [gtf_path, bam_path]:
        if not os.path.exists(path):
            raise Exception("File not exists!: {}.".format(path))

    # NOTE: GENCODE data format;
    # https://www.gencodegenes.org/pages/data_format.html
    columns_selected = ['seqname', 'strand', 'start', 'end',
                        'gene_id', 'transcript_id', 'exon_id', 'exon_number']
    annotations = read_gtf(gtf_path).query("feature == 'exon'").filter(columns_selected).sort_values(
        by=['seqname', 'transcript_id', 'exon_number'])

    annotations['offset'] = offsets(annotations)

    alignments = pysam.AlignmentFile(bam_path, "rb")
    try:
        alignments.check_index()
    except Exception:
        pysam.index(bam_path)
        sys.stderr.write("Index is not found.")

    # NOTE: pysam.AlignedSegment;
    # https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment
    stats = {
        'mapped': {
            'paired': 0,
            'unpaired': 0,
            'tx_unmapped': 0
        },
        'unmapped': 0}
    ids = {
        'mapped': {
            'paired': [],
            'unpaired': [],
            'tx_unmapped': []
        },
        'unmapped': []}

    chromosome = None
    for i, alignment in enumerate(alignments):
        if alignment.is_unmapped:
            stats['unmapped'] += 1
            ids['unmapped'].append(i)
            continue

        reference_name = alignment.reference_name
        if chromosome != reference_name:
            chromosome = reference_name
            annotations_chr = annotations.query("seqname == '{}'".format(chromosome))

        positions = alignment.get_reference_positions()
        strand = '-' if alignment.is_reverse else '+'

        exons = exons_overlapped(annotations_chr, strand, positions)
        if len(exons) < 1:
            stats['mapped']['tx_unmapped'] += 1
            ids['mapped']['tx_unmapped'].append(i)
            continue

        if derived_from(alignment.query_name) in set(exons['transcript_id']):
            paired_or_unpaired = 'paired' if alignment.is_paired else 'unpaired'
            stats['mapped'][paired_or_unpaired] += 1
            ids['mapped'][paired_or_unpaired].append(i)
        else:
            stats['mapped']['tx_unmapped'] += 1
            ids['mapped']['tx_unmapped'].append(i)

        # NOTE: Obsoleted exact match
        # for transcript_id, regions in transcript_relregions(
        #         list_relpositions(positions, exons)):
        #     print("\t".join([i, derived_from(alignment.query_name), transcript_id, str(regions)]))

    alignments.close()

    print('-' * 20)
    print(yaml.dump(stats))
    print(yaml.dump(ids))
    print('-' * 20)


if __name__ == '__main__':
    main()
