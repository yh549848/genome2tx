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


# @profile
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
    list_ = []
    transcript_id = exon_previous = None

    for _, exon in annotations.iterrows():
        if transcript_id != exon['transcript_id']:
            transcript_id = exon['transcript_id']
            offset = ~-exon['start']
        else:
            offset = offset + ~-exon['start'] - exon_previous['end']

        exon_previous = exon
        list_.append(offset)

    return list_


# @profile
def exons_overlapped(annotations, strand: str, positions: List[int]):
    for r in to_ranges(positions):
        q = "strand == '{0}' & "\
            "(start >= {1} & start <= {2} | end >= {1} & end <= {2} | start <= {1} & end >= {2})".format(
                strand, r.start, ~-r.stop)
        exons = annotations.query(q)
        yield exons


# @profile
def list_relpositions(positions: List[int], exons):
    list_ = []
    for _, exon in exons.iterrows():
        exon_positions = range(exon['start'], -~exon['end'])
        positions_overlapped = sorted(set(positions).intersection(set(exon_positions)))
        relpositions = [p - exon['offset'] for p in positions_overlapped]
        list_.append((exon['transcript_id'], relpositions))

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
            raise Exception("File not exists: {}.".format(path))

    # NOTE: GENCODE data format;
    # https://www.gencodegenes.org/pages/data_format.html
    columns_selected = ['seqname', 'strand', 'start', 'end', 'gene_id', 'transcript_id', 'exon_id',
                        'exon_number']
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

        if chromosome != alignment.reference_name:
            chromosome = alignment.reference_name
            annotations_chr = annotations.query("seqname == '{}'".format(chromosome))

        positions = alignment.get_reference_positions()
        strand = '-' if alignment.is_reverse else '+'

        is_assigned = False
        for exons in exons_overlapped(annotations_chr, strand, positions):
            if len(exons) < 1:
                continue

            if derived_from(alignment.query_name) in set(exons['transcript_id']):
                paired_or_unpaired = 'paired' if alignment.is_paired else 'unpaired'
                stats['mapped'][paired_or_unpaired] += 1
                ids['mapped'][paired_or_unpaired].append(i)
                is_assigned = True
                break

            # NOTE: Obsoleted exact match
            # for transcript_id, regions in transcript_relregions(
            #         list_relpositions(positions, exons)):
            #     print("\t".join([i, derived_from(alignment.query_name), transcript_id, str(regions)]))

        if not is_assigned:
            stats['mapped']['tx_unmapped'] += 1
            ids['mapped']['tx_unmapped'].append(i)

    alignments.close()

    print('-' * 20)
    print(yaml.dump(stats))
    print(yaml.dump(ids))
    print('-' * 20)


if __name__ == '__main__':
    main()
