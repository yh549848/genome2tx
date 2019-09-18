"""
This is test for genome2tx
"""
import sys
import unittest
import copy

from genome2tx import main


class TestUtils(unittest.TestCase):
    def test_to_ranges(self):
        set_ = {1, 5, 6, 8, 3, 2}

        expected = [range(1, 4), range(5, 7), range(8, 9)]
        actual = list(main.to_ranges(set_))

        self.assertEqual(expected, actual)

    def strand_str(self):
        pass

    def test_offsets(self):
        pass

    def create_db_engine(self):
        pass

    def overlaps(self):
        pass

    def test_derived_from(self):
        # NEAT1_2
        str_ = 'read1/ENST00000501122.2'

        expected = 'ENST00000501122.2'
        actual = main.transcript_id_derived_from(str_)

        self.assertEqual(expected, actual)

    def test_build_query(self):
        # NEAT1_2
        strand = '+'
        # 5'
        positions = list(range(65422800, 65423212))
        # 3'
        positions.extend(list(range(65445540 - 100, 65445540)))
        table = 'annotations'
        columns = ['seqname', 'start', 'end', 'transcript_id', 'transcript_name']

        actual = main.build_query(strand, positions, table, columns)
        print(sys._getframe().f_code.co_name)
        print(actual)

    def test_init_stats(self):
        expected = {
            '0':
                {
                    '0': [],
                    '1': []
                },
            '1': []
        }

        keys = [{'0': ['0', '1']}, '1']
        actual = main.init_stats(keys)

        self.assertEqual(expected, actual)

    def test_count_up_dict(self):
        keys = [{'0': ['0', '1']}, '1']
        actual = main.init_stats(keys)
        expected = copy.deepcopy(actual)

        actual['0']['0'].append(1)
        for i in range(10):
            actual['1'].append(1)

        main.counted_up_dict(actual)

        expected['0']['0'] = 1
        expected['0']['1'] = 0
        expected['1'] = 10

        self.assertEqual(expected, actual)


if __name__ == '__main__':
    unittest.main()
