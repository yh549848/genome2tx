"""
This is test for geneome2tx
"""

import unittest

from genome2tx import main


class TestUtils(unittest.TestCase):
    def test_to_ranges(self):
        set_ = {1, 5, 6, 8, 3, 2}

        expected = [range(1, 4), range(5, 7), range(8, 9)]
        actual = list(main.to_ranges(set_))

        self.assertEqual(expected, actual)

    def test_offsets(self):
        pass

    def test_exons_overlapped(self):
        pass

    def test_derived_from(self):
        pass


if __name__ == '__main__':
    unittest.main()
