"""Tests for `taeper` package."""
import unittest
import logging
from bam2fast5 import bam2fast5

logging.disable(logging.CRITICAL)


class TestReadIdExtraction(unittest.TestCase):
    """Test the read id extracxtion functions for bam, sam and fastq."""

    def test_BamReadIdExtractionTB_SixReadIds(self):
        """Test read id extracxtion from bam with 6 mapped reads"""
        zulu_time = "2018-01-03T16:45:30Z"
        result = taeper._zulu_to_epoch_time(zulu_time)
        expected = 1514997930.0
        self.assertEqual(result, expected)