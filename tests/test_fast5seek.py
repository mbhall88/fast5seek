"""Tests for fast5seek package."""
import unittest
import logging
import io
from contextlib import redirect_stdout
from fast5seek import fast5seek

logging.disable(logging.CRITICAL)


class TestBamReadIdExtraction(unittest.TestCase):
    """Test the read id extracxtion functions for bam and sam"""

    def test_BamReadIdExtractionTBMapped_SixReadIds(self):
        """Test read id extracxtion from bam with 6 mapped reads"""
        mapped = True
        bam = 'tests/data/bam/tb.bam'
        result = fast5seek.get_sam_read_ids(bam, mapped)
        expected = {
            '57d4cd63-3189-4006-93ec-bf3c8bfb2ce1',
            'bbd563e9-1bf8-4268-92d5-45ccb8e3da72',
            '8aecf428-af00-4791-b065-5d4abd798a29',
            'd707ff64-6ade-477a-8b68-0b3c394ef3b1',
            'c967d421-3da4-4e11-accf-0d3fb0155840',
            '28acaa47-1cec-4a91-9abf-c780d27e6cc4'
        }
        self.assertSetEqual(result, expected)

    def test_BamReadIdExtractionTBUnmapped_EightReadIds(self):
        """Test read id extracxtion from bam with 6 mapped reads"""
        mapped = False
        bam = 'tests/data/bam/tb.bam'
        result = fast5seek.get_sam_read_ids(bam, mapped)
        expected = {
            '57d4cd63-3189-4006-93ec-bf3c8bfb2ce1',
            'bbd563e9-1bf8-4268-92d5-45ccb8e3da72',
            '8aecf428-af00-4791-b065-5d4abd798a29',
            'd707ff64-6ade-477a-8b68-0b3c394ef3b1',
            'c967d421-3da4-4e11-accf-0d3fb0155840',
            '28acaa47-1cec-4a91-9abf-c780d27e6cc4',
            '6cf511b6-1724-46bd-b5a4-59c18bb57343',
            '6c26d9b5-d892-4fc6-b035-abe575895c88'
        }
        self.assertSetEqual(result, expected)

    def test_SamReadIdExtractionTBMapped_SixReadIds(self):
        """Test read id extracxtion from bam with 6 mapped reads"""
        mapped = True
        sam = 'tests/data/sam/tb.sam'
        result = fast5seek.get_sam_read_ids(sam, mapped)
        expected = {
            '57d4cd63-3189-4006-93ec-bf3c8bfb2ce1',
            'bbd563e9-1bf8-4268-92d5-45ccb8e3da72',
            '8aecf428-af00-4791-b065-5d4abd798a29',
            'd707ff64-6ade-477a-8b68-0b3c394ef3b1',
            'c967d421-3da4-4e11-accf-0d3fb0155840',
            '28acaa47-1cec-4a91-9abf-c780d27e6cc4'
        }
        self.assertSetEqual(result, expected)

    def test_SamReadIdExtractionEcoliMapped_TwoReadIds(self):
        """Test read id extracxtion from bam with 6 mapped reads"""
        mapped = True
        sam = 'tests/data/sam/ecoli.sam'
        result = fast5seek.get_sam_read_ids(sam, mapped)
        expected = {
            '6cf511b6-1724-46bd-b5a4-59c18bb57343',
            '6c26d9b5-d892-4fc6-b035-abe575895c88'
        }
        self.assertSetEqual(result, expected)


class TestFastqReadIdExtraction(unittest.TestCase):
    """Test the read id extracxtion functions for fastq"""

    def test_FastqReadIdExtractionTB_SixReadIds(self):
        """Test read id extracxtion from gzipped fastq with 6 mapped reads"""
        fastq = 'tests/data/fastq/tb_mapped.fastq'
        result = fast5seek.get_fastq_read_ids(fastq)
        expected = {
            '57d4cd63-3189-4006-93ec-bf3c8bfb2ce1',
            'bbd563e9-1bf8-4268-92d5-45ccb8e3da72',
            '8aecf428-af00-4791-b065-5d4abd798a29',
            'd707ff64-6ade-477a-8b68-0b3c394ef3b1',
            'c967d421-3da4-4e11-accf-0d3fb0155840',
            '28acaa47-1cec-4a91-9abf-c780d27e6cc4'
        }
        self.assertSetEqual(result, expected)

    def test_FastqGzipReadIdExtractionEcoli_TwoReadIds(self):
        """Test read id extracxtion from fastq with 2 mapped reads"""
        fastq = 'tests/data/fastq/ecoli_mapped.fastq.gz'
        result = fast5seek.get_fastq_read_ids(fastq)
        expected = {
            '6cf511b6-1724-46bd-b5a4-59c18bb57343',
            '6c26d9b5-d892-4fc6-b035-abe575895c88'
        }
        self.assertSetEqual(result, expected)

    def test_FastqGzipReadIdExtractionAll_EightReadIds(self):
        """Test read id extracxtion from fastq with 8 mapped reads"""
        fastq = 'tests/data/fastq/basecalled.fastq.gz'
        result = fast5seek.get_fastq_read_ids(fastq)
        expected = {
            '57d4cd63-3189-4006-93ec-bf3c8bfb2ce1',
            'bbd563e9-1bf8-4268-92d5-45ccb8e3da72',
            '8aecf428-af00-4791-b065-5d4abd798a29',
            'd707ff64-6ade-477a-8b68-0b3c394ef3b1',
            'c967d421-3da4-4e11-accf-0d3fb0155840',
            '28acaa47-1cec-4a91-9abf-c780d27e6cc4',
            '6cf511b6-1724-46bd-b5a4-59c18bb57343',
            '6c26d9b5-d892-4fc6-b035-abe575895c88'
        }
        self.assertSetEqual(result, expected)


class TestReadIdExtraction(unittest.TestCase):
    """Test the read id extracxtion functions for fastq, bam and sam"""

    def test_ReadIdExtractionWithFast5_EmptySet(self):
        """Test read id extracxtion from gzipped fastq with 6 mapped reads"""
        mapped = False
        files = ['tests/data/fast5/tb1.fast5']
        result = fast5seek.extract_read_ids(files, mapped)
        expected = set()
        self.assertSetEqual(result, expected)

    def test_ReadIdExtractionWithFast5Fastq_TwoIds(self):
        """Test read id extracxtion from gzipped fastq with 6 mapped reads"""
        mapped = False
        files = [
            'tests/data/fast5/tb1.fast5',
            'tests/data/fastq/ecoli_mapped.fastq.gz'
        ]
        result = fast5seek.extract_read_ids(files, mapped)
        expected = {
            '6cf511b6-1724-46bd-b5a4-59c18bb57343',
            '6c26d9b5-d892-4fc6-b035-abe575895c88'
        }
        self.assertSetEqual(result, expected)

    def test_ReadIdExtractionWithFast5FastqBam_TwoIds(self):
        """Test read id extracxtion from gzipped fastq with 6 mapped reads"""
        mapped = True
        files = [
            'tests/data/fast5/tb1.fast5',
            'tests/data/fastq/ecoli_mapped.fastq.gz',
            'tests/data/bam/ecoli.bam'
        ]
        result = fast5seek.extract_read_ids(files, mapped)
        expected = {
            '6cf511b6-1724-46bd-b5a4-59c18bb57343',
            '6c26d9b5-d892-4fc6-b035-abe575895c88'
        }
        self.assertSetEqual(result, expected)

    def test_ReadIdExtractionWithFast5FastqBam_EightIds(self):
        """Test read id extracxtion from gzipped fastq with 6 mapped reads"""
        mapped = True
        files = [
            'tests/data/fast5/tb1.fast5',
            'tests/data/fastq/ecoli_mapped.fastq.gz',
            'tests/data/bam/tb.bam'
        ]
        result = fast5seek.extract_read_ids(files, mapped)
        expected = {
            '57d4cd63-3189-4006-93ec-bf3c8bfb2ce1',
            'bbd563e9-1bf8-4268-92d5-45ccb8e3da72',
            '8aecf428-af00-4791-b065-5d4abd798a29',
            'd707ff64-6ade-477a-8b68-0b3c394ef3b1',
            'c967d421-3da4-4e11-accf-0d3fb0155840',
            '28acaa47-1cec-4a91-9abf-c780d27e6cc4',
            '6cf511b6-1724-46bd-b5a4-59c18bb57343',
            '6c26d9b5-d892-4fc6-b035-abe575895c88'
        }
        self.assertSetEqual(result, expected)


class TestFastqRunIdExtraction(unittest.TestCase):
    """Test the read id extracxtion functions for fastq"""

    def test_FastqRundIdExtractionTB_EmptySet(self):
        """Test run id extracxtion from non-albacore fastq."""
        fastq = 'tests/data/fastq/tb_mapped.fastq'
        result = fast5seek.get_fastq_run_ids([fastq])
        expected = set()
        self.assertSetEqual(result, expected)

    def test_FastqRundIdExtractionAlbacore_TwoRunIds(self):
        """Test run id extracxtion from non-albacore fastq."""
        fastq = 'tests/data/fastq/basecalled.fastq.gz'
        result = fast5seek.get_fastq_run_ids([fastq])
        expected = {
            'bfa81348704ecd62c348b404e974a37daf030951',
            'dc6ee09815f8baff16d92e7189e3a46d855f02b4'
        }
        self.assertSetEqual(result, expected)


class TestGetFast5ReadAndRunId(unittest.TestCase):
    """Test read id extraction from fast5 file.
    TODO: Find a fast5 with more than one read id to test handling
    """
    def test_TBFast5File_OneReadAndRunID(self):
        filepath = 'tests/data/fast5/tb1.fast5'
        result = fast5seek.get_read_and_run_id(filepath)
        expected = ('d707ff64-6ade-477a-8b68-0b3c394ef3b1',
                    'dc6ee09815f8baff16d92e7189e3a46d855f02b4')
        self.assertTupleEqual(result, expected)

    def test_EcoliFast5File_OneReadAndRunID(self):
        filepath = 'tests/data/fast5/ecoli1.fast5'
        result = fast5seek.get_read_and_run_id(filepath)
        expected = ('6cf511b6-1724-46bd-b5a4-59c18bb57343',
                    'bfa81348704ecd62c348b404e974a37daf030951')
        self.assertTupleEqual(result, expected)

    def test_EmptyFast5File_EmptyTuple(self):
        filepath = 'tests/data/fast5/empty.fast5'
        result = fast5seek.get_read_and_run_id(filepath)
        expected = ('', '')
        self.assertTupleEqual(result, expected)


class TestCollectAllFast5Filepaths(unittest.TestCase):
    """Test function that collects all the unique fast5 filepaths."""
    def test_OneFast5Directory_EightFilepaths(self):
        fast5_dir = ['tests/data/fast5']
        result = fast5seek.collect_all_fast5_filepaths(fast5_dir)
        expected = {
            'tests/data/fast5/ecoli1.fast5',
            'tests/data/fast5/ecoli2.fast5',
            'tests/data/fast5/empty.fast5',
            'tests/data/fast5/tb1.fast5',
            'tests/data/fast5/tb2.fast5',
            'tests/data/fast5/tb3.fast5',
            'tests/data/fast5/tb4.fast5',
            'tests/data/fast5/tb5.fast5',
            'tests/data/fast5/tb6.fast5',
        }
        self.assertSetEqual(result, expected)

    def test_TwoFast5Directory_EightFilepaths(self):
        fast5_dir = ['tests/data/fast5', 'tests/data']
        result = fast5seek.collect_all_fast5_filepaths(fast5_dir)
        expected = {
            'tests/data/fast5/ecoli1.fast5',
            'tests/data/fast5/ecoli2.fast5',
            'tests/data/fast5/empty.fast5',
            'tests/data/fast5/tb1.fast5',
            'tests/data/fast5/tb2.fast5',
            'tests/data/fast5/tb3.fast5',
            'tests/data/fast5/tb4.fast5',
            'tests/data/fast5/tb5.fast5',
            'tests/data/fast5/tb6.fast5',
        }
        self.assertSetEqual(result, expected)


class TestCollectPresentFast5Filepaths(unittest.TestCase):
    """Test function that collects all the present fast5 filepaths."""
    def test_OneFast5FileOneReadIdOneRunId_Present(self):
        filepath = {'tests/data/fast5/tb1.fast5'}
        read_ids = {'d707ff64-6ade-477a-8b68-0b3c394ef3b1'}
        run_ids = {'dc6ee09815f8baff16d92e7189e3a46d855f02b4'}
        result = fast5seek.collect_present_fast5_filepaths(filepath, read_ids,
                                                           run_ids, None)
        expected = ['tests/data/fast5/tb1.fast5']
        self.assertListEqual(result, expected)

    def test_OneFast5FileOneReadIdNoRunId_Present(self):
        filepath = {'tests/data/fast5/tb1.fast5'}
        read_ids = {'d707ff64-6ade-477a-8b68-0b3c394ef3b1'}
        run_ids = set()
        result = fast5seek.collect_present_fast5_filepaths(filepath, read_ids,
                                                           run_ids, None)
        expected = ['tests/data/fast5/tb1.fast5']
        self.assertListEqual(result, expected)

    def test_OneFast5FileNoReadIdNoRunId_EmptyList(self):
        filepath = {'tests/data/fast5/tb1.fast5'}
        read_ids = set()
        run_ids = set()
        result = fast5seek.collect_present_fast5_filepaths(filepath, read_ids,
                                                           run_ids, None)
        expected = []
        self.assertListEqual(result, expected)

    def test_TwoFast5FileOneReadIdNoRunId_OneFast5Path(self):
        filepath = {'tests/data/fast5/tb1.fast5', 'tests/data/fast5/tb3.fast5'}
        read_ids = {'d707ff64-6ade-477a-8b68-0b3c394ef3b1'}
        run_ids = set()
        result = fast5seek.collect_present_fast5_filepaths(filepath, read_ids,
                                                           run_ids, None)
        expected = ['tests/data/fast5/tb1.fast5']
        self.assertListEqual(result, expected)

    def test_TwoFast5FileOneReadIdNoRunId_NoFast5Path(self):
        filepath = {'tests/data/fast5/tb1.fast5', 'tests/data/fast5/tb3.fast5'}
        read_ids = {'not a real read id'}
        run_ids = set()
        result = fast5seek.collect_present_fast5_filepaths(filepath, read_ids,
                                                           run_ids, None)
        expected = []
        self.assertListEqual(result, expected)

    def test_TwoFast5FileOneReadIdNoRunId_OneFast5PathNoError(self):
        filepath = {'tests/data/fast5/tb1.fast5',
                    'tests/data/fast5/empty.fast5'}
        read_ids = {'d707ff64-6ade-477a-8b68-0b3c394ef3b1'}
        run_ids = set()
        result = fast5seek.collect_present_fast5_filepaths(filepath, read_ids,
                                                           run_ids, None)
        expected = ['tests/data/fast5/tb1.fast5']
        self.assertListEqual(result, expected)

    def test_TwoFast5FileTwoReadIdOneRunId_TwoFast5Paths(self):
        filepath = {'tests/data/fast5/tb1.fast5',
                    'tests/data/fast5/ecoli1.fast5'}
        read_ids = {'d707ff64-6ade-477a-8b68-0b3c394ef3b1',
                    '6cf511b6-1724-46bd-b5a4-59c18bb57343'}
        run_ids = set('bfa81348704ecd62c348b404e974a37daf030951')
        result = fast5seek.collect_present_fast5_filepaths(filepath, read_ids,
                                                           run_ids, None)
        expected = ['tests/data/fast5/ecoli1.fast5',
                    'tests/data/fast5/tb1.fast5']
        self.assertCountEqual(result, expected)


class TestArgs(object):
    """Class to create args object like that from argparse."""
    def __init__(self, fast5_dir, reference, mapped, output=None):
        self.output = output
        self.fast5_dir = fast5_dir
        self.reference = reference
        self.mapped = mapped
        self.log_level = 4
        self.no_progress_bar = True


class TestMain(unittest.TestCase):
    """Test the main function of the program."""
    def test_TBTestFast5DataFastq_SixFilepaths(self):
        fast5_dir = ['tests/data/fast5']
        reference = ['tests/data/fastq/tb_mapped.fastq.gz']
        args = TestArgs(fast5_dir, reference, True)

        f = io.StringIO()
        with redirect_stdout(f):
            fast5seek.main(args)
        stdout = f.getvalue()
        result = [x for x in stdout.split('\n') if x]
        expected = ['tests/data/fast5/tb5.fast5', 'tests/data/fast5/tb6.fast5',
                    'tests/data/fast5/tb4.fast5', 'tests/data/fast5/tb3.fast5',
                    'tests/data/fast5/tb2.fast5', 'tests/data/fast5/tb1.fast5']
        self.assertCountEqual(result, expected)

    def test_TBTestFast5DataBam_SixFilepaths(self):
        fast5_dir = ['tests/data/fast5']
        reference = ['tests/data/bam/tb.bam']
        args = TestArgs(fast5_dir, reference, True)

        f = io.StringIO()
        with redirect_stdout(f):
            fast5seek.main(args)
        stdout = f.getvalue()
        result = [x for x in stdout.split('\n') if x]
        expected = ['tests/data/fast5/tb5.fast5', 'tests/data/fast5/tb6.fast5',
                    'tests/data/fast5/tb4.fast5', 'tests/data/fast5/tb3.fast5',
                    'tests/data/fast5/tb2.fast5', 'tests/data/fast5/tb1.fast5']
        self.assertCountEqual(result, expected)

    def test_TBTestFast5DataSam_SixFilepaths(self):
        fast5_dir = ['tests/data/fast5']
        reference = ['tests/data/sam/tb.sam']
        args = TestArgs(fast5_dir, reference, True)

        f = io.StringIO()
        with redirect_stdout(f):
            fast5seek.main(args)
        stdout = f.getvalue()
        result = [x for x in stdout.split('\n') if x]
        expected = ['tests/data/fast5/tb5.fast5', 'tests/data/fast5/tb6.fast5',
                    'tests/data/fast5/tb4.fast5', 'tests/data/fast5/tb3.fast5',
                    'tests/data/fast5/tb2.fast5', 'tests/data/fast5/tb1.fast5']
        self.assertCountEqual(result, expected)

    def test_TBTestFast5DataSamAndFastq_SixFilepaths(self):
        fast5_dir = ['tests/data/fast5']
        reference = ['tests/data/sam/tb.sam',
                     'tests/data/fastq/tb_mapped.fastq.gz']
        args = TestArgs(fast5_dir, reference, True)

        f = io.StringIO()
        with redirect_stdout(f):
            fast5seek.main(args)
        stdout = f.getvalue()
        result = [x for x in stdout.split('\n') if x]
        expected = ['tests/data/fast5/tb5.fast5', 'tests/data/fast5/tb6.fast5',
                    'tests/data/fast5/tb4.fast5', 'tests/data/fast5/tb3.fast5',
                    'tests/data/fast5/tb2.fast5', 'tests/data/fast5/tb1.fast5']
        self.assertCountEqual(result, expected)

    def test_EcoliTestFast5DataSam_TwoFilepaths(self):
        fast5_dir = ['tests/data/fast5']
        reference = ['tests/data/sam/ecoli.sam']
        args = TestArgs(fast5_dir, reference, True)

        f = io.StringIO()
        with redirect_stdout(f):
            fast5seek.main(args)
        stdout = f.getvalue()
        result = [x for x in stdout.split('\n') if x]
        expected = ['tests/data/fast5/ecoli1.fast5',
                    'tests/data/fast5/ecoli2.fast5']
        self.assertCountEqual(result, expected)

    def test_TBTestFast5DataBamUnmapped_EightFilepaths(self):
        fast5_dir = ['tests/data/fast5']
        reference = ['tests/data/bam/tb.bam']
        args = TestArgs(fast5_dir, reference, False)

        f = io.StringIO()
        with redirect_stdout(f):
            fast5seek.main(args)
        stdout = f.getvalue()
        result = [x for x in stdout.split('\n') if x]
        expected = ['tests/data/fast5/ecoli1.fast5',
                    'tests/data/fast5/ecoli2.fast5',
                    'tests/data/fast5/tb5.fast5',
                    'tests/data/fast5/tb6.fast5',
                    'tests/data/fast5/tb4.fast5',
                    'tests/data/fast5/tb3.fast5',
                    'tests/data/fast5/tb2.fast5',
                    'tests/data/fast5/tb1.fast5']
        self.assertCountEqual(result, expected)

    def test_TBTestFast5DataAlbacoreFastq_EightFilepaths(self):
        fast5_dir = ['tests/data/fast5']
        reference = ['tests/data/fastq/basecalled.fastq.gz']
        args = TestArgs(fast5_dir, reference, False)

        f = io.StringIO()
        with redirect_stdout(f):
            fast5seek.main(args)
        stdout = f.getvalue()
        result = [x for x in stdout.split('\n') if x]
        expected = ['tests/data/fast5/ecoli1.fast5',
                    'tests/data/fast5/ecoli2.fast5',
                    'tests/data/fast5/tb5.fast5',
                    'tests/data/fast5/tb6.fast5',
                    'tests/data/fast5/tb4.fast5',
                    'tests/data/fast5/tb3.fast5',
                    'tests/data/fast5/tb2.fast5',
                    'tests/data/fast5/tb1.fast5']
        self.assertCountEqual(result, expected)
