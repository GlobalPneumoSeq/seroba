import unittest
import os
import filecmp
from seroba import summary

modules_dir = os.path.dirname(os.path.abspath(summary.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestSummary(unittest.TestCase):
    def test__summary(self):
        expected = os.path.join(data_dir,'exp_summary.csv')
        got=summary.Summary.summarise(os.path.join(data_dir,'summ'))
        self.assertTrue(filecmp.cmp('summary.csv',expected), 'files are not equal')
        os.remove('summary.csv')
