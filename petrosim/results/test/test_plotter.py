import csv
from io import StringIO
import os
import unittest

from petrosim.models.ecafc import plotter


EXAMPLE_CSVFILE = os.path.join(os.path.dirname(__file__), 'test.csv')


class TestTopLevelFuncs(unittest.TestCase):

    def test_combine_header_rows(self):
        test_csv_one_line_header = """Tm,Ta\n1.0,2.0"""
        f_single = StringIO(test_csv_one_line_header)
        reader = csv.reader(f_single)
        header_single = plotter.combine_header_rows(reader, 1)
        assert header_single[1] == 'Ta'

        test_csv_two_line_header = """,Sr\nTm,Cm\n1.0,2.0"""
        f_double = StringIO(test_csv_two_line_header)
        reader = csv.reader(f_double)
        header_double = plotter.combine_header_rows(reader, 2)
        assert header_double[1] == 'Sr Cm'

    def test_csv_reader(self):
        reader = list(plotter.csv_reader(EXAMPLE_CSVFILE))
        assert reader[0]['Tm,norm'] == 1.0

    def test_get_data(self):
        reader = list(plotter.csv_reader(EXAMPLE_CSVFILE))
        data = plotter.get_data(reader, 'Tm')
        assert len(data) == 175

    def test_plot_cols(self):
        reader = list(plotter.csv_reader(EXAMPLE_CSVFILE))
        plotter.plot_cols(reader, 'Tm', 'Ta', show=False)
        expected_file = 'Tm-Ta.png'
        assert os.path.exists(expected_file)
        if os.path.exists(expected_file):
            os.remove(expected_file)

    def test_pretty_label(self):
        """
        This test is less useful because the function only affects the display
        of the string in matplotlib.
        """

        string = '143Nd/144Nd em'
        expected_str = r'$^1$$^4$$^3$Nd/$^1$$^4$$^4$Nd $ε_m$'
        assert plotter.pretty_label(string) == expected_str

    def test_isotopize(self):
        """
        This test is less useful because the function only affects the display
        of the string in matplotlib.
        """

        string = '143Nd'
        expected_str = r'$^1$$^4$$^3$Nd'
        assert plotter.isotopize(string) == expected_str

    def test_construct_filename(self):
        filename = plotter.construct_filename('Tm', 'Sr Cm')
        expected_name = 'Tm-Sr_Cm'
        assert filename == expected_name