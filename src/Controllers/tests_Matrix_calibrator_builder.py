__author__ = 'francisco'

import unittest
import os
import glob

import numpy as np

import src.Model.Antenna as Antenna
import src.Controllers.RFDNCreator as RFDNCreator
import src.Controllers.Matrix_Calibrator_builder as MatCalBuilder


class MyTestCase(unittest.TestCase):
    def setUp(self):
        self.filename = "test"
        self.separation = 1
        self.__quantity_rows = 1
        self.__quantity_columns = 2
        self.__create_antenna(self.__quantity_rows, self.__quantity_columns, self.separation)
        self.equations = dict([
            ((0, 0), 20), ((0, 1), 19),
            ((1, 0), 19), ((1, 1), 20),
            ((0, None), 10)])

    def test_linear_equations_are_not_executable(self):
        calibrator = self.__initialize_calibrator(MatCalBuilder.LinearBuilder(), self.equations)

        [a, b_trans] = calibrator.get_tx_matrix()
        self.assertEqual(a.size, 0)
        self.assertEqual(b_trans.size, 0)

        [a, b_rec] = calibrator.get_rx_matrix()
        self.assertEqual(a.size, 0)
        self.assertEqual(b_rec.size, 0)

    def test_cross_equations_are_not_executable(self):
        calibrator = self.__initialize_calibrator(MatCalBuilder.CrossBuilder(), self.equations)

        [a, b_trans] = calibrator.get_tx_matrix()
        self.assertEqual(a.size, 0)
        self.assertEqual(b_trans.size, 0)

        [a, b_rec] = calibrator.get_rx_matrix()
        self.assertEqual(a.size, 0)
        self.assertEqual(b_rec.size, 0)

    def test_the_linear_equations_are_correct(self):
        quantity_rows = 1
        quantity_columns = 3
        equations = dict([
            ((0, 0), 20), ((0, 1), 18), ((0, 2), 18),
            ((1, 0), 19), ((1, 1), 20), ((1, 2), 18),
            ((2, 0), 18), ((2, 1), 19), ((2, 2), 20),
            ((0, None), 10)])
        self.__create_antenna(quantity_rows, quantity_columns, self.separation)
        calibrator = self.__initialize_calibrator(MatCalBuilder.LinearBuilder(), equations)

        [a, b_trans] = calibrator.get_tx_matrix()
        self.assertTrue((a == np.matrix([[1, 0, -1]])).all())
        self.assertTrue((b_trans == np.array([-1])).all())

        [a, b_rec] = calibrator.get_rx_matrix()
        self.assertTrue((a == np.matrix([[1, 0, -1]])).all())
        self.assertTrue((b_rec == np.array([1])).all())

    def test_the_cross_equations_are_correct(self):
        quantity_rows = 2
        quantity_columns = 2
        equations = dict([
            ((0, 0), 20), ((0, 1), 19), ((0, 2), 18), ((0, 3), 17),
            ((1, 0), 19), ((1, 1), 20), ((1, 2), 19), ((1, 3), 18),
            ((2, 0), 18), ((2, 1), 19), ((2, 2), 20), ((2, 3), 19),
            ((3, 0), 17), ((3, 1), 18), ((3, 2), 19), ((3, 3), 20),
            ((0, None), 10)])
        self.__create_antenna(quantity_rows, quantity_columns, self.separation)
        calibrator = self.__initialize_calibrator(MatCalBuilder.CrossBuilder(), equations)
        [a, b_trans] = calibrator.get_tx_matrix()
        self.assertTrue((a == np.matrix([[0, 1, -1, 0], [-1, 0, 0, 1], [1, 0, 0, -1], [0, -1, 1, 0]])).all())
        self.assertTrue((b_trans == np.array([1, 1, 1, 1])).all())

        [a, b_rec] = calibrator.get_tx_matrix()
        self.assertTrue((a == np.matrix([[0, 1, -1, 0], [-1, 0, 0, 1], [1, 0, 0, -1], [0, -1, 1, 0]])).all())
        self.assertTrue((b_rec == np.array([1, 1, 1, 1])).all())

    def test_the_tiny_equations_are_correct(self):
        quantity_rows = 2
        quantity_columns = 2
        equations = dict([
            ((0, 0), 20), ((0, 1), 19), ((0, 2), 18), ((0, 3), 17),
            ((1, 0), 19), ((1, 1), 20), ((1, 2), 19), ((1, 3), 18),
            ((2, 0), 18), ((2, 1), 19), ((2, 2), 20), ((2, 3), 19),
            ((3, 0), 17), ((3, 1), 18), ((3, 2), 19), ((3, 3), 20),
            ((0, None), 10)])
        self.__create_antenna(quantity_rows, quantity_columns, self.separation)
        calibrator = self.__initialize_calibrator(MatCalBuilder.TinyBuilder(), equations)
        [a, b_trans] = calibrator.get_tx_matrix()
        self.assertTrue((a == np.matrix([[2, 0, -2, 0], [0, 2, 0, -2], [2, -2, 0, 0], [0, 0, 2, -2]])).all())
        self.assertTrue((b_trans == np.array([0, 0, 0, 0])).all())

        [a, b_rec] = calibrator.get_tx_matrix()
        self.assertTrue((a == np.matrix([[2, 0, -2, 0], [0, 2, 0, -2], [2, -2, 0, 0], [0, 0, 2, -2]])).all())
        self.assertTrue((b_rec == np.array([0, 0, 0, 0])).all())

    def test_the_default_equations_are_correct(self):
        quantity_rows = 2
        quantity_columns = 2
        equations = dict([
            ((0, 0), 20), ((0, 1), 19), ((0, 2), 18), ((0, 3), 17),
            ((1, 0), 19), ((1, 1), 20), ((1, 2), 19), ((1, 3), 18),
            ((2, 0), 18), ((2, 1), 19), ((2, 2), 20), ((2, 3), 19),
            ((3, 0), 17), ((3, 1), 18), ((3, 2), 19), ((3, 3), 20),
            ((0, None), 10)])
        self.__create_antenna(quantity_rows, quantity_columns, self.separation)
        calibrator = self.__initialize_calibrator(MatCalBuilder.DefaultBuilder(), equations)
        [a, b_trans] = calibrator.get_tx_matrix()
        self.assertTrue((a == np.matrix([[1, 0, 0, 0]])).all())
        self.assertTrue((b_trans == np.array([10])).all())

        [a, b_rec] = calibrator.get_tx_matrix()
        self.assertTrue((a == np.matrix([[1, 0, 0, 0]])).all())
        self.assertTrue((b_rec == np.array([10])).all())

    def __initialize_calibrator(self, calibrator, equations):
        calibrator.initialize_matrix_builder(self.antenna, equations)
        calibrator.build_matrix()
        return calibrator

    def __create_antenna(self, quantity_rows, quantity_columns, separation):
        rms = quantity_columns * quantity_rows
        sequence_items = ["cable", "PSC1{0}".format(rms), "cable", "TRM", "circulator", "cable", "RM"]
        creator = RFDNCreator.AntennaCreator(quantity_rows, separation, separation)
        creator.create_structure(self.filename, sequence_items, 0.)
        self.antenna = Antenna.Antenna()
        self.antenna.initialize(separation, separation, self.filename)

    def tearDown(self):
        for filename in glob.glob(self.filename + "_*"):
            os.remove(filename)
            pass


if __name__ == '__main__':
    unittest.main()
