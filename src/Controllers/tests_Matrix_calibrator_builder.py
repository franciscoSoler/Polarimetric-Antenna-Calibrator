__author__ = 'francisco'

import unittest
import os
import glob

import numpy as np

import src.Model.Antenna as Antenna
import src.Utilities.Antenna_Common as AntennaCommon
import src.Controllers.Antenna_Creator as RFDNCreator
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

        a, tx_gain, tx_phase = calibrator.get_tx_matrix()
        self.assertEqual(a.size, 0)
        self.assertEqual(tx_gain.size, 0)
        self.assertEqual(tx_phase.size, 0)

        a, rx_gain, rx_phase = calibrator.get_rx_matrix()
        self.assertEqual(a.size, 0)
        self.assertEqual(rx_gain.size, 0)
        self.assertEqual(rx_phase.size, 0)

    def test_cross_equations_are_not_executable(self):
        calibrator = self.__initialize_calibrator(MatCalBuilder.CrossBuilder(), self.equations)

        a, tx_gain, tx_phase = calibrator.get_tx_matrix()
        self.assertEqual(a.size, 0)
        self.assertEqual(tx_gain.size, 0)
        self.assertEqual(tx_phase.size, 0)

        a, rx_gain, rx_phase = calibrator.get_rx_matrix()
        self.assertEqual(a.size, 0)
        self.assertEqual(rx_gain.size, 0)
        self.assertEqual(rx_phase.size, 0)

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

        a, tx_gain, tx_phase = calibrator.get_tx_matrix()
        np.testing.assert_equal(a, [[1, 0, -1]])
        np.testing.assert_equal(tx_gain, AntennaCommon.v2db(18) - AntennaCommon.v2db(19))
        np.testing.assert_equal(tx_phase, [0])

        a, rx_gain, rx_phase = calibrator.get_rx_matrix()
        np.testing.assert_equal(a, [[1, 0, -1]])
        np.testing.assert_equal(rx_gain, AntennaCommon.v2db(19) - AntennaCommon.v2db(18))
        np.testing.assert_equal(rx_phase, [0])

    def test_the_cross_equations_are_correct(self):
        f = lambda x: AntennaCommon.v2db(x)
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

        a, tx_gain, tx_phase = calibrator.get_tx_matrix()
        np.testing.assert_equal(a, [[0, 1, -1, 0], [-1, 0, 0, 1], [1, 0, 0, -1], [0, -1, 1, 0]])
        np.testing.assert_equal(tx_gain, [f(19) - f(18), f(19) - f(18), f(19) - f(18), f(19) - f(18)])
        np.testing.assert_equal(tx_phase, [0, 0, 0, 0])

        a, rx_gain, rx_phase = calibrator.get_tx_matrix()
        np.testing.assert_equal(rx_gain, [f(19) - f(18), f(19) - f(18), f(19) - f(18), f(19) - f(18)])
        np.testing.assert_equal(rx_phase, [0, 0, 0, 0])

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
        a, tx_gain, tx_phase = calibrator.get_tx_matrix()
        np.testing.assert_equal(a, [[2, 0, -2, 0], [0, 2, 0, -2], [2, -2, 0, 0], [0, 0, 2, -2]])
        np.testing.assert_equal(tx_gain, [0, 0, 0, 0])
        np.testing.assert_equal(tx_phase, [0, 0, 0, 0])

        a, rx_gain, rx_phase = calibrator.get_tx_matrix()
        np.testing.assert_equal(a, [[2, 0, -2, 0], [0, 2, 0, -2], [2, -2, 0, 0], [0, 0, 2, -2]])
        np.testing.assert_equal(rx_gain, [0, 0, 0, 0])
        np.testing.assert_equal(rx_phase, [0, 0, 0, 0])

    def test_the_default_equations_are_correct(self):
        f = lambda x: AntennaCommon.v2db(x)

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
        a, tx_gain, tx_phase = calibrator.get_tx_matrix()
        np.testing.assert_equal(a, [[1, 0, 0, 0]])
        np.testing.assert_equal(tx_gain, [f(10)])
        np.testing.assert_equal(tx_phase, [0])

        a, rx_gain, rx_phase = calibrator.get_tx_matrix()
        np.testing.assert_equal(a, [[1, 0, 0, 0]])
        np.testing.assert_equal(rx_gain, [f(10)])
        np.testing.assert_equal(rx_phase, [0])

    def __initialize_calibrator(self, calibrator, equations):
        calibrator.initialize_matrix_builder(self.antenna, equations)
        calibrator.build_matrix()
        return calibrator

    def __create_antenna(self, quantity_rows, quantity_columns, separation):
        att = 0.5           # [neper/m]
        c = 299792458       # [m/seg]
        f = 1275000000      # [Hz]
        wavelenght = c/f    # [m]

        length1 = 0.45      # [m]
        length2 = 8         # [m]
        length3 = 0.5       # [m]

        trm_gain = 0       # []
        trm_ph_shift = 0   # [deg]

        psc_out_ports = quantity_columns * quantity_rows

        cable1 = [AntennaCommon.Cable, [att, wavelenght, length1]]
        psc = ["{0}1{1}".format(AntennaCommon.Psc, psc_out_ports), [psc_out_ports]]
        cable2 = [AntennaCommon.Cable, [att, wavelenght, length2]]
        trm = [AntennaCommon.Trm, [trm_gain, trm_ph_shift]]
        circulator = [AntennaCommon.Circulator, []]
        cable3 = [AntennaCommon.Cable, [att, wavelenght, length3]]
        rm = [AntennaCommon.Rm, []]

        sequence_items = [cable1, psc, cable2, trm, circulator, cable3, rm]
        creator = RFDNCreator.AntennaCreator(quantity_rows, separation, separation)
        creator.create_structure(self.filename, sequence_items)
        self.antenna = Antenna.Antenna()
        self.antenna.initialize(separation, separation, self.filename)

    def tearDown(self):
        for filename in glob.glob(self.filename + "_*"):
            os.remove(filename)
            pass


if __name__ == '__main__':
    unittest.main()
