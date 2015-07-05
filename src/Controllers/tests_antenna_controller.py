__author__ = 'francisco'

import unittest
import glob
import os

import numpy as np

import src.Controllers.Antenna_Calibrator as AntennaCalibrator
import src.Controllers.RFDNCreator as RFDNCreator


class MyTestCase(unittest.TestCase):
    def setUp(self):
        self.filename = "test"
        self.power = 20
        separation = 1
        self.__create_antenna(2, 2, separation)
        self.calibrator = AntennaCalibrator.MutualCalibrator(self.power, separation, separation, self.filename)

    def test_the_antenna_retrieve_the_correct_cal_paths(self):
        print(self.calibrator.generate_cal_paths(AntennaCalibrator.every_one_to_one_path_strategy))
        # todo correct this strategies
        # print(self.calibrator.generate_cal_paths(AntennaCalibrator.strategy_2))
        # print(self.calibrator.generate_cal_paths(AntennaCalibrator.strategy_3))

    def test_antenna_retrieve_the_correct_powers(self):
        self.calibrator.generate_cal_paths(AntennaCalibrator.every_one_to_one_path_strategy)
        np.testing.assert_almost_equal(self.calibrator.get_transmission_power(), [[20, 20], [20, 20]])
        np.testing.assert_almost_equal(self.calibrator.get_reception_power(), [[0, 0], [0, 0]])

    def test_signal_shifts_are_correct(self):
        self.calibrator.generate_cal_paths(AntennaCalibrator.every_one_to_one_path_strategy)
        np.testing.assert_almost_equal(self.calibrator.get_rx_signal_shifts(0), [[0, 0], [0, 0]])
        np.testing.assert_almost_equal(self.calibrator.get_tx_signal_shifts(20), [[0, 0], [0, 0]])

    def __create_antenna(self, quantity_rows, quantity_columns, separation):
        rms = quantity_columns * quantity_rows
        sequence_items = ["cable", "PSC1{0}".format(rms), "cable", "TRM", "circulator", "cable", "RM"]
        creator = RFDNCreator.AntennaCreator(quantity_rows, separation, separation)
        creator.create_structure(self.filename, sequence_items)

    def tearDown(self):
        for filename in glob.glob(self.filename + "_*"):
            os.remove(filename)


if __name__ == '__main__':
    unittest.main()
