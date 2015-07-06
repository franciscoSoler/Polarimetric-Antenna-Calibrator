__author__ = 'francisco'

import unittest
import glob
import os

import numpy as np

import src.Controllers.Antenna_Calibrator as AntennaCalibrator
import src.Model.Antenna as Antenna
import src.Utilities.Antenna_Common as AntennaCommon
import src.Controllers.Antenna_Creator as RFDNCreator


class MyTestCase(unittest.TestCase):
    def setUp(self):
        self.filename = "test"
        self.power = 20
        self.phase = 0
        self.separation = 1
        self.__create_antenna(2, 2, self.separation)
        self.calibrator = AntennaCalibrator.MutualCalibrator(self.power, self.phase, self.separation, self.separation,
                                                             self.filename)

    def test_the_antenna_retrieve_the_correct_cal_paths(self):
        print(self.calibrator.generate_cal_paths(AntennaCalibrator.every_one_to_one_path_strategy))
        # todo correct this strategies
        # print(self.calibrator.generate_cal_paths(AntennaCalibrator.strategy_2))
        # print(self.calibrator.generate_cal_paths(AntennaCalibrator.strategy_3))

    def test_antenna_retrieve_the_correct_powers(self):
        self.calibrator.generate_cal_paths(AntennaCalibrator.every_one_to_one_path_strategy)
        tx_power, tx_phase = self.calibrator.get_transmission_power()
        np.testing.assert_almost_equal(tx_power, [[40.784171182715156, 40.784171182715156],
                                                  [40.784171182715156, 40.784171182715156]])
        np.testing.assert_almost_equal(tx_phase, [[32.979812787692396, 32.979812787692396],
                                                  [32.979812787692396, 32.979812787692396]])

        rx_power, rx_phase = self.calibrator.get_reception_power()
        np.testing.assert_almost_equal(rx_power, [[20.784171182715156, 20.784171182715156],
                                                  [20.784171182715156, 20.784171182715156]])
        np.testing.assert_almost_equal(rx_phase, [[32.979812787692396, 32.979812787692396],
                                                  [32.979812787692396, 32.979812787692396]])

    def test_signal_shifts_are_correct(self):
        self.calibrator.generate_cal_paths(AntennaCalibrator.every_one_to_one_path_strategy)
        rx_power_shift, rx_phase_shift = self.calibrator.get_rx_signal_shifts(0)
        np.testing.assert_almost_equal(rx_power_shift, [[-20.784171182715156, -20.784171182715156],
                                                        [-20.784171182715156, -20.784171182715156]])
        np.testing.assert_almost_equal(rx_phase_shift, [[-32.979812787692396, -32.979812787692396],
                                                        [-32.979812787692396, -32.979812787692396]])

        tx_power_shift, tx_phase_shift = self.calibrator.get_tx_signal_shifts(20)
        np.testing.assert_almost_equal(tx_power_shift, [[-20.784171182715156, -20.784171182715156],
                                                        [-20.784171182715156, -20.784171182715156]])
        np.testing.assert_almost_equal(tx_phase_shift, [[-32.979812787692396, -32.979812787692396],
                                                        [-32.979812787692396, -32.979812787692396]])

    def test_every_one_to_one_strategy(self):
        with open("../base_rfdn") as f:
            with open(self.filename + "_rfdn", "w") as g:
                g.write(f.read())

        with open("../base_panel") as f:
            with open(self.filename + "_panel", "w") as g:
                g.write(f.read())
        antenna = Antenna.Antenna()
        antenna.initialize(self.separation, self.separation, self.filename)
        self.calibrator = AntennaCalibrator.MutualCalibrator(self.power, self.phase, self.separation, self.separation,
                                                             self.filename)
        self.calibrator.generate_cal_paths(AntennaCalibrator.every_one_to_one_path_strategy)

        f = lambda x: np.matrix([list(map(lambda z: AntennaCommon.v2db(abs(z.item(1, 0))), y)) for y in x])
        np.testing.assert_almost_equal(f(antenna.get_gain_paths("TxH")[0]) + self.power,
                                       self.calibrator.get_transmission_power()[0])
        np.testing.assert_almost_equal(f(antenna.get_gain_paths("RxV")[0]), self.calibrator.get_reception_power()[0])

    def __create_antenna(self, quantity_rows, quantity_columns, separation):
        att = 0.1           # [neper/m]
        c = 299792458       # [m/seg]
        f = 1275000000      # [Hz]
        wavelenght = c/f    # [m]

        length1 = 0.45      # [m]
        length2 = 8         # [m]
        length3 = 0.5       # [m]

        trm_gain = 10       # []
        trm_ph_shift = 10   # [deg]

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

    def tearDown(self):
        for filename in glob.glob(self.filename + "_*"):
            os.remove(filename)


if __name__ == '__main__':
    unittest.main()
