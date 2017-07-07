__author__ = 'francisco'

import unittest
import glob
import os
import sys

import numpy as np

import Controllers.Mutual_Calibrator as AntennaCalibrator
import Controllers.Antenna_Calibrator as AntennaCal
import Model.Antenna as Antenna
import Utilities.Antenna_Common as common
import Controllers.Antenna_Creator as RFDNCreator


class MyTestCase(unittest.TestCase):
    def setUp(self):
        self.filename = "test"
        self.power = 20
        self.phase = 0
        self.separation = 1
        self.row_steering = 30
        self.column_steering = 0
        self.rows = 3
        self.cols = 2
        self.__create_antenna(self.rows, self.cols, self.separation)
        self.calibrator = AntennaCalibrator.MutualCalibrator(self.power, self.phase, self.column_steering, self.row_steering, self.separation, self.separation,
                                                             self.filename)

    def test_the_antenna_retrieve_the_correct_cal_paths(self):
        antenna = Antenna.Antenna()
        antenna.initialize(self.separation, self.separation, self.filename)

        [tx_network, rx_network] = antenna.get_gain_paths('TxH-RxV')
        coupling = antenna.get_mutual_coupling_front_panel()

        cal = self.calibrator.generate_cal_paths(AntennaCal.every_one_to_one_path_strategy)
        # sys.stdout.write('\ncal paths {}\n'.format(list(cal.keys())))
        # sys.stdout.write('\nvalues {}\n'.format(list(cal.values())))
        # sys.stdout.write('\nphase paths {}\n'.format(np.angle(list(cal.values()), deg=True)))

        for path in cal.keys():
            if path[0] is None:
                rx_index = antenna.index_to_row_col(path[1])
                np.testing.assert_almost_equal(cal[path], rx_network[rx_index[0], rx_index[1], 1, 0] * common.db2v(self.power))
            elif path[1] is None:
                tx_index = antenna.index_to_row_col(path[0])
                np.testing.assert_almost_equal(cal[path], tx_network[tx_index[0], tx_index[1], 1, 0] * common.db2v(self.power))
            else:
                tx_index = antenna.index_to_row_col(path[0])
                rx_index = antenna.index_to_row_col(path[1])
                tx_t_matrix = common.s2t_parameters(tx_network[tx_index[0], tx_index[1]])
                rx_t_matrix = common.s2t_parameters(rx_network[rx_index[0], rx_index[1]])
                coupling_matrix = common.s2t_parameters(AntennaCal.obtain_submatrix(coupling, *path))
                np.testing.assert_almost_equal(cal[path], common.t2s_parameters(tx_t_matrix * coupling_matrix * rx_t_matrix)[1, 0] * common.db2v(self.power))

    def test_antenna_calibrator_retrieve_the_correct_gain(self):
        antenna = Antenna.Antenna()
        antenna.initialize(self.separation, self.separation, self.filename)

        [tx_network, rx_network] = antenna.get_gain_paths('TxH-RxV')

        self.calibrator.generate_cal_paths(AntennaCal.every_one_to_one_path_strategy)
        tx_power, tx_phase = self.calibrator.get_transmission_power()

        np.testing.assert_almost_equal(tx_power, np.abs(tx_network[:, :, 1, 0]))
        np.testing.assert_almost_equal(tx_phase, np.angle(tx_network[:, :, 1, 0], deg=True))

        rx_power, rx_phase = self.calibrator.get_reception_power()
        
        np.testing.assert_almost_equal(rx_power, np.abs(rx_network[:, :, 1, 0]))
        np.testing.assert_almost_equal(rx_phase, np.angle(rx_network[:, :, 1, 0], deg=True))
        
    def test_signal_shifts_are_correct(self):
        antenna = Antenna.Antenna()
        antenna.initialize(self.separation, self.separation, self.filename)
        [tx_network, rx_network] = antenna.get_gain_paths('TxH-RxV')

        self.calibrator.generate_cal_paths(AntennaCal.every_one_to_one_path_strategy)

        desired_power = 0
        desired_shift_power, desired_shift_phase = self.__get_desired_shifts(desired_power, rx_network)

        rx_power_shift, rx_phase_shift = self.calibrator.get_rx_signal_shifts(desired_power)
        np.testing.assert_almost_equal(rx_power_shift, desired_shift_power)
        np.testing.assert_almost_equal(rx_phase_shift, desired_shift_phase)

        desired_power = 20
        desired_shift_power, desired_shift_phase = self.__get_desired_shifts(desired_power, tx_network)

        tx_power_shift, tx_phase_shift = self.calibrator.get_tx_signal_shifts(desired_power)
        np.testing.assert_almost_equal(tx_power_shift, desired_shift_power)
        np.testing.assert_almost_equal(tx_phase_shift, desired_shift_phase)

    def test_every_one_to_one_strategy(self):
        antenna = Antenna.Antenna()
        antenna.initialize(self.separation, self.separation, self.filename)

        [tx_network, rx_network] = antenna.get_gain_paths('TxH-RxV')
        coupling = antenna.get_mutual_coupling_front_panel()

        gains, paths = AntennaCal.every_one_to_one_path_strategy(antenna, tx_network, coupling, rx_network)

        for i, path in enumerate(paths):

            if path[0] is None:
                rx_index = antenna.index_to_row_col(path[1])
                np.testing.assert_almost_equal(gains[i], rx_network[rx_index[0], rx_index[1]])
            elif path[1] is None:
                tx_index = antenna.index_to_row_col(path[0])
                np.testing.assert_almost_equal(gains[i], tx_network[tx_index[0], tx_index[1]])
            else:
                tx_index = antenna.index_to_row_col(path[0])
                rx_index = antenna.index_to_row_col(path[1])
                tx_t_matrix = common.s2t_parameters(tx_network[tx_index[0], tx_index[1]])
                rx_t_matrix = common.s2t_parameters(rx_network[rx_index[0], rx_index[1]])
                coupling_matrix = common.s2t_parameters(AntennaCal.obtain_submatrix(coupling, *path))
                np.testing.assert_almost_equal(gains[i], common.t2s_parameters(tx_t_matrix * coupling_matrix * rx_t_matrix))

    def __get_desired_shifts(self, desired_gain, network):
        desired_shift_power = np.abs(desired_gain) - np.abs(network[:, :, 1, 0])
        desired_shift_phase = np.angle(desired_gain, deg=True) - np.angle(network[:, :, 1, 0], deg=True)
        return desired_shift_power, desired_shift_phase

    def __s2t(self, param):
        return common.s2t_parameters(param)

    def __gain_paths_to_db(self, paths):
        return np.array([list(map(lambda z: common.v2db(abs(z.item(1, 0))), y)) for y in paths])

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

        # row_steering = 0
        # column_steering = 0

        psc_out_ports = quantity_columns * quantity_rows

        cable1 = [common.Cable, [att, wavelenght, length1]]
        psc = ["{0}1{1}".format(common.Psc, psc_out_ports), [psc_out_ports]]
        cable2 = [common.Cable, [att, wavelenght, length2]]
        trm = [common.Trm, [trm_gain, trm_ph_shift]]
        circulator = [common.Circulator, []]
        cable3 = [common.Cable, [att, wavelenght, length3]]
        rm = [common.Rm, []]

        sequence_items = [cable1, psc, cable2, trm, circulator, cable3, rm]
        creator = RFDNCreator.AntennaCreator(quantity_rows, quantity_columns, separation, separation)

        creator.create_structure(self.filename, sequence_items, self.row_steering, self.column_steering)

    def tearDown(self):
        pass
        # for filename in glob.glob(self.filename + "_*"):
        #     os.remove(filename)


if __name__ == '__main__':
    unittest.main()
