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
        self.row_steering = 0
        self.column_steering = 0
        self.rows = 3
        self.cols = 2
        self.tx_transmission = self.__calculate_tx_power()
        self.__create_antenna(self.rows, self.cols, self.separation)
        self.calibrator = AntennaCalibrator.MutualCalibrator(self.power, self.phase, self.column_steering, self.row_steering, self.separation, self.separation,
                                                             self.filename)

    def __calculate_tx_power(self):
        cable1 = np.matrix([["0", "(0.896401138946-0.539109712212j)"], ["(0.896401138946-0.539109712212j)", "0"]]).astype(complex)
        psc = np.matrix([["0.0", "0.31622776601683794"], ["0.31622776601683794", "0.0"]]).astype(complex)
        trm = np.matrix([["0", "0"],["(-4.19918462144-9.07561835442j)", "0"]]).astype(complex)
        circulator = np.matrix([["0", "0"],["1", "0"]]).astype(complex)
        cable2 = np.matrix([["0", "(0.736458202447+0.750200129382j)"], ["(0.736458202447+0.750200129382j)", "0"]]).astype(complex)

        return common.t2s_parameters(self.__s2t(cable1) * self.__s2t(psc) * self.__s2t(trm) * self.__s2t(circulator) * self.__s2t(cable2))
        

    def test_the_antenna_retrieve_the_correct_cal_paths(self):
        cal = self.calibrator.generate_cal_paths(AntennaCal.every_one_to_one_path_strategy)
        sys.stdout.write('\ncal paths {}\n'.format(list(cal.keys())))
        sys.stdout.write('\nvalues {}\n'.format(list(cal.values())))
        sys.stdout.write('\nphase paths {}\n'.format(np.angle(list(cal.values()), deg=True)))
        # todo correct this strategies
        # print(self.calibrator.generate_cal_paths(AntennaCalibrator.strategy_2))
        # print(self.calibrator.generate_cal_paths(AntennaCalibrator.strategy_3))

    def test_antenna_retrieve_the_correct_powers(self):
        self.calibrator.generate_cal_paths(AntennaCal.every_one_to_one_path_strategy)
        tx_power, tx_phase = self.calibrator.get_transmission_power()

        np.testing.assert_almost_equal(tx_power, [[np.abs(self.tx_transmission)[1, 0]]*self.cols]*self.rows)
        np.testing.assert_almost_equal(tx_phase, [[np.angle(self.tx_transmission)[1, 0]]*self.cols]*self.rows)

        rx_power, rx_phase = self.calibrator.get_reception_power()
        
        np.testing.assert_almost_equal(rx_power, [[np.abs(self.tx_transmission)[1, 0]]*self.cols]*self.rows)
        np.testing.assert_almost_equal(rx_phase, [[np.angle(self.tx_transmission)[1, 0]]*self.cols]*self.rows)
        
    def test_signal_shifts_are_correct(self):
        self.calibrator.generate_cal_paths(AntennaCal.every_one_to_one_path_strategy)
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
        antenna = Antenna.Antenna()
        antenna.initialize(self.separation, self.separation, self.filename)

        [tx_network, rx_network] = antenna.get_gain_paths('TxH-RxV')
        coupling = antenna.get_mutual_coupling_front_panel()

        gains, paths = AntennaCal.every_one_to_one_path_strategy(antenna, tx_network, coupling, rx_network)

        for i, path in enumerate(paths):

            if path[0] is None:
                rx_index = antenna.index_to_row_col(path[1])
                np.testing.assert_almost_equal(gains[i], rx_network[rx_index[0]][rx_index[1]])
            elif path[1] is None:
                tx_index = antenna.index_to_row_col(path[0])
                np.testing.assert_almost_equal(gains[i], tx_network[tx_index[0]][tx_index[1]])
            else:
                tx_index = antenna.index_to_row_col(path[0])
                rx_index = antenna.index_to_row_col(path[1])
                tx_t_matrix = common.s2t_parameters(tx_network[tx_index[0]][tx_index[1]])
                rx_t_matrix = common.s2t_parameters(rx_network[rx_index[0]][rx_index[1]])
                coupling_matrix = common.s2t_parameters(AntennaCal.obtain_submatrix(coupling, *path))
                np.testing.assert_almost_equal(gains[i], common.t2s_parameters(tx_t_matrix * coupling_matrix * rx_t_matrix))

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
