__author__ = 'francisco'

import unittest
import glob
import os

import numpy as np

import src.Controllers.Antenna_Calibrator as AntennaCalibrator
import src.Model.Antenna as Antenna
import src.Utilities.Antenna_Common as AntennaCommon
import src.Controllers.RFDNCreator as RFDNCreator


class MyTestCase(unittest.TestCase):
    def setUp(self):
        self.filename = "test"
        self.power = 20
        self.separation = 1
        self.__create_antenna(2, 2, self.separation)
        self.calibrator = AntennaCalibrator.AntennaCalibrator(self.power, self.separation, self.separation,
                                                              self.filename)

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

    def test_every_one_to_one_strategy(self):
        with open("../base_rfdn") as f:
            with open(self.filename + "_rfdn", "w") as g:
                g.write(f.read())

        with open("../base_panel") as f:
            with open(self.filename + "_panel", "w") as g:
                g.write(f.read())
        antenna = Antenna.Antenna()
        antenna.initialize(self.separation, self.separation, self.filename)
        self.calibrator = AntennaCalibrator.AntennaCalibrator(self.power, self.separation, self.separation,
                                                              self.filename)
        self.calibrator.generate_cal_paths(AntennaCalibrator.every_one_to_one_path_strategy)
        antenna.get_gain_paths("TxH")
        antenna.get_gain_paths("RxV")
        f = lambda x: 20*np.log10(10*np.array(x))
        g = lambda x: [[mtx.item(1, 0) for mtx in row]for row in antenna.get_gain_paths("TxH")[0]]
        print("antena", f(g(antenna.get_gain_paths("TxH"))))
        print("calibrated", self.calibrator.get_transmission_power())
        np.testing.assert_almost_equal(f(g(antenna.get_gain_paths("TxH"))), self.calibrator.get_transmission_power())
        np.testing.assert_almost_equal(f(g(antenna.get_gain_paths("RxH"))), self.calibrator.get_reception_power())
        """
        tx_network = [[np.matrix([[-0.06039018+0.j,  0.96957050+0.j],
                                  [0.98155780+0.j, -0.05896126-0.j]]),
                       np.matrix([[-0.20441523+0.j,  1.30932119+0.j],
                                  [1.21843198+0.j, -0.47589587-0.j]])],
                      [np.matrix([[-0.06468088+0.j,  1.10330958+0.j],
                                  [1.09883445+0.j,  0.23394620+0.j]]),
                       np.matrix([[-0.04692100+0.j,  1.19473155+0.j],
                                  [1.00275439+0.j,  0.23362583+0.j]])]]


        same_coup = np.matrix([[0., 1.],
                               [1., 0.]])

        neighbour_coup = np.matrix([[0., 3.],
                                    [3., 0.]])

        cross_coup = np.matrix([[0., 3.82842712],
                                [3.82842712, 0.]])

        coupling_network = {(0, 1): {(0, 1): np.matrix([[0., 1.],
                                                        [1., 0.]]),
                                     (1, 0): np.matrix([[0., 3.82842712],
                                                        [3.82842712, 0.]]),
                                     (0, 0): np.matrix([[0., 3.],
                                                        [3., 0.]]),
                                     (1, 1): np.matrix([[0., 3.],
                                                        [3., 0.]])},
                            (1, 0): {(0, 1): np.matrix([[0., 3.82842712],
                                                        [3.82842712, 0.]]),
                                     (1, 0): np.matrix([[0., 1.],
                                                        [1., 0.]]),
                                     (0, 0): np.matrix([[0., 3.],
                                                        [3., 0.]]),
                                     (1, 1): np.matrix([[0., 3.],
                                                        [3., 0.]])},
                            (0, 0): {(0, 1): np.matrix([[0., 3.],
                                                        [3., 0.]]),
                                     (1, 0): np.matrix([[0., 3.],
                                                        [3., 0.]]),
                                     (0, 0): np.matrix([[0., 1.],
                                                        [1., 0.]]),
                                     (1, 1): np.matrix([[0., 3.82842712],
                                                        [3.82842712, 0.]])},
                            (1, 1): {(0, 1): np.matrix([[0., 3.],
                                                        [3., 0.]]),
                                     (1, 0): np.matrix([[0., 3.],
                                                        [3., 0.]]),
                                     (0, 0): np.matrix([[0., 3.82842712],
                                                        [3.82842712, 0.]]),
                                     (1, 1): np.matrix([[0., 1.],
                                                        [1., 0.]])}}

        rx_network = [[np.matrix([[-0.28491053-0.j,  1.28334637+0.j],
                                  [1.22235657+0.j, -0.50988760+0.j]]),
                       np.matrix([[0.03935427+0.j,  1.31926821+0.j],
                                  [1.13068183+0.j, -0.51933755+0.j]])],
                      [np.matrix([[0.13253064+0.j,  0.92211135+0.j],
                                  [0.75821942+0.j,  0.09271972+0.j]]),
                       np.matrix([[-0.13835146-0.j,  1.17647078+0.j],
                                  [1.43692311+0.j, -0.53852356+0.j]])]]
        f = lambda x: AntennaCommon.s2t_parameters(x)
        g = lambda x: AntennaCommon.t2s_parameters(x)
        print(g(f(tx_network[0][0])*f(same_coup)*f(rx_network[0][0])))
        print(g(f(tx_network[0][1])*f(neighbour_coup)*f(rx_network[1][1])))
        print(g(f(tx_network[1][0])*f(neighbour_coup)*f(rx_network[1][1])))
        print(g(f(tx_network[0][0])*f(cross_coup)*f(rx_network[1][1])))


        antenna = Antenna.Antenna()
        antenna.initialize(1, 1, self.filename)
        print("que ida?")
        f = lambda x: [[AntennaCommon.s2t_parameters(x[row][col]) for col in range(2)] for row in range(2)]
        print(AntennaCalibrator.every_one_to_one_path_strategy(antenna, f(tx_network), coupling_network, f(rx_network)))
        """

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
