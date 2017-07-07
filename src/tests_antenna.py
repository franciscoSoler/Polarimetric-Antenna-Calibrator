__author__ = 'fsoler'

import unittest
import glob
import os

import Model.Antenna as Antenna
import Utilities.Antenna_Common as common
import Controllers.Antenna_Creator as RFDNCreator
import numpy as np


class TestAntenna(unittest.TestCase):
    def setUp(self):
        self.filename = "test"
        self.antenna = Antenna.Antenna()
        self.quantity_rows = 3
        self.quantity_columns = 3
        self.separation = 1
        self.__create_antenna(self.quantity_rows, self.quantity_columns, self.separation)

    def test_get_qtty_antennas(self):
        self.assertEqual(self.antenna.get_qtty_antennas(), self.quantity_rows * self.quantity_columns)

    def test_antenna_raises_exception_if_a_pol_mode_is_incorrect(self):
        self.assertRaises(Exception, self.antenna.get_gain_paths, "TxT")

    def test_the_antenna_retrieve_the_proper_row_col_from_index(self):
        self.assertEqual(self.antenna.index_to_row_col(0), (0, 0))
        self.assertEqual(self.antenna.index_to_row_col(1), (0, 1))
        self.assertEqual(self.antenna.index_to_row_col(8), (2, 2))

    def test_the_antenna_retrieve_the_proper_index_from_row_col(self):
        self.assertEqual(self.antenna.row_col_to_index(0, 0), 0)
        self.assertEqual(self.antenna.row_col_to_index(0, 1), 1)
        self.assertEqual(self.antenna.row_col_to_index(2, 2), 8)

    def test_the_antenna_retrieve_the_gain_paths_correctly(self):
        cable1 = np.matrix([["0", "(-1.01356527725-0.144343237369j)"], ["(-1.01356527725-0.144343237369j)", "0"]]).astype(complex)
        psc = np.matrix([["0.0", "0.31622776601683794"], ["0.31622776601683794", "0.0"]]).astype(complex)
        cable2 = np.matrix([["0", "(1.00638378932+0.187992579698j)"], ["(1.00638378932+0.187992579698j)", "0"]]).astype(complex)
        trm = np.matrix([["0", "0"],["(9.84807753012+1.73648177667j)", "0"]]).astype(complex)
        circulator = np.matrix([["0", "0"],["1", "0"]]).astype(complex)
        cable3 = np.matrix([["0", "(-1.00597283446+0.190179383088j)"], ["(-1.00597283446+0.190179383088j)", "0"]]).astype(complex)

        tx_transmission = common.t2s_parameters(self.__s2t(cable1) * self.__s2t(psc) * self.__s2t(cable2) * self.__s2t(trm) * self.__s2t(circulator) * self.__s2t(cable3))

        self.assertEqual(len(self.antenna.get_gain_paths("TxV-RxH")), 2)
        self.assertEqual(len(self.antenna.get_gain_paths("TxH-RxV")), 2)
        self.assertEqual(len(self.antenna.get_gain_paths("TxV")), 1)
        self.assertEqual(len(self.antenna.get_gain_paths("TxH")), 1)
        self.assertEqual(len(self.antenna.get_gain_paths("RxV")), 1)
        self.assertEqual(len(self.antenna.get_gain_paths("RxH")), 1)

        np.testing.assert_almost_equal(self.antenna.get_gain_paths("TxV")[0], [[tx_transmission]*self.quantity_columns]*self.quantity_rows)

    def test_the_antenna_retrieve_the_correct_mutual_coupling(self):
        base_coupling = np.array([
        [
            "0j",
            "(-0.02042955017362778+1.1049820775197547j)",
            "(-1.2205680251195763-0.04514857358729857j)",
            "(-0.02042955017362778+1.1049820775197547j)",
            "(1.1470873443358467+0.10529513573702565j)",
            "(-1.2481762086656558-0.07748850273686755j)",
            "(-1.2205680251195763-0.04514857358729857j)",
            "(-1.2481762086656558-0.07748850273686755j)",
            "(1.3047223099255865+0.24156543524813448j)"
        ],
        [
            "(-0.02042955017362778+1.1049820775197547j)",
            "0j",
            "(-0.02042955017362778+1.1049820775197547j)",
            "(1.1470873443358467+0.10529513573702565j)",
            "(-0.02042955017362778+1.1049820775197547j)",
            "(1.1470873443358467+0.10529513573702565j)",
            "(-1.2481762086656558-0.07748850273686755j)",
            "(-1.2205680251195763-0.04514857358729857j)",
            "(-1.2481762086656558-0.07748850273686755j)"
        ],
        [
            "(-1.2205680251195763-0.04514857358729857j)",
            "(-0.02042955017362778+1.1049820775197547j)",
            "0j",
            "(-1.2481762086656558-0.07748850273686755j)",
            "(1.1470873443358467+0.10529513573702565j)",
            "(-0.02042955017362778+1.1049820775197547j)",
            "(1.3047223099255865+0.24156543524813448j)",
            "(-1.2481762086656558-0.07748850273686755j)",
            "(-1.2205680251195763-0.04514857358729857j)"
        ],
        [
            "(-0.02042955017362778+1.1049820775197547j)",
            "(1.1470873443358467+0.10529513573702565j)",
            "(-1.2481762086656558-0.07748850273686755j)",
            "0j",
            "(-0.02042955017362778+1.1049820775197547j)",
            "(-1.2205680251195763-0.04514857358729857j)",
            "(-0.02042955017362778+1.1049820775197547j)",
            "(1.1470873443358467+0.10529513573702565j)",
            "(-1.2481762086656558-0.07748850273686755j)"
        ],
        [
            "(1.1470873443358467+0.10529513573702565j)",
            "(-0.02042955017362778+1.1049820775197547j)",
            "(1.1470873443358467+0.10529513573702565j)",
            "(-0.02042955017362778+1.1049820775197547j)",
            "0j",
            "(-0.02042955017362778+1.1049820775197547j)",
            "(1.1470873443358467+0.10529513573702565j)",
            "(-0.02042955017362778+1.1049820775197547j)",
            "(1.1470873443358467+0.10529513573702565j)"
        ],
        [
            "(-1.2481762086656558-0.07748850273686755j)",
            "(1.1470873443358467+0.10529513573702565j)",
            "(-0.02042955017362778+1.1049820775197547j)",
            "(-1.2205680251195763-0.04514857358729857j)",
            "(-0.02042955017362778+1.1049820775197547j)",
            "0j",
            "(-1.2481762086656558-0.07748850273686755j)",
            "(1.1470873443358467+0.10529513573702565j)",
            "(-0.02042955017362778+1.1049820775197547j)"
        ],
        [
            "(-1.2205680251195763-0.04514857358729857j)",
            "(-1.2481762086656558-0.07748850273686755j)",
            "(1.3047223099255865+0.24156543524813448j)",
            "(-0.02042955017362778+1.1049820775197547j)",
            "(1.1470873443358467+0.10529513573702565j)",
            "(-1.2481762086656558-0.07748850273686755j)",
            "0j",
            "(-0.02042955017362778+1.1049820775197547j)",
            "(-1.2205680251195763-0.04514857358729857j)"
        ],
        [
            "(-1.2481762086656558-0.07748850273686755j)",
            "(-1.2205680251195763-0.04514857358729857j)",
            "(-1.2481762086656558-0.07748850273686755j)",
            "(1.1470873443358467+0.10529513573702565j)",
            "(-0.02042955017362778+1.1049820775197547j)",
            "(1.1470873443358467+0.10529513573702565j)",
            "(-0.02042955017362778+1.1049820775197547j)",
            "0j",
            "(-0.02042955017362778+1.1049820775197547j)"
        ],
        [
            "(1.3047223099255865+0.24156543524813448j)",
            "(-1.2481762086656558-0.07748850273686755j)",
            "(-1.2205680251195763-0.04514857358729857j)",
            "(-1.2481762086656558-0.07748850273686755j)",
            "(1.1470873443358467+0.10529513573702565j)",
            "(-0.02042955017362778+1.1049820775197547j)",
            "(-1.2205680251195763-0.04514857358729857j)",
            "(-0.02042955017362778+1.1049820775197547j)",
            "0j"
        ]]).astype(complex)
        sha = (self.quantity_rows * self.quantity_columns, self.quantity_rows * self.quantity_columns)
        self.assertEqual(self.antenna.get_mutual_coupling_front_panel().shape, base_coupling.shape)
        np.testing.assert_array_equal(self.antenna.get_mutual_coupling_front_panel(), base_coupling)

    def test_gain_change(self):
        rx_h = self.antenna.get_gain_paths("RxH")
        rx_v = self.antenna.get_gain_paths("RxV")
        tx_h = self.antenna.get_gain_paths("TxH")
        tx_v = self.antenna.get_gain_paths("TxV")
        shift = 2
        att_shifts = np.matrix([[shift]* self.quantity_columns] * self.quantity_rows)

        self.antenna.change_trm_tx_params(att_shifts, "vPolarization")
        self.antenna.change_trm_rx_params(att_shifts, "hPolarization")

        new_tx_v = self.antenna.get_gain_paths("TxV")
        new_rx_h = self.antenna.get_gain_paths("RxH")

        np.testing.assert_equal(new_tx_v[0][:, :, 1, 0], tx_v[0][:, :, 1, 0] * shift)
        np.testing.assert_equal(new_rx_h[0][:, :, 1, 0], rx_h[0][:, :, 1, 0] * shift)

        shift = 4
        att_shifts = np.matrix([[shift]* self.quantity_columns] * self.quantity_rows)

        self.antenna.change_trm_tx_params(att_shifts, "hPolarization")
        self.antenna.change_trm_rx_params(att_shifts, "vPolarization")

        new_tx_h = self.antenna.get_gain_paths("TxH")
        new_rx_v = self.antenna.get_gain_paths("RxV")

        np.testing.assert_equal(new_tx_h[0][:, :, 1, 0], tx_h[0][:, :, 1, 0] * shift)
        np.testing.assert_equal(new_rx_v[0][:, :, 1, 0], rx_v[0][:, :, 1, 0] * shift)

    def test_the_cross_gain_pol_doesnt_change_when_a_gain_change_is_performed(self):
        rx_h = self.antenna.get_gain_paths("RxH")
        rx_v = self.antenna.get_gain_paths("RxV")
        tx_v = self.antenna.get_gain_paths("TxV")
        att_shifts = np.matrix([[1]* self.quantity_columns] * self.quantity_rows)

        self.antenna.change_trm_tx_params(att_shifts, "hPolarization")

        rx_h_post_cal = self.antenna.get_gain_paths("RxH")
        rx_v_post_cal = self.antenna.get_gain_paths("RxV")
        tx_v_post_cal = self.antenna.get_gain_paths("TxV")

        np.testing.assert_equal(rx_h, rx_h_post_cal)
        np.testing.assert_equal(rx_v, rx_v_post_cal)
        np.testing.assert_equal(tx_v, tx_v_post_cal)

        rx_h_post_cal = self.antenna.get_gain_paths("RxH")
        tx_h_post_cal = self.antenna.get_gain_paths("TxH")
        tx_v_post_cal = self.antenna.get_gain_paths("TxV")

        self.antenna.change_trm_rx_params(att_shifts, "vPolarization")

        rx_h_post_second_cal = self.antenna.get_gain_paths("RxH")
        tx_h_post_second_cal = self.antenna.get_gain_paths("TxH")
        tx_v_post_second_cal = self.antenna.get_gain_paths("TxV")

        np.testing.assert_equal(rx_h_post_cal, rx_h_post_second_cal)
        np.testing.assert_equal(tx_h_post_cal, tx_h_post_second_cal)
        np.testing.assert_equal(tx_v_post_cal, tx_v_post_second_cal)

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

        row_steering = 0
        column_steering = 0

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
        creator.create_structure(self.filename, sequence_items, row_steering, column_steering)
        self.antenna.initialize(separation, separation, self.filename)

    def __s2t(self, param):
        return common.s2t_parameters(param)

    def tearDown(self):
        pass
        # for filename in glob.glob(self.filename + "_*"):
        #     os.remove(filename)


if __name__ == '__main__':
    unittest.main()
