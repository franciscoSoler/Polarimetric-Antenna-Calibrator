__author__ = 'fsoler'

import unittest
import glob
import os

import src.Model.Antenna as Antenna
import src.Utilities.Antenna_Common as AntennaCommon
import src.Controllers.Antenna_Creator as RFDNCreator
import numpy as np


class TestAntenna(unittest.TestCase):
    def setUp(self):
        self.filename = "test"
        self.antenna = Antenna.Antenna()
        self.__create_antenna(1, 2, 1)

    def test_get_qtty_antennas(self):
        self.assertEqual(self.antenna.get_qtty_antennas(), 2)

    def test_antenna_raises_exception_if_a_pol_mode_is_incorrect(self):
        self.assertRaises(Exception, self.antenna.get_gain_paths, "TxT")

    def test_the_antenna_retrieve_the_proper_row_col_from_index(self):
        self.assertEqual(self.antenna.index_to_row_col(0), (0, 0))
        self.assertEqual(self.antenna.index_to_row_col(1), (0, 1))

    def test_the_antenna_retrieve_the_proper_index_from_row_col(self):
        self.assertEqual(self.antenna.row_col_to_index(0, 0), 0)
        self.assertEqual(self.antenna.row_col_to_index(0, 1), 1)

    def test_the_antenna_retrieve_the_gain_paths_correctly(self):
        self.assertEqual(len(self.antenna.get_gain_paths("TxV-RxH")), 2)
        self.assertEqual(len(self.antenna.get_gain_paths("TxH-RxV")), 2)
        self.assertEqual(len(self.antenna.get_gain_paths("TxV")), 1)
        self.assertEqual(len(self.antenna.get_gain_paths("TxH")), 1)
        self.assertEqual(len(self.antenna.get_gain_paths("RxV")), 1)
        self.assertEqual(len(self.antenna.get_gain_paths("RxH")), 1)
        print(self.antenna.get_gain_paths("TxH-RxV"))
        # TODO: now i must check if the equations are well performed for every polarization and mode

    def test_the_antenna_retrieve_the_correct_mutual_coupling(self):
        self.assertEqual(len(self.antenna.get_mutual_coupling_front_panel()), 2)

    def test_the_cross_gain_pol_doesnt_change_when_a_gain_change_is_performed(self):
        rx_h = self.antenna.get_gain_paths("RxH")
        rx_v = self.antenna.get_gain_paths("RxV")
        tx_v = self.antenna.get_gain_paths("TxV")
        att_shifts = np.matrix([[1, 2]])

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

        cable1 = [AntennaCommon.Cable, [att, wavelenght, length1]]
        psc = ["{0}1{1}".format(AntennaCommon.Psc, psc_out_ports), [psc_out_ports]]
        cable2 = [AntennaCommon.Cable, [att, wavelenght, length2]]
        trm = [AntennaCommon.Trm, [trm_gain, trm_ph_shift]]
        circulator = [AntennaCommon.Circulator, []]
        cable3 = [AntennaCommon.Cable, [att, wavelenght, length3]]
        rm = [AntennaCommon.Rm, []]

        sequence_items = [cable1, psc, cable2, trm, circulator, cable3, rm]
        creator = RFDNCreator.AntennaCreator(quantity_columns, separation, separation)
        creator.create_structure(self.filename, sequence_items, row_steering, column_steering)
        self.antenna.initialize(separation, separation, self.filename)

    def tearDown(self):
        for filename in glob.glob(self.filename + "_*"):
            os.remove(filename)


if __name__ == '__main__':
    unittest.main()
