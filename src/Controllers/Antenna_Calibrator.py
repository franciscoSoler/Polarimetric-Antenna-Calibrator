__author__ = 'fsoler'

import Model.Antenna as Antenna
import Utilities.Antenna_Common as common

import numpy as np
from abc import ABCMeta, abstractmethod


class AntennaCalibrator(object):
    __metaclass__ = ABCMeta
    _Available_calibration_modes = ("TxH-RxV", "TxV-RxH")

    def __init__(self, input_power, input_phase, dist_rows, dist_columns,
                 filename="antenna"):
        self._antenna = Antenna.Antenna()
        self._antenna.initialize(dist_rows, dist_columns, filename)

        self._input_power = input_power
        self._input_phase = input_phase

        self._input_delta_phase = 0
        self._input_delta_power = 0

        self._pol_mode = self._Available_calibration_modes[0]

        self._power_calculated = False
        self._tx_power = None
        self._rx_power = None

        self._tx_phase = None
        self._rx_phase = None

    def _set_pol_mode(self, pol_mode):
        if not [pol_mode for cal_mode in self._Available_calibration_modes if pol_mode == cal_mode]:
            raise Exception("the pol_mode is not valid: ", pol_mode)
        self._power_calculated = False
        self._pol_mode = pol_mode

    def _calibrate_antenna(self, desired_tx_power, desired_tx_phase, desired_rx_power, desired_rx_phase):
        """
        :param desired_tx_power: desired tx power in dB
        :param desired_tx_phase: desired tx phase in degrees
        :param desired_rx_power: desired rx power in dB
        :param desired_rx_phase: desired rx phase in degrees
        :return:
        """
        if not self._power_calculated:
            self._obtain_tx_rx_power()

        tx_shift = common.pol2rec(common.db2v(desired_tx_power - self._tx_power),
                                         desired_tx_phase - self._tx_phase)
        rx_shift = common.pol2rec(common.db2v(desired_rx_power - self._rx_power),
                                         desired_rx_phase - self._rx_phase)

        modes = common.parse_polarization_mode(self._pol_mode)
        self._antenna.change_trm_tx_params(tx_shift, modes[0][1])
        self._antenna.change_trm_rx_params(rx_shift, modes[1][1])

    @abstractmethod
    def _add_calibration_errors(self, errors):
        pass

    @abstractmethod
    def _obtain_tx_rx_power(self):
        pass

    @abstractmethod
    def calibrate_antenna(self, desired_tx_power, desired_tx_phase, desired_rx_power, desired_rx_phase):
        pass

    def add_calibration_errors(self, errors):
        if not isinstance(errors, list) or len(errors) == 0 or [True for error in errors if len(error) != 2]:
            raise Exception('errors are not well created')
        if [True for error in errors if error[0] == common.Inter_pulse_gain_err]:
            self._input_delta_power, self._input_delta_phase = [err[1] for err in errors
                                                                if err[0] == common.Inter_pulse_gain_err].pop()
        self._add_calibration_errors(errors)

    def get_antenna_gain_paths(self, complete_calibration=True):
        f = lambda x: [list(map(lambda z: common.v2db(abs(z.item(1, 0))), y)) for y in x]
        g = lambda x: [list(map(lambda z: np.angle(z.item(1, 0), deg=True), y)) for y in x]
        tx_signal, rx_signal = self._antenna.get_gain_paths(self._pol_mode, complete_calibration)
        return f(tx_signal), g(tx_signal), f(rx_signal), g(rx_signal)

    def get_reception_power(self):
        if not self._power_calculated:
            self._obtain_tx_rx_power()
        return self._rx_power.tolist(), self._rx_phase.tolist()

    def get_transmission_power(self):
        if not self._power_calculated:
            self._obtain_tx_rx_power()
        return self._tx_power.tolist(), self._tx_phase.tolist()

    def get_tx_signal_shifts(self, desired_power):
        if not self._power_calculated:
            self._obtain_tx_rx_power()
        return (abs(desired_power) - self._tx_power).tolist(), (np.angle(desired_power) - self._tx_phase).tolist()

    def get_rx_signal_shifts(self, desired_power):
        if not self._power_calculated:
            self._obtain_tx_rx_power()
        return (abs(desired_power) - self._rx_power).tolist(), (np.angle(desired_power) - self._rx_phase).tolist()

    @property
    def input_power(self):
        return self._input_power

    @input_power.setter
    def input_power(self, power):
        self._input_power = power

    @property
    def input_phase(self):
        return self._input_phase

    @input_phase.setter
    def input_phase(self, phase):
        self._input_phase = phase


def obtain_submatrix(matrix, ix1, ix2):
    it11 = matrix.item(ix1, ix1)
    it12 = matrix.item(ix1, ix2)
    it21 = matrix.item(ix2, ix1)
    it22 = matrix.item(ix2, ix2)

    return np.matrix([[it11, it12],[it21, it22]])


def every_one_to_one_path_strategy(antenna, tx_network, rm_coupling, rx_network):
    """
    This strategy calculates every path that unify every RM. The configuration is use one in transmission and one in
    reception mode at time.
    :return:
    """
    rows, columns = antenna.shape
    rowrange = range(rows)
    colrange = range(columns)

    s2t = lambda x: common.s2t_parameters(x)
    t2s = lambda x: common.t2s_parameters(x)

    g = lambda tx_row, tx_col, rx_row, rx_col: t2s(tx_network[tx_row][tx_col] *
                                                   s2t(obtain_submatrix(rm_coupling, 
                                                                        antenna.row_col_to_index(tx_row, tx_col),
                                                                        antenna.row_col_to_index(rx_row, rx_col))) *
                                                   rx_network[rx_row][rx_col])
    out = [[g(tx_row, tx_col, rx_row, rx_col),
            (antenna.row_col_to_index(tx_row, tx_col),
             antenna.row_col_to_index(rx_row, rx_col))] for rx_row in rowrange for rx_col in colrange
           for tx_col in colrange for tx_row in rowrange 
           if antenna.row_col_to_index(tx_row, tx_col) != antenna.row_col_to_index(rx_row, rx_col)]

    half_row = rows // 2
    half_col= columns // 2
    unique_path = [[t2s(tx_network[half_row][half_col]), (antenna.row_col_to_index(half_row, half_col), None)],
                   [t2s(rx_network[half_row][half_col]), (None, antenna.row_col_to_index(half_row, half_col))]]

    return list(zip(*(out + unique_path)))


def strategy_2(antenna, tx_network, rm_coupling, rx_network):
    """
    I don't know what strategy I intend to perform, this makes the same as the first one
    :param antenna:
    :param tx_network:
    :param rm_coupling:
    :param rx_network:
    :return:
    """
    (last_row, last_col) = max(rm_coupling.keys())
    rows = range(last_row + 1)
    columns = range(last_col + 1)

    s2t = lambda x: common.s2t_parameters(x)
    t2s = lambda x: common.t2s_parameters(x)
    g = lambda tx_row, tx_col, rx_row, rx_col: t2s(tx_network[tx_row][tx_col] *
                                                   s2t(rm_coupling[tx_row, tx_col][rx_row, rx_col]) *
                                                   rx_network[rx_row][rx_col])
    out = [[g(tx_row, tx_col, rx_row, rx_col),
            (antenna.row_col_to_index(tx_row, tx_col),
             antenna.get_index_from_rm_separation(abs(tx_row - rx_row), abs(tx_col - rx_col)),
             antenna.row_col_to_index(rx_row, rx_col))] for rx_row in rows for rx_col in columns
           for tx_col in columns for tx_row in rows]

    return list(zip(*out))


def strategy_3(antenna, tx_network, rm_coupling, rx_network):
    """
    calculates every faced RM and every vertical/horizontal reception neighbour
    :param antenna:
    :param tx_network:
    :param rm_coupling:
    :param rx_network:
    :return:
    """
    (last_row, last_col) = max(rm_coupling.keys())
    rows = range(last_row + 1)
    columns = range(last_col + 1)

    s2t = lambda x: common.s2t_parameters(x)
    t2s = lambda x: common.t2s_parameters(x)
    g = lambda tx_row, tx_col, rx_row, rx_col: t2s(tx_network[tx_row][tx_col] *
                                                   s2t(rm_coupling[tx_row, tx_col][rx_row, rx_col]) *
                                                   rx_network[rx_row][rx_col])
    h = lambda x, y: x + y
    row_fix = lambda row: row + -2*int((row + 1)/len(rm_coupling)) + 1 if len(rm_coupling) > 1 else row
    col_fix = lambda col: col + -2*int((col + 1)/len(rm_coupling[0])) + 1 if len(rm_coupling[0]) > 1 else col
    out_double = [[h(g(tx_row, tx_col, row_fix(tx_row), tx_col), g(tx_row, tx_col, tx_row, col_fix(tx_col))),
                   (antenna.row_col_to_index(tx_row, tx_col),
                    (antenna.get_index_from_rm_separation(abs(row_fix(tx_row) - tx_row), 0),
                     antenna.get_index_from_rm_separation(0, abs(col_fix(tx_col) - tx_col))),
                    (antenna.row_col_to_index(row_fix(tx_row), tx_col),
                     antenna.row_col_to_index(tx_row, col_fix(tx_col))))] for tx_col in columns for tx_row in rows]

    out_simple = [[g(tx_row, tx_col, row_fix(tx_row), tx_col),
                   (antenna.row_col_to_index(tx_row, tx_col),
                    antenna.get_index_from_rm_separation(abs(row_fix(tx_row) - tx_row), 0),
                    antenna.row_col_to_index(row_fix(tx_row), tx_col))] for tx_col in columns for tx_row in rows]

    return list(zip(*(out_simple + out_double)))
