__author__ = 'fsoler'

import src.Model.Antenna as Antenna
import src.Utilities.Antenna_Common as AntennaCommon
import src.Controllers.Matrix_Calibrator_builder as MatrixBuilder
import numpy as np
from abc import ABCMeta
import random


class AntennaCalibrator(object):
    __metaclass__ = ABCMeta
    Available_errors = []
    errors = []

    def __init__(self, input_power, input_phase, dist_rows, dist_columns, filename="antenna"):
        self._antenna = Antenna.Antenna()
        self._antenna.initialize(dist_rows, dist_columns, filename)

        self._input_power = input_power
        self._input_phase = input_phase

        self._input_delta_phase = 0
        self._input_delta_power = 0

    def add_calibration_errors(self, errors):
        if not isinstance(errors, list) or len(errors) == 0 or [True for error in errors if len(error) != 2]:
            raise Exception('errors are not well created')
        self._input_delta_power = [error[1] for error in errors if error[0] == AntennaCommon.Inter_pulse_power_err][0]
        self._input_delta_phase = [error[1] for error in errors if error[0] == AntennaCommon.Inter_pulse_phase_err][0]

    def add_errors(self, errors):
        [self.errors.append(error) for error in errors if error not in self.errors]

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


class ClassicCalibrator(AntennaCalibrator):

    def __init__(self, input_power, input_phase, dist_rows, dist_columns, filename):
        super(ClassicCalibrator, self).__init__(input_power, input_phase, dist_rows, dist_columns, filename)

    def __obtain_walsh_matrix(self):
        pass

class MutualCalibrator(AntennaCalibrator):
    Available_calibration_modes = ("TxH-RxV", "TxV-RxH")
    Available_power_modes = ("TxH", "TxV", "RxH", "RxV")

    def __init__(self, input_power, input_phase, dist_rows, dist_columns, filename):
        super(MutualCalibrator, self).__init__(input_power, input_phase, dist_rows, dist_columns, filename)

        self.__matrix_builder = MatrixBuilder.LinearBuilder()
        cross = MatrixBuilder.CrossBuilder()
        double = MatrixBuilder.DoubleBuilder()
        tiny = MatrixBuilder.TinyBuilder()
        tiny.set_successor(MatrixBuilder.DefaultBuilder())
        double.set_successor(tiny)
        cross.set_successor(double)
        self.__matrix_builder.set_successor(cross)

        self.__pol_mode = None

        self.__rm_coupling = self._antenna.get_mutual_coupling_front_panel()
        self.__tx_network = None
        self.__rx_network = None

        self.__power_calculated = False
        self.__tx_power = None
        self.__rx_power = None

        self.__tx_phase = None
        self.__rx_phase = None

        self.__last_cal_paths = None

        self.__equations = None

    def __fix_rx_power_and_phase(self):
        """
        this method is a patch, a biiiiig patch, how it works:
            measured_gain - tx_network_cal - rx_network_cal = delta_rx_cal, then rx_cal = rx_network_cal + delta_rx_cal
            measured_phase - tx_phase - rx_phase_cal = delta_rx_ph, then rx_ph_cal = rx_phase_cal + delta_rx_ph
        :return:
        """
        power_shift = AntennaCommon.v2db(abs(self.__equations[(0, 0)])) - self.__tx_power.item(0, 0) - \
            self.__rx_power.item(0, 0)
        phase_shift = np.angle(self.__equations[(0, 0)], deg=True) - self.__tx_phase.item(0, 0) - \
            self.__rx_phase.item(0, 0)

        self.__rx_power += power_shift
        self.__rx_phase += phase_shift

    def generate_cal_paths(self, strategy, pol_mode="TxH-RxV"):
        self.__power_calculated = False
        self.__last_cal_paths = [strategy, pol_mode]

        if not [pol_mode for cal_mode in self.Available_calibration_modes if pol_mode == cal_mode]:
            raise Exception("the pol_mode is not valid: ", pol_mode)

        self.__pol_mode = pol_mode
        [self.__tx_network, self.__rx_network] = self._antenna.get_gain_paths(pol_mode)
        (last_row, last_col) = max(self.__rm_coupling.keys())
        rows = range(last_row + 1)
        columns = range(last_col + 1)
        f = lambda x: [[AntennaCommon.s2t_parameters(x[row][col]) for col in columns] for row in rows]

        tx_network = f(self.__tx_network)
        rx_network = f(self.__rx_network)

        [b, a] = strategy(self._antenna, tx_network, self.__rm_coupling, rx_network)

        rand = lambda x, y: x + random.uniform(-y, y)
        f = lambda x: x * AntennaCommon.pol2rec(AntennaCommon.db2v(rand(self._input_power, self._input_delta_power)),
                                                rand(self._input_phase, self._input_delta_phase))
        # f = lambda x: x * AntennaCommon.pol2rec(input_voltage + random.uniform(0, self._input_delta_power), self._input_phase)
        self.__equations = dict([a[i], f(b[i].item(1, 0))] for i in range(len(a)))

        self.__matrix_builder.initialize_matrix_builder(self._antenna, self.__equations)
        self.__matrix_builder.build_matrix()
        return self.__equations

    def __obtain_tx_rx_power(self):
        """
        This method build tx_power and rx_power. each power value is in db plus deg format.
        :return:
        """
        self.__power_calculated = True
        format_phase = lambda x: (x + 180) % 360 - 180
        least_squares = lambda a_mx, b: np.matrix(np.linalg.lstsq(a_mx, b)[0]).reshape(self._antenna.quantity_rows,
                                                                                       self._antenna.quantity_columns)

        a, tx_gain, tx_phase = self.__matrix_builder.get_tx_matrix()
        self.__tx_power = least_squares(a, tx_gain)
        self.__tx_phase = least_squares(a, format_phase(tx_phase))

        a, rx_gain, rx_phase = self.__matrix_builder.get_rx_matrix()
        self.__rx_power = least_squares(a, rx_gain)
        self.__rx_phase = least_squares(a, format_phase(rx_phase))

        self.__fix_rx_power_and_phase()

    def calibrate_antenna(self, desired_tx_power, desired_tx_phase, desired_rx_power, desired_rx_phase):
        """

        :param desired_tx_power: desired tx power in dB
        :param desired_tx_phase: desired tx phase in degrees
        :param desired_rx_power: desired rx power in dB
        :param desired_rx_phase: desired rx phase in degrees
        :return:
        """
        if not self.__power_calculated:
            self.__obtain_tx_rx_power()

        # pol2rec = lambda mod, ang: np.multiply(mod, np.exp(1j*AntennaCommon.deg2rad(ang)))
        tx_shift = AntennaCommon.pol2rec(AntennaCommon.db2v(desired_tx_power - self.__tx_power),
                                         desired_tx_phase - self.__tx_phase)
        rx_shift = AntennaCommon.pol2rec(AntennaCommon.db2v(desired_rx_power - self.__rx_power),
                                         desired_rx_phase - self.__rx_phase)

        f = lambda x: [list(map(lambda z: AntennaCommon.v2db(abs(z.item(1, 0))), y)) for y in x]
        g = lambda x: [list(map(lambda z: np.angle(z.item(1, 0), deg=True), y)) for y in x]

        # print("Tx_power non-cal inside of calibrate_antenna", f(self.__antenna.get_gain_paths("TxH")[0]))
        print("Tx_phase non-cal inside of calibrate_antenna", g(self._antenna.get_gain_paths("TxH")[0]))
        # print("rx_power non-cal inside of calibrate_antenna", f(self.__antenna.get_gain_paths("RxV")[0]))
        print("rx_phase non-cal inside of calibrate_antenna", g(self._antenna.get_gain_paths("RxV")[0]))

        modes = AntennaCommon.parse_polarization_mode(self.__pol_mode)
        self._antenna.change_trm_tx_params(tx_shift, modes[0][1])
        self._antenna.change_trm_rx_params(rx_shift, modes[1][1])

        # print("Tx_power cal inside of calibrate_antenna", f(self.__antenna.get_gain_paths("TxH")[0]))
        print("Tx_phase cal inside of calibrate_antenna", g(self._antenna.get_gain_paths("TxH")[0]))
        # print("rx_power cal inside of calibrate_antenna", f(self.__antenna.get_gain_paths("RxV")[0]))
        print("rx_phase cal inside of calibrate_antenna", g(self._antenna.get_gain_paths("RxV")[0]))
        self.generate_cal_paths(*self.__last_cal_paths)

    def get_reception_power(self):
        if not self.__power_calculated:
            self.__obtain_tx_rx_power()
        return self.__rx_power.tolist(), self.__rx_phase.tolist()

    def get_transmission_power(self):
        if not self.__power_calculated:
            self.__obtain_tx_rx_power()
        return self.__tx_power.tolist(), self.__tx_phase.tolist()

    def get_tx_signal_shifts(self, desired_power):
        if not self.__power_calculated:
            self.__obtain_tx_rx_power()
        return (abs(desired_power) - self.__tx_power).tolist(), (np.angle(desired_power) - self.__tx_phase).tolist()

    def get_rx_signal_shifts(self, desired_power):
        if not self.__power_calculated:
            self.__obtain_tx_rx_power()
        return (abs(desired_power) - self.__rx_power).tolist(), (np.angle(desired_power) - self.__rx_phase).tolist()


def every_one_to_one_path_strategy(antenna, tx_network, rm_coupling, rx_network):
    """
    This strategy calculates every path that unify every RM. The configuration is use one in transmission and one in
    reception mode at time.
    :return:
    """

    (last_row, last_col) = max(rm_coupling.keys())
    rows = range(last_row + 1)
    columns = range(last_col + 1)

    s2t = lambda x: AntennaCommon.s2t_parameters(x)
    t2s = lambda x: AntennaCommon.t2s_parameters(x)

    g = lambda tx_row, tx_col, rx_row, rx_col: t2s(tx_network[tx_row][tx_col] *
                                                   s2t(rm_coupling[tx_row, tx_col][rx_row, rx_col]) *
                                                   rx_network[rx_row][rx_col])
    out = [[g(tx_row, tx_col, rx_row, rx_col),
            (antenna.row_col_to_index(tx_row, tx_col),
             antenna.row_col_to_index(rx_row, rx_col))] for rx_row in rows for rx_col in columns
           for tx_col in columns for tx_row in rows]
    unique_path = [[t2s(tx_network[0][0]), (antenna.row_col_to_index(0, 0), None)]]
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

    s2t = lambda x: AntennaCommon.s2t_parameters(x)
    t2s = lambda x: AntennaCommon.t2s_parameters(x)
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

    s2t = lambda x: AntennaCommon.s2t_parameters(x)
    t2s = lambda x: AntennaCommon.t2s_parameters(x)
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
