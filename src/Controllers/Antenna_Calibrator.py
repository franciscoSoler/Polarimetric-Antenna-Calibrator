__author__ = 'fsoler'

import Model.Antenna as Antenna
import Utilities.Antenna_Common as AntennaCommon
import Controllers.Matrix_Calibrator_builder as MatrixBuilder
import Utilities.Walsh_mtx_creator as WalshCreator
import Utilities.Chirp_creator as ChirpCreator
import numpy as np
from abc import ABCMeta, abstractmethod
import random


class AntennaCalibrator(object):
    __metaclass__ = ABCMeta
    _Available_calibration_modes = ("TxH-RxV", "TxV-RxH")

    def __init__(self, calibrationType, input_power, input_phase, dist_rows, dist_columns,
                 filename="antenna"):
        self.__complete_calibration = False if calibrationType == AntennaCommon.CalibratorType.Classical else True
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

        tx_shift = AntennaCommon.pol2rec(AntennaCommon.db2v(desired_tx_power - self._tx_power),
                                         desired_tx_phase - self._tx_phase)
        rx_shift = AntennaCommon.pol2rec(AntennaCommon.db2v(desired_rx_power - self._rx_power),
                                         desired_rx_phase - self._rx_phase)

        modes = AntennaCommon.parse_polarization_mode(self._pol_mode)
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
        self._input_delta_power, self._input_delta_phase = [err[1] for err in errors if err[0] == AntennaCommon.Inter_pulse_gain_err].pop()
        self._add_calibration_errors(errors)

    def get_antenna_gain_paths(self):
        f = lambda x: [list(map(lambda z: AntennaCommon.v2db(abs(z.item(1, 0))), y)) for y in x]
        g = lambda x: [list(map(lambda z: np.angle(z.item(1, 0), deg=True), y)) for y in x]
        tx_signal, rx_signal = self._antenna.get_gain_paths(self._pol_mode, self.__complete_calibration)
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


class ClassicCalibrator(AntennaCalibrator):
    __Dec = 11
    __Swl = 2*AntennaCommon.tp  # sampling window length [sec]

    def __init__(self, input_power, input_phase, dist_rows, dist_columns, filename):
        super(ClassicCalibrator, self).__init__(AntennaCommon.CalibratorType.Classical, input_power, input_phase, dist_rows,
                                                dist_columns, filename)
        self.__walsh_creator = WalshCreator.WalshMatrixCreator()
        self.__chirp_creator = ChirpCreator.ChirpCreator(self._input_power, self._input_phase)
        self.__quantity_elements = self._antenna.get_qtty_antennas()

    def _add_calibration_errors(self, errors):
        self.__walsh_creator.add_walsh_errors(errors)
        self.__chirp_creator.add_chirp_errors(errors)

    def __obtain_estimated_signal(self, power_signal, phase_signal):
        att_db = np.reshape(power_signal, (-1, 1))
        ph_deg = np.reshape(phase_signal, (-1, 1))
        sequences = 2**np.ceil(np.log2(self.__quantity_elements))  # quantity of sequences (and mode pulses)

        # built of walsh sequence matrix, which is a sequence of phase shift per element (rows) in each pulse (cols).
        # in rad
        walsh_phi_m = self.__walsh_creator.create_ideal_phase_walsh(self.__quantity_elements)
        walsh_phi_m_err = self.__walsh_creator.create_phase_walsh_matrix(self.__quantity_elements)

        # built of ICAL LONG LOOP chirp and ICAL SHORT LOOP chirp
        chirp_parameters = [AntennaCommon.fs, AntennaCommon.fc, AntennaCommon.bw, AntennaCommon.tp, self.__Swl]
        chirp = [self.__chirp_creator.create_chirp(*chirp_parameters) for _ in range(self.__quantity_elements)]
        chirp_rep = np.matrix(self.__chirp_creator.create_chirp_replica(*chirp_parameters))

        amp = AntennaCommon.db2v(-att_db)
        ph_rad = AntennaCommon.deg2rad(ph_deg)

        """
        RAW DATA CODING
        """
        # added phase per loop (real setting + walsh coding, with phase shift errors)
        phi0 = np.tile(ph_rad, sequences) + walsh_phi_m_err[:self.__quantity_elements, :]
        # Built of every loop signal and summed among them
        acq = np.multiply(np.dot(amp.T, np.exp(1j * phi0)).T, chirp)

        """
        RAW DATA DECODING
        """
        signal = np.tile((acq * chirp_rep.H).T, (self.__quantity_elements, 1))
        integral = np.multiply(signal, np.exp(-1j * walsh_phi_m[:self.__quantity_elements, :]))  # integral of every term
        signal_est = np.dot(integral, np.ones(sequences)) / (sequences * chirp_rep * chirp_rep.H)
        """
            the divisor is composed by:
              sequences: is ||cj||². This is not correct, the calculation is the correlation of generation matrices,
                one with errors with the other, without errors.
              chirpRep*chirpRep' is the module of both chirps in which the signal was multiplied.
        """

        estimated_phase = np.around(np.angle(signal_est, deg=True), decimals=self.__Dec)
        estimated_power = np.around(-AntennaCommon.v2db(abs(signal_est)), decimals=self.__Dec)
        f = lambda x: np.resize(x, (self._antenna.quantity_rows, self._antenna.quantity_columns))
        return f(estimated_power), f(estimated_phase)

    def _obtain_tx_rx_power(self):
        self._power_calculated = True
        # tx_signal, rx_signal = self._antenna.get_gain_paths(self._pol_mode)
        tx_power, tx_phase, rx_power, rx_phase = self.get_antenna_gain_paths()
        self._tx_power, self._tx_phase = self.__obtain_estimated_signal(tx_power, tx_phase)
        self._rx_power, self._rx_phase = self.__obtain_estimated_signal(rx_power, rx_phase)

    def set_pol_mode(self, pol_mode):
        self._set_pol_mode(pol_mode)

    def calibrate_antenna(self, desired_tx_power, desired_tx_phase, desired_rx_power, desired_rx_phase):
        self._calibrate_antenna(desired_tx_power, desired_tx_phase, desired_rx_power, desired_rx_phase)
        self._power_calculated = False


class MutualCalibrator(AntennaCalibrator):
    def __init__(self, input_power, input_phase, row_steering, column_steering, dist_rows, dist_columns, filename):
        super(MutualCalibrator, self).__init__(AntennaCommon.CalibratorType.Mutual, input_power, input_phase,
                                               dist_rows, dist_columns, filename)

        f = lambda row, col: (row * row_steering + col * column_steering + 41 + 180) % 360 - 180
        self.__trm_setting = np.array([f(row, col) for row in range(self._antenna.quantity_rows)
                                       for col in range(self._antenna.quantity_columns)])

        self.__matrix_builder = MatrixBuilder.LinearBuilder()
        cross = MatrixBuilder.CrossBuilder()
        double = MatrixBuilder.DoubleBuilder()
        tiny = MatrixBuilder.TinyBuilder()
        tiny.set_successor(MatrixBuilder.DefaultBuilder())
        double.set_successor(tiny)
        cross.set_successor(double)
        self.__matrix_builder.set_successor(cross)

        self.__rm_coupling = self._antenna.get_mutual_coupling_front_panel()
        self.__tx_network = None
        self.__rx_network = None

        self.__last_cal_paths = None

        self.__equations = None

    def _add_calibration_errors(self, errors):
        pass

    def __fix_rx_power_and_phase(self):
        """
        this method is a patch, a biiiiig patch, how it works:
            measured_gain - tx_network_cal - rx_network_cal = delta_rx_cal, then rx_cal = rx_network_cal + delta_rx_cal
            measured_phase - tx_phase - rx_phase_cal = delta_rx_ph, then rx_ph_cal = rx_phase_cal + delta_rx_ph
        :return:
        """
        power_shift = AntennaCommon.v2db(abs(self.__equations[(0, 0)])) - self._tx_power.item(0, 0) - \
            self._rx_power.item(0, 0)
        phase_shift = np.mod(np.angle(self.__equations[(0, 0)], deg=True) - self._tx_phase.item(0, 0) -
                             self._rx_phase.item(0, 0) + 180, 360) - 180
        self._rx_power += power_shift
        self._rx_phase += phase_shift

    def generate_cal_paths(self, strategy, pol_mode="TxH-RxV"):
        self._set_pol_mode(pol_mode)

        self.__last_cal_paths = [strategy, pol_mode]

        [self.__tx_network, self.__rx_network] = self._antenna.get_gain_paths(pol_mode)
        (last_row, last_col) = max(self.__rm_coupling.keys())
        rows = range(last_row + 1)
        columns = range(last_col + 1)
        f = lambda x: [[AntennaCommon.s2t_parameters(x[row][col]) for col in columns] for row in rows]

        tx_network = f(self.__tx_network)
        rx_network = f(self.__rx_network)

        [b, a] = strategy(self._antenna, tx_network, self.__rm_coupling, rx_network)

        rand = lambda x, y: np.random.normal(x, y) if y else x
        f = lambda x: x * AntennaCommon.pol2rec(AntennaCommon.db2v(rand(self._input_power, self._input_delta_power)),
                                                rand(self._input_phase, self._input_delta_phase))
        self.__equations = dict([a[i], f(b[i].item(1, 0))] for i in range(len(a)))

        self.__matrix_builder.initialize_matrix_builder(self._antenna, self.__equations)
        self.__matrix_builder.build_matrix()

        return self.__equations

    def _obtain_tx_rx_power(self):
        """
        This method build tx_power and rx_power. each power value is in db plus deg format.
        :return:
        """
        self._power_calculated = True
        cut = 180
        format_phase = lambda x: (x + cut) % 360 - cut
        # format_phase = lambda x: np.mod(AntennaCommon.deg2rad(x) + np.pi, 2*np.pi) - np.pi
        # format_phase = lambda x: np.array([-np.sign(val)*(360-abs(val)) if 130 < val < 149 or 159 < val < 168 or 205 < val < 260 else val for val in x.tolist()])
        least_squares = lambda a_mx, b: np.matrix(np.linalg.lstsq(a_mx, b)[0]).reshape(self._antenna.quantity_rows,
                                                                                       self._antenna.quantity_columns)

        a, tx_gain, tx_phase = self.__matrix_builder.get_tx_matrix()
        self._tx_power = least_squares(a, tx_gain)
        #self._tx_phase = least_squares(a, self.__correct_phase(a, tx_phase))
        self._tx_phase = format_phase(least_squares(a, self.__correct_phase(a, tx_phase)))
        # self._tx_phase = least_squares(a, self.__correct_phase(a, format_phase(tx_phase)))
        # self._tx_phase = least_squares(a, format_phase(tx_phase))
        # self._tx_phase = format_phase(least_squares(a, format_phase(tx_phase)))
        # self._tx_phase = AntennaCommon.rad2deg(least_squares(a, format_phase(tx_phase)))
        """
        print(a)
        print(tx_phase)
        print(format_phase(tx_phase))
        print(AntennaCommon.rad2deg(least_squares(a, format_phase(tx_phase))))
        """
        a, rx_gain, rx_phase = self.__matrix_builder.get_rx_matrix()
        self._rx_power = least_squares(a, rx_gain)
        #self._rx_phase = least_squares(a, self.__correct_phase(a, rx_phase))
        self._rx_phase = format_phase(least_squares(a, self.__correct_phase(a, rx_phase)))
        # self._rx_phase = least_squares(a, self.__correct_phase(a, format_phase(rx_phase)))
        # self._rx_phase = least_squares(a, format_phase(rx_phase))
        # self._rx_phase = format_phase(least_squares(a, format_phase(rx_phase)))
        # self._rx_phase = AntennaCommon.rad2deg(least_squares(a, format_phase(rx_phase)))

        #self.__fix_rx_power_and_phase()

    def calibrate_antenna(self, desired_tx_power, desired_tx_phase, desired_rx_power, desired_rx_phase):
        self._calibrate_antenna(desired_tx_power, desired_tx_phase, desired_rx_power, desired_rx_phase)
        self.generate_cal_paths(*self.__last_cal_paths)

    def __change_values(self, phase, indexes, increment):
        for idx in indexes:
            phase[idx] += increment

    def __correct_phase(self, a, phase):
        ideal_phase = np.array(np.dot(a, self.__trm_setting)).reshape(-1)

        increment = 360
        new_phase = np.round((ideal_phase - phase) / increment) * increment + phase
        """
        print("old", [phase])
        print("ideal", [ideal_phase])
        print(new_phase)
        """
        # exit()
        return new_phase


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
    unique_path = [[t2s(tx_network[0][0]), (antenna.row_col_to_index(0, 0), None)], [t2s(rx_network[0][0]), (None, antenna.row_col_to_index(0, 0))]]
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
