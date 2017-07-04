import numpy as np
import Controllers.Antenna_Calibrator as calibrator
import Utilities.Antenna_Common as common
import Controllers.Matrix_Calibrator_builder as MatrixBuilder
import itertools
import logging


class MutualCalibrator(calibrator.AntennaCalibrator):
    def __init__(self, input_power, input_phase, row_steering, column_steering, dist_rows, dist_columns, filename):
        super(MutualCalibrator, self).__init__(input_power, input_phase, dist_rows, dist_columns, filename)
        self.__logger = logging.getLogger('Mutual')

        self.__trm_setting = common.obtain_shift_phases(column_steering, row_steering,
                                                               self._antenna.quantity_columns,
                                                               self._antenna.quantity_rows, dist_columns, dist_rows,
                                                               common.f)

        self.__trm_setting = list(itertools.chain.from_iterable(self.__trm_setting))
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
        power_shift = common.v2db(abs(self.__equations[(0, 0)])) - self._tx_power.item(0, 0) - \
            self._rx_power.item(0, 0)
        phase_shift = np.mod(np.angle(self.__equations[(0, 0)], deg=True) - self._tx_phase.item(0, 0) -
                             self._rx_phase.item(0, 0) + 180, 360) - 180
        self._rx_power += power_shift
        self._rx_phase += phase_shift

    def generate_cal_paths(self, strategy, pol_mode="TxH-RxV"):
        self.__logger.debug('Polarization mode used: %s', pol_mode)
        self._set_pol_mode(pol_mode)

        self.__last_cal_paths = [strategy, pol_mode]

        [self.__tx_network, self.__rx_network] = self._antenna.get_gain_paths(pol_mode)

        rows, columns = self._antenna.shape
        self.__logger.debug('Antenna shape: %s x %s', rows, columns)

        f = lambda x: [[common.s2t_parameters(x[row][col]) for col in range(columns)] for row in range(rows)]

        tx_network = f(self.__tx_network)
        rx_network = f(self.__rx_network)

        self.__logger.debug('Applying strategy')
        [b, a] = strategy(self._antenna, tx_network, self.__rm_coupling, rx_network)

        rand = lambda x, y: np.random.normal(x, y) if y else x
        f = lambda x: x * common.pol2rec(common.db2v(rand(self._input_power, self._input_delta_power)),
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
        # format_phase = lambda x: np.mod(common.deg2rad(x) + np.pi, 2*np.pi) - np.pi
        # format_phase = lambda x: np.array([-np.sign(val)*(360-abs(val)) if 130 < val < 149 or 159 < val < 168 or 205 < val < 260 else val for val in x.tolist()])
        least_squares = lambda a_mx, b: np.matrix(np.linalg.lstsq(a_mx, b)[0]).reshape(self._antenna.quantity_rows,
                                                                                       self._antenna.quantity_columns)

        a, tx_gain, tx_phase = self.__matrix_builder.get_tx_matrix()
        self._tx_power = least_squares(a, tx_gain)
        self._tx_phase = format_phase(least_squares(a, self.__correct_phase(a, tx_phase)))

        """
        print(a)
        print(tx_phase)
        print(format_phase(tx_phase))
        print(common.rad2deg(least_squares(a, format_phase(tx_phase))))
        """
        a, rx_gain, rx_phase = self.__matrix_builder.get_rx_matrix()
        self._rx_power = least_squares(a, rx_gain)
        self._rx_phase = format_phase(least_squares(a, self.__correct_phase(a, rx_phase)))

        #self.__fix_rx_power_and_phase()

    def calibrate_antenna(self, desired_tx_power, desired_tx_phase, desired_rx_power, desired_rx_phase):
        self.__logger.debug('Calibrating antenna')

        self._calibrate_antenna(desired_tx_power, desired_tx_phase, desired_rx_power, desired_rx_phase)
        self.__logger.info('Generating calibration paths')
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
