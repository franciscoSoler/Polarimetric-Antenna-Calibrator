import numpy as np
import Controllers.Antenna_Calibrator as calibrator
import Utilities.Antenna_Common as common
import Utilities.Walsh_mtx_creator as WalshCreator
import Utilities.Chirp_creator as ChirpCreator


class ClassicCalibrator(calibrator.AntennaCalibrator):
    __Dec = 11
    __Swl = 2*common.tp  # sampling window length [sec]

    def __init__(self, input_power, input_phase, dist_rows, dist_columns, filename):
        super(ClassicCalibrator, self).__init__(input_power, input_phase, dist_rows, dist_columns, filename)
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
        chirp_parameters = [common.fs, common.fc, common.bw, common.tp, self.__Swl]
        # chirp = [self.__chirp_creator.create_chirp(*chirp_parameters) for _ in range(self.__quantity_elements)]
        chirp = [self.__chirp_creator.create_chirp(*chirp_parameters) for _ in range(int(sequences))]
        chirp_rep = np.matrix(self.__chirp_creator.create_chirp_replica(*chirp_parameters))

        amp = common.db2v(att_db)
        ph_rad = common.deg2rad(ph_deg)

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
        estimated_power = np.around(common.v2db(abs(signal_est)), decimals=self.__Dec)
        f = lambda x: np.resize(x, (self._antenna.quantity_rows, self._antenna.quantity_columns))
        return f(estimated_power), f(estimated_phase)

    def _obtain_tx_rx_power(self):
        self._power_calculated = True
        # tx_signal, rx_signal = self._antenna.get_gain_paths(self._pol_mode)
        tx_power, tx_phase, rx_power, rx_phase = self.get_antenna_gain_paths(complete_calibration=False)
        self._tx_power, self._tx_phase = self.__obtain_estimated_signal(tx_power, tx_phase)
        self._rx_power, self._rx_phase = self.__obtain_estimated_signal(rx_power, rx_phase)

    def set_pol_mode(self, pol_mode):
        self._set_pol_mode(pol_mode)

    def calibrate_antenna(self, desired_tx_power, desired_tx_phase, desired_rx_power, desired_rx_phase):
        self._calibrate_antenna(desired_tx_power, desired_tx_phase, desired_rx_power, desired_rx_phase)
        self._power_calculated = False
