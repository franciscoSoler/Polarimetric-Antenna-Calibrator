__author__ = 'fsoler'

import Utilities.Antenna_Common as AntennaCommon
import numpy as np


class ChirpCreator:
    @staticmethod
    def __create_unitary_chirp(fs, fc, bw, tp, swl):
        kr = bw/tp                      # Chirp rate in range
        t = np.arange(0, swl, 1/fs)     # Time for each sample
        f0 = fc-bw/2                    # Initial frequency

        return np.array(np.exp(1j*2*np.pi*(f0*t + kr/2*t**2)))

    def __init__(self, power, phase):
        """
        :param power: in dB
        :param phase: in deg
        :return:
        """
        self.__gain = AntennaCommon.db2v(power)
        self.__phase = phase
        self.__std_gain = 0
        self.__std_phase = 0
        self.__delta_chirp_rep_phase = 0
        self.__delta_chirp_rep_gain = 0

    def __add_normal_errors(self, initial_gain, initial_phase):
        f = lambda x, std: np.random.normal(x, std) if std else x
        gain = f(initial_gain, self.__std_gain)
        phase = f(initial_phase, self.__std_phase)
        return AntennaCommon.pol2rec(gain, phase)


    def __set_std_errors(self, errors):
        self.__std_gain = AntennaCommon.db2v(errors[0])
        self.__std_phase = errors[1]

    def __set_rep_std_errors(self, errors):
        self.__delta_chirp_rep_gain = AntennaCommon.db2v(errors[0])
        self.__delta_chirp_rep_phase = errors[1]

    def add_chirp_errors(self, errors):
        err = {AntennaCommon.Inter_pulse_gain_err: self.__set_std_errors,
               AntennaCommon.Chirp_rep_err: self.__set_rep_std_errors}
        [err[error[0]](error[1]) for error in errors if error[0] in err]

    def create_ideal_chirp(self, fs, fc, bw, tp, swl):
        return AntennaCommon.pol2rec(self.__gain, self.__phase) * self.__create_unitary_chirp(fs, fc, bw, tp, swl)

    def create_chirp_replica(self, fs, fc, bw, tp, swl):
        amp = self.__add_normal_errors(self.__gain + self.__delta_chirp_rep_gain,
                                               self.__phase + self.__delta_chirp_rep_phase)

        return amp * self.__create_unitary_chirp(fs, fc, bw, tp, swl)

    def create_chirp(self, fs, fc, bw, tp, swl):
        """
        :param fs: frequency sampling
        :param fc: central frequency
        :param bw: bandwidth
        :param tp: pulse duration
        :param swl: sampling windows length
        :return:
        """
        amp = self.__add_normal_errors(self.__gain, self.__phase)
        return amp * self.__create_unitary_chirp(fs, fc, bw, tp, swl)
