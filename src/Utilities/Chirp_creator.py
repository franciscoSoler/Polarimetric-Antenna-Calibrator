__author__ = 'fsoler'

import src.Utilities.Antenna_Common as AntennaCommon
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
        self.__add_gain_error = False
        self.__add_phase_error = False
        self.__std_gain = 0
        self.__std_phase = 0

    def add_chirp_errors(self, errors):
        if [True for error in errors if error[0] == AntennaCommon.Phase_chirp_rep_err]:
            self.__add_phase_error = True
            self.__std_phase = [error[1] for error in errors if error[0] == AntennaCommon.Phase_chirp_rep_err].pop()

        if [True for error in errors if error[0] == AntennaCommon.Gain_chirp_rep_err]:
            self.__add_gain_error = True
            self.__std_gain = AntennaCommon.db2v([error[1] for error in errors
                                                  if error[0] == AntennaCommon.Gain_chirp_rep_err].pop())

    def create_ideal_chirp(self, fs, fc, bw, tp, swl):
        return self.__gain * np.exp(1j*self.__phase) * self.__create_unitary_chirp(fs, fc, bw, tp, swl)

    def create_chirp(self, fs, fc, bw, tp, swl):
        """
        :param fs: frequency sampling
        :param fc: central frequency
        :param bw: bandwidth
        :param tp: pulse duration
        :param swl: sampling windows length
        :return:
        """
        f = lambda x, y, z: x + np.random.normal(scale=y) if z else x
        gain = f(self.__gain, self.__std_gain, self.__add_gain_error)
        phase = AntennaCommon.deg2rad(f(self.__phase, self.__std_phase, self.__add_phase_error))

        return gain * np.exp(1j*phase) * self.__create_unitary_chirp(fs, fc, bw, tp, swl)
