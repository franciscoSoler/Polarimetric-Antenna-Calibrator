__author__ = 'fsoler'

import src.Utilities.Antenna_Common as AntennaCommon
import numpy as np


class ChirpCreator:

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
        self.__add_phase_error, self.__std_phase = [True if isinstance(i, str) else i for err in errors for i in err
                                                    if err[0] == AntennaCommon.Phase_chirp_rep_err].pop()
        self.__add_gain_error, self.__std_gain = [True if isinstance(i, str) else i for err in errors for i in err
                                                  if err[0] == AntennaCommon.Gain_chirp_rep_err].pop()

    def create_chirp(self, fs, fc, bw, tp, swl):
        """
        :param fs: frequency sampling
        :param fc: central frequency
        :param bw: bandwidth
        :param tp: pulse duration
        :param swl: sampling windows length
        :return:
        """

        kr = bw/tp                      # Chirp rate in range
        t = np.arange(0, swl, 1/fs)     # Time for each sample
        f0 = fc-bw/2                    # Initial frequency

        f = lambda x, y, z: x + np.random.normal(sigma=y) if z else x
        gain = f(self.__gain, self.__std_gain, self.__add_gain_error)
        phase = AntennaCommon.deg2rad(f(self.__phase, self.__std_phase, self.__add_phase_error))

        return gain * np.exp(1j*phase) * np.exp(1j*2*np.pi*(f0*t + kr/2*t**2))