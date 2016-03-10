__author__ = 'fsoler'
import numpy as np
import Utilities.Antenna_Common as common
import Utilities.pattern as pat


class PatternGenerator:
    # speed of light [m/sec]
    C_0 = 299792458

    def __init__(self, freq, d_x, d_y):
        self.__freq = freq
        self.__d_x = d_x
        self.__d_y = d_y

    def generate_pattern(self, output_power, output_phase, start_stop_angle=[], phi=0):
        """
        This function generates the beam pattern from the antenna matrix.
        :param output_power: is the matrix of gain-phase transmitted pairs. The format is, the global list represent the
         rows the inner list contain one gain-phase pair for each column of the antenna.
        :param start_stop_angle: list containing the start and stop angle values in degrees.
        :param phi: represents the cut in which the pattern is generated (default = 0)
        :return:
            pattern -- list of angle-power pairs in watts
        """
        amn = common.pol2rec(common.db2v(output_power), output_phase)
        return pat.Pattern(self.__d_x, self.__d_y, self.__freq, start_stop_angle, phi, amn)
