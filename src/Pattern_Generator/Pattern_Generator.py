__author__ = 'fsoler'
import math
import numpy as np


def deg_to_rad(deg):
    return deg * math.pi / 180


def rad_to_deg(rad):
    return rad * 180 / math.pi


class PatternGenerator:
    # speed of light [m/sec]
    C_0 = 299792458

    def __init__(self, freq, d_x, d_y):
        self.__freq = freq
        self.__d_x = d_x
        self.__d_y = d_y

    def __calculate_directivity_pattern(self, amn, theta, phi=0):
        """
        This function generates the beam pattern from the antenna matrix.
        the directivity should be calculated as = sum(sum(amn*exp(j[ku(mdx+ deltan) + kvndy + kwzmn) if contains
        deformations
        k = 2pi/lambda

        :param amn: the power/phase pairs of each element, this is a list of lists, the first list has the rows of the
                   antenna, the second one has the pairs per each element.
        :param theta: angles in z plane
        :param phi: angles in xy plane (default = 0)
        :return:
        the directivity pattern
        """

        k = 2 * math.pi * self.__freq / PatternGenerator.C_0

        # n is the quantity of rows.
        # m is the quantity of columns.
        f = lambda x: np.array(range(int(-x/2), math.ceil(x/2)))

        n = f(len(amn))
        m = f(len(amn[0]))

        # u is the length in x direction
        # v is the length in y direction
        # w is the length in z direction
        u = np.matrix([[math.sin(deg_to_rad(ang)) * math.cos(deg_to_rad(phi))] for ang in theta])
        v = np.matrix([[math.sin(deg_to_rad(ang)) * math.sin(deg_to_rad(phi))] for ang in theta])
        # w = np.array(math.cos(theta))
        angle_x = u * k * m * self.__d_x
        angle_y = v * k * n * self.__d_y

        # angle is a list of matrix
        g = lambda x, y, z: np.repeat(x, len(y), axis=z)
        h = lambda x: np.exp(1j * x).tolist()

        angle = [h(g(angle_x[i], n, 0) + g(angle_y[i].T, m, 1)) for i in range(len(theta))]

        rows = range(len(amn))
        cols = range(len(amn[0]))
        return [sum([amn[row][col] * ang[row][col] for row in rows for col in cols]) for ang in angle]

    def generate_pattern(self, output_power, start_stop_angle, phi=0):
        """
        This function generates the beam pattern from the antenna matrix.
        :param output_power: is the matrix of gain-phase transmitted pairs. The format is, the global list represent the
         rows the inner list contain one gain-phase pair for each column of the antenna.
        :param start_stop_angle: list containing the start and stop angle values in degrees.
        :param phi: represents the cut in which the pattern is generated (default = 0)
        :return:
            pattern -- list of angle-power pairs
        """

        eps = 0.000001
        step = 0.1
        decimals = 2
        angle_range = [round(x, decimals) for x in np.arange(start_stop_angle[0], start_stop_angle[1] + eps,
                                                             step).tolist()]

        b = self.__calculate_directivity_pattern(output_power, angle_range, phi)
        return angle_range, b
