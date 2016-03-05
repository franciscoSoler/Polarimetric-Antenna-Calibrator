__author__ = 'fsoler'
import numpy as np
import Utilities.Antenna_Common as common


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

        k = 2 * np.pi * self.__freq / PatternGenerator.C_0

        # n is the quantity of rows.
        # m is the quantity of columns.
        f = lambda x: np.array(range(int(-x/2), int(np.ceil(x/2))))

        n = f(len(amn))
        m = f(len(amn[0]))

        # u is the length in x direction
        # v is the length in y direction
        # w is the length in z direction
        u = np.matrix(np.sin(common.deg2rad(theta)) * np.cos(common.deg2rad(phi)))
        v = np.matrix(np.sin(common.deg2rad(theta)) * np.sin(common.deg2rad(phi)))
        # w = np.array(np.cos(theta))
        angle_x = u.T * k * m * self.__d_x
        angle_y = v.T * k * n * self.__d_y

        # angle is a list of matrix
        g = lambda x, y, z: np.repeat(x, len(y), axis=z)
        h = lambda x: np.exp(1j * x).tolist()

        angle = [h(g(angle_x[i], n, 0) + g(angle_y[i].T, m, 1)) for i in range(len(theta))]

        rows = range(len(amn))
        cols = range(len(amn[0]))
        return [sum([amn[row][col] * ang[row][col] for row in rows for col in cols]) for ang in angle][::-1]

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
        angle_range = [round(x, decimals) for x in np.arange(start_stop_angle[0], start_stop_angle[1] + eps, step)]

        power = common.db2p(abs(output_power)/6)
        b = self.__calculate_directivity_pattern(common.pol2rec(power, np.angle(output_power, deg=True)), angle_range, phi)
        return angle_range, b
        return angle_range, common.pol2rec(common.p2db(abs(b)), np.angle(b, deg=True))
