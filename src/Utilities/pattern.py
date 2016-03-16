import Utilities.Antenna_Common as common
import numpy as np


class Pattern():

    def __init__(self, d_x, d_y, freq, start_stop_angle, phi, amn):
        eps = 0.000001
        self.__step = 0.1
        decimals = 2
        self.__theta = [round(x, decimals) for x in np.arange(start_stop_angle[0], start_stop_angle[1] + eps, self.__step)]
        self.__pattern = self.__calculate_directivity_pattern(amn, phi, d_x, d_y, freq)

    def __calculate_directivity_pattern(self, amn, phi, d_x, d_y, freq):
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
        k = 2 * np.pi * freq / common.C

        # n is the quantity of rows.
        # m is the quantity of columns.
        f = lambda x: np.array(range(int(-x/2), int(np.ceil(x/2))))

        n = f(len(amn))
        m = f(len(amn[0]))

        # u is the length in x direction
        # v is the length in y direction
        # w is the length in z direction
        u = np.matrix(np.sin(common.deg2rad(self.__theta)) * np.cos(common.deg2rad(phi)))
        v = np.matrix(np.sin(common.deg2rad(self.__theta)) * np.sin(common.deg2rad(phi)))
        # w = np.array(np.cos(self.__theta))
        angle_x = u.T * k * m * d_x
        angle_y = v.T * k * n * d_y

        # angle is a list of matrix
        g = lambda x, y, z: np.repeat(x, len(y), axis=z)
        h = lambda x: np.exp(1j * x).tolist()

        angle = [h(g(angle_x[i], n, 0) + g(angle_y[i].T, m, 1)) for i in range(len(self.__theta))]

        rows = range(len(amn))
        cols = range(len(amn[0]))
        return np.array([sum([amn[row][col] * ang[row][col] for row in rows for col in cols]) for ang in angle][::-1])

    @staticmethod
    def __search_next_max_idx(pattern, idx_initial, delta_idx):
        decreasing = True

        for idx in range(idx_initial + delta_idx, len(pattern) if delta_idx > 0 else 0, delta_idx):
            if decreasing and pattern[idx] > pattern[idx - delta_idx]:
                decreasing = False
            if not decreasing and pattern[idx] < pattern[idx - delta_idx]:
                return pattern[idx - delta_idx]
        return 0

    def __get_width_lobe(self, pattern, idx_max):
        for idx in range(idx_max, len(pattern)):
            if pattern[idx_max] - pattern[idx] > 3:
                return 2*(idx - idx_max)*self.__step
        return 0

    def get_pattern_properties(self):
        """
        this function returns a list with the peak power, the sidelobesPower and the width of the main beam.
        """
        pat = common.v2db(abs(self.__pattern))
        idx_max = np.argmax(pat)
        peak = np.max(pat)

        return [peak - self.__search_next_max_idx(pat, idx_max, -1), peak,
                peak - self.__search_next_max_idx(pat, idx_max, 1), self.__get_width_lobe(pat, idx_max)]

    def get_db(self):
        return common.v2db(abs(self.__pattern))

    @property
    def pattern(self):
        return self.__pattern

    @property
    def angles(self):
        return self.__theta
