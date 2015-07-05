import re
import numpy as np
import math
import cmath

Transmission = 'T'
Reception = 'R'
H_pol = 'H'
V_pol = 'V'

# For RFDN
Rfdn_h_pol = 'hPolarization'
Rfdn_v_pol = 'vPolarization'
Extreme = 'extremeAttached'
SParams = 'sParameters'

Trm_gain_shift = 10
Cable_length = 1
Cable = 'cable'
Circulator = 'circulator'
Trm = 'TRM'
Rm = 'RM'
Psc = 'PSC'
Ph_shifter = 'phase_shifter'
Amplifier = 'amplifier'

c = 299792458       # [m/seg]
f = 1275000000      # [Hz]


def is_cable(component):
    return re.match(Cable, component) if isinstance(component, str) else False


def is_circulator(component):
    return re.match(Circulator, component) if isinstance(component, str) else False


def is_trm(component):
    return re.match(Trm, component) if isinstance(component, str) else False


def is_rm(component):
    return re.match(Rm, component) if isinstance(component, str) else False


def is_psc(component):
    return re.match(Psc + "1.*", component) if isinstance(component, str) else False


def get_qtty_output_ports(component):
    return int(re.match(Psc + "1(.*)", component).group(1)) if is_psc(component) else 1


def get_rm_position(component):
    if not is_rm(component):
        raise Exception("the component should be a RM")
    return eval(re.match(Rm + "(.*)", component).group(1))

"""
def get_qtty_output_ports(psc):
    return int(re.match("PSC1(.*)", psc).group(1))
"""


def get_qtty_ports(component):
    if re.match("{0}|{1}|{2}|{3}|{4}".format(Psc, Cable, Rm, Trm, Circulator), component) is None:
        raise Exception("the input component({0}) is not from the network", component)

    match = re.match("{0}1(.*)|{1}|{2}".format(Psc, Cable, Rm), component)

    # three ports because component is a circulator or a TRM, the second is a PSC, the third corresponds to a cable and
    # the last one a RM
    return 3 if match is None else int(match.group(1)) + 1 if match.group(1) is not None else 2 if re.match(
        Cable, component) else 1


def s2t_parameters(s_matrix):
    s11 = s_matrix.item(0, 0)
    s12 = s_matrix.item(0, 1)
    s21 = s_matrix.item(1, 0)
    s22 = s_matrix.item(1, 1)
    return 1/s21 * np.matrix([[s12*s21 - s11*s22, s11], [-s22, 1]])
    # return 1/s21 * np.matrix([[1, -s22], [s11, s12*s21 - s11*s22]])


def t2s_parameters(t_matrix):
    t11 = t_matrix.item(0, 0)
    t12 = t_matrix.item(0, 1)
    t21 = t_matrix.item(1, 0)
    t22 = t_matrix.item(1, 1)
    return 1/t22 * np.matrix([[t12, t11*t22 - t12*t21], [1, -t21]])
    # return 1/t11 * np.matrix([[t21, t11*t22 - t12*t21], [1, -t12]])


def get_s2p(component, sxp_matrix, mode, idx):
    if is_cable(component):
        fir_p = 0 if mode == Transmission else 1
        sec_p = 1 if mode == Transmission else 0
    elif mode == Reception:
        fir_p = 1 if is_circulator(component) else 2 if is_trm(component) else idx + 1
        sec_p = 2 if is_circulator(component) else 0
    else:
        fir_p = 0
        sec_p = idx + 1

    return np.matrix([[sxp_matrix[fir_p][fir_p], sxp_matrix[fir_p][sec_p]], [sxp_matrix[sec_p][fir_p],
                                                                             sxp_matrix[sec_p][sec_p]]])


def db2v(decibel):
    return np.power(10, decibel/20)


def v2db(voltage):
    if voltage < 0:
        raise Exception("the voltage delivered was negative:", voltage)
    return 20*np.log10(voltage)


def rad2deg(rad):
    return rad*180/cmath.pi


def deg2rad(deg):
    return deg*cmath.pi/180


def pol2rec(module, angle):
    """
    :param module: voltage
    :param angle: degrees
    :return:
    """
    return np.multiply(module, np.exp(1j*deg2rad(angle)))

def parse_polarization_mode(mode):
    """
    :keyword parameters:
    mode -- the required mode, contains the polarization and the transmission reception mode.

    :return:
    modes -- list of S parameter (that determines if is transmission or reception) and polarization pairs
    """
    modes = re.findall("(\w)x(\w)", mode)
    # f = lambda x: "S12" if x == "T" else "S21"
    g = lambda x: Rfdn_h_pol if x == H_pol else Rfdn_v_pol
    return [[mode[0], g(mode[1])]for mode in modes]


def calculate_distances_between_rms(column_length, row_length, dist_columns, dist_rows, row_shift):
    """
    This method calculates all the distances of every radiant module against the one positioned in the upper left
    antenna corner.

    :return matrix_distances: the key value is a tuple containing the distance in (row, column) format and the value
    the position in which the distance can be retrieved
    :return dist: list containing only all the distances to every rm to the one positioned in the corner of the
    antenna. The position of each distance is the value of the matrix_distances dictionary
    """
    get_shift_length = lambda x: (x % 2)/2 if row_shift else 0

    if dist_columns <= dist_rows:
        d_min = dist_columns
        d_max = dist_rows
        antennas_in_max_dir = column_length
        antennas_in_min_dir = row_length

        distance_calculator = lambda x, y: math.sqrt(((get_shift_length(y) + x)*d_min)**2 + (y*d_max)**2)
        position_calculator = lambda x, y: (y, x)
    else:
        d_min = dist_rows
        d_max = dist_columns
        antennas_in_max_dir = row_length
        antennas_in_min_dir = column_length

        distance_calculator = lambda x, y: math.sqrt((x*d_min)**2 + ((get_shift_length(x) + y)*d_max)**2)
        position_calculator = lambda x, y: (x, y)

    rm_used = [0] * antennas_in_max_dir
    dist = list()
    matrix_distances = {}

    for j in range(antennas_in_max_dir):
        for i in range(rm_used[j], antennas_in_min_dir):
            __append_next_distance(j, j + 1, rm_used, matrix_distances, dist, position_calculator, distance_calculator)

    return matrix_distances, dist


def __append_next_distance(j, j2, rm_used, matrix_distances, distances, pos_calculator, dist_calculator):
    # i must compare with the value of the last column in order to cut the
    # unnecessary comparisons
    if j2 == len(rm_used) or rm_used[-1] == rm_used[j2 - 1]:
        matrix_distances[pos_calculator(rm_used[j], j)] = len(distances)
        distances.append(dist_calculator(rm_used[j], j))
        rm_used[j] += 1
        return

    dist_1 = dist_calculator(rm_used[j], j)
    dist_2 = dist_calculator(rm_used[j2], j2)

    if dist_1 > dist_2:
        # if this occurs, meaning that the compared distance is lower than
        # the actual
            # I should get out of here after append all the numbers, j2 must
            # not be greater than i, because of the matrix is not square
            # this must check if the next value is greater than the next
        __append_next_distance(j2, j2+1, rm_used, matrix_distances, distances, pos_calculator, dist_calculator)

        # now i must recheck if the distance with the current i and j is lower
        # than the rest until this value is added
        __append_next_distance(j, j2, rm_used, matrix_distances, distances, pos_calculator, dist_calculator)

    elif dist_1 == dist_2:
        # i know that the last value added is equal to this, then the
        # distance is the same, i can add the distance again and increment
        # the value of the RM
        __append_next_distance(j2, j2+1, rm_used, matrix_distances, distances, pos_calculator, dist_calculator)
        matrix_distances[pos_calculator(rm_used[j], j)] = len(distances) - 1
        matrix_distances[pos_calculator(rm_used[j2] - 1, j2)] = len(distances) - 1
        rm_used[j] += 1
    else:
        # instead of adding the distance, I must check if this value is equal
        # to another distance.
        __append_next_distance(j, j2+1, rm_used, matrix_distances, distances, pos_calculator, dist_calculator)