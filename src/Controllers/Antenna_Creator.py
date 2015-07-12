#!/usr/bin/python3.3

import json
import collections
import functools
import math
import numpy as np

import src.Utilities.Antenna_Common as AntennaCommon
import src.Utilities.Scattering_Parameters as ScatteringParameters


class AntennaCreator:

    def __init__(self, cols, dist_rows, dist_columns, row_shift=False):
        """

        :param cols: quantity of rm in a row
        :param dist_rows: physical distance between rows of rm
        :param dist_columns: physical distance between columns of rm
        :param row_shift: is the shift of odd antenna rows (default False)
        :return:
        """
        self.__scattering_handler = ScatteringParameters.CableScatteringParameters()
        rm_handler = ScatteringParameters.RmScatteringParameters()
        trm_handler = ScatteringParameters.TrmScatteringParameters()
        circulator_handler = ScatteringParameters.CirculatorScatteringParameters()
        psc_handler = ScatteringParameters.PscScatteringParameters()

        self.__scattering_handler.set_successor(psc_handler)
        psc_handler.set_successor(trm_handler)
        trm_handler.set_successor(circulator_handler)
        circulator_handler.set_successor(rm_handler)

        self.__scattering_handler.initialize()
        rm_handler.initialize()
        trm_handler.initialize()
        circulator_handler.initialize()
        psc_handler.initialize()

        self.__quantity_cols = cols
        self.__quantity_rows = 0
        self.__dist_rows = dist_rows
        self.__dist_columns = dist_columns
        self.__row_shift = row_shift

    def __row_col_to_index(self, row, col):
        return row * self.__quantity_cols + col

    def add_errors(self, errors):
        if not isinstance(errors, list) or len(errors) == 0 or [True for error in errors if len(error) != 2]:
            raise Exception('errors are not well created')

        rm_handler = ScatteringParameters.RmScatteringParameters()
        trm_handler = ScatteringParameters.TrmScatteringParameters()
        circulator_handler = ScatteringParameters.CirculatorScatteringParameters()
        psc_handler = ScatteringParameters.PscScatteringParameters()

        self.__scattering_handler.set_successor(psc_handler)
        psc_handler.set_successor(trm_handler)
        trm_handler.set_successor(circulator_handler)
        circulator_handler.set_successor(rm_handler)

        f = lambda x: [True if isinstance(idx, str) else idx for error in errors for idx in error if error[0] == x]
        self.__scattering_handler.initialize(*f(AntennaCommon.Cable_error))
        rm_handler.initialize(*f(AntennaCommon.Rm_error))
        trm_handler.initialize(*f(AntennaCommon.Trm_error))
        circulator_handler.initialize(*f(AntennaCommon.Circulator_error))
        psc_handler.initialize(*f(AntennaCommon.Psc_error))

    def __build_front_panel_structure(self, filename):
        parameters = [self.__quantity_rows, self.__quantity_cols, self.__dist_columns, self.__dist_rows, self.__row_shift]
        (matrix_distances, distances) = AntennaCommon.calculate_distances_between_rms(*parameters)

        att = 0.1
        wavelength = AntennaCommon.c / AntennaCommon.f

        dispersion_params = list(map(lambda x: self.__scattering_handler.get_scattering_matrix("cable",
                                                                                               [att, wavelength,
                                                                                                x]), distances))

        keys = [(row, col) for col in range(self.__quantity_cols) for row in range(self.__quantity_rows)]

        front_panel = []
        f = lambda x: AntennaCommon.Rm + " " + str(x)
        g = lambda x: list(map(lambda y: list(map(str, y)), x))
        for key in keys:

            parameters = [[f(new_key), g(dispersion_params[matrix_distances[
                tuple(map(lambda x, y: abs(x-y), key, new_key))]])] for new_key in keys]
            front_panel.append([f(key), parameters])

        with open(filename + "_panel", "w") as f:
            f.write(json.dumps(front_panel, sort_keys=False, indent=4, separators=(',', ': ')))

    def __calculate_distances_between_rms(self):
        """
        This method calculates all the distances of every radiant module to the one positioned in the upper left antenna
        corner.

        :return matrix_distances: the key value is a tuple containing the distance in (row, column) format and the value
        the position in which the distance can be retrieved
        :return dist: list containing only all the distances to every rm to the one positioned in the corner of the
        antenna. The position of each distance is the value of the matrix_distances dictionary
        """
        f = lambda x: (x % 2)/2 if self.__row_shift else 0

        if self.__dist_columns <= self.__dist_rows:
            d_min = self.__dist_columns
            d_max = self.__dist_rows
            antennas_in_max_dir = self.__quantity_rows
            antennas_in_min_dir = self.__quantity_cols

            distance_calculator = lambda x, y: math.sqrt(((f(y) + x)*d_min)**2 + (y*d_max)**2)
            position_calculator = lambda x, y: (y, x)
        else:
            d_min = self.__dist_rows
            d_max = self.__dist_columns
            antennas_in_max_dir = self.__quantity_cols
            antennas_in_min_dir = self.__quantity_rows

            distance_calculator = lambda x, y: math.sqrt((x*d_min)**2 + ((f(x) + y)*d_max)**2)
            position_calculator = lambda x, y: (x, y)

        rm_used = [0] * antennas_in_max_dir
        dist = list()
        matrix_distances = {}

        for j in range(antennas_in_max_dir):
            for i in range(rm_used[j], antennas_in_min_dir):
                self.__append_next_distance(j, j + 1, rm_used, matrix_distances, dist, position_calculator,
                                            distance_calculator)

        return matrix_distances, dist

    def __append_next_distance(self, j, j2, rm_used, matrix_distances, distances, pos_calculator, dist_calculator):
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
            self.__append_next_distance(j2, j2+1, rm_used, matrix_distances, distances, pos_calculator, dist_calculator)

            # now i must recheck if the distance with the current i and j is lower
            # than the rest until this value is added
            self.__append_next_distance(j, j2, rm_used, matrix_distances, distances, pos_calculator, dist_calculator)

        elif dist_1 == dist_2:
            # i know that the last value added is equal to this, then the
            # i know that the last value added is equal to this, then the
            # distance is the same, i can add the distance again and increment
            # the value of the RM
            self.__append_next_distance(j2, j2+1, rm_used, matrix_distances, distances, pos_calculator, dist_calculator)
            matrix_distances[rm_used[j], j] = len(distances) - 1
            rm_used[j] += 1
        else:
            # instead of adding the distance, I must check if this value is equal
            # to another distance.
            self.__append_next_distance(j, j2+1, rm_used, matrix_distances, distances, pos_calculator, dist_calculator)

    def __build_rfdn_structure(self, sequence, rm_iterator, trm_state_iterator):
        # el TRM se comporta igual que el cable, no tengo que distinguirlos realmente,
        # los diferentes son los PSC y los RM
        if AntennaCommon.is_rm(sequence[0][0]):
            return sequence[0][0] + next(rm_iterator)

        [component, parameters] = sequence[0]

        structure = collections.OrderedDict()

        if AntennaCommon.is_trm(component):
            is_dead, steering_shift = next(trm_state_iterator)
            parameters = list(parameters)
            parameters[1] = np.mod(parameters[1] + steering_shift, 360)
            parameters.append(is_dead)
            structure[AntennaCommon.Dead] = str(is_dead)

        structure[AntennaCommon.SParams] = [list(map(str, si)) for si in
                                            self.__scattering_handler.get_scattering_matrix(component, parameters)]

        if AntennaCommon.is_psc(component):
            list_cables = []
            for _ in range(AntennaCommon.get_qtty_output_ports(component)):
                semi_structure = self.__build_rfdn_structure(sequence[1:], rm_iterator, trm_state_iterator)
                list_cables.append(semi_structure)
            extreme = list_cables
        else:
            extreme = self.__build_rfdn_structure(sequence[1:], rm_iterator, trm_state_iterator)

        structure[AntennaCommon.Extreme] = extreme
        return {component: structure}

    def create_structure(self, filename, sequence, row_steering, column_steering, dead_trms=()):
        quantity_signal_splitters = [AntennaCommon.get_qtty_output_ports(component[0]) for component in sequence if
                                     AntennaCommon.is_psc(component[0])]

        quantity_rms = functools.reduce(lambda x, y: x*y, quantity_signal_splitters)
        if quantity_rms % self.__quantity_cols != 0:
            raise Exception("quantity of rms: {0} is not multiple of row length: {1}".format(quantity_rms,
                                                                                             self.__quantity_cols))
        self.__quantity_rows = int(quantity_rms / self.__quantity_cols)

        trms_dead = [False] * quantity_rms
        if dead_trms:
            if max(dead_trms)[0] >= self.__quantity_rows or max(dead_trms)[1] >= self.__quantity_cols:
                raise Exception("ERROR: TRM dead has an index bigger than the antenna size, TRM: {}".format(max(dead_trms)))

            for idx in [self.__row_col_to_index(*pair) for pair in dead_trms]:
                trms_dead[idx] = True

        # i must create an iterator with pair true (if dead or not) and angle for the TRM shift

        f = lambda row, col: np.mod(row * row_steering + col * column_steering, 360)
        steering_angle = [f(row, col) for col in range(self.__quantity_cols) for row in range(self.__quantity_rows)]
        trm_state = list(zip(trms_dead, steering_angle))
        structure = collections.OrderedDict()

        g = lambda: [" "+str((row, col)) for col in range(self.__quantity_cols) for row in range(self.__quantity_rows)]
        rm_iterator = iter(g())
        steering_iterator = iter(trm_state)
        structure[AntennaCommon.Rfdn_v_pol] = self.__build_rfdn_structure(sequence, rm_iterator, steering_iterator)
        rm_iterator = iter(g())
        steering_iterator = iter(trm_state)
        structure[AntennaCommon.Rfdn_h_pol] = self.__build_rfdn_structure(sequence, rm_iterator, steering_iterator)

        with open(filename + "_rfdn", "w") as f:
            f.write(json.dumps(structure, sort_keys=False, indent=4, separators=(',', ': ')))

        # this is the other file, the one that have all distances between RMs
        self.__build_front_panel_structure(filename)
