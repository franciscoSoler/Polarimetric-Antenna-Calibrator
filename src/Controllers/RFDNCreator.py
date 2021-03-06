#!/usr/bin/python3.3

import sys
import random
import json
import collections
import functools
import math

import src.Utilities.Antenna_Common as AntennaCommon
import src.Utilities.Scattering_Parameters as ScatteringParameters


class AntennaCreator:

    def __init__(self, row_length, dist_rows, dist_columns, row_shift=False):
        """
        :param row_length: quantity of rm in a row
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

        self.__scattering_handler.initialize(delta=0, add_errors=False)
        rm_handler.initialize(delta=0.1, add_errors=True)
        trm_handler.initialize(delta=1.9, add_errors=True)
        circulator_handler.initialize(delta=0.5, add_errors=True)
        psc_handler.initialize(delta=0.5, add_errors=True)

        self.__row_length = row_length
        self.__column_length = 0
        self.__dist_rows = dist_rows
        self.__dist_columns = dist_columns
        self.__row_shift = row_shift

    def __build_front_panel_structure(self, filename):
        (matrix_distances, distances) = AntennaCommon.calculate_distances_between_rms(self.__column_length,
                                                                                      self.__row_length,
                                                                                      self.__dist_columns,
                                                                                      self.__dist_rows,
                                                                                      self.__row_shift)

        att = 0.1
        wavelength = AntennaCommon.c / AntennaCommon.f

        dispersion_params = list(map(lambda x: self.__scattering_handler.get_scattering_matrix("cable",
                                                                                               [att, wavelength,
                                                                                                x]), distances))

        keys = [(row, col) for col in range(self.__row_length) for row in range(self.__column_length)]

        front_panel = []
        f = lambda x: "RM " + str(x)
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
            antennas_in_max_dir = self.__column_length
            antennas_in_min_dir = self.__row_length

            distance_calculator = lambda x, y: math.sqrt(((f(y) + x)*d_min)**2 + (y*d_max)**2)
            position_calculator = lambda x, y: (y, x)
        else:
            d_min = self.__dist_rows
            d_max = self.__dist_columns
            antennas_in_max_dir = self.__row_length
            antennas_in_min_dir = self.__column_length

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

    def __build_rfdn_structure(self, sequence, rm_iterator):
        # el TRM se comporta igual que el cable, no tengo que distinguirlos realmente,
        # los diferentes son los PSC y los RM
        if AntennaCommon.is_rm(sequence[0][0]):
            return sequence[0][0] + next(rm_iterator)

        [component, parameters] = sequence[0]
        structure = collections.OrderedDict()
        structure["sParameters"] = [list(map(str, si)) for si in self.__scattering_handler.get_scattering_matrix(component, parameters)]

        if AntennaCommon.is_psc(component):
            list_cables = []
            for _ in range(AntennaCommon.get_qtty_output_ports(component)):
                semi_structure = self.__build_rfdn_structure(sequence[1:], rm_iterator)
                list_cables.append(semi_structure)
            extreme = list_cables
        else:
            extreme = self.__build_rfdn_structure(sequence[1:], rm_iterator)

        structure["extremeAttached"] = extreme
        return {component: structure}

    def create_structure(self, filename, sequence):
        quantity_signal_splitters = [AntennaCommon.get_qtty_output_ports(component[0]) for component in sequence if
                                     AntennaCommon.is_psc(component[0])]

        quantity_rms = functools.reduce(lambda x, y: x*y, quantity_signal_splitters)
        if quantity_rms % self.__row_length != 0:
            raise Exception("quantity of rms: {0} is not multiple of row length: {1}".format(quantity_rms,
                                                                                             self.__row_length))

        self.__column_length = int(quantity_rms / self.__row_length)

        structure = collections.OrderedDict()

        g = lambda: [" " + str((col, row)) for row in range(self.__row_length) for col in range(self.__column_length)]
        rm_iterator = iter(g())
        structure["vPolarization"] = self.__build_rfdn_structure(sequence, rm_iterator)
        rm_iterator = iter(g())
        structure["hPolarization"] = self.__build_rfdn_structure(sequence, rm_iterator)

        with open(filename + "_rfdn", "w") as f:
            f.write(json.dumps(structure, sort_keys=False, indent=4, separators=(',', ': ')))

        # this is the other file, the one that have all distances between RMs
        self.__build_front_panel_structure(filename)

    @staticmethod
    def __create_parameters(large):
        return [[1 for _ in range(large)] for _ in range(large)]

    def __add_component_behaviour(self, sequence):
        return [[component[0], self.__scattering_handler.get_scattering_matrix(component[0], component[1])]
                for component in sequence]
