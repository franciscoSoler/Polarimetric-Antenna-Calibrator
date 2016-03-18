#!/usr/bin/python3.3

import json
import collections
import itertools
import functools
import math
import numpy as np

import Utilities.Antenna_Common as AntennaCommon
import Utilities.Scattering_Parameters as ScatteringParameters


class AntennaCreator:

    def __init__(self, rows, cols, dist_rows, dist_columns, row_shift=False):
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
        self.__quantity_rows = rows
        self.__dist_rows = dist_rows
        self.__dist_columns = dist_columns
        self.__row_shift = row_shift

    def __row_col_to_index(self, row, col):
        return row * self.__quantity_cols + col

    def add_errors(self, errors):
        if not isinstance(errors, list):
            raise Exception('errors are not well created. The variable is not a list')

        if len(errors) == 0:
            return
        if [True for error in errors if len(error) != 2]:
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
        parameters = [self.__quantity_rows, self.__quantity_cols, self.__dist_rows, self.__dist_columns, self.__row_shift]
        (matrix_distances, distances) = AntennaCommon.calculate_distances_between_rms(*parameters)

        att = 0.1
        wavelength = AntennaCommon.C / AntennaCommon.f

        dispersion_params = list(map(lambda x: self.__scattering_handler.get_scattering_matrix("cable",
                                                                                               [att, x, wavelength]),
                                     distances))
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
            parameters[1] = np.mod(parameters[1] + steering_shift + 180, 360) - 180
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
        if quantity_rms != self.__quantity_cols * self.__quantity_rows:
            raise Exception("quantity of rms from RFDN: {0} is different from the front of the antenna: {1}".format(
                                                            quantity_rms, self.__quantity_cols*self.__quantity_rows))

        trms_dead = [False] * quantity_rms
        if dead_trms:
            if max(dead_trms)[0] >= self.__quantity_rows or max(dead_trms)[1] >= self.__quantity_cols:
                raise Exception("ERROR: TRM dead has an index bigger than the antenna size, TRM: {}".format(max(dead_trms)))

            for idx in [self.__row_col_to_index(*pair) for pair in dead_trms]:
                trms_dead[idx] = True

        steering_angle = AntennaCommon.obtain_shift_phases(column_steering, row_steering, self.__quantity_cols,
                                                           self.__quantity_rows, self.__dist_columns, self.__dist_rows,
                                                           AntennaCommon.f)

        trm_state = list(zip(trms_dead, list(itertools.chain.from_iterable(steering_angle))))
        structure = collections.OrderedDict()

        g = lambda: [" "+str((row, col)) for row in range(self.__quantity_rows) for col in range(self.__quantity_cols)]
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
