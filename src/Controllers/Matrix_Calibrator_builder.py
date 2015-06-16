__author__ = 'fsoler'

from abc import ABCMeta, abstractmethod
import numpy as np
import math
import src.Utilities.Antenna_Common as Common


class MatrixCalibratorBuilder(object):
    __metaclass__ = ABCMeta
    _tx_a = []
    _rx_a = []
    _tx_gain = []
    _tx_phase = []
    _rx_gain = []
    _rx_phase = []

    def __init__(self):
        self._antenna = None
        self._equations = None
        self.__successor = None

    def initialize_matrix_builder(self, antenna, equations):
        self._antenna = antenna
        self._equations = equations
        self.__initialize_output_matrix()
        if self.__successor is not None:
            self.__successor.initialize_matrix_builder(antenna, equations)

    def __initialize_output_matrix(self):
        self._tx_a.clear()
        self._rx_a.clear()
        self._tx_gain.clear()
        self._rx_gain.clear()
        self._tx_phase.clear()
        self._rx_phase.clear()

    @abstractmethod
    def _add_equations(self):
        pass

    @abstractmethod
    def _is_executable(self):
        pass

    def build_matrix(self):
        if self._is_executable():
            self._add_equations()
        if self.__successor is not None:
            self.__successor.build_matrix()

    def get_tx_matrix(self):
        return [np.matrix(self._tx_a), np.array(self._tx_gain), np.array(self._tx_phase)]

    def get_rx_matrix(self):
        return [np.matrix(self._rx_a), np.array(self._rx_gain), np.array(self._rx_phase)]

    def set_successor(self, successor):
        self.__successor = successor

    def _add_trans_rec_file(self, first_rm, common_rm, second_rm):
        equation = [0] * self._antenna.quantity_columns * self._antenna.quantity_rows
        equation[first_rm] = 1
        equation[second_rm] = -1

        self._tx_a.append(equation)
        self._rx_a.append(equation)

        f = lambda x, y: Common.v2db(abs(x)) - Common.v2db(abs(y))
        self._tx_gain.append(f(self._equations[(first_rm, common_rm)], self._equations[(second_rm, common_rm)]))
        self._rx_gain.append(f(self._equations[(common_rm, first_rm)], self._equations[(common_rm, second_rm)]))

        g = lambda x, y: np.angle(x, deg=True) - np.angle(y, deg=True)
        self._tx_phase.append(g(self._equations[(first_rm, common_rm)], self._equations[(second_rm, common_rm)]))
        self._rx_phase.append(g(self._equations[(common_rm, first_rm)], self._equations[(common_rm, second_rm)]))


class LinearBuilder(MatrixCalibratorBuilder):

    __Minimum_length = 3

    def __init__(self):
        super(LinearBuilder, self).__init__()

    def _is_executable(self):
        return False if max(self._antenna.quantity_columns, self._antenna.quantity_rows) < self.__Minimum_length \
            else True

    def __add_aligned_rms(self, common_rm, row, column, vertical_dist, horizontal_dist):
        first_rm = self._antenna.row_col_to_index(row - horizontal_dist, column - vertical_dist)
        second_rm = self._antenna.row_col_to_index(row + horizontal_dist, column + vertical_dist)
        self._add_trans_rec_file(first_rm, common_rm, second_rm)

        if vertical_dist != 0 and horizontal_dist != 0:
            first_rm = self._antenna.row_col_to_index(row - horizontal_dist, column + vertical_dist)
            second_rm = self._antenna.row_col_to_index(row + horizontal_dist, column - vertical_dist)
            self._add_trans_rec_file(first_rm, common_rm, second_rm)

    def _add_equations(self):
        for column in range(self._antenna.quantity_columns):
            for row in range(self._antenna.quantity_rows):
                common_rm = self._antenna.row_col_to_index(row, column)
                horizontal_dist = 1
                while row - horizontal_dist >= 0 and row + horizontal_dist < self._antenna.quantity_rows:
                    vertical_dist = 0
                    while(vertical_dist <= horizontal_dist and column - vertical_dist >= 0 and
                            column + vertical_dist < self._antenna.quantity_columns):

                        self.__add_aligned_rms(common_rm, row, column, vertical_dist, horizontal_dist)
                        vertical_dist += 1
                    horizontal_dist += 1

                vertical_dist = 1
                while column - vertical_dist >= 0 and column + vertical_dist < self._antenna.quantity_columns:
                    horizontal_dist = 0
                    while(horizontal_dist < vertical_dist and row - horizontal_dist >= 0 and
                            row + horizontal_dist < self._antenna.quantity_rows):

                        self.__add_aligned_rms(common_rm, row, column, vertical_dist, horizontal_dist)
                        horizontal_dist += 1
                    vertical_dist += 1


class CrossBuilder(MatrixCalibratorBuilder):

    __Minimum_length = 2

    def __init__(self):
        super(CrossBuilder, self).__init__()

    def _is_executable(self):
        return False if max(self._antenna.quantity_columns, self._antenna.quantity_rows) < self.__Minimum_length \
            else True

    def __add_cross_h_rms(self, common_rm, row, column, rm_dist):
        f = lambda x: int(math.copysign(x, rm_dist))
        dist = abs(rm_dist)
        if column + dist < self._antenna.quantity_columns:
            rotation = 0
            while column - rotation >= 0 and rotation <= dist:
                first_rm = self._antenna.row_col_to_index(row + f(rotation), column + dist)
                second_rm = self._antenna.row_col_to_index(row + rm_dist, column - rotation)
                self._add_trans_rec_file(first_rm, common_rm, second_rm)
                rotation += 1

        if column - dist >= 0:
            rotation = 0
            while column + rotation < self._antenna.quantity_columns and rotation < dist:
                first_rm = self._antenna.row_col_to_index(row + f(rotation), column - dist)
                second_rm = self._antenna.row_col_to_index(row + rm_dist, column + rotation)
                self._add_trans_rec_file(first_rm, common_rm, second_rm)
                rotation += 1

    def __add_cross_v_rms(self, common_rm, row, column, rm_dist):
        f = lambda x: int(math.copysign(x, rm_dist))
        dist = abs(rm_dist)
        if row + dist < self._antenna.quantity_rows:
            rotation = 0
            while row - rotation >= 0 and rotation <= dist:
                first_rm = self._antenna.row_col_to_index(row + dist, column + f(rotation))
                second_rm = self._antenna.row_col_to_index(row - rotation, column + rm_dist)
                self._add_trans_rec_file(first_rm, common_rm, second_rm)
                rotation += 1

        if row - dist >= 0:
            rotation = 0
            while row + rotation < self._antenna.quantity_rows and rotation < dist:
                first_rm = self._antenna.row_col_to_index(row - dist, column + f(rotation))
                second_rm = self._antenna.row_col_to_index(row + rotation, column + rm_dist)
                self._add_trans_rec_file(first_rm, common_rm, second_rm)
                rotation += 1

    def _add_equations(self):
        """
        this is only for the RMs of the edge
        """
        if self._antenna.quantity_rows != 2 or self._antenna.quantity_columns == 2:
            for column in range(self._antenna.quantity_columns):
                row = 0
                common_rm = self._antenna.row_col_to_index(row, column)
                rm_dist = 1
                while row + rm_dist < self._antenna.quantity_rows:
                    self.__add_cross_h_rms(common_rm, row, column, rm_dist)
                    rm_dist += 1

                row = self._antenna.quantity_rows - 1
                common_rm = self._antenna.row_col_to_index(row, column)
                rm_dist = -1
                while row + rm_dist >= 0:
                    self.__add_cross_h_rms(common_rm, row, column, rm_dist)
                    rm_dist -= 1

        if self._antenna.quantity_columns != 2:
            for row in range(self._antenna.quantity_rows):
                column = 0
                common_rm = self._antenna.row_col_to_index(row, column)
                rm_dist = 1
                while column + rm_dist < self._antenna.quantity_columns:
                    self.__add_cross_v_rms(common_rm, row, column, rm_dist)
                    rm_dist += 1

                column = self._antenna.quantity_columns - 1
                common_rm = self._antenna.row_col_to_index(row, column)
                rm_dist = -1
                while column + rm_dist >= 0:
                    self.__add_cross_v_rms(common_rm, row, column, rm_dist)
                    rm_dist -= 1


class AsymmetricalBuilder(MatrixCalibratorBuilder):
    """
    i think this class must be deprecated, I don't know if it is necessary to add these other equations
    """
    def __init__(self):
        super(AsymmetricalBuilder, self).__init__()

    def _is_executable(self):
        pass

    def _add_equations(self):
        pass


class TinyBuilder(MatrixCalibratorBuilder):

    def __init__(self):
        super(TinyBuilder, self).__init__()

    def _is_executable(self):
        return True if self._antenna.quantity_rows == 2 or self._antenna.quantity_columns == 2 else False

    def __add_file(self, first_rm, second_rm):
        equation = [0] * self._antenna.quantity_columns * self._antenna.quantity_rows
        equation[first_rm] = 2
        equation[second_rm] = -2
        self._tx_a.append(equation)
        self._rx_a.append(equation)

        f = lambda x, y, z, w: Common.v2db(abs(x)) + Common.v2db(abs(y)) - Common.v2db(abs(z)) - Common.v2db(abs(w))
        self._tx_gain.append(f(self._equations[(first_rm, first_rm)], self._equations[(first_rm, second_rm)],
                               self._equations[(second_rm, first_rm)], self._equations[(second_rm, second_rm)]))
        self._rx_gain.append(f(self._equations[(first_rm, first_rm)], self._equations[(second_rm, first_rm)],
                               self._equations[(first_rm, second_rm)], self._equations[(second_rm, second_rm)]))

        g = lambda x, y, z, w: np.angle(x, deg=True) + np.angle(y, deg=True) - np.angle(z, deg=True) - \
            np.angle(w, deg=True)
        self._tx_phase.append(g(self._equations[(first_rm, first_rm)], self._equations[(first_rm, second_rm)],
                                self._equations[(second_rm, first_rm)], self._equations[(second_rm, second_rm)]))
        self._rx_phase.append(g(self._equations[(first_rm, first_rm)], self._equations[(second_rm, first_rm)],
                                self._equations[(first_rm, second_rm)], self._equations[(second_rm, second_rm)]))

    def _add_equations(self):
        if self._antenna.quantity_rows == 2:
            row = 0
            for column in range(self._antenna.quantity_columns):
                first_rm = self._antenna.row_col_to_index(row, column)
                second_rm = self._antenna.row_col_to_index(row + 1, column)
                self.__add_file(first_rm, second_rm)
        if self._antenna.quantity_columns == 2:
            column = 0
            for row in range(self._antenna.quantity_rows):
                first_rm = self._antenna.row_col_to_index(row, column)
                second_rm = self._antenna.row_col_to_index(row, column + 1)
                self.__add_file(first_rm, second_rm)


class DoubleBuilder(MatrixCalibratorBuilder):

    def __init__(self):
        super(DoubleBuilder, self).__init__()

    def _is_executable(self):
        return True if self._antenna.quantity_rows != 2 and self._antenna.quantity_columns != 2 else False

    def __add_file(self, first_rm, second_rm, first_int_rm=0, second_int_rm=0):
        # print("Double builder, rms", first_rm, second_rm)
        equation = [0] * self._antenna.quantity_columns * self._antenna.quantity_rows
        equation[first_rm] = 2
        equation[second_rm] = -2
        self._tx_a.append(equation)
        self._rx_a.append(equation)

        f = lambda x, y, z, w: Common.v2db(abs(x)) + Common.v2db(abs(y)) - Common.v2db(abs(z)) - Common.v2db(abs(w))
        self._tx_gain.append(f(self._equations[(first_rm, first_rm)], self._equations[(first_rm, second_rm)],
                               self._equations[(second_rm, first_rm)], self._equations[(second_rm, second_rm)]))
        self._rx_gain.append(f(self._equations[(first_rm, first_rm)], self._equations[(second_rm, first_rm)],
                               self._equations[(first_rm, second_rm)], self._equations[(second_rm, second_rm)]))

        g = lambda x, y, z, w: np.angle(x, deg=True) + np.angle(y, deg=True) - np.angle(z, deg=True) - \
            np.angle(w, deg=True)
        self._tx_phase.append(g(self._equations[(first_rm, first_rm)], self._equations[(first_rm, second_rm)],
                                self._equations[(second_rm, first_rm)], self._equations[(second_rm, second_rm)]))
        self._rx_phase.append(g(self._equations[(first_rm, first_rm)], self._equations[(second_rm, first_rm)],
                                self._equations[(first_rm, second_rm)], self._equations[(second_rm, second_rm)]))

        if first_int_rm != 0:
            # print("Double builder, int rms", first_int_rm, second_int_rm)

            self._tx_a.append(equation)
            rx_equation = [0] * self._antenna.quantity_columns * self._antenna.quantity_rows
            rx_equation[first_int_rm] = 2
            rx_equation[second_int_rm] = -2
            self._rx_a.append(equation)

            self._tx_gain.append(f(self._equations[(first_rm, first_int_rm)],
                                   self._equations[(first_rm, second_int_rm)],
                                   self._equations[(second_rm, first_int_rm)],
                                   self._equations[(second_rm, second_int_rm)]))
            self._rx_gain.append(f(self._equations[(first_rm, first_int_rm)],
                                   self._equations[(second_rm, first_int_rm)],
                                   self._equations[(first_rm, second_int_rm)],
                                   self._equations[(second_rm, second_int_rm)]))

            g = lambda x, y, z, w: np.angle(x, deg=True) + np.angle(y, deg=True) - np.angle(z, deg=True) - \
                np.angle(w, deg=True)
            self._tx_phase.append(g(self._equations[(first_rm, first_int_rm)],
                                    self._equations[(first_rm, second_int_rm)],
                                    self._equations[(second_rm, first_int_rm)],
                                    self._equations[(second_rm, second_int_rm)]))
            self._rx_phase.append(g(self._equations[(first_rm, first_int_rm)],
                                    self._equations[(second_rm, first_int_rm)],
                                    self._equations[(first_rm, second_int_rm)],
                                    self._equations[(second_rm, second_int_rm)]))

    def _add_equations(self):
        if self._antenna.quantity_rows == 1:
            row = 0
            for first_col in range(self._antenna.quantity_columns - 1):
                for sec_col in range(first_col + 1, self._antenna.quantity_columns):
                    first_rm = self._antenna.row_col_to_index(row, first_col)
                    second_rm = self._antenna.row_col_to_index(row, sec_col)
                    self.__add_file(first_rm, second_rm)

                    if ((sec_col - first_col) % 2) != 0 and sec_col - first_col > 1:
                        for delta in range(1, int((sec_col - first_col + 1)/2)):
                            first_int_rm = self._antenna.row_col_to_index(row, first_col + delta)
                            second_int_rm = self._antenna.row_col_to_index(row, sec_col - delta)
                            self.__add_file(first_rm, second_rm, first_int_rm, second_int_rm)

        elif self._antenna.quantity_columns == 1:
            column = 0
            for first_row in range(self._antenna.quantity_rows - 1):
                for sec_row in range(first_row + 1, self._antenna.quantity_rows):
                    first_rm = self._antenna.row_col_to_index(first_row, column)
                    second_rm = self._antenna.row_col_to_index(sec_row, column)
                    self.__add_file(first_rm, second_rm)

                    if ((sec_row - first_row) % 2) != 0 and sec_row - first_row > 1:
                        for delta in range(1, int((sec_row - first_row + 1)/2)):
                            first_int_rm = self._antenna.row_col_to_index(first_row + delta, column)
                            second_int_rm = self._antenna.row_col_to_index(sec_row - delta, column)
                            self.__add_file(first_rm, second_rm, first_int_rm, second_int_rm)

        else:
            for column in range(self._antenna.quantity_columns):
                # for row in range(1, min(self._antenna.quantity_rows, column + 2)):
                for first_row in range(1, self._antenna.quantity_rows - 1):
                    for second_row in range(first_row + 1, self._antenna.quantity_rows):
                        first_rm = self._antenna.row_col_to_index(first_row, column)
                        second_rm = self._antenna.row_col_to_index(second_row, column)
                        self.__add_file(first_rm, second_rm)

            for row in range(self._antenna.quantity_rows):
                # for column in range(1, min(self._antenna.quantity_columns, row + 2)):
                for first_column in range(1, self._antenna.quantity_columns - 1):
                    for second_column in range(first_column + 1, self._antenna.quantity_columns):
                        first_rm = self._antenna.row_col_to_index(row, first_column)
                        second_rm = self._antenna.row_col_to_index(row, second_column)
                        self.__add_file(first_rm, second_rm)


class DefaultBuilder(MatrixCalibratorBuilder):

    def __init__(self):
        super(DefaultBuilder, self).__init__()

    def _is_executable(self):
        return True

    def __add_file(self, rm, a, gain, phase, key):
        equation = [0] * self._antenna.quantity_columns * self._antenna.quantity_rows
        equation[rm] = 1
        a.append(equation)

        gain.append(Common.v2db(abs(self._equations[key])))

        phase.append(np.angle(self._equations[key], deg=True))

    def _add_equations(self):
        tx_equations = lambda x: [x[0], self._tx_a, self._tx_gain, self._tx_phase, x]
        rx_equations = lambda x: [x[1], self._rx_a, self._rx_gain, self._rx_phase, x]

        [self.__add_file(*tx_equations(key)) for key in self._equations if key[1] is None]
        [self.__add_file(*rx_equations(key)) for key in self._equations if key[0] is None]