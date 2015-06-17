__author__ = 'francisco'
from abc import ABCMeta, abstractmethod
import src.Utilities.Antenna_Common as AntennaCommon
import numpy as np
import random


class ScatteringParametersHandler(object):
    __metaclass__ = ABCMeta

    def __init__(self, check_component):
        self._successor = None
        self._wavelength = 0
        self._delta = 0
        self.__errors = False
        self.__check_component = check_component

    def initialize(self, wavelength, delta=0, add_errors=True):
        random.seed(None)
        self._wavelength = wavelength
        self._delta = delta
        self.__errors = add_errors

    def _add_error(self, param):
        return param * (1 + random.uniform(0, self._delta)) if self.__errors else param

    @abstractmethod
    def _get_scattering_matrix(self, length):
        pass

    def get_scattering_matrix(self, component, length=None):
        """

        :param component:
        :return:
            the scattering parameters matrix
        """
        if self.__check_component(component):
            return self._get_scattering_matrix(length)
        elif self._successor is not None:
            return self._successor.get_scattering_matrix(component, length)
        else:
            raise Exception("the component is not valid: ", component)

    def set_successor(self, successor):
        self._successor = successor


class CableScatteringParameters(ScatteringParametersHandler):

    def __init__(self):
        super(CableScatteringParameters, self).__init__(AntennaCommon.is_cable)

    def _get_scattering_matrix(self, length):
        alpha = self._add_error(self._delta)
        sij = lambda x: np.exp((alpha + 2j*np.pi/self._wavelength) * x)
        s_parameters = lambda x: [[0, sij(x)], [sij(x), 0]]
        return [s_parameters(len_i) for len_i in length] if isinstance(length, list) else s_parameters(length)


class PscScatteringParameters(ScatteringParametersHandler):

    def __init__(self):
        super(PscScatteringParameters, self).__init__(AntennaCommon.is_psc)

    def _get_scattering_matrix(self, length):
        val = np.sqrt(1/length)
        add_errors = lambda x: list(map(lambda y: self._add_error(y), x))
        s = np.matrix([[val] * length for _ in range(length)])
        s[range(length), range(length)] = 0
        return [add_errors(row) for row in s.tolist()]


class CirculatorScatteringParameters(ScatteringParametersHandler):

    def __init__(self):
        super(CirculatorScatteringParameters, self).__init__(AntennaCommon.is_circulator)

    def _get_scattering_matrix(self, length):
        add_errors = lambda x: list(map(lambda y: self._add_error(y), x))
        return [add_errors([0, 0, 1]), add_errors([1, 0, 1]), add_errors([0, 1, 0])]


class TrmScatteringParameters(ScatteringParametersHandler):

    def __init__(self):
        super(TrmScatteringParameters, self).__init__(AntennaCommon.is_trm)

    def _get_scattering_matrix(self, gain_shift):
        add_errors = lambda x: list(map(lambda y: self._add_error(y), x))
        return [add_errors([0, 0, gain_shift]), add_errors([gain_shift, 0, 0]), add_errors([gain_shift, 0, 0])]


class RmScatteringParameters(ScatteringParametersHandler):

    def __init__(self):
        super(RmScatteringParameters, self).__init__(AntennaCommon.is_rm)

    def _get_scattering_matrix(self, length):
        return [[self._add_error(0)]]
