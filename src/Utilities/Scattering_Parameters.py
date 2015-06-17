__author__ = 'francisco'
from abc import ABCMeta, abstractmethod
import src.Utilities.Antenna_Common as AntennaCommon
import numpy as np
import random


class ScatteringParametersHandler(object):
    __metaclass__ = ABCMeta

    def __init__(self, check_component):
        self._successor = None
        self._delta = 0
        self.__errors = False
        self.__check_component = check_component

    def initialize(self, delta=0, add_errors=True):
        random.seed(None)
        self._delta = delta
        self.__errors = add_errors

    def _add_error(self, param):
        return param * (1 + random.uniform(0, self._delta)) if self.__errors else param

    @abstractmethod
    def _get_scattering_matrix(self, attributes):
        pass

    def get_scattering_matrix(self, component, attributes=None):
        """

        :param component:
        :return:
            the scattering parameters matrix
        """
        if self.__check_component(component):
            return self._get_scattering_matrix(attributes)
        elif self._successor is not None:
            return self._successor.get_scattering_matrix(component, attributes)
        else:
            raise Exception("the component is not valid: ", component)

    def set_successor(self, successor):
        self._successor = successor


class CableScatteringParameters(ScatteringParametersHandler):

    def __init__(self):
        super(CableScatteringParameters, self).__init__(AntennaCommon.is_cable)

    def __get_scattering_matrix(self, attenuation, wavelength, length):
        """
        alpha = self._add_error(self._delta)
        sij = lambda x: np.exp((alpha + 2j*np.pi/self._wavelength) * x)
        s_parameters = lambda x: [[0, sij(x)], [sij(x), 0]]
        return [s_parameters(len_i) for len_i in attributes] if isinstance(attributes, list) else s_parameters(
        attributes)
        """
        sij = np.exp((attenuation + 2j*np.pi/wavelength) * length)
        return [[0, self._add_error(sij)], [self._add_error(sij), 0]]

    def _get_scattering_matrix(self, attributes):
        return self.__get_scattering_matrix(*attributes)


class PscScatteringParameters(ScatteringParametersHandler):

    def __init__(self):
        super(PscScatteringParameters, self).__init__(AntennaCommon.is_psc)

    def __get_scattering_matrix(self, output_ports):
        ports = output_ports + 1
        val = np.sqrt(1/ports)
        add_errors = lambda x: list(map(self._add_error, x))
        s = np.matrix([[val] * ports for _ in range(ports)])
        s[range(ports), range(ports)] = 0
        return [add_errors(row) for row in s.tolist()]

    def _get_scattering_matrix(self, attributes):
        return self.__get_scattering_matrix(*attributes)


class CirculatorScatteringParameters(ScatteringParametersHandler):

    def __init__(self):
        super(CirculatorScatteringParameters, self).__init__(AntennaCommon.is_circulator)

    def _get_scattering_matrix(self, attributes):
        add_errors = lambda x: list(map(self._add_error, x))
        return [add_errors([0, 0, 1]), add_errors([1, 0, 0]), add_errors([0, 1, 0])]


class TrmScatteringParameters(ScatteringParametersHandler):

    def __init__(self):
        super(TrmScatteringParameters, self).__init__(AntennaCommon.is_trm)

    def __get_scattering_matrix(self, gain, phase_shift):
        # add_errors = lambda x: list(map(self._add_error, x))
        # return [add_errors([0, 0, gain_shift]), add_errors([gain_shift, 0, 0]), add_errors([0, 0, 0])]
        sij = gain * np.exp(1j*AntennaCommon.deg2rad(phase_shift))
        add_errors = lambda x: list(map(self._add_error, x))
        param = [add_errors([0, 0, sij]), add_errors([sij, 0, 0]), add_errors([0, 0, 0])]
        print("param", param)
        return param

    def _get_scattering_matrix(self, attributes):
        return self.__get_scattering_matrix(*attributes)


class RmScatteringParameters(ScatteringParametersHandler):

    def __init__(self):
        super(RmScatteringParameters, self).__init__(AntennaCommon.is_rm)

    def _get_scattering_matrix(self, attributes):
        return [[self._add_error(0)]]
