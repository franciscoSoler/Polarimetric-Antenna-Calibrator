__author__ = 'francisco'
from abc import ABCMeta, abstractmethod
import Utilities.Antenna_Common as AntennaCommon
import numpy as np


class ScatteringParametersHandler(object):
    __metaclass__ = ABCMeta

    def __init__(self, check_component):
        self._successor = None
        self.__errors = False
        self.__delta_gain = 0
        self.__delta_phase = 0
        self.__check_component = check_component

    def initialize(self, add_errors=False, delta=[0,0]):

        self.__errors = add_errors
        self.__delta_gain = delta[0]
        self.__delta_phase = delta[1]

    def _add_error(self, param):
        if param == 0:
            return param
        module = np.random.normal(abs(param), self.__delta_gain) if self.__delta_gain else abs(param)
        angle = AntennaCommon.deg2rad(np.random.normal(np.angle(param, deg=True), self.__delta_phase)) \
                    if self.__delta_phase else np.angle(param)
        return module * np.exp(1j * angle) if self.__errors else param

    @abstractmethod
    def _get_scattering_matrix(self, attributes):
        pass

    def get_scattering_matrix(self, component, attributes=None):
        if self.__check_component(component):
            return self._get_scattering_matrix(attributes)
        elif self._successor is not None:
            return self._successor.get_scattering_matrix(component, attributes)
        else:
            raise Exception("the component is not valid: {}".format(component))

    def set_successor(self, successor):
        self._successor = successor


class CableScatteringParameters(ScatteringParametersHandler):

    def __init__(self):
        super(CableScatteringParameters, self).__init__(AntennaCommon.is_cable)

    def __get_scattering_matrix(self, attenuation, length, wavelength):
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
    __Dead_gain = 0.00000001

    def __init__(self):
        super(TrmScatteringParameters, self).__init__(AntennaCommon.is_trm)

    def __get_scattering_matrix(self, gain, phase_shift, is_dead):
        sij = self.__Dead_gain if is_dead else gain * np.exp(1j*AntennaCommon.deg2rad(phase_shift))
        add_errors = lambda x: list(map(self._add_error, x))
        return [add_errors([0, 0, sij]), add_errors([sij, 0, 0]), add_errors([0, 0, 0])]

    def _get_scattering_matrix(self, attributes):
        return self.__get_scattering_matrix(*attributes)


class RmScatteringParameters(ScatteringParametersHandler):

    def __init__(self):
        super(RmScatteringParameters, self).__init__(AntennaCommon.is_rm)

    def _get_scattering_matrix(self, attributes):
        return [[self._add_error(0)]]
