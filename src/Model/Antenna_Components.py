__author__ = 'francisco'

import random
import numpy as np
from abc import ABCMeta, abstractmethod
import src.Utilities.Antenna_Common as AntennaCommon


class Component(object):
    __metaclass__ = ABCMeta

    def __init__(self, name, matrix, delta, add_errors):
        self._component = name
        self.__add_errors = add_errors
        self._delta = delta
        self._s_parameters = np.matrix([list(map(self.__add_error, row)) for row in matrix])

    def __add_error(self, param):
        """
        could be different
        :param param:
        :return:
        """
        return param * (1 + random.uniform(0, self._delta)) if self.__add_errors else param


class ExternalComponent(Component):
    __Available_modes = [AntennaCommon.Trans, AntennaCommon.Rec]

    def __init__(self, name, matrix, delta, add_errors):
        random.seed(None)
        self.__errors = add_errors
        super(ExternalComponent, self).__init__(name, matrix, delta, add_errors)

    def get_scattering_matrix(self, mode, leaf_index):
        if mode not in self.__Available_modes:
            raise Exception("the mode is not valid: ", mode)
        return self._get_tx_scattering_mtx(leaf_index) if mode == AntennaCommon.Trans else self._get_rx_scattering_mtx(
            leaf_index)

    def _add_error(self, param):
        """
        could be different
        :param param:
        :return:
        """
        return param * (1 + random.uniform(0, self._delta)) if self.__errors else param

    @abstractmethod
    def _get_tx_scattering_mtx(self, out_port_number):
        pass

    @abstractmethod
    def _get_rx_scattering_mtx(self, out_port_number):
        pass


class InternalComponent(Component):

    def __init__(self, name, matrix, delta, add_errors):
        super(InternalComponent, self).__init__(name, matrix, delta, add_errors)

    @property
    def s_parameters(self):
        return self._s_parameters


class Cable(ExternalComponent):

    def __init__(self, attenuation=0, length=1, wavelength=1, delta=0, add_errors=False):
        self.__attenuation = attenuation
        self.__length = length
        self.__wavelength = wavelength

        sij = np.exp((attenuation + 2j*np.pi/self.__wavelength) * self.__length)
        super(Cable, self).__init__(AntennaCommon.Cable, [[0, sij], [sij, 0]], delta, add_errors)

        self.__rx_items = [[[1, 1], [1, 0]], [[0, 0], [1, 0]]]

    def _get_tx_scattering_mtx(self, out_port_number):
        return self._s_parameters

    def _get_rx_scattering_mtx(self, out_port_number):
        return self._s_parameters[list(zip(*self.__rx_items))]


class PowerSplitterCombiner(ExternalComponent):

    def __init__(self, quantity_ports, delta=0, add_errors=False):
        ports = quantity_ports - 1
        val = np.sqrt(1/ports)
        s = np.matrix([[val] * quantity_ports for _ in range(quantity_ports)])
        s[range(quantity_ports), range(quantity_ports)] = 0
        super(PowerSplitterCombiner, self).__init__("{0}1-{1}".format(AntennaCommon.Psc, ports), s.tolist(), delta,
                                                    add_errors)

        self.__tx_items = lambda x: [[[0, 0], [0, x]], [[x, x], [0, x]]]
        self.__rx_items = lambda x: [[[x, x], [x, 0]], [[0, 0], [x, 0]]]

    def _get_tx_scattering_mtx(self, out_port_number):
        return self._s_parameters[list(zip(*self.__tx_items(out_port_number + 1)))]

    def _get_rx_scattering_mtx(self, out_port_number):
        return self._s_parameters[list(zip(*self.__rx_items(out_port_number + 1)))]


class TransmitterReceiverModule(ExternalComponent):

    __Available_modes = [AntennaCommon.Trans, AntennaCommon.Rec]

    def __init__(self, tx_gain, rx_gain, tx_phase_shift, rx_phase_shift, delta=0, add_errors=False):
        self.__tx_amplifier = Amplifier(tx_gain, delta, add_errors)
        self.__rx_amplifier = Amplifier(rx_gain, delta, add_errors)
        self.__tx_phase_shifter = PhaseShifter(tx_phase_shift, delta, add_errors)
        self.__rx_phase_shifter = PhaseShifter(rx_phase_shift, delta, add_errors)
        super(TransmitterReceiverModule, self).__init__(AntennaCommon.Trm, [[]], delta, add_errors)

    def _get_tx_scattering_mtx(self, out_port_number):
        s2t = AntennaCommon.s2t_parameters
        t2s = AntennaCommon.t2s_parameters
        return t2s(s2t(self.__tx_amplifier) * s2t(self.__tx_phase_shifter))

    def _get_rx_scattering_mtx(self, out_port_number):
        s2t = AntennaCommon.s2t_parameters
        t2s = AntennaCommon.t2s_parameters
        return t2s(s2t(self.__rx_amplifier) * s2t(self.__rx_phase_shifter))

    def change_gain(self, delta_gain, mode):
        if mode not in self.__Available_modes:
            raise Exception("the mode is not valid:", mode)

        phase = np.angle(delta_gain)
        gain = abs(delta_gain)
        if mode == AntennaCommon.Trans:
            self.__tx_amplifier.change_gain(gain)
            self.__tx_phase_shifter.change_phase(phase)
        else:
            self.__rx_amplifier.change_gain(gain)
            self.__rx_phase_shifter.change_phase(phase)


class Circulator(ExternalComponent):

    def __init__(self, delta=0, add_errors=False):
        super(Circulator, self).__init__(AntennaCommon.Circulator, [[0, 0, 1], [1, 0, 0], [0, 1, 0]], delta, add_errors)

        self.__tx_items = [[[0, 0], [0, 1]], [[1, 1], [0, 1]]]
        self.__rx_items = [[[1, 1], [1, 2]], [[2, 2], [1, 2]]]

    def _get_tx_scattering_mtx(self, out_port_number):
        return self._s_parameters[list(zip(*self.__tx_items))]

    def _get_rx_scattering_mtx(self, out_port_number):
        return self._s_parameters[list(zip(*self.__rx_items))]


class RadiantModule(ExternalComponent):

    def __init__(self, row, col, delta=0, add_errors=False):
        super(RadiantModule, self).__init__(AntennaCommon.Rm, [[0]], delta, add_errors)
        self.__row = row
        self.__col = col

    def _get_tx_scattering_mtx(self, out_port_number):
        return self._s_parameters

    def _get_rx_scattering_mtx(self, out_port_number):
        return self._s_parameters

    @property
    def row(self):
        return self.__row

    @property
    def col(self):
        return self.__col

    def get_coordinates(self):
        return self.__row, self.__col


class Amplifier(InternalComponent):

    def __init__(self, gain, delta, add_errors):
        self.__gain = gain
        super(Amplifier, self).__init__(AntennaCommon.Amplifier, [[0, 0], [self.__gain, 0]], delta, add_errors)

    @property
    def gain(self):
        return self.__gain

    def change_gain(self, gain_shift):
        self.__gain *= gain_shift


class PhaseShifter(InternalComponent):

    def __init__(self, phase, delta, add_errors):
        self.__phase = phase
        sij = np.exp(-1j*AntennaCommon.deg2rad(self.__phase))
        super(PhaseShifter, self).__init__(AntennaCommon.Ph_shifter, [[0, sij], [sij, 0]], delta, add_errors)

    @property
    def phase(self):
        return self.__phase

    def change_phase(self, phase_shift):
        self.__phase += phase_shift
        self._s_parameters[1, 0] *= np.exp(1j*AntennaCommon.deg2rad(phase_shift))
        self._s_parameters[0, 1] *= np.exp(1j*AntennaCommon.deg2rad(phase_shift))