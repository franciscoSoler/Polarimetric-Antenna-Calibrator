__author__ = 'fsoler'
import numpy as np
import matplotlib.pyplot as plt
import itertools
import cmath


def rad2deg(rad):
    return rad*180/cmath.pi


def deg2rad(deg):
    return deg*cmath.pi/180


class VisualComparator:

    def __init__(self):
        self.__figure_number = 0
        self.__att_delta = 5 / 2
        self.__ph_delta = 5.625 / 2

        # limits
        self.__upper_att_limit = None
        self.__upper_ph_limit = None
        self.__lower_att_limit = None
        self.__lower_ph_limit = None

    def __get_figure_number(self):
        self.__figure_number += 1
        return self.__figure_number

    @staticmethod
    def __rad2deg(radians):
        return list(map(lambda x: rad2deg(x), radians))

    @staticmethod
    def __deg2rad(degrees):
        return list(map(lambda x: deg2rad(x), degrees))

    @staticmethod
    def __format_signal(signal):
        decimals = 6
        f = lambda x: list(map(lambda y: eval(repr(round(y, decimals))), x))
        # f = lambda x: list(map(lambda y: eval(repr(round(y.real, decimals) + round(y.imag, decimals) * 1j)), x))
        return f(itertools.chain.from_iterable(signal))

    @staticmethod
    def __separate_signal_2_power_phase(signal):
        """
        deprecated
        :param signal:
        :return:
        """
        return zip(*map(lambda x: cmath.polar(x), signal))

    @staticmethod
    def __set_plot_environment(title, y_label, x_label, plots, legends):
        plt.title(title)
        plt.ylabel(y_label)
        plt.xlabel(x_label)
        plt.legend(plots, legends, bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, borderaxespad=0.,
                   prop={'size': 8})
        plt.grid(True)

    def __set_limits(self, ideal_att, ideal_ph):
        self.__upper_att_limit = np.array(ideal_att) + self.__att_delta
        self.__upper_ph_limit = np.array(ideal_ph) + self.__ph_delta
        self.__lower_att_limit = self.__upper_att_limit - 2 * self.__att_delta
        self.__lower_ph_limit = self.__upper_ph_limit - 2 * self.__ph_delta

    def compare_signals_against_ideal(self, power_one, phase_one, power_two, phase_two, id_power, id_phase, title=""):

        power = self.__format_signal(power_one)
        cal_power = self.__format_signal(power_two)
        ideal_power = self.__format_signal(id_power)
        phase = self.__format_signal(phase_one)
        cal_phase = self.__format_signal(phase_two)
        ideal_phase = self.__format_signal(id_phase)

        antennas = range(len(ideal_phase))
        f = lambda x: len(x) != len(ideal_phase)
        if f(power) or f(cal_power) or f(ideal_power) or f(phase) or f(cal_phase):
            raise Exception("both signals must contain the same length")

        self.__set_limits(ideal_power, ideal_phase)

        plt.figure(self.__get_figure_number())
        plt.subplot(211)
        p1, = plt.plot(antennas, power, "bo")
        plt.plot(antennas, power, "b")
        p2, = plt.plot(antennas, cal_power, "g^")
        plt.plot(antennas, cal_power, "g")
        p3, = plt.plot(antennas, ideal_power, "k")
        p4, = plt.plot(antennas, self.__upper_att_limit, "r--")
        plt.plot(antennas, self.__lower_att_limit, "r--")
        legends = ["non_cal_att", "cal_att", "ideal_att", "upper_lower_limit"]
        self.__set_plot_environment(title, "Attenuation [dB]", "RMs", [p1, p2, p3, p4], legends)

        plt.subplot(212)
        print("reception shits", phase)
        p1, = plt.plot(antennas, phase, "bo")
        plt.plot(antennas, phase, "b")
        p2, = plt.plot(antennas, cal_phase, "g^")
        plt.plot(antennas, cal_phase, "g")
        p3, = plt.plot(antennas, ideal_phase, "k")
        p4, = plt.plot(antennas, self.__upper_ph_limit, "r--")
        plt.plot(antennas, self.__lower_ph_limit, "r--")
        legends = ["non_cal_ph", "cal_ph", "ideal_ph", "upper_lower_limit"]
        self.__set_plot_environment(title, "Phase [deg]", "RMs", [p1, p2, p3, p4], legends)

    def compare_signals(self, signal_one, signal_two, title):
        """
        This function compares two signals, the first one should be the ideal
        :param signal_one:
        :param signal_two:
        :param title:
        :return:
        """

        first_signal = self.__format_signal(signal_one)
        second_signal = self.__format_signal(signal_two)

        if len(first_signal) != len(second_signal):
            raise Exception("both signals must contain the same length")

        antennas = range(len(first_signal))
        ideal_att, ideal_ph = zip(*map(lambda x: cmath.polar(x), first_signal))
        cal_att, cal_ph = zip(*map(lambda x: cmath.polar(x), second_signal))

        self.__set_limits(ideal_att, ideal_ph)

        plt.figure(self.__get_figure_number())
        plt.subplot(211)
        p2, = plt.plot(antennas, cal_att, "g^")
        plt.plot(antennas, cal_att, "g")
        p3, = plt.plot(antennas, ideal_att, "k")
        p4, = plt.plot(antennas, self.__upper_att_limit, "r--")
        plt.plot(antennas, self.__lower_att_limit, "r--")
        legends = ["s1_att", "s2_att", "upper_lower_limit"]
        self.__set_plot_environment(title, "Attenuation [dB]", "RMs", [p2, p3, p4], legends)

        plt.subplot(212)
        p2, = plt.plot(antennas, self.__rad2deg(cal_ph), "g^")
        # plt.plot(antennas, cal_ph, "g")
        p3, = plt.plot(antennas, self.__rad2deg(ideal_ph), "k")
        p4, = plt.plot(antennas, self.__upper_ph_limit, "r--")
        plt.plot(antennas, self.__lower_ph_limit, "r--")
        legends = ["s1_ph", "s2_ph", "upper_lower_limit"]
        self.__set_plot_environment(title, "Phase [deg]", "RMs", [p2, p3, p4], legends)

    def compare_patterns(self, angles, pattern_one, pattern_two, title):
        """
        plots both patterns in order to get a visual comparison.
        :param angles:
        :param pattern_one:
        :param pattern_two:
        :param title:
        :return:
        """

        f = lambda x: len(x) != len(angles)
        if f(pattern_one) or f(pattern_two):
            raise Exception("both patterns must have the same length")

        first_pattern = map(lambda x: abs(x), pattern_one)
        second_pattern = map(lambda x: abs(x), pattern_two)

        plt.figure(self.__get_figure_number())
        p2, = plt.plot(angles, first_pattern, "g")
        p3, = plt.plot(angles, second_pattern, "b")
        plt.title(title)
        plt.ylabel("Gain [dB]")
        plt.xlabel("Angle [deg]")
        plt.legend([p2, p3], ["pat1", "pat2"],
                   bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, borderaxespad=0., prop={'size': 8})
        plt.grid(True)

    def compare_patterns_against_ideal(self, angles, non_calibrated, calibrated, ideal, title):
        """
        plots the three signals in order to get a visual comparison.
        :param angles:
        :param non_calibrated:
        :param calibrated:
        :param ideal:
        :param title:
        :return:
        """

        f = lambda x: len(x) != len(angles)
        if f(non_calibrated) or f(calibrated) or f(ideal):
            raise Exception("the signals must have the same length")

        ideal_pattern = map(lambda x: abs(x), ideal)
        cal_pattern = map(lambda x: abs(x), calibrated)
        non_cal_pattern = map(lambda x: abs(x), non_calibrated)

        plt.figure(self.__get_figure_number())
        p1, = plt.plot(angles, ideal_pattern, "k")
        p2, = plt.plot(angles, cal_pattern, "g")
        p3, = plt.plot(angles, non_cal_pattern, "b")
        plt.title(title)
        plt.ylabel("Gain [dB]")
        plt.xlabel("Angle [deg]")
        plt.legend([p1, p2, p3], ["ideal pat", "cal pat", "non-cal pat"],
                   bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, borderaxespad=0., prop={'size': 8})
        plt.grid(True)

    @staticmethod
    def show_graphics():
        plt.show()
