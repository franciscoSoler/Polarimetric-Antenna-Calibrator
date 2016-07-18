__author__ = 'fsoler'
import numpy as np
import matplotlib.pyplot as plt
import xlsxwriter
import itertools
import cmath
import os


def rad2deg(rad):
    return rad*180/cmath.pi


def deg2rad(deg):
    return deg*cmath.pi/180


class VisualComparator:

    def __init__(self, save_files, graphics_to_csv):
        plt.close('all')
        self.__save_files = save_files
        self.__to_csv = graphics_to_csv
        self.__figure_number = 0
        self.__att_delta = 5 / 2
        self.__ph_delta = 5.625 / 2

        # limits
        self.__upper_att_limit = None
        self.__upper_ph_limit = None
        self.__lower_att_limit = None
        self.__lower_ph_limit = None
        self.__path_to_save = '/media/francisco/Datos/francisco/Documents/mutual/written/thesis/gfx'
        self.__path_to_save_csv = '/media/francisco/Datos/francisco/Documents/mutual/written/thesis/gfx/csv'

    def __get_figure_number(self):
        self.__figure_number += 1
        return self.__figure_number

    def __save_plots(self, filename, plt):
        plt.tight_layout()
        if self.__save_files:
            plt.savefig(os.path.join(self.__path_to_save, filename + ".png"), bbox_inches='tight')

    def __save_into_csv(self, filename, data):
        if self.__to_csv:
            workbook = xlsxwriter.Workbook(os.path.join(self.__path_to_save_csv, filename + '.xlsx'))
            worksheet = workbook.add_worksheet()
            worksheet.write_row(0, 0, ['ideal', 'pre_cal', 'cal'])
            [worksheet.write_column(1, i, d) for i, d in enumerate(data)]
            #np.savetxt(os.path.join(self.__path_to_save_csv, filename + '.csv'), data, delimiter=',', fmt='%f', header='Created by Numpy')

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
    def __set_plot_environment(title, y_label, x_label, locc=None):
        plt.title(title)
        plt.ylabel(y_label)
        plt.xlabel(x_label)

        p1, = plt.plot([], [], label="non-cal", color="b", marker="o")
        p2, = plt.plot([], [], label="cal", color="g", marker="^", markersize=8)
        p3, = plt.plot([], [], label="ideal", color="k", marker="*", linewidth=2, markersize=15)

        if locc is None:
            plt.legend(handles=[p1, p2, p3], bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, borderaxespad=0.,
                   prop={'size': 8})
        else:
            plt.legend(handles=[p1, p2, p3], loc=locc)
        plt.grid(True)

    def __set_limits(self, ideal_att, ideal_ph):
        self.__upper_att_limit = np.array(ideal_att) + self.__att_delta
        self.__upper_ph_limit = np.array(ideal_ph) + self.__ph_delta
        self.__lower_att_limit = self.__upper_att_limit - 2 * self.__att_delta
        self.__lower_ph_limit = self.__upper_ph_limit - 2 * self.__ph_delta

    def compare_signals_against_ideal(self, power_one, phase_one, power_two, phase_two, id_power, id_phase, title="", filename=""):

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
        plt.plot(antennas, power, "b", linewidth=2)
        p2, = plt.plot(antennas[::2], cal_power[::2], "g^", markersize=8)
        p3, = plt.plot(antennas[::4], ideal_power[::4], color="k", marker="*", linewidth=2, markersize=15)
        plt.plot(antennas, cal_power, "g", linewidth=2)

        self.__set_plot_environment(title, "Power [dBm]", "ERs", 4)

        plt.subplot(212)
        plt.plot(antennas, phase, "bo")
        plt.plot(antennas, phase, "b", linewidth=2)
        plt.plot(antennas, ideal_phase, "k", linewidth=2)
        plt.plot(antennas[::2], cal_phase[::2], "g^", markersize=8)
        plt.plot(antennas[::5], ideal_phase[::5], "k*", markersize=15)
        plt.plot(antennas, cal_phase, "g", linewidth=2)

        self.__set_plot_environment(title, "Phase [deg]", "ERs", 4)
        self.__save_plots(filename, plt)
        self.__save_into_csv(filename + 'Power', [ideal_power, power, cal_power])
        self.__save_into_csv(filename + 'Phase', [ideal_phase, phase, cal_phase])

    def compare_signals(self, power_one, phase_one, power_two, phase_two, title=""):
        """
        This function compares two signals, the first one should be the ideal
        :param signal_one:
        :param signal_two:
        :param title:
        :return:
        """
        power = self.__format_signal(power_one)
        cal_power = self.__format_signal(power_two)

        phase = self.__format_signal(phase_one)
        cal_phase = self.__format_signal(phase_two)

        antennas = range(len(cal_phase))
        f = lambda x: len(x) != len(cal_phase)
        if f(power) or f(cal_power) or f(phase):
            raise Exception("both signals must contain the same length")

        plt.figure(self.__get_figure_number())
        plt.subplot(211)
        p1, = plt.plot(antennas, power, "bo")
        plt.plot(antennas, power, "b")
        p2, = plt.plot(antennas, cal_power, "g^")
        plt.plot(antennas, cal_power, "g")
        legends = ["non_cal_power", "cal_power"]
        self.__set_plot_environment(title, "Power [dB]", "ERs", locc=4)

        plt.subplot(212)
        p1, = plt.plot(antennas, phase, "bo")
        plt.plot(antennas, phase, "b")
        p2, = plt.plot(antennas, cal_phase, "g^")
        plt.plot(antennas, cal_phase, "g")
        legends = ["non_cal_ph", "cal_ph"]
        self.__set_plot_environment("", "Phase [deg]", "ERs", locc=4)

    def compare_patterns(self, angles, pattern_one, pattern_two, title, filename):
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

        plt.figure(self.__get_figure_number())
        p2, = plt.plot(angles, pattern_one, "b", linewidth=2)
        p3, = plt.plot(angles, pattern_two, "k", linewidth=2)
        plt.plot(angles, pattern_two+2.5, "r--", linewidth=1)
        plt.plot(angles, pattern_two-2.5, "r--", linewidth=1)
        plt.title(title)
        plt.ylabel("Gain [dB]")
        plt.xlabel("Angle [deg]")
        plt.grid(True)

    def compare_patterns_against_ideal(self, angles, non_calibrated, calibrated, ideal, title, filename):
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

        cut = -10
        p_ideal = [p if p > cut else cut for p in ideal]
        p_calibrated = [p if p > cut else cut for p in calibrated]
        p_non_calibrated = [p if p > cut else cut for p in non_calibrated]

        delta = 1
        sup = [p if p > cut else cut for p in ideal + delta]
        inf = [p if p > cut else cut for p in ideal - delta]

        plt.figure(self.__get_figure_number())
        plt.plot(angles, p_ideal, "k", label="ideal", linewidth=2)
        plt.plot(angles, p_calibrated, "g", label="cal", linewidth=2)
        plt.plot(angles, p_non_calibrated, "b", label="non-cal", linewidth=2)

        plt.plot(angles[::125], p_ideal[::125], "k*", markersize=13)
        plt.plot(angles[::150], p_calibrated[::150], "g^", markersize=10)
        plt.plot(angles[::175], p_non_calibrated[::175], "bo", markersize=8)

        plt.plot(angles, sup, "r--", linewidth=1)
        plt.plot(angles, inf, "r--", linewidth=1)

        self.__set_plot_environment(title, "Gain [dB]", "Angle [deg]", 4)
        self.__save_plots(filename, plt)
        self.__save_into_csv(filename, [p_ideal, p_non_calibrated, p_calibrated])

    @staticmethod
    def show_graphics():
        plt.tight_layout()
        plt.show()
