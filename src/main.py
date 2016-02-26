#!/usr/bin/python3.4
__author__ = 'fsoler'

import sys
import Visual_Comparator.Visual_Comparator as VisualComparator
import Pattern_Generator.Pattern_Generator as PatternGenerator
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import Controllers.Antenna_Creator as AntennaCreator
import Controllers.Antenna_Calibrator as AntennaCalibrator
import Utilities.Antenna_Common as Common
import json


def generate_pattern():
    """
    l_band_freq = 1275MHz
    s_band_freq = 2GHz
    x_band_freq = 8 GHz
    :return:
    """
    generator = PatternGenerator.PatternGenerator(1275000000, 0.127, 0.110)
    weights = [[1 for _ in range(20)]for _ in range(1)]

    angles, power = generator.generate_pattern(weights, [-30, 30], 0)
    plt.plot(angles, 20*np.log10([abs(p) for p in power]))
    plt.show()


def create_antenna2(filename):
    with open("base_rfdn") as f:
        with open(filename + "_rfdn", "w") as g:
            g.write(f.read())

    with open("base_panel") as f:
        with open(filename + "_panel", "w") as g:
            g.write(f.read())


def append_signal_into_signals(signals, power_signal, phase_signal):
    signals.append(power_signal)
    signals.append(phase_signal)


def compare_estimated_gains_against_real(calibrator, visual_comparator, title):
    tx_signals = []
    rx_signals = []

    tx_power, tx_phase = calibrator.get_transmission_power()
    rx_power, rx_phase = calibrator.get_reception_power()
    append_signal_into_signals(tx_signals, tx_power, tx_phase)
    append_signal_into_signals(rx_signals, rx_power, rx_phase)

    tx_ant_power, tx_ant_phase, rx_ant_power, rx_ant_phase = calibrator.get_antenna_gain_paths()
    append_signal_into_signals(tx_signals, tx_ant_power, tx_ant_phase)
    append_signal_into_signals(rx_signals, rx_ant_power, rx_ant_phase)

    visual_comparator.compare_signals(*tx_signals, title="{}: Tx chain".format(title))
    visual_comparator.compare_signals(*rx_signals, title="{}: Rx chain".format(title))
    """
    print("{}:".format(title), "Tx antenna power    ", tx_ant_power)
    print("{}:".format(title), "Tx estimated power  ", tx_power)
    print("{}:".format(title), "Tx abs power error  ", (np.array(tx_ant_power) - np.array(tx_power)).tolist())
    """
    print("{}:".format(title), "Tx antenna phase    ", tx_ant_phase)
    print("{}:".format(title), "Tx estimated phase  ", tx_phase)
    print("{}:".format(title), "Tx abs phase error  ", (np.array(tx_ant_phase) - np.array(tx_phase)).tolist())
    print("")
    """
    print("{}:".format(title), "Rx antenna power    ", rx_ant_power)
    print("{}:".format(title), "Rx estimated power  ", rx_power)
    print("{}:".format(title), "Rx abs power error  ", (np.array(rx_ant_power) - np.array(rx_power)).tolist())
    """
    print("{}:".format(title), "Rx antenna phase    ", rx_ant_phase)
    print("{}:".format(title), "Rx estimated phase  ", rx_phase)
    print("{}:".format(title), "Rx abs phase error  ", (np.array(rx_ant_phase) - np.array(rx_phase)).tolist())
    print("")


class Simulator:
    def __init__(self):
        self.__config = ""
        with open("configurationFile") as f:
            self.__config = json.load(f)

    def __build_sequence(self, wavelength):
        seq = []
        sequence = self.__config[Common.Conf_ant][Common.Conf_comp_seq]
        for component in sequence:
            params = [l[1] for l in sorted(self.__config[Common.Conf_component][component].items())]
            component_type = params.pop()
            if component_type == Common.Cable:
                params.append(wavelength)
            seq.append([component_type, params])
        return seq

    def __build_errors(self, sequence):
        seq = []
        for component in sequence:
            std_deviation = self.__config[Common.Conf_std_err][component]
            seq.append([component, std_deviation])
        return seq

    def __build_component_errors(self):
        return self.__build_errors(self.__config[Common.Conf_ant][Common.Conf_errors])

    def __build_calibration_errors(self):
        return self.__build_errors(self.__config[Common.Conf_cal_param][Common.Conf_errors])

    def __get_desired_phases(self):
        quantity_rows = self.__config[Common.Conf_ant][Common.Conf_qtty_rows]
        quantity_columns = self.__config[Common.Conf_ant][Common.Conf_qtty_cols]

        row_steering = self.__config[Common.Conf_in_param][Common.Conf_row_steer]
        column_steering = self.__config[Common.Conf_in_param][Common.Conf_col_steer]

        f = lambda row, col: np.mod(row * row_steering + col * column_steering + 180, 360) - 180
        return [[f(row, col) for col in range(quantity_columns)] for row in range(quantity_rows)]

    def __create_antenna(self):
        quantity_columns = self.__config[Common.Conf_ant][Common.Conf_qtty_cols]
        quantity_rows = self.__config[Common.Conf_ant][Common.Conf_qtty_rows]
        row_separation = self.__config[Common.Conf_ant][Common.Conf_vert_sep]
        col_separation = self.__config[Common.Conf_ant][Common.Conf_horiz_sep]
        row_steering = self.__config[Common.Conf_in_param][Common.Conf_row_steer]
        column_steering = self.__config[Common.Conf_in_param][Common.Conf_col_steer]
        filename = self.__config[Common.Conf_ant][Common.Conf_filename]
        dead_trms = self.__config[Common.Conf_ant][Common.Conf_dead_trm]
        c = 299792458       # [m/seg]
        f = self.__config[Common.Conf_in_param][Common.Conf_freq]      # [Hz]
        wavelength = c/f    # [m]

        sequence_items = self.__build_sequence(wavelength)
        component_errors = self.__build_component_errors()

        creator = AntennaCreator.AntennaCreator(quantity_rows, quantity_columns, row_separation, col_separation)
        creator.add_errors(component_errors)
        creator.create_structure(filename, sequence_items, row_steering, column_steering, dead_trms)

    def __create_calibrator(self):
        in_power = self.__config[Common.Conf_in_param][Common.Conf_power]
        in_phase = self.__config[Common.Conf_in_param][Common.Conf_phase]
        row_steering = self.__config[Common.Conf_in_param][Common.Conf_row_steer]
        column_steering = self.__config[Common.Conf_in_param][Common.Conf_col_steer]
        row_separation = self.__config[Common.Conf_ant][Common.Conf_vert_sep]
        col_separation = self.__config[Common.Conf_ant][Common.Conf_horiz_sep]
        filename = self.__config[Common.Conf_ant][Common.Conf_filename]

        calibration_errors = self.__build_calibration_errors()

        calibrator = AntennaCalibrator.MutualCalibrator(in_power, in_phase, row_steering, column_steering,
                                                        row_separation, col_separation, filename)
        # calibrator.add_calibration_errors(calibration_errors)
        calibrator.generate_cal_paths(AntennaCalibrator.every_one_to_one_path_strategy)
        """
        calibrator = AntennaCalibrator.ClassicCalibrator(in_power, in_phase, row_separation, col_separation, filename)
        # calibrator.add_calibration_errors(calibration_errors)
        """
        return calibrator

    def __calibrate_antenna(self, calibrator):
        desired_tx_power = self.__config[Common.Conf_cal_param][Common.Conf_id_tx_power]
        desired_rx_power = self.__config[Common.Conf_cal_param][Common.Conf_id_rx_power]

        desired_phase = self.__get_desired_phases()

        desired_signals = [desired_tx_power, desired_phase, desired_rx_power, desired_phase]
        calibrator.calibrate_antenna(*desired_signals)

    def __compare_estimated_gains_against_ideal(self, calibrator, visual_comparator, title):
        quantity_rows = self.__config[Common.Conf_ant][Common.Conf_qtty_rows]
        quantity_columns = self.__config[Common.Conf_ant][Common.Conf_qtty_cols]

        desired_tx_power = self.__config[Common.Conf_cal_param][Common.Conf_id_tx_power]
        desired_rx_power = self.__config[Common.Conf_cal_param][Common.Conf_id_rx_power]

        tx_signals = []
        rx_signals = []

        tx_ant_power, tx_ant_phase, rx_ant_power, rx_ant_phase = calibrator.get_antenna_gain_paths()
        append_signal_into_signals(tx_signals, tx_ant_power, tx_ant_phase)
        append_signal_into_signals(rx_signals, rx_ant_power, rx_ant_phase)

        tx_power, tx_phase = calibrator.get_transmission_power()
        rx_power, rx_phase = calibrator.get_reception_power()
        append_signal_into_signals(tx_signals, tx_power, tx_phase)
        append_signal_into_signals(rx_signals, rx_power, rx_phase)

        tx_ideal_power = [[desired_tx_power] * quantity_columns] * quantity_rows
        rx_ideal_power = [[desired_rx_power] * quantity_columns] * quantity_rows
        ideal_phase = self.__get_desired_phases()
        """
        tx_ideal_phase = [[desired_tx_phase] * quantity_columns] * quantity_rows
        rx_ideal_phase = [[desired_rx_phase] * quantity_columns] * quantity_rows

        append_signal_into_signals(tx_signals, tx_ideal_power, tx_ideal_phase)
        append_signal_into_signals(rx_signals, rx_ideal_power, rx_ideal_phase)
        """
        append_signal_into_signals(tx_signals, tx_ideal_power, ideal_phase)
        append_signal_into_signals(rx_signals, rx_ideal_power, ideal_phase)

        visual_comparator.compare_signals_against_ideal(*tx_signals, title="{}: Tx chain".format(title))
        visual_comparator.compare_signals_against_ideal(*rx_signals, title="{}: Rx chain".format(title))

    def __remove_antenna(self):
        filename = self.__config[Common.Conf_ant][Common.Conf_filename]
        for filename in glob.glob(filename + "_*"):
            os.remove(filename)

    def run(self):
        visual_comparator = VisualComparator.VisualComparator()

        self.__create_antenna()

        calibrator = self.__create_calibrator()
        compare_estimated_gains_against_real(calibrator, visual_comparator, "BEFORE CALIBRATION")

        self.__calibrate_antenna(calibrator)

        compare_estimated_gains_against_real(calibrator, visual_comparator, "AFTER CALIBRATION")
        self.__compare_estimated_gains_against_ideal(calibrator, visual_comparator, "")

        visual_comparator.show_graphics()

        #self.__remove_antenna()

        return 0


if __name__ == "__main__":
    sys.exit(Simulator().run())
    # sys.exit(compare_signals())
    # sys.exit(generate_pattern())
