#!/usr/bin/python3.4
__author__ = 'fsoler'

import sys
import Visual_Comparator.Visual_Comparator as VisualComparator
import Pattern_Generator.Pattern_Generator as PatternGenerator
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import re
import Controllers.Antenna_Creator as AntennaCreator
import Controllers.Antenna_Calibrator as AntennaCalibrator
import Utilities.Antenna_Common as Common
import json


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
        self.__calibrator = ""
        self.__properties = "../written/thesis/capitulos/capitulo06.tex"
        #self.__properties = "testFindReplace"
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

        row_separation = self.__config[Common.Conf_ant][Common.Conf_vert_sep]
        col_separation = self.__config[Common.Conf_ant][Common.Conf_horiz_sep]

        row_steering = self.__config[Common.Conf_in_param][Common.Conf_row_steer]
        column_steering = self.__config[Common.Conf_in_param][Common.Conf_col_steer]

        freq = self.__config[Common.Conf_in_param][Common.Conf_freq]      # [Hz]

        return Common.obtain_shift_phases(column_steering, row_steering, quantity_columns, quantity_rows,
                                          col_separation, row_separation, freq)

    def create_antenna(self):
        quantity_columns = self.__config[Common.Conf_ant][Common.Conf_qtty_cols]
        quantity_rows = self.__config[Common.Conf_ant][Common.Conf_qtty_rows]
        row_separation = self.__config[Common.Conf_ant][Common.Conf_vert_sep]
        col_separation = self.__config[Common.Conf_ant][Common.Conf_horiz_sep]
        row_steering = self.__config[Common.Conf_in_param][Common.Conf_row_steer]
        column_steering = self.__config[Common.Conf_in_param][Common.Conf_col_steer]
        filename = self.__config[Common.Conf_ant][Common.Conf_filename]
        dead_trms = self.__config[Common.Conf_ant][Common.Conf_dead_trm]
        f = self.__config[Common.Conf_in_param][Common.Conf_freq]      # [Hz]
        wavelength = Common.C/f    # [m]

        sequence_items = self.__build_sequence(wavelength)
        component_errors = self.__build_component_errors()

        creator = AntennaCreator.AntennaCreator(quantity_rows, quantity_columns, row_separation, col_separation)
        creator.add_errors(component_errors)
        creator.create_structure(filename, sequence_items, row_steering, column_steering, dead_trms)

    def __create_calibrator(self, calibr):
        self.__calibrator = calibr

        in_power = self.__config[Common.Conf_in_param][Common.Conf_power]
        in_phase = self.__config[Common.Conf_in_param][Common.Conf_phase]
        row_steering = self.__config[Common.Conf_in_param][Common.Conf_row_steer]
        column_steering = self.__config[Common.Conf_in_param][Common.Conf_col_steer]
        row_separation = self.__config[Common.Conf_ant][Common.Conf_vert_sep]
        col_separation = self.__config[Common.Conf_ant][Common.Conf_horiz_sep]
        filename = self.__config[Common.Conf_ant][Common.Conf_filename]

        calibration_errors = self.__build_calibration_errors()

        if calibr == Common.MutualCal:
            calibrator = AntennaCalibrator.MutualCalibrator(in_power, in_phase, row_steering, column_steering,
                                                            row_separation, col_separation, filename)
        else:
            calibrator = AntennaCalibrator.ClassicCalibrator(in_power, in_phase, row_separation, col_separation, filename)

        if calibration_errors:
            calibrator.add_calibration_errors(calibration_errors)

        if calibr == Common.MutualCal:
            calibrator.generate_cal_paths(AntennaCalibrator.every_one_to_one_path_strategy)

        return calibrator

    def __calibrate_antenna(self, calibrator):
        desired_tx_power = self.__config[Common.Conf_cal_param][Common.Conf_id_tx_power]
        desired_rx_power = self.__config[Common.Conf_cal_param][Common.Conf_id_rx_power]

        desired_phase = self.__get_desired_phases()

        desired_signals = [desired_tx_power, desired_phase, desired_rx_power, desired_phase]
        calibrator.calibrate_antenna(*desired_signals)

    @staticmethod
    def __get_pattern_properties(non_cal_pat, cal_pat, ideal_pat):
        return [non_cal_pat.get_pattern_properties(), cal_pat.get_pattern_properties(), ideal_pat.get_pattern_properties()]

    def clear_configuration(self):
        os.remove(self.__properties)

    def __save_properties(self, text, label):
        patt = r"\\begin{table}.*((?:\n.+(?!begin{table}))+(?:\n.*\\label{%s}))" %label
        with open(self.__properties, "r+") as f:
            tex = f.read()
            f.seek(0)
            f.write(re.sub(patt, repr(text[:-1])[1:-1], tex, flags=re.M))

    def __save_pattern_properties(self, non_cal_pat0, cal_pat0, ideal_pat0, non_cal_pat90, cal_pat90,
                                  ideal_pat90, name):
        decimals = 2
        namee = "properties"
        with open("Utilities/" + namee) as f:
            text = f.read()

        t = text.replace("@0", "{} & {} & {} & {}".format(*non_cal_pat0))
        t = t.replace("@1", "{} & {} & {} & {}".format(*non_cal_pat90))
        t = t.replace("@2", "{} & {} & {} & {}".format(*cal_pat0))
        t = t.replace("@3", "{} & {} & {} & {}".format(*cal_pat90))
        t = t.replace("@4", "{} & {} & {} & {}".format(*ideal_pat0))
        t = t.replace("@5", "{} & {} & {} & {}".format(*ideal_pat90))
        t = t.replace("@6", "{} & {} & {} & {}".format(*map(lambda x,y: round(x-y, decimals), non_cal_pat0, ideal_pat0)))
        t = t.replace("@7", "{} & {} & {} & {}".format(*map(lambda x,y: round(x-y, decimals), non_cal_pat90, ideal_pat90)))
        t = t.replace("@8", "{} & {} & {} & {}".format(*map(lambda x,y: round(x-y, decimals), cal_pat0, ideal_pat0)))
        t = t.replace("@9", "{} & {} & {} & {}".format(*map(lambda x,y: round(x-y, decimals), cal_pat90, ideal_pat90)))
        t = t.replace("@a", name)
        self.__save_properties(t, "tab:{}".format(name))


    def __compare_final_pattern_against_initial(self, calibrator, visual_comparator, title, filename):
        row_separation = self.__config[Common.Conf_ant][Common.Conf_vert_sep]
        col_separation = self.__config[Common.Conf_ant][Common.Conf_horiz_sep]
        freq = self.__config[Common.Conf_in_param][Common.Conf_freq]
        limits = [-100, 100]
        # limits = [-60, 60]
        phi = 0


        generator = PatternGenerator.PatternGenerator(freq, col_separation, row_separation)
        tx_power, tx_phase, _, _ = calibrator.get_antenna_gain_paths()

        non_cal_pattern = generator.generate_pattern(self.__tx_ini_ant_power, self.__tx_ini_ant_phase, limits, phi)
        ideal_pattern = generator.generate_pattern(*self.__create_ideal_output_power(),
                                                           start_stop_angle=limits, phi=phi)
        cal_pattern = generator.generate_pattern(tx_power, tx_phase, limits, phi)
        angles = cal_pattern.angles
        name = filename + "AzCut"
        visual_comparator.compare_patterns_against_ideal(angles, non_cal_pattern.get_db(), cal_pattern.get_db(),
                                                         ideal_pattern.get_db(), title + "Corte horizontal", name)

        # visual_comparator.compare_patterns(angles, non_cal_pattern.get_db(), ideal_pattern.get_db(), title + "Corte vertical", name)
        # visual_comparator.compare_patterns(angles, cal_pattern.get_db(), ideal_pattern.get_db(), title + "Corte vertical", name)

        props = self.__get_pattern_properties(non_cal_pattern, cal_pattern, ideal_pattern)

        phi = 90
        non_cal_pattern = generator.generate_pattern(self.__tx_ini_ant_power, self.__tx_ini_ant_phase, limits, phi)
        ideal_pattern = generator.generate_pattern(*self.__create_ideal_output_power(),
                                                           start_stop_angle=limits, phi=phi)
        cal_pattern = generator.generate_pattern(tx_power, tx_phase, limits, phi)
        name = filename + "ElCut"
        visual_comparator.compare_patterns_against_ideal(angles, non_cal_pattern.get_db(), cal_pattern.get_db(),
                                                         ideal_pattern.get_db(), title + "Corte vertical", name)

        # visual_comparator.compare_patterns(angles, non_cal_pattern.get_db(), ideal_pattern.get_db(), title + "Corte vertical", name)
        # visual_comparator.compare_patterns(angles, cal_pattern.get_db(), ideal_pattern.get_db(), title + "Corte vertical", name)

        props.extend(self.__get_pattern_properties(non_cal_pattern, cal_pattern, ideal_pattern))
        decimals = 2
        props = list(map(lambda x: list(map(lambda y: round(y, decimals), x)), props))
        props.append(filename)
        self.__save_pattern_properties(*props)

    def __compare_final_gain_against_initial(self, calibrator, visual_comparator, title, filename):
        quantity_rows = self.__config[Common.Conf_ant][Common.Conf_qtty_rows]
        quantity_columns = self.__config[Common.Conf_ant][Common.Conf_qtty_cols]

        desired_tx_power = self.__config[Common.Conf_cal_param][Common.Conf_id_tx_power]

        tx_signals = []
        append_signal_into_signals(tx_signals, self.__tx_ini_ant_power, self.__tx_ini_ant_phase)

        tx_ant_power, tx_ant_phase, _, _ = calibrator.get_antenna_gain_paths()
        append_signal_into_signals(tx_signals, tx_ant_power, tx_ant_phase)

        tx_ideal_power = [[desired_tx_power] * quantity_columns] * quantity_rows
        ideal_phase = self.__get_desired_phases()

        append_signal_into_signals(tx_signals, tx_ideal_power, ideal_phase)

        visual_comparator.compare_signals_against_ideal(*tx_signals, title="{}: Tx chain".format(title),
                                                        filename=filename)

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

        append_signal_into_signals(tx_signals, tx_ideal_power, ideal_phase)
        append_signal_into_signals(rx_signals, rx_ideal_power, ideal_phase)

        visual_comparator.compare_signals_against_ideal(*tx_signals, title="{}: Tx chain".format(title))
        visual_comparator.compare_signals_against_ideal(*rx_signals, title="{}: Rx chain".format(title))

    def __remove_antenna(self):
        filename = self.__config[Common.Conf_ant][Common.Conf_filename]
        for filename in glob.glob(filename + "_*"):
            os.remove(filename)

    def __create_ideal_output_power(self):
        quantity_rows = self.__config[Common.Conf_ant][Common.Conf_qtty_rows]
        quantity_columns = self.__config[Common.Conf_ant][Common.Conf_qtty_cols]

        ideal_power = self.__config[Common.Conf_cal_param][Common.Conf_id_tx_power]
        ideal_phase = self.__get_desired_phases()
        return ideal_power, ideal_phase

    def __get_error_name(self):
        names = {Common.Inter_pulse_gain_err: "chirpErr", Common.Chirp_rep_err: "chirpRepErr",
                 Common.Walsh_phase_err: "wallErr"}
        cal_errors = self.__config[Common.Conf_cal_param][Common.Conf_errors]

        if len(cal_errors) == 1:
            return names[cal_errors[0]]
        elif self.__config[Common.Conf_ant][Common.Conf_dead_trm]:
            return "deadTRMs"
        elif self.__config[Common.Conf_ant][Common.Conf_errors]:
            return "compErr"
        else:
            return "nonErr"

    def __get_angle(self):
        row_steering = self.__config[Common.Conf_in_param][Common.Conf_row_steer]
        column_steering = self.__config[Common.Conf_in_param][Common.Conf_col_steer]


        def f(x, y=''):
            return str(x) + "deg" + y

        return f(column_steering, "Col") if column_steering else f(row_steering, "Row") if row_steering else f(0)

    def run(self, calibr, , save_files=True):
        visual_comparator = VisualComparator.VisualComparator(save_files)

        # self.create_antenna()
        calibrator = self.__create_calibrator(calibr)
        prefix = self.__get_error_name() + self.__calibrator + self.__get_angle()

        self.__tx_ini_ant_power, self.__tx_ini_ant_phase, _, _ = calibrator.get_antenna_gain_paths()
        # compare_estimated_gains_against_real(calibrator, visual_comparator, "BEFORE CALIBRATION")

        self.__calibrate_antenna(calibrator)

        # compare_estimated_gains_against_real(calibrator, visual_comparator, "AFTER CALIBRATION")
        # self.__compare_estimated_gains_against_ideal(calibrator, visual_comparator, "")
        self.__compare_final_gain_against_initial(calibrator, visual_comparator, "RESULTS", prefix)

        self.__compare_final_pattern_against_initial(calibrator, visual_comparator, "patterns ", prefix)

        visual_comparator.show_graphics()

        #self.__remove_antenna()

        return 0


if __name__ == "__main__":
    sim = Simulator()
    sim.create_antenna()
    #sys.exit(sim.run(Common.ClassicCal))
    sys.exit(sim.run(Common.MutualCal, save_files=False))
