__author__ = 'fsoler'

import sys
import src.Visual_Comparator.Visual_Comparator as VisualComparator
import src.Pattern_Generator.Pattern_Generator as PatternGenerator
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import src.Controllers.Antenna_Creator as AntennaCreator
import src.Controllers.Antenna_Calibrator as AntennaCalibrator
import src.Utilities.Antenna_Common as Common
import json


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

    print("{}:".format(title), "Tx antenna power    ", tx_ant_power)
    print("{}:".format(title), "Tx estimated power  ", tx_power)
    print("{}:".format(title), "Tx abs power error  ", (np.array(tx_ant_power) - np.array(tx_power)).tolist())
    print("{}:".format(title), "Tx antenna phase    ", tx_ant_phase)
    print("{}:".format(title), "Tx estimated phase  ", tx_phase)
    print("{}:".format(title), "Tx abs phase error  ", (np.array(tx_ant_phase) - np.array(tx_phase)).tolist())
    print("")

    print("{}:".format(title), "Rx antenna power    ", rx_ant_power)
    print("{}:".format(title), "Rx estimated power  ", rx_power)
    print("{}:".format(title), "Rx abs power error  ", (np.array(rx_ant_power) - np.array(rx_power)).tolist())
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

    def __create_antenna(self):
        # quantity_rows = self.__config[Common.Conf_ant][Common.Conf_qtty_rows]
        quantity_columns = self.__config[Common.Conf_ant][Common.Conf_qtty_cols]
        row_separation = self.__config[Common.Conf_ant][Common.Conf_vert_sep]
        col_separation = self.__config[Common.Conf_ant][Common.Conf_horiz_sep]
        row_steering = self.__config[Common.Conf_in_param][Common.Conf_row_steer]
        column_steering = self.__config[Common.Conf_in_param][Common.Conf_col_steer]
        filename = self.__config[Common.Conf_ant][Common.Conf_filename]

        # att = 0.1           # [neper/m]
        c = 299792458       # [m/seg]
        f = self.__config[Common.Conf_in_param][Common.Conf_freq]      # [Hz]
        wavelength = c/f    # [m]
        """
        length1 = 0.45      # [m]
        length2 = 8         # [m]
        length3 = 0.5       # [m]

        trm_gain = 10       # []
        trm_ph_shift = 10   # [deg]

        psc_out_ports = quantity_columns * quantity_rows

        cable1 = [Common.Cable, [att, wavelenght, length1]]
        psc = ["{0}1{1}".format(Common.Psc, psc_out_ports), [psc_out_ports]]
        cable2 = [Common.Cable, [att, wavelenght, length2]]
        trm = [Common.Trm, [trm_gain, trm_ph_shift]]
        circulator = [Common.Circulator, []]
        cable3 = [Common.Cable, [att, wavelenght, length3]]
        rm = [Common.Rm, []]
        """
        sequence_items = self.__build_sequence(wavelength)
        component_errors = self.__build_component_errors()
        """
        sequence_items = [cable1, psc, cable2, trm, circulator, cable3, rm]
        component_errors = [[Common.Rm_error, 0.1], [Common.Trm_error, 1.9],
                            [Common.Circulator_error, 0.5], [Common.Psc_error, 0.5]]
        """
        creator = AntennaCreator.AntennaCreator(quantity_columns, row_separation, col_separation)
        # creator.add_errors(component_errors)
        creator.create_structure(filename, sequence_items, row_steering, column_steering, [[0, 0], [1, 1]])
        """
        create_antenna2(filename)
        """

    def __create_calibrator(self):
        in_power = self.__config[Common.Conf_in_param][Common.Conf_power]
        in_phase = self.__config[Common.Conf_in_param][Common.Conf_phase]
        row_separation = self.__config[Common.Conf_ant][Common.Conf_vert_sep]
        col_separation = self.__config[Common.Conf_ant][Common.Conf_horiz_sep]
        filename = self.__config[Common.Conf_ant][Common.Conf_filename]

        calibration_errors = self.__build_calibration_errors()
        """
        calibration_errors = [[Common.Inter_pulse_power_err, 0.5], [Common.Inter_pulse_phase_err, 5],
                              [Common.Gain_chirp_rep_err, 1], [Common.Phase_chirp_rep_err, 5]]
        """

        calibrator = AntennaCalibrator.MutualCalibrator(in_power, in_phase, row_separation, col_separation, filename)
        # calibrator.add_calibration_errors(calibration_errors)
        calibrator.generate_cal_paths(AntennaCalibrator.every_one_to_one_path_strategy)
        """
        calibrator = AntennaCalibrator.ClassicCalibrator(input_power, input_phase, separation, separation, filename)
        # calibrator.add_calibration_errors(calibration_errors)
        """
        return calibrator

    def __calibrate_antenna(self, calibrator):
        desired_tx_power = self.__config[Common.Conf_cal_param][Common.Conf_id_tx_power]
        desired_rx_power = self.__config[Common.Conf_cal_param][Common.Conf_id_rx_power]

        desired_tx_phase = self.__config[Common.Conf_cal_param][Common.Conf_id_tx_phase]
        desired_rx_phase = self.__config[Common.Conf_cal_param][Common.Conf_id_rx_phase]

        desired_signals = [desired_tx_power, desired_tx_phase, desired_rx_power, desired_rx_phase]
        calibrator.calibrate_antenna(*desired_signals)

    def __compare_estimated_gains_against_ideal(self, calibrator, visual_comparator, title):
        quantity_rows = self.__config[Common.Conf_ant][Common.Conf_qtty_rows]
        quantity_columns = self.__config[Common.Conf_ant][Common.Conf_qtty_cols]

        desired_tx_power = self.__config[Common.Conf_cal_param][Common.Conf_id_tx_power]
        desired_rx_power = self.__config[Common.Conf_cal_param][Common.Conf_id_rx_power]

        desired_tx_phase = self.__config[Common.Conf_cal_param][Common.Conf_id_tx_phase]
        desired_rx_phase = self.__config[Common.Conf_cal_param][Common.Conf_id_rx_phase]

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
        tx_ideal_phase = [[desired_tx_phase] * quantity_columns] * quantity_rows
        rx_ideal_phase = [[desired_rx_phase] * quantity_columns] * quantity_rows

        append_signal_into_signals(tx_signals, tx_ideal_power, tx_ideal_phase)
        append_signal_into_signals(rx_signals, rx_ideal_power, rx_ideal_phase)

        visual_comparator.compare_signals_against_ideal(*tx_signals, title="{}: Tx chain".format(title))
        visual_comparator.compare_signals_against_ideal(*rx_signals, title="{}: Rx chain".format(title))


    def __remove_antenna(self):
        filename = self.__config[Common.Conf_ant][Common.Conf_filename]
        for filename in glob.glob(filename + "_*"):
            os.remove(filename)

    def run(self):
        visual_comparator = VisualComparator.VisualComparator()
        """
        filename = "test"
        separation = 1
        quantity_columns = 2
        quantity_rows = 2
        row_steering = 10
        column_steering = 10
        """
        self.__create_antenna()
        """
        input_power = 10
        input_phase = 0
        """
        calibrator = self.__create_calibrator()
        compare_estimated_gains_against_real(calibrator, visual_comparator, "BEFORE CALIBRATION")

        self.__calibrate_antenna(calibrator)

        compare_estimated_gains_against_real(calibrator, visual_comparator, "AFTER CALIBRATION")
        self.__compare_estimated_gains_against_ideal(calibrator, visual_comparator, "")
        """
        quantity_rows = self.__config[Common.Conf_ant][Common.Conf_qtty_rows]
        quantity_columns = self.__config[Common.Conf_ant][Common.Conf_qtty_cols]
        params = [calibrator, visual_comparator, quantity_rows, quantity_columns]
        params.extend(desired_signals)
        params.append("")
        compare_estimated_gains_against_ideal(*params)
        """
        visual_comparator.show_graphics()

        self.__remove_antenna()

        return 0


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


def create_antenna(quantity_columns, quantity_rows, separation, row_steering, column_steering, filename):
    att = 0.1           # [neper/m]
    c = 299792458       # [m/seg]
    f = 1275000000      # [Hz]
    wavelenght = c/f    # [m]

    length1 = 0.45      # [m]
    length2 = 8         # [m]
    length3 = 0.5       # [m]

    trm_gain = 10       # []
    trm_ph_shift = 10   # [deg]

    psc_out_ports = quantity_columns * quantity_rows

    cable1 = [Common.Cable, [att, wavelenght, length1]]
    psc = ["{0}1{1}".format(Common.Psc, psc_out_ports), [psc_out_ports]]
    cable2 = [Common.Cable, [att, wavelenght, length2]]
    trm = [Common.Trm, [trm_gain, trm_ph_shift]]
    circulator = [Common.Circulator, []]
    cable3 = [Common.Cable, [att, wavelenght, length3]]
    rm = [Common.Rm, []]

    sequence_items = [cable1, psc, cable2, trm, circulator, cable3, rm]
    component_errors = [[Common.Rm_error, 0.1], [Common.Trm_error, 1.9],
                        [Common.Circulator_error, 0.5], [Common.Psc_error, 0.5]]
    """
    rms = quantity_columns * quantity_rows
    sequence_items = ["cable", "PSC1{0}".format(rms), "cable", "TRM", "circulator", "cable", "RM"]
    """
    creator = AntennaCreator.AntennaCreator(quantity_columns, separation, separation)
    # creator.add_errors(component_errors)
    creator.create_structure(filename, sequence_items, row_steering, column_steering, [[0, 0], [1, 1]])
    """
    create_antenna2(filename)
    """


def create_calibrator(input_power, input_phase, separation, filename):

    calibration_errors = [[Common.Inter_pulse_power_err, 0.5], [Common.Inter_pulse_phase_err, 5],
                          [Common.Gain_chirp_rep_err, 1], [Common.Phase_chirp_rep_err, 5]]

    calibrator = AntennaCalibrator.MutualCalibrator(input_power, input_phase, separation, separation, filename)
    # calibrator.add_calibration_errors(calibration_errors)
    calibrator.generate_cal_paths(AntennaCalibrator.every_one_to_one_path_strategy)
    """
    calibrator = AntennaCalibrator.ClassicCalibrator(input_power, input_phase, separation, separation, filename)
    # calibrator.add_calibration_errors(calibration_errors)
    """
    return calibrator


def compare_estimated_gains_against_ideal(calibrator, visual_comparator, quantity_rows, quantity_columns,
                                          desired_tx_power, desired_tx_phase, desired_rx_power, desired_rx_phase,
                                          title):
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
    tx_ideal_phase = [[desired_tx_phase] * quantity_columns] * quantity_rows
    rx_ideal_phase = [[desired_rx_phase] * quantity_columns] * quantity_rows

    append_signal_into_signals(tx_signals, tx_ideal_power, tx_ideal_phase)
    append_signal_into_signals(rx_signals, rx_ideal_power, rx_ideal_phase)

    visual_comparator.compare_signals_against_ideal(*tx_signals, title="{}: Tx chain".format(title))
    visual_comparator.compare_signals_against_ideal(*rx_signals, title="{}: Rx chain".format(title))


def remove_antenna(filename):
    for filename in glob.glob(filename + "_*"):
        os.remove(filename)


def main():
    visual_comparator = VisualComparator.VisualComparator()
    filename = "test"
    separation = 1
    quantity_columns = 2
    quantity_rows = 2
    row_steering = 10
    column_steering = 10

    create_antenna(quantity_columns, quantity_rows, separation, row_steering, column_steering, filename)

    input_power = 10
    input_phase = 0

    calibrator = create_calibrator(input_power, input_phase, separation, filename)
    compare_estimated_gains_against_real(calibrator, visual_comparator, "BEFORE CALIBRATION")

    desired_tx_power = 20
    desired_rx_power = 0

    desired_tx_phase = 0
    desired_rx_phase = 0

    desired_signals = [desired_tx_power, desired_tx_phase, desired_rx_power, desired_rx_phase]
    calibrator.calibrate_antenna(*desired_signals)
    compare_estimated_gains_against_real(calibrator, visual_comparator, "AFTER CALIBRATION")

    params = [calibrator, visual_comparator, quantity_rows, quantity_columns]
    params.extend(desired_signals)
    params.append("")
    compare_estimated_gains_against_ideal(*params)

    visual_comparator.show_graphics()

    remove_antenna(filename)

    return 0


if __name__ == "__main__":
    sys.exit(Simulator().run())
    # sys.exit(compare_signals())
    # sys.exit(generate_pattern())
