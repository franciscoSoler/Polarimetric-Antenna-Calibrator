__author__ = 'fsoler'

import sys
import src.Visual_Comparator.Visual_Comparator as VisualComparator
import src.Pattern_Generator.Pattern_Generator as Pattern_Generator
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import src.Controllers.RFDNCreator as RFDNCreator
import src.Controllers.Antenna_Calibrator as AntennaCalibrator
import src.Utilities.Antenna_Common as AntennaCommon


def generate_pattern():
    """
    l_band_freq = 1275MHz
    s_band_freq = 2GHz
    x_band_freq = 8 GHz
    :return:
    """
    generator = Pattern_Generator.PatternGenerator(1275000000, 0.127, 0.110)
    weights = [[1 for _ in range(20)]for _ in range(1)]

    angles, power = generator.generate_pattern(weights, [-30, 30], 0)
    plt.plot(angles, 20*np.log10([abs(p) for p in power]))
    plt.show()


def create_antenna(filename):
    with open("base_rfdn") as f:
        with open(filename + "_rfdn", "w") as g:
            g.write(f.read())

    with open("base_panel") as f:
        with open(filename + "_panel", "w") as g:
            g.write(f.read())


def append_signal_into_signals(signals, power_signal, phase_signal):
    signals.append(power_signal)
    signals.append(phase_signal)


def main():
    filename = "test"
    power = 20
    separation = 1
    quantity_columns = 2
    quantity_rows = 2

    desired_tx_power = 20
    desired_rx_power = 0

    desired_tx_phase = 0
    desired_rx_phase = 0

    desired_signals = [desired_tx_power, desired_tx_phase, desired_rx_power, desired_rx_phase]
    tx_signals = []
    rx_signals = []

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

    cable1 = [AntennaCommon.Cable, [att, wavelenght, length1]]
    psc = ["{0}1{1}".format(AntennaCommon.Psc, psc_out_ports), [psc_out_ports]]
    cable2 = [AntennaCommon.Cable, [att, wavelenght, length2]]
    trm = [AntennaCommon.Trm, [trm_gain, trm_ph_shift]]
    circulator = [AntennaCommon.Circulator, []]
    cable3 = [AntennaCommon.Cable, [att, wavelenght, length3]]
    rm = [AntennaCommon.Rm, []]

    sequence_items = [cable1, psc, cable2, trm, circulator, cable3, rm]

    """
    rms = quantity_columns * quantity_rows
    sequence_items = ["cable", "PSC1{0}".format(rms), "cable", "TRM", "circulator", "cable", "RM"]
    """
    creator = RFDNCreator.AntennaCreator(quantity_columns, separation, separation)
    creator.create_structure(filename, sequence_items, 0.1)
    """
    create_antenna(filename)
    """

    calibrator = AntennaCalibrator.AntennaCalibrator(power, separation, separation, filename)
    print("equations", calibrator.generate_cal_paths(AntennaCalibrator.every_one_to_one_path_strategy))

    tx_power, tx_phase = calibrator.get_transmission_power()
    rx_power, rx_phase = calibrator.get_reception_power()
    tx_power_shift, tx_phase_shift = calibrator.get_tx_signal_shifts(desired_tx_power)
    rx_power_shift, rx_phase_shift = calibrator.get_rx_signal_shifts(desired_rx_power)

    append_signal_into_signals(tx_signals, tx_power, tx_phase)
    append_signal_into_signals(rx_signals, rx_power, rx_phase)

    calibrator.calibrate_antenna(*desired_signals)

    tx_cal_power, tx_cal_phase = calibrator.get_transmission_power()
    rx_cal_power, rx_cal_phase = calibrator.get_reception_power()
    tx_cal_power_shift, tx_cal_phase_shift = calibrator.get_tx_signal_shifts(desired_tx_power)
    rx_cal_power_shift, rx_cal_phase_shift = calibrator.get_rx_signal_shifts(desired_rx_power)

    append_signal_into_signals(tx_signals, tx_cal_power, tx_cal_phase)
    append_signal_into_signals(rx_signals, rx_cal_power, rx_cal_phase)
    """
    print("Tx power         ", tx_power)
    print("Tx shift         ", tx_power_shift)
    print("Tx power + shift ", (np.array(tx_power) + np.array(tx_power_shift)).tolist())
    print("Tx cal power     ", tx_cal_power)
    """
    print("Tx phase         ", tx_phase)
    print("Tx shift         ", tx_phase_shift)
    print("Tx phase + shift ", (np.array(tx_phase) + np.array(tx_phase_shift)).tolist())
    print("Tx cal phase     ", tx_cal_phase)
    print("")
    """
    print("Rx power         ", rx_power)
    print("Rx shift         ", rx_power_shift)
    print("Rx power + shift ", (np.array(rx_power) + np.array(rx_power_shift)).tolist())
    print("Rx cal power     ", rx_cal_power)
    """
    print("Rx phase         ", rx_phase)
    print("Rx shift         ", rx_phase_shift)
    print("Rx phase + shift ", (np.array(rx_phase) + np.array(rx_phase_shift)).tolist())
    print("Rx cal phase     ", rx_cal_phase)

    visual_comparator = VisualComparator.VisualComparator()

    tx_ideal_power = [[desired_tx_power] * quantity_columns] * quantity_rows
    rx_ideal_power = [[desired_rx_power] * quantity_columns] * quantity_rows
    tx_ideal_phase = [[desired_tx_phase] * quantity_columns] * quantity_rows
    rx_ideal_phase = [[desired_rx_phase] * quantity_columns] * quantity_rows

    append_signal_into_signals(tx_signals, tx_ideal_power, tx_ideal_phase)
    append_signal_into_signals(rx_signals, rx_ideal_power, rx_ideal_phase)

    visual_comparator.compare_signals_against_ideal(*tx_signals, title="Tx power")
    visual_comparator.compare_signals_against_ideal(*rx_signals, title="Rx power")

    visual_comparator.show_graphics()
    """
    for filename in glob.glob(filename + "_*"):
        os.remove(filename)
    """
    return 0


if __name__ == "__main__":
    sys.exit(main())
    # sys.exit(compare_signals())
    # sys.exit(generate_pattern())