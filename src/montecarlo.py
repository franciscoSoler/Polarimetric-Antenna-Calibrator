#!/usr/bin/python3.4
import Utilities.Antenna_Common as common
import simulator
import numpy as np
import json
import sys


def run_simulator(simm, mut_gain, mut_phase, clas_gain, clas_phase):
    simm.create_antenna()
    m_gain, m_phase = simm.run(common.MutualCal, save_files=False, show_graph=False)
    c_gain, c_phase = simm.run(common.ClassicCal, save_files=False, show_graph=False)

    mut_gain.append(m_gain)
    mut_phase.append(m_phase)
    clas_gain.append(c_gain)
    clas_phase.append(c_phase)


def save_file(filename, structure):
    with open(filename, "w") as f:
        f.write(json.dumps(structure, sort_keys=False, indent=4, separators=(',', ': ')))


def load_file(filename):
    with open(filename) as f:
        config = json.loads(f.read())
    return config


def create_config_file(errors):
    config_base = "configurationFileBase"
    config_filename = "configurationFile"

    config = load_file(config_base)
    config[common.Conf_in_param][common.Conf_row_steer] = 0
    config[common.Conf_in_param][common.Conf_col_steer] = 0

    for error in errors:
        config[error[0][0]][error[0][1]] = error[1]
    save_file(config_filename, config)

    mut_gain = []
    mut_phase = []
    clas_gain = []
    clas_phase = []
    sim = simulator.Simulator()

    [run_simulator(sim, mut_gain, mut_phase, clas_gain, clas_phase) for _ in range(100)]
    print("mutual", np.std(mut_gain, axis=0), np.std(mut_phase, axis=0))
    print("classical", np.std(clas_gain, axis=0), np.std(clas_phase, axis=0))


def main():
    errors = [[((common.Conf_cal_param), (common.Conf_errors)), [common.Inter_pulse_gain_err]],
              [((common.Conf_cal_param), (common.Conf_errors)), [common.Chirp_rep_err]],
              [((common.Conf_cal_param), (common.Conf_errors)), [common.Walsh_phase_err]],
              [((common.Conf_ant), (common.Conf_errors)), [common.Circulator_error, common.Trm_error, common.Psc_error]],
              [((common.Conf_ant), (common.Conf_dead_trm)), [[1, 0], [4, 3], [5, 5]]]]
    create_config_file(errors)
    return 0

if __name__ == "__main__":
    sys.exit(main())
