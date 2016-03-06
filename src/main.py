#!/usr/bin/python3.4
import Utilities.Antenna_Common as common
import simulator
import json
import sys
from pprint import pprint


def run_simulator(simm):
    simm.create_antenna()
    simm.run(common.MutualCal)
    simm.run(common.ClassicCal)


def save_file(filename, structure):
    with open(filename, "w") as f:
        f.write(json.dumps(structure, sort_keys=False, indent=4, separators=(',', ': ')))


def load_file(filename):
    with open(filename) as f:
        config = json.loads(f.read())
    return config


def create_config_file(row_steering, column_steering, key, value):
    config_base = "configurationFileBase"
    config_filename = "configurationFile"

    config = load_file(config_base)
    config[common.Conf_in_param][common.Conf_row_steer] = row_steering
    config[common.Conf_in_param][common.Conf_col_steer] = column_steering

    if key is not None:
        config[key[0]][key[1]] = value
    save_file(config_filename, config)

    simm = simulator.Simulator()
    run_simulator(simm)


def add_errors(row_steering, column_steering):
    if row_steering and column_steering:
        return
    errors = [[((common.Conf_cal_param), (common.Conf_errors)), [common.Inter_pulse_gain_err]],
              [((common.Conf_cal_param), (common.Conf_errors)), [common.Chirp_rep_err]],
              [((common.Conf_cal_param), (common.Conf_errors)), [common.Walsh_phase_err]],
              [((common.Conf_ant), (common.Conf_errors)), [common.Circulator_error, common.Trm_error, common.Psc_error]],
              [((common.Conf_ant), (common.Conf_dead_trm)), [[1, 0], [3, 4], [5, 5]]]]
    [create_config_file(row_steering, column_steering, item[0], item[1]) for item in errors]
    create_config_file(row_steering, column_steering, None, "")


def main():
    row_steering = [0, 10]
    column_steering = [0, 10]

    [add_errors(r_steer, c_steer) for r_steer in row_steering for c_steer in column_steering]


if __name__ == "__main__":
    sys.exit(main())
