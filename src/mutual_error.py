#!/usr/bin/python3.4
import Utilities.Antenna_Common as common
import matplotlib.pyplot as plt
import simulator
import numpy as np
import json
import sys
import os


def plot_results(phases, max_phase, title, filename):
    plt.plot(phases, max_phase)
    plt.xlabel("Elementos")
    plt.ylabel("Max. Diferencia Fase [deg]")
    plt.savefig(os.path.join('/media/francisco/Datos/francisco/Documents/mutual/written/thesis/gfx', filename + ".png"), bbox_inches='tight')


def run_simulator(simm, phase):
    simm.create_antenna()
    phase.append(simm.run(common.MutualCal, mutual_error=True))


def save_file(filename, structure):
    with open(filename, "w") as f:
        f.write(json.dumps(structure, sort_keys=False, indent=4, separators=(',', ': ')))


def load_file(filename):
    with open(filename) as f:
        config = json.loads(f.read())
    return config


def create_config_file(sim, trm_config, max_phases):
    config_base = "configurationFileBase"
    config_filename = "configurationFile"

    config = load_file(config_base)
    config[common.Conf_in_param][common.Conf_row_steer] = 50
    config[common.Conf_in_param][common.Conf_col_steer] = 50
    config[common.Conf_component][common.Trm.lower()][common.Phase_shift] = trm_config

    save_file(config_filename, config)

    phase = []
    [run_simulator(sim, phase) for _ in range(100)]
    max_phases.append(np.max(phase))


def main():
    sim = simulator.Simulator()
    phases = range(0, 100, 10)
    max_phases = []
    [create_config_file(sim, trm_phase, max_phases) for trm_phase in phases]
    plot_results(phases, max_phases, "error fases", "errorMutualEstimationPhase")


if __name__ == "__main__":
    sys.exit(main())
