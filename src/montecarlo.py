#!/usr/bin/python3.4
import Utilities.Antenna_Common as common
import matplotlib.pyplot as plt
import simulator
import numpy as np
import xlsxwriter
import json
import sys
import os


def save_into_csv(filename, data):
    path = '/media/francisco/Datos/francisco/Documents/mutual/written/thesis/gfx/csv'
    workbook = xlsxwriter.Workbook(os.path.join(path, filename + '.xlsx'))
    worksheet = workbook.add_worksheet()
    # [worksheet.write_row(i, 0, d) for i, d in enumerate(data)]
    worksheet.write_row(0, 0, ['Elementos Radiantes', 'Classical', 'ClassicalStd', 'Mutual', 'MutualStd'])
    [worksheet.write_column(1, i, d) for i, d in enumerate(data)]


def plot_results(elements, gain, phase, title, fig, filename):
    mean_gain = np.mean(gain, axis=0).flatten()
    mean_phase = np.mean(phase, axis=0).flatten()

    plt.figure(fig)
    plt.title(title)
    plt.subplot(211)
    plt.errorbar(elements, mean_gain, np.std(gain, axis=0).flatten(), fmt='ok', lw=3)
    plt.errorbar(elements, mean_gain, [mean_gain - np.min(gain, 0).flatten(), np.max(gain, 0).flatten() - mean_gain],
                 fmt='.k', ecolor='gray', lw=1)
    plt.xlabel("Elementos")
    plt.ylabel("Ganancia [dB]")

    plt.subplot(212)
    plt.errorbar(elements, mean_phase, np.std(phase, axis=0).flatten(), fmt='ok', lw=3)
    plt.errorbar(elements, mean_phase, [mean_phase - np.min(phase, 0).flatten(), np.max(phase, 0).flatten() - mean_phase],
                 fmt='.k', ecolor='gray', lw=1)
    plt.xlabel("Elementos")
    plt.ylabel("Fase [deg]")
    plt.savefig(os.path.join('/media/francisco/Datos/francisco/Documents/mutual/written/thesis/gfx', filename + ".png"), bbox_inches='tight')


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


def get_mean_and_std(data):
    return [np.mean(data, axis=0), np.std(data, axis=0)]


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

    [run_simulator(sim, mut_gain, mut_phase, clas_gain, clas_phase) for _ in range(500)]

    plt.close('all')
    elements = range(len(mut_gain[0])*len(mut_gain[0][0]))

    # plot_results(elements, mut_gain, mut_phase, "Calibración con acoplamientos mutuos", 1, "mutualMontecarlo")
    # plot_results(elements, clas_gain, clas_phase, "Calibración clásica", 2, "classicalMontecarlo")

    els = list(elements)

    dataClas = [[it for l in values for it in l] for values in clas_gain]
    dataMut = [[it for l in values for it in l] for values in mut_gain]
    save_into_csv('MontecarloPower', [els] + get_mean_and_std(dataClas) + get_mean_and_std(dataMut))

    dataClas = [[it for l in values for it in l] for values in clas_phase]
    dataMut = [[it for l in values for it in l] for values in mut_phase]
    save_into_csv('MontecarloPhase', [els] + get_mean_and_std(dataClas) + get_mean_and_std(dataMut))

    plt.show()


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
