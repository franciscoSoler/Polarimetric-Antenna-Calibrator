#!/usr/bin/python3.4
import Pattern_Generator.Pattern_Generator as pat_gen
import Utilities.Antenna_Common as common
import numpy as np
import json
import sys
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


class Simulator:

    def __init__(self):
        self.__start = 0.1
        self.__stop = 1
        self.__step = 0.1
        self.__config = ""
        self.__read_config()
        self.__quantity_rows = self.__config[common.Conf_ant][common.Conf_qtty_rows]
        self.__quantity_columns = self.__config[common.Conf_ant][common.Conf_qtty_cols]

    def __read_config(self):
        with open("configurationFile") as f:
            self.__config = json.load(f)

    def __compute_realization(self, generator, sigma_gain, sigma_phase):
        power = np.random.normal(0, sigma_gain, (self.__quantity_rows, self.__quantity_columns))
        phase = np.random.normal(0, sigma_phase, (self.__quantity_rows, self.__quantity_columns))
        return generator.generate_pattern(power, phase, [-0.5, 0.5]).get_pattern_peak()

    def __obtain_sigma(self, generator, sigma_gain, sigma_phase):
        gain = [self.__compute_realization(generator, sigma_gain, sigma_phase) for _ in range(1000)]
        return np.std(gain)

    def __plot_results(self, x, y ,z):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot_trisurf(x, y, z, cmap=cm.jet, linewidth=0.2)
        ax.set_xlabel('sigma fase')
        ax.set_ylabel('sigma ganancia')
        ax.set_zlabel('Ganancia')

        plt.show()

    def run(self):
        d_y = self.__config[common.Conf_ant][common.Conf_vert_sep]
        d_x = self.__config[common.Conf_ant][common.Conf_horiz_sep]
        generator = pat_gen.PatternGenerator(common.f, d_x, d_y)

        quantity = (self.__stop - self.__start)/self.__step
        x = np.linspace(self.__start, self.__stop*20, quantity) # phase
        y = np.linspace(self.__start, self.__stop, quantity) # gain

        z = []
        [z.append(self.__obtain_sigma(generator, i, j)) for j in x for i in y]

        y = np.tile(y, quantity)
        x = np.repeat(x, quantity)
        self.__plot_results(x, y, z)


if __name__ == "__main__":
    sim = Simulator()
    sys.exit(sim.run())
