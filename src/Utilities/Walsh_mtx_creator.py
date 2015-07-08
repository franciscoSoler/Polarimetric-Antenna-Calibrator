__author__ = 'francisco'

import numpy as np
import src.Utilities.Antenna_Common as AntennaCommon


class WalshMatrixCreator:
    @staticmethod
    def __get_order(amount_elements):
        return np.ceil(np.log2(amount_elements))

    def __init__(self):
        np.random.seed(seed=None)
        self.__add_errors = False
        self.__std_phase = 0
        self.__walkm = 1

    def add_walsh_errors(self, errors):
        if [True for error in errors if error[0] == AntennaCommon.Walsh_phase_err]:
            self.__add_errors = True
            self.__std_phase = [error[1] for error in errors if error[0] == AntennaCommon.Walsh_phase_err].pop()

    def __initialize_walsh(self, order):
        walkm = 1
        for i in range(int(order)):
            walkm = np.vstack((np.hstack((walkm, walkm)), np.hstack((walkm, -walkm))))
        self.__walkm = walkm

    def create_ideal_phase_walsh(self, amount_elements, ph_shift=np.pi/2):
        order = self.__get_order(amount_elements)
        self.__initialize_walsh(order)
        return self.__walkm * ph_shift if ph_shift != 0 else self.__walkm

    def create_phase_walsh_matrix(self, amount_elements, ph_shift=np.pi/2):
        order = self.__get_order(amount_elements)
        walkm = self.create_ideal_phase_walsh(amount_elements, ph_shift)
        return walkm + np.random.normal(sigma=self.__std_phase, size=(order, order)) if self.__add_errors else walkm
