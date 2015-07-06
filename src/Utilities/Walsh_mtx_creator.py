__author__ = 'francisco'

import numpy as np


class WalshMatrixCreator:

    def __init__(self):
        np.random.seed(seed=None)
        self.__add_errors = False
        self.__std_phase = 0
        self.__walkm = 1

    def add_walsh_errors(self, errors):
        self.__add_errors, self.__std_phase = [True if isinstance(i, str) else i for error in errors for i in error
                                                 if error[0] == 'WalPhaseErrors']

    def __initialize_walsh(self, order):
        walkm = 1
        for i in range(order):
            self.__walkm = np.vstack((np.hstack((walkm, walkm)), np.hstack((walkm, -walkm))))

    def create_phase_walsh_matrix(self, amount_elements, ph_shift=np.pi/2):
        order = np.ceil(np.log2(amount_elements))
        self.__initialize_walsh(order)

        walkm = self.__walkm * ph_shift if ph_shift != 0 else self.__walkm
        return walkm + np.random.normal(sigma=self.__add_errors, size=(order, order)) if self.__add_errors else walkm
