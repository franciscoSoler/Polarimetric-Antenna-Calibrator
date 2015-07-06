__author__ = 'fsoler'
import json
import numpy as np
import re
import src.Utilities.Antenna_Common as AntennaCommon


class Antenna:
    __AvailableModes = ["TxH", "TxV", "RxH", "RxV", "TxH-RxV", "TxV-RxH"]
    __Polarizations = [AntennaCommon.Rfdn_h_pol, AntennaCommon.Rfdn_v_pol]

    def __init__(self):
        self.__json_rfdn = None
        self.__json_panel = None
        self.__quantity_rows = 0
        self.__quantity_columns = 0
        self.__matrix_distances = None
        pass

    @property
    def quantity_rows(self):
        return self.__quantity_rows

    @property
    def quantity_columns(self):
        return self.__quantity_columns

    def initialize(self, dist_rows, dist_columns, filename="antenna"):
        """
        :param dist_rows: physical distance between rows of rm
        :param dist_columns: physical distance between columns of rm
        :param filename: filename where the antenna model is written
        :return:
        """
        with open(filename + "_rfdn") as f:
            self.__json_rfdn = json.load(f)
        with open(filename + "_panel") as f:
            self.__json_panel = json.load(f)
        self.__quantity_rows = int(re.match("^RM \(([\d]+).*", self.__json_panel[-1][0]).group(1)) + 1
        self.__quantity_columns = int(re.match("^RM \([\d]+, ([\d]+).*", self.__json_panel[-1][0]).group(1)) + 1
        (self.__matrix_distances, _) = AntennaCommon.calculate_distances_between_rms(self.__quantity_rows,
                                                                                     self.__quantity_columns, dist_rows,
                                                                                     dist_columns, False)

    def __get_qtty_antennas(self, structure):
        if type(structure) is not dict:
            return 1

        component, value = list(structure.items())[0]

        if AntennaCommon.is_psc(component):
            qtty_ports = AntennaCommon.get_qtty_output_ports(component)
            acumulated_ports = self.__get_qtty_antennas(value["extremeAttached"][0])
            return qtty_ports * acumulated_ports
        else:
            return self.__get_qtty_antennas(value["extremeAttached"])

    def get_qtty_antennas(self):
        return self.__get_qtty_antennas(next(iter(self.__json_rfdn.values())))

    def __get_attenuation_paths(self, structure, mode):
        """

        :param structure:
        :param mode: is Transmission or Reception, used to determine which sij must use for TRM
        :return:
        """
        if type(structure) is not dict:
            return [[structure, np.identity(2)]]

        component, value = list(structure.items())[0]

        children = []
        if AntennaCommon.is_psc(component):
            [children.extend(self.__get_attenuation_paths(extreme, mode)) for extreme in value[AntennaCommon.Extreme]]
        else:
            children.extend(self.__get_attenuation_paths(value[AntennaCommon.Extreme], mode))

        # section to obtain the correct PSC port number
        out_ports = AntennaCommon.get_qtty_output_ports(component)
        idx = lambda x: self.row_col_to_index(*AntennaCommon.get_rm_position(x))
        f = lambda x: int((idx(x) * out_ports / len(children)) % out_ports)

        # convert the S parameters of the component in complex values
        parameters = [list(map(complex, si)) for si in value[AntennaCommon.SParams]]
        build_cascade = lambda x, y: x * y if mode == AntennaCommon.Transmission else y * x
        return [[rm_pos, build_cascade(AntennaCommon.s2t_parameters(AntennaCommon.get_s2p(component, parameters, mode,
                                                                                          f(rm_pos))), child_param)]
                for rm_pos, child_param in children]

    def get_gain_paths(self, pol_mode):
        """
        :keyword parameters:
        pol_mode -- must be one of the available modes: TxH, TxV, RxH, RxV, TxH-RxV or TxV-RxH

        :return:
        att_path -- list of S parameter of every path of the antenna depending on the polarization required.
        """
        if not [mode for mode in Antenna.__AvailableModes if pol_mode == mode]:
            raise Exception("polarization {0} is not a valid mode", pol_mode)

        modes = AntennaCommon.parse_polarization_mode(pol_mode)

        f = lambda x, y: [AntennaCommon.t2s_parameters(matrix[1]) for matrix in sorted(x, key=lambda z: z[0])]
        format_list = lambda x: list(map(lambda y: x[self.__quantity_columns * y: self.__quantity_columns * (y+1)],
                                         range(self.__quantity_rows)))
        return [format_list(f(self.__get_attenuation_paths(self.__json_rfdn[mode[1]], mode[0]),
                              mode[0])) for mode in modes]

    def get_mutual_coupling_front_panel(self):
        """
        calculates the s parameters between all RMs
        :return:
        s_parameters -- dict with its keys are every RMs id in transmission mode and its values are a pair of RMid-S
        parameter in reception mode
        """
        str2cplx = lambda x: list(map(lambda y: list(map(complex, y)), x))
        f = lambda x: eval(re.match(".*(\(.*\))", x).group(1))
        create_dict = lambda x: dict([[f(bla[0]), np.matrix(str2cplx(bla[1]))] for bla in x])
        return {f(self.__json_panel[idx][0]): create_dict(self.__json_panel[idx][1])
                for idx in range(self.__quantity_columns * self.__quantity_rows)}

    def __change_trm_param(self, structure, power_shifts, f):
        if AntennaCommon.is_rm(structure):
            return tuple(map(lambda x: int(x), re.match(".*([\d]+), ([\d]+).*", structure).groups()))

        component, value = list(structure.items())[0]

        if AntennaCommon.is_psc(component):
            [self.__change_trm_param(extreme, power_shifts, f) for extreme in value["extremeAttached"]]
        else:
            rm_position = self.__change_trm_param(value["extremeAttached"], power_shifts, f)
            if AntennaCommon.is_trm(component):
                parameters = [list(map(complex, si)) for si in value["sParameters"]]
                parameters[f[0]][f[1]] *= power_shifts.item(rm_position)
                value["sParameters"] = [list(map(str, si)) for si in parameters]
            return rm_position

    def change_trm_rx_params(self, power_shift, polarization):
        if not [mode for mode in Antenna.__Polarizations if polarization == mode]:
            raise Exception("polarization {0} is not a valid mode", polarization)

        self.__change_trm_param(self.__json_rfdn[polarization], power_shift, f=[0, 2])

    def change_trm_tx_params(self, power_shift, polarization):
        """
        changes attenuation and phase of every Trans-Rec Module.
        :param power_shift: list of lists containing every attenuation and phase shift for every RM in dB
        :param polarization: hPolarization or vPolarization
        :return:
        """
        if not [mode for mode in Antenna.__Polarizations if polarization == mode]:
            raise Exception("polarization", polarization, "is not a valid mode")

        self.__change_trm_param(self.__json_rfdn[polarization], power_shift, f=[1, 0])

    def index_to_row_col(self, index):
        return divmod(index, self.__quantity_columns)

    def row_col_to_index(self, row, col):
        return row * self.__quantity_columns + col

    def get_index_from_rm_separation(self, row, col):
        return self.__matrix_distances[(row, col)]
