__author__ = 'francisco'

import unittest
import numpy as np
import src.Model.Antenna as Antenna
import src.Utilities.Antenna_Common as AntennaCommon


class NetworkTests(unittest.TestCase):
    def setUp(self):
        filename = "test"
        separation = 1
        self.__antenna = Antenna.Antenna()
        self.__create_antenna(filename, separation)

    def test_the_antenna_should_retrieve_a_correct_tx_s_parameters(self):
        # RM (0, 0)
        # cable
        s1 = np.matrix([[0.9958289498424986, 1.062822311415281],
                        [0.967442386293221, 0.0497725333530283]])

        # psc 14
        s2_a = np.matrix([[0.0296899042791072, 1.077112813624679],
                          [1.0182915320992103, 0.9553381008957869]])

        # cable
        s3 = np.matrix([[0.9831867162483208, 0.9313943667381208],
                        [0.9473643489647968, 0.0027731018616812]])

        # trm
        s4 = np.matrix([[0.015776122528384, 0.9140463321701637],
                        [1.0612094040146005, 0.9520974976976104]])

        # circulator
        s5 = np.matrix([[0.9600394972242867, 1.0984531613627018],
                        [1.054266537774246, 0.0504451320941863]])

        # cable
        s6 = np.matrix([[0.9245209536045264, 1.0527664905425866],
                        [1.0927176618768548, 0.0801894569638044]])

        f = lambda x: AntennaCommon.s2t_parameters(x)
        tx_1_matrix = AntennaCommon.t2s_parameters(f(s1)*f(s2_a)*f(s3)*f(s4)*f(s5)*f(s6))
        np.testing.assert_almost_equal(tx_1_matrix, self.__antenna.get_gain_paths("TxH")[0][0][0], decimal=14)

        # RM (0, 1)
        s2_b = np.matrix([[0.0296899042791072, 0.9307609512829832],
                          [0.9461543143665295, 0.9063045734770023]])
        # cable
        s3 = np.matrix([[0.9500599390422649, 1.05715064287121],
                        [1.0709140163535578, 0.9347302224777585]])
        # trm
        s4 = np.matrix([[0.0016234697821107, 1.029668886262384],
                        [1.0578039615665382, 0.9981837638934482]])
        # circulator
        s5 = np.matrix([[0.905739137163666, 1.0773654965268085],
                        [0.9299092320738062, 0.062223190694619]])
        # cable
        s6 = np.matrix([[0.9731243674865884, 0.9727099535342559],
                        [1.0890488566605971, 0.9010411205247378]])

        tx_2_matrix = AntennaCommon.t2s_parameters(f(s1)*f(s2_b)*f(s3)*f(s4)*f(s5)*f(s6))
        np.testing.assert_almost_equal(tx_2_matrix, self.__antenna.get_gain_paths("TxH")[0][0][1], decimal=14)

        # RM (1, 0)
        s2_c = np.matrix([[0.0296899042791072, 0.9404586608683853],
                          [0.9252378009469628, 0.9763915177125125]])
        # cable
        s3 = np.matrix([[0.9603450472744973, 1.0121026938948903],
                        [0.9051313326294861, 0.0069072925968607]])
        # trm
        s4 = np.matrix([[0.0583251036799202, 0.935980963233763],
                        [1.089640970311811, 0.0551623486549793]])
        # circulator
        s5 = np.matrix([[0.0761793545361649, 0.9374654398210459],
                        [0.9220476712267959, 0.9388495932799226]])
        # cable
        s6 = np.matrix([[0.9638591058440866, 0.9797275159698008],
                        [1.0640463319333375, 0.9891252748649797]])

        tx_3_matrix = AntennaCommon.t2s_parameters(f(s1)*f(s2_c)*f(s3)*f(s4)*f(s5)*f(s6))
        np.testing.assert_almost_equal(tx_3_matrix, self.__antenna.get_gain_paths("TxH")[0][1][0], decimal=14)

        # RM (1, 1)
        s2_d = np.matrix([[0.0296899042791072, 0.921582732133075],
                          [1.0122259069560258, 0.9597850800907802]])
        # cable
        s3 = np.matrix([[0.004530276705505, 1.0494231218650965],
                        [1.0710941360699917, 0.0174968479907784]])
        # trm
        s4 = np.matrix([[0.9778433462072884, 1.0643881323060533],
                        [0.9821270128376145, 0.9257936880594658]])
        # circulator
        s5 = np.matrix([[0.072893886225584, 0.9312125897788912],
                        [0.9020375873057322, 0.0040086768451177]])
        # cable
        s6 = np.matrix([[0.0905861642768713, 0.9817369357682371],
                        [ 0.9034234523476304, 0.0637755241693794]])

        tx_4_matrix = AntennaCommon.t2s_parameters(f(s1)*f(s2_d)*f(s3)*f(s4)*f(s5)*f(s6))
        np.testing.assert_almost_equal(tx_4_matrix, self.__antenna.get_gain_paths("TxH")[0][1][1], decimal=14)
        print("TxH", self.__antenna.get_gain_paths("TxH")[0])

    def test_the_antenna_should_retrieve_the_correct_rx_s_parameters(self):
        # RM (0, 0)
        # cable
        s1 = np.matrix([[0.0801894569638044, 1.0927176618768548],
                        [1.0527664905425866, 0.9245209536045264]])

        # circulator
        s2 = np.matrix([[0.9998757319224704, 1.0324552889581091],
                        [1.0508948937213702, 0.0504451320941863]])

        # trm
        s3 = np.matrix([[0.9284530034533963, 0.90397432852924],
                        [0.9837935705346328, 0.015776122528384]])

        # cable
        s4 = np.matrix([[0.0027731018616812, 0.9473643489647968],
                        [0.9313943667381208, 0.9831867162483208]])

        # psc 14
        s5_a = np.matrix([[0.9553381008957869, 1.0182915320992103],
                          [1.077112813624679, 0.0296899042791072]])

        # common cable
        s6 = np.matrix([[0.0497725333530283, 0.967442386293221],
                        [1.062822311415281, 0.9958289498424986]])

        f = lambda x: AntennaCommon.s2t_parameters(x)
        rx_1_matrix = AntennaCommon.t2s_parameters(f(s1)*f(s2)*f(s3)*f(s4)*f(s5_a)*f(s6))
        np.testing.assert_almost_equal(rx_1_matrix, self.__antenna.get_gain_paths("RxH")[0][0][0], decimal=14)

        # RM (0, 1)
        # cable
        s1 = np.matrix([[0.9010411205247378, 1.0890488566605971],
                        [0.9727099535342559, 0.9731243674865884]])

        # circulator
        s2 = np.matrix([[0.0661970346133671, 1.0855382724817848],
                        [1.0444426083417648, 0.062223190694619]])

        # trm
        s3 = np.matrix([[0.9414686420665633, 0.9269847346376758],
                        [0.9860255672592524, 0.0016234697821107]])

        # cable
        s4 = np.matrix([[0.9347302224777585, 1.0709140163535578],
                        [1.05715064287121, 0.9500599390422649]])

        # psc 14
        s5_b = np.matrix([[0.9063045734770023, 0.9461543143665295],
                          [0.9307609512829832, 0.0296899042791072]])

        f = lambda x: AntennaCommon.s2t_parameters(x)
        rx_2_matrix = AntennaCommon.t2s_parameters(f(s1)*f(s2)*f(s3)*f(s4)*f(s5_b)*f(s6))
        np.testing.assert_almost_equal(rx_2_matrix, self.__antenna.get_gain_paths("RxH")[0][0][1], decimal=14)

        # RM (1, 0)
        # cable
        s1 = np.matrix([[0.9891252748649797, 1.0640463319333375],
                        [0.9797275159698008, 0.9638591058440866]])

        # circulator
        s2 = np.matrix([[0.0596211949125214, 0.9311591940805913],
                        [1.0677349122779725, 0.9388495932799226]])

        # trm
        s3 = np.matrix([[0.0554018186960843, 0.9121776642410232],
                        [0.9953315626985285, 0.0583251036799202]])

        # cable
        s4 = np.matrix([[0.0069072925968607, 0.9051313326294861],
                        [1.0121026938948903, 0.9603450472744973]])

        # psc 14
        s5_c = np.matrix([[0.9763915177125125, 0.9252378009469628],
                          [0.9404586608683853, 0.0296899042791072]])

        f = lambda x: AntennaCommon.s2t_parameters(x)
        rx_3_matrix = AntennaCommon.t2s_parameters(f(s1)*f(s2)*f(s3)*f(s4)*f(s5_c)*f(s6))
        np.testing.assert_almost_equal(rx_3_matrix, self.__antenna.get_gain_paths("RxH")[0][1][0], decimal=14)

        # RM (1, 1)
        # cable
        s1 = np.matrix([[0.0637755241693794, 0.9034234523476304],
                        [0.9817369357682371, 0.0905861642768713]])

        # circulator
        s2 = np.matrix([[0.0849035304525922, 0.9385528577728693],
                        [1.0146230756127568, 0.0040086768451177]])

        # trm
        s3 = np.matrix([[0.9734516304631368, 1.010306751816013],
                        [0.9613481965348015, 0.9778433462072884]])

        # cable
        s4 = np.matrix([[0.0174968479907784, 1.0710941360699917],
                        [1.0494231218650965, 0.004530276705505]])

        # psc 14
        s5_d = np.matrix([[0.9597850800907802, 1.0122259069560258],
                          [0.921582732133075, 0.0296899042791072]])

        f = lambda x: AntennaCommon.s2t_parameters(x)
        rx_4_matrix = AntennaCommon.t2s_parameters(f(s1)*f(s2)*f(s3)*f(s4)*f(s5_d)*f(s6))
        np.testing.assert_almost_equal(rx_4_matrix, self.__antenna.get_gain_paths("RxH")[0][1][1], decimal=14)
        print("RxH", self.__antenna.get_gain_paths("RxH")[0])

    def __create_antenna(self, filename, separation):
        with open("base_rfdn") as f:
            with open(filename + "_rfdn", "w") as g:
                g.write(f.read())

        with open("base_panel") as f:
            with open(filename + "_panel", "w") as g:
                g.write(f.read())
        self.__antenna.initialize(separation, separation, filename)

if __name__ == '__main__':
    unittest.main()
