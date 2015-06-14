__author__ = 'francisco'


class CalSignalGenerator():
    """
    generates the signal paths
    :keyword attributes:
    Tx_Rx_modes -- the first word means the quantity of TRMs emitting at the same time.
                       The second word means the the same but in reception
    Tx_separation -- the quantity of RMs are between the two emitting/receiving RMs
    Tx_Rx_separation -- the RM shift that are between the higher emitter RM and the higher receiver RM
                        e.g. if Tx_Rx_modes="two-to-one" and Tx_Rx_separation=0, the TxRm, which has the higher index,
                        and RxRM are aligned (face to face).
    """
    Tx_Rx_modes = ["one-to-one", "two-to-one", "one-to-two"]
    Tx_separation = list(range(20))
    Tx_Rx_separation = list(range(20))

    def __init__(self):
        pass

    def generate_cal_signals(self, trans_rec_mode, rm_separation=0, trans_to_rec_separation=None):
        """
        generates the signal calibration paths
        :keyword parameters:
        trans_rec_mode -- the first word means the quantity of TRMs emitting at the same time.
                          The second word means the the same but in reception.
        rm_separation -- the quantity of RMs are between the two emitting or receiving RMs (default = 0)
        trans_to_rec_separation -- the RM shift that are between the higher emitter RM and the higher receiver RM
                                   (default = None). If None, all transmitter/receivers attached to "one mode" are used
        :returns:
        cal_signals -- this list of lists depends on the mode used, the format is with the same emitter, all receivers
                       are stored consecutively as att-phShift pairs.-
        """
        if not [mode for mode in CalSignalGenerator.Tx_Rx_modes if mode == trans_rec_mode]:
            raise Exception("the mode {0} is not valid", trans_rec_mode)
        if not [separation for separation in CalSignalGenerator.Tx_separation if separation == rm_separation]:
            raise Exception("the separation between RMs {0} is not valid", trans_rec_mode)
        if trans_to_rec_separation is not None and not [separation for separation in CalSignalGenerator.Tx_Rx_separation
                                                        if separation == trans_to_rec_separation]:
            raise Exception("the transmission to reception separation {0} is not valid", trans_to_rec_separation)
        pass