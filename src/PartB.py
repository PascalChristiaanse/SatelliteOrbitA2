from Constants import Constants


class PartB:
    def __init__(self):
        pass

    def Q1(self):
        return """
        The PRN or pseudo random noise codes are part of the C/A or coarse/acquisition code, which is a 1023 bit sequence that is repeated every millisecond. 
        The code encrypts several things, but one of the most important parts is the "time" at which the code was sent.
        This time called a chip, and resolves the time within the millisecond window and is sent around every microsecond, giving around that precision. 
        This is not an absolute time but a code that can be aligned with a clock using a special algorithm, thus finding the time it was sent.
        The navigational message contains the local time at which it was sent and needs to be resolved using a receiver clock and using data from multiple satellites.
        By comparing the time the code was sent with the local time, the range to the satellite can be found.
        If the code was static or non-changing, it couldnt contain a time signal and this method of ranging would not be possible.
        
        """
