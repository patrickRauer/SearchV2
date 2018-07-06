import armapy
from .data import Data
# TODO: implement SED


class SED:
    data = None
    seds = None

    def __init__(self, data):
        if type(data) is not Data:
            data = Data(data)
        self.data = data

    def __download_properties__(self):
        self.data.head.download_SVO_properties()

    def __calc_sed__(self):
        pass
