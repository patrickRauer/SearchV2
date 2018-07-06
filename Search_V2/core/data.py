import os
import armapy
import numpy as np
from astropy.io import fits
from astropy.table import Table
from numpy.lib.recfunctions import append_fields
import pandas
from Search_V2.core.log import Log

ID_NAME = 'row_id'


class SurveyItem:

    name = ''
    cols = None
    effective_wavelength = None
    flux_zeropoint = None

    def __init__(self, name, cols):
        self.name = name
        self.cols = cols

    def get_magnitude_column(self, mag_name):
        """
        Returns the magnitude column name in the default format

        :param mag_name: The name of the magnitude
        :type mag_name: str
        :return: The name of the magnitude column in the right format
        :rtype: str
        """
        if mag_name in self.cols:
            pass
        else:
            return None

    def __download_properties__(self):
        for c in self.cols:
            properties = armapy.get_filter_information(self.name, self.name, c)


class Head:
    mag_cols = ['g', 'r', 'i', 'z', 'y', 'J', 'H', 'K']
    mag_str = '{}mag'

    backup_path = './temp/'

    def __init__(self, mag_cols=None):
        """
        Contains the meta information of a Data object.
        This includes for example column names.
        """
        if mag_cols is not None:
            self.mag_cols = mag_cols

    def get_mag_name(self, c):
        """
        Returns the column name for a given band
        :param c: The name of the band
        :type c: str
        :return: The magnitude column name based on mag_str layout
        :rtype: str
        """
        return self.mag_str.format(c)

    def __add__(self, other):
        """
        Add the information of this head object to a second one and return the results

        :param other: The second head object
        :type other: Search_V2.core.data.Head
        :return: The combined head object
        :rtype: Search_V2.core.data.Head
        """
        mag_cols = self.mag_cols
        mag_cols.extend(other.mag_cols)
        h = Head(mag_cols=mag_cols)
        return h


class Data:
    path = None
    head = None
    data = None
    data_original = None
    log = None
    backup = False
    backup_count = 0

    # random id of the
    data_id = np.random.randint(1, 10000, 1)[0]

    def __init__(self, data=None, path=None, head=None, backup=False):
        """
        This class is the basic data storage. All search and analyse actions will use this class.

        :param path: Path to the data
        :type path: str
        :param head: Meta information object of the data. Default is None which standard Vizier meta data will be used
        :type head: Search_V2.core.data.Head
        :param backup:
            If True it will make a backup (saves the current data) before every new process step.
            The place can be specified in the Head. Default is False and no backups will be made.
        :type backup: bool
        """
        self.log = Log()
        self.path = path

        if path is not None:
            self.data, head = load_data(self.path)
        elif data is not None:
            self.data = data
        else:
            raise AttributeError('No data set!')
        # self.data = pandas.DataFrame(self.data)
        # add an ID field which is in the beginning just the row number
        # self.data = append_fields(self.data, ID_NAME,
        #                           np.linspace(1, len(self.data), num=len(self.data)))
        self.data[ID_NAME] = np.linspace(1, len(self.data),
                                         num=len(self.data),
                                         dtype=np.int32)
        # if head is None, use the default Header
        if head is None:
            head = Head()
        self.head = head

        self.log.add_log('load data')
        self.__calc_colors__()

        self.data_original = self.data.copy()
        self.backup = backup

    def __calc_colors__(self):
        """
        Computes all possible colors
        :return:
        """
        for i, c1 in enumerate(self.head.mag_cols):
            for j, c2 in enumerate(self.head.mag_cols):
                if i < j:
                    m1 = self.head.get_mag_name(c1)
                    m2 = self.head.get_mag_name(c2)
                    self.data['{}{}'.format(c1, c2)] = self.data[m1]-self.data[m2]
                    
        self.log.add_log('colors computed')

    def set_data(self, data, method, comment=''):
        """
        Replaces the old data and add a new log item to the log.
        If backup is set to True, it will save a backup before it
        replaces the old data.

        :param data: The new data-set
        :type data: numpy.ndarray
        :param method: The performed method to generate the new data set
        :type method: str
        :param comment: Additional comment
        :type comment: str
        :return:
        """
        if self.backup:
            self.__make_backup__()
        self.data = data
        self.log.add_log(method, comment)

    def __make_backup__(self):
        """
        Creates the backup path name and calls the writing method
        to save the current data set

        :return:
        """
        path = os.path.join(self.head.backup_path,
                            'backup_{}_{}.fits'.format(self.data_id,
                                                       self.backup_count))
        self.write(path, True)
        self.backup_count += 1

    def write(self, path, backup=False):
        """
        Writes the current data to file
        :param path: The path to the file
        :type path: str
        :param backup:
            True if it should be a back up (only the id's will be saved) or False otherwise.
            default is False
        :type backup: bool
        :return:
        """
        # check if the path to the directory exists
        if not os.path.exists(os.path.split(path)[0]):
            os.makedirs(os.path.split(path)[0])

        if '.fits' in path:
            # TODO: include header
            # for backups the ids are the only needed items
            if backup:
                hdu = Table.from_pandas(self.data[[ID_NAME]])
                hdu.write(path, overwrite=True)
            else:
                hdu = fits.BinTableHDU(np.array(self.data))
                hdu.writeto(path, overwrite=True)


def load_data(path):
    """
    Loads the data from a file and convert them to a numpy.structuredarray.

    Current supported file types are fits and csv. For larger files it is highly recommented to use
    fits instead of csv files.

    :param path: Path to the file
    :type path: str
    :return: Loaded data
    :type: numpy.ndarray
    """
    if 'fits' in os.path.split(path)[-1]:
        data = Table.read(path).to_pandas()
        head = None
    elif 'csv' in os.path.split(path)[-1]:
        # TODO: implement csv loads
        data = None
        head = None
    else:
        raise AttributeError('Unknown or unsupported file format!')
    return data, head


SDSS_2MASS_HEAD = Head(['u', 'g', 'r', 'i', 'z',
                        'J', 'H', 'K'])
