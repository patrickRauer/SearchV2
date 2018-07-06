from astroquery.vizier import Vizier as Viz
from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack
from astropy import units as u
from threading import Thread
from Search_V2.core.match import match_catalogs
import numpy as np
Viz.ROW_LIMIT = -1


class Vizier:
    columns = ['*']
    catalog = ''
    radius = 2*u.arcsec
    vizier = Viz
    name = ''
    ra_name = 'RAJ2000'
    dec_name = 'DEJ2000'
    cache = None

    def __init__(self, catalog):
        """
        Basic class for an interface between the Vizier class of the astroquery package
        :param catalog: The name of the catalog in the Vizier-database
        :type catalog: str
        """
        self.catalog = catalog

    def __query__(self, coordinates):
        """
        Perform the query with a split of the coordinates in 200 item
        blocks to reduce a problem with the connection

        :param coordinates: The list of coordinates of the targets
        :type coordinates: SkyCoord
        :return: A array or a list of arrays with the results of the queries
        :rtype: list, numpy.ndarray
        """

        if self.name != '' and self.name not in self.ra_name:
            self.ra_name = '{}_{}'.format(self.name, self.ra_name)
            self.dec_name = '{}_{}'.format(self.name, self.dec_name)
        steps = 200
        i = 0
        out = []

        if type(self.catalog) == list:
            for i in range(len(self.catalog)):
                out.append([])
        else:
            out = [[]]

        # perform the query with maximal 200 targets at the same time
        # this is necessary because otherwise it happens from time to time that a timeout appears
        while i*steps < len(coordinates):
            coords = coordinates[i*steps: (i+1)*steps]
            rs = self.vizier.query_region(coords, self.radius,
                                          catalog=self.catalog)
            for j, r in enumerate(rs):
                out[j].append(r)
            i += 1

        # stack the results
        for i in range(len(out)):
            if len(out[i]) != 0:
                o = vstack(out[i])

                # add the name of a prefix to the column names
                if self.name != '':
                    for c in o.colnames:
                        o.rename_column(c, '{}_{}'.format(self.name, c))
                out[i] = np.array(o)

        # if only one catalog query is performed
        if len(out) == 1:
            out = out[0]
            self.cache = out
            return out
        else:
            self.cache = out
            return out

    def query(self, data, ra_name='ra', dec_name='dec'):
        """
        Start a query around the input positions in the Vizier database

        :param data: The input data
        :type data: numpy.ndarray, astropy.table.Table
        :param ra_name: The name of the RA column
        :type ra_name: str
        :param dec_name: The name of the Dec column
        :type dec_name: str
        :return: A array with the results of the query
        :rtype: numpy.ndarray
        """
        if type(data) == Table:
            # if no units are add to the coordinate columns
            if data[ra_name].unit is None:
                coordinates = SkyCoord(data[ra_name]*u.deg,
                                       data[dec_name]*u.deg)
            # if units are add to the coordinate columns
            else:
                coordinates = SkyCoord(data[ra_name],
                                       data[dec_name])
        # if the data are for example numpy.ndarray
        else:
            coordinates = SkyCoord(data[ra_name]*u.deg,
                                   data[dec_name]*u.deg)

        return self.__query__(coordinates)


class Gaia(Vizier):
    columns = ['RA_ICRS', 'DE_ICRS',
               'Plx', 'e_Plx',
               'pmRA', 'e_pmRA',
               'pmDE', 'e_pmDE',
               'Gmag', 'e_Gmag',
               'BPmag', 'e_BPmag',
               'RPmag', 'e_RPmag',
               'epsi', 'sepsi']
    name = 'GAIA'
    ra_name = 'GAIA_RA_ICRS'
    dec_name = 'GAIA_DE_ICRS'

    def __init__(self):
        """
        Child class to query GAIA DR2 data
        """
        Vizier.__init__(self, 'I/345/gaia2')
        self.vizier = Viz(columns=self.columns)
        self.vizier.ROW_LIMIT = -1


class TwoMass(Vizier):

    columns = ['RAJ2000', 'DEJ2000',
               'Jmag', 'e_Jmag',
               'Hmag', 'e_Hmag',
               'Kmag', 'e_Kmag']
    name = '2MASS'

    def __init__(self):
        """
        Child class to query 2MASS data
        """
        Vizier.__init__(self, 'II/246/out')
        self.vizier = Viz(columns=self.columns)
        self.vizier.ROW_LIMIT = -1


class Wise(Vizier):
    columns = ['RAJ2000', 'DEJ2000',
               'W1mag', 'e_W1mag',
               'W2mag', 'e_W2mag',
               'W3mag', 'e_W3mag',
               'W4mag', 'e_W4mag']
    name = 'WISE'

    def __init__(self):
        """
        Child class to query AllWise data
        """
        Vizier.__init__(self, 'II/328/allwise')
        self.vizier = Viz(columns=self.columns)
        self.vizier.ROW_LIMIT = -1


class PanStarrs(Vizier):
    columns = ['RAJ2000', 'DEJ2000',
               'gmag', 'e_gmag',
               'rmag', 'e_rmag',
               'imag', 'e_imag',
               'zmag', 'e_zmag',
               'ymag', 'e_ymag']
    name = 'PS'

    def __init__(self):
        """
        Child class to query PanSTARRS DR1 data
        """
        Vizier.__init__(self, 'II/349/ps1')
        self.vizier = Viz(columns=self.columns)
        self.vizier.ROW_LIMIT = -1


class SDSS(Vizier):
    columns = ['RAJ2000', 'DEJ2000',
               'umag', 'e_umag',
               'gmag', 'e_gmag',
               'rmag', 'e_rmag',
               'imag', 'e_imag',
               'zmag', 'e_zmag']
    name = 'SDSS'

    def __init__(self):
        """
        Child class to query SDSS DR12 data
        """
        Vizier.__init__(self, 'V/147/sdss12')
        self.vizier = Viz(columns=self.columns)
        self.vizier.ROW_LIMIT = -1


class Apass(Vizier):
    columns = ['RAJ2000', 'DEJ2000',
               'Vmag', 'e_Vmag',
               'Bmag', 'e_Bmag',
               'gpmag', 'e_gpmag',
               'rpmag', 'e_rpmag',
               'ipmag', 'e_ipmag']
    name = 'APASS'

    def __init__(self):
        """
        Child class to query APASS DR9 data
        """
        Vizier.__init__(self, 'II/336/apass9')
        self.vizier = Viz(columns=self.columns)
        self.vizier.ROW_LIMIT = -1


class Galex(Vizier):
    columns = ['RAJ2000', 'DEJ2000',
               'NUV', 'e_NUV',
               'FUV', 'e_FUV']
    name = 'GALEX'

    def __init__(self):
        """
        Child class to query GALEX DR5 data
        """
        Vizier.__init__(self, 'II/312/ais')
        self.vizier = Viz(columns=self.columns)
        self.vizier.ROW_LIMIT = -1


class MultiSurvey:
    surveys = []

    def __init__(self, sdss=True, panstarrs=True,
                 apass=True, two_mass=True,
                 all_wise=True, gaia=True,
                 galex=True):
        """
        This class provides a simple way to query multiple surveys at the same time.
        It also merges all the resulting data into one file with outer join (no-matches will get empty cells).

        :param sdss: True if the query should include SDSS DR12, else False
        :type sdss: bool
        :param panstarrs: True if the query should include PanSTARRS DR1, else False
        :type panstarrs: bool
        :param apass: True if the query should include APASS DR9, else False
        :type apass: bool
        :param two_mass: True if the query should include 2MASS, else False
        :type two_mass: bool
        :param all_wise: True if the query should include AllWISE, else False
        :type all_wise: bool
        :param gaia: True if the query should include GAIA DR2, else False
        :type gaia: bool
        :param galex: True if the query should include GALEX DR5, else False
        :type galex: bool
        """
        if sdss:
            self.surveys.append(SDSS())

        if panstarrs:
            self.surveys.append(PanStarrs())

        if apass:
            self.surveys.append(Apass())

        if two_mass:
            self.surveys.append(TwoMass())

        if all_wise:
            self.surveys.append(Wise())

        if gaia:
            self.surveys.append(Gaia())

        if galex:
            self.surveys.append(Galex())

    def query(self, data, ra_name='ra', dec_name='dec'):
        """
        Queries the chosen surveys

        :param data: The input data
        :type data: numpy.ndarray, astropy.table.Table
        :param ra_name: The name of the RA column
        :type ra_name: str
        :param dec_name: The name of the Dec column
        :type dec_name: str
        :return: A array with the results of the query
        :rtype: numpy.ndarray
        """
        th = []
        for s in self.surveys:
            t = Thread(target=s.query,
                       args=(data, ra_name, dec_name))
            t.start()
            th.append(t)

        # join the threads
        for t in th:
            t.join()

        data = data.copy()
        # match all the collected data in one table
        for s in self.surveys:
            if len(s.cache) > 0:
                data = match_catalogs(data, s.cache,
                                      ra_name, dec_name,
                                      s.ra_name, s.dec_name)
        return data
