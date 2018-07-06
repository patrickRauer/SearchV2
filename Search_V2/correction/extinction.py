# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 13:19:30 2018

This script provides the extinction calculation based one the IRSA dust map.

@author: Jean Patrick Rauer
"""

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table, Column
from astropy.units.core import UnitTypeError
from astropy.wcs import WCS

# TODO: replace the path
path = '/Users/patrickr/Documents/GitHub/SimulationSearch/SimulationSearch/databases/data/'
a = [0.32999, -0.77530, 0.01979, 0.72085, -0.02427, -0.50447, 0.17699, 1]
a = np.poly1d(a)
b = [-2.09002, 5.30260, -0.62251, -5.38434, 1.07233, 2.28305, 1.41338, 0]
b = np.poly1d(b)

R_v = 3.1


def nir_al_av(wavelength):
    """
    NIR extinction computation

    .. math:
        E(\lambda ) = (0.574 - 0.527/R_v) \lambda^1.61

    :param wavelength: Central wavelength of the required bands
    :type wavelength: list
    :return: NIR extinction
    """
    # TODO: check the exponents
    return 0.574 * np.power(wavelength, 1.61) - 0.527 * np.power(wavelength, 1.61) / R_v


def opt_al_av(wavelength):
    """
    Optical extinction computation

    :param wavelength: Central wavelength of the required bands
    :type wavelength: list
    :return: Optical extinction
    """
    return a(wavelength) + b(wavelength) / R_v


def a_lambda_a_v(wavelength):
    """
    Calculation of the extinctions values

    :param wavelength: Central wavelength of the required bands
    :type wavelength: list
    :return: Optical extinction
    """
    wavelength = np.array(wavelength)
    al_av = np.zeros(wavelength.shape)

    nir = wavelength <= 1.1
    opt = wavelength > 1.1

    al_av[nir] = nir_al_av(wavelength[nir])
    al_av[opt] = opt_al_av(wavelength[opt] - 1.82)
    return al_av


def extract_bv_values(s, file_name):
    with fits.open(path + file_name) as fi:
        wcs = WCS(fi[0].header)
        x, y = wcs.all_world2pix(s.galactic.l, s.galactic.b, 0)
        e_bv = fi[0].data[(np.int32(y), np.int32(x))]
        return e_bv


def extictions_direct_multi(s, wavelength):
    """
    Computes the extinction in multiple bands at a specific coordinate

    :param s: The coordinates of the required extinction
    :type s: astropy.coordinates.SkyCoord
    :param wavelength: The wavelength of required extinctions
    :type wavelength: numpy.ndarray, list
    :return: The extinctions at the input wavelengths
    :rtype: numpy.ndarray
    """
    # split the coordinates in two samples
    # one above the galactic plane
    north = s[s.galactic.b >= 0]
    # and one below the galactic plane
    south = s[s.galactic.b < 0]
    e_bv = np.zeros(len(s))

    # if there are targets above the galactic plane
    if len(north) > 0:
        file_name = 'SFD_dust_4096_ngp.fits'
        e_bv[s.galactic.b >= 0] = extract_bv_values(north, file_name)

    # if there are targets below the galactic plane
    if len(south) > 0:
        file_name = 'SFD_dust_4096_sgp.fits'
        e_bv[s.galactic.b < 0] = extract_bv_values(south, file_name)

    al_av = a_lambda_a_v(wavelength)
    av = R_v * e_bv
    al_av = np.tile(al_av, (len(av), 1))
    av = np.tile(av, (al_av.shape[1], 1))
    al = np.transpose(av) * al_av
    return al


def get_extinctions(ra, dec):
    """
    Computes the extinction for the coordinates

    :param ra: The RA coordinate in degree
    :type ra: float
    :param dec: The Dec coordinate in degree
    :type dec: float
    :return: An table with the extinctions in the different bands
    :rtype: astropy.table.Table
    """
    # TODO: replace try-except with an if-statement
    try:
        s = SkyCoord(ra * u.deg, dec * u.deg)
    except UnitTypeError:
        s = SkyCoord(ra, dec)

    wavelengths = [0.3734, 0.4309, 0.5517, 0.6520, 0.8007,
                   0.4621, 0.6546, 0.8111, 0.3587, 0.4717,
                   0.6165, 0.7476, 0.8923, 1.248, 1.659,
                   2.190, 1.23, 1.64, 2.16, 3.52, 4.46,
                   5.66, 7.68, 3.32, 4.57]
    wavelengths = np.array(wavelengths)
    filter_names = ['CTIO U', 'CTIO B', 'CTIO V', 'CTIO R',
                    'CTIO I', 'DSS-II g', 'DSS-II r', 'DSS-II i',
                    'SDSS u', 'SDSS g', 'SDSS r', 'SDSS i',
                    'SDSS z', 'UKIRT J', 'UKIRT H', 'UKIRT K',
                    '2MASS J', '2MASS H', '2MASS Ks', 'IRAC-1',
                    'IRAC-2', 'IRAC-3', 'IRAC-4', 'WISE-1',
                    'WISE-2']
    extinc = extictions_direct_multi(s, 1. / wavelengths)
    return Table(rows=extinc, names=filter_names)


def add_extinction(tab, ra_name, dec_name):
    """
    Adds extinction columns to a table

    :param tab: The table to add the extinction columns
    :type tab: astropy.table.Table
    :param ra_name: The name of the RA column
    :type ra_name: str
    :param dec_name: The name of the Dec column
    :type dec_name: str
    :return: The input table with the addtional extinctions columns
    :rtype: astropy.table.Table
    """
    extinction = get_extinctions(tab[ra_name], tab[dec_name])

    for c in extinction.colnames:
        print(c)
        tab[c] = np.round(extinction[c], 2)
    return tab


def extinctions_direct(ra, dec, wavelength):
    """
    Computes the extinction for a single or a set of coordinates and for a set of wavelengths

    :param ra: The RA coordinate
    :type ra: float, list
    :param dec: The Dec coordinate
    :type dec: float, list
    :param wavelength: The wavelengths for the extinction
    :type wavelength: float, list
    :return:
    """
    s = SkyCoord(ra * u.deg, dec * u.deg)
    if type(ra) == list or type(ra) == np.ndarray or type(ra) == Column:
        return extictions_direct_multi(s, wavelength)
    else:

        # decide which file is needed
        # northern part of the galactic plane
        if s.galactic.b >= 0:
            file_name = 'SFD_dust_4096_ngp.fits'
        # southern part of the galactic plane
        else:
            file_name = 'SFD_dust_4096_sgp.fits'

        with fits.open(path + file_name) as fi:
            wcs = WCS(fi[0].header)
            x, y = wcs.all_world2pix(s.galactic.l, s.galactic.b, 0)
            e_bv = fi[0].data[int(y), int(x)]
            al_av = a_lambda_a_v(wavelength)

            av = R_v * e_bv

            return al_av * av
