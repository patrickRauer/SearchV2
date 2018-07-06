import time
import urllib
from threading import Thread

import numpy as np
import sewpy
from astropy.io import fits
from astropy.table import Table, vstack

# Sextractor columns
params = ["X_IMAGE", "Y_IMAGE", 'X_WORLD', 'Y_WORLD',
          "FLUX_RADIUS(3)", "FLAGS", 'FWHM_IMAGE', 'FWHM_WORLD',
          'KRON_RADIUS', 'PETRO_RADIUS', 'THETA_IMAGE', 'A_IMAGE',
          'B_IMAGE', 'A_WORLD', 'B_WORLD', 'ELLIPTICITY', 'CLASS_STAR']

# Sextractor short cut
sextractor_name = 'sextractor'
TEMP_PANSTARRS_PATH = './temp/'

try:
    sewpy.SEW(params=params,
              config={"DETECT_MINAREA": 10, "PHOT_FLUXFRAC": "0.3, 0.5, 0.8",
                      'SATUR_LEVEL': 1111060, 'SEEING_FWHM': 1.31, 'GAIN': 1.030748,
                      'PIXEL_SCALE': 0.25})
except RuntimeError:
    sextractor_name = 'sex'

# SExtractor properties for the different bands
se = {'g': sewpy.SEW(params=params,
                     config={"DETECT_MINAREA": 10, "PHOT_FLUXFRAC": "0.3, 0.5, 0.8",
                             'SATUR_LEVEL': 1111060, 'SEEING_FWHM': 1.31, 'GAIN': 1.030748,
                             'PIXEL_SCALE': 0.25}),

      'r': sewpy.SEW(params=params,
                     config={"DETECT_MINAREA": 10, "PHOT_FLUXFRAC": "0.3, 0.5, 0.8",
                             'SATUR_LEVEL': 1582662, 'SEEING_FWHM': 1.19, 'GAIN': 1.100708,
                             'PIXEL_SCALE': 0.25}),
      'i': sewpy.SEW(params=params,
                     config={"DETECT_MINAREA": 10, "PHOT_FLUXFRAC": "0.3, 0.5, 0.8",
                             'SATUR_LEVEL': 1994015, 'SEEING_FWHM': 1.11, 'GAIN': 1.057824,
                             'PIXEL_SCALE': 0.25}),
      'z': sewpy.SEW(params=params,
                     config={"DETECT_MINAREA": 10, "PHOT_FLUXFRAC": "0.3, 0.5, 0.8",
                             'SATUR_LEVEL': 1089845, 'SEEING_FWHM': 1.07, 'GAIN': 1.030766,
                             'PIXEL_SCALE': 0.25}),
      'y': sewpy.SEW(params=params,
                     config={"DETECT_MINAREA": 10, "PHOT_FLUXFRAC": "0.3, 0.5, 0.8",
                             'SATUR_LEVEL': 2514071, 'SEEING_FWHM': 1.02, 'GAIN': 1.066191,
                             'PIXEL_SCALE': 0.25}),
      }


def get_all_image_urls(ra, dec, size):
    """
    Extract the image urls from the PAN-STARRS image web-site

    :param ra: The RA component of the coordinate of the target
    :type ra: float
    :param dec: The DEC component of the coordinate of the target
    :type dec: float
    :param size: The size of the image in acrsec
    :type size: float
    :return: Dict with where the keys are the name of the bands and the items are the urls to images
    :rtype: dict
    """
    bands = ['g', 'r', 'i', 'z', 'y']
    user_agent = ''.join(['Mozilla/5.0 (Macintosh; U; Intel Mac OS X 10_6_4; en-US) ',
                          'AppleWebKit/534.3 (KHTML, like Gecko) Chrome/6.0.472.63 Safari/534.3'])
    headers = {'User-Agent': user_agent}

    # create the url of the panstarrs image website
    url = '{}{}{}{}{}{}{}{}'.format('http://ps1images.stsci.edu/cgi-bin/ps1cutouts?pos=',
                                    ra,
                                    '+%09',
                                    dec,
                                    '&filter=color&filter=g&filter=r&filter=i&filter=z&filter=y',
                                    '&filetypes=stack&auxiliary=data&auxiliary=mask&size=',
                                    int(size * 3600 / 0.25),
                                    '&output_size=0&verbose=0&autoscale=99.500000&catlist=')

    # call the website and download the source code
    req = urllib.request.Request(url, None, headers)
    response = urllib.request.urlopen(req)
    page = response.read()
    response.close()
    page = str(page)

    # extract the url of the cutout from the source of the image website
    page = page.split('\\n')
    out = {}
    for band in bands:
        img_url = ''
        bband = ''.join(['.', band, '.'])
        for n in page:
            if 'FITS-cutout' in n and bband in n:
                img_url = n
                break
        img_url = img_url.split('Download FITS cutout')[-1]
        img_url = img_url.split('="')[-1]
        img_url = img_url.split('">')[0]
        out.update({band: img_url})
    return out


def get_all_panstarrs_images(ra, dec, size):
    """
    Downloads all 5 PAN-STARRS images and read the data in.

    :param ra: The center of the image in RA direction in degree
    :type ra: float
    :param dec: The center of the image in DEC direction in degree
    :type dec: float
    :param size: The size of the image in degrees
    :type size: float
    :return: A dict where the keys are the name of the different filters and the items are the images itself
    :rtype: dict
    """
    paths = "id_{}_band_{}.fits"
    download_all_bands(ra, dec, size, paths, 0)
    bands = ['g', 'r', 'i', 'z', 'y']
    out = {}
    for b in bands:
        with fits.open(paths.format(0, b)) as fi:
            out.update({b: {'header': fi[0].header, 'data': fi[0].data}})
    return out


def download_all_bands(ra, dec, size, save_path, extract_id):
    """
    Downloads all five filter-images

    :param ra: The ra component of the coordinate in deg
    :type ra: float
    :param dec: The dec component of the coordinate in deg
    :type dec: float
    :param size: The size of the image in degrees
    :type size: float
    :param save_path: The path to the save place with two parameters
    :type save_path: str
    :param extract_id:
        Unique key for the run to distinguish between different extraction runs
    :type extract_id: int
    """
    # get the image urls
    img_urls = get_all_image_urls(ra, dec, size)
    ths = []
    # go through all available bands
    for b in img_urls.keys():
        # start a new thread to download every bands
        th = Thread(target=urllib.request.urlretrieve,
                    args=('http:' + img_urls[b],
                          save_path.format(extract_id, b),))
        th.start()
        ths.append(th)

    # wait until all bands are finished or until a timeout of 30sec is reached
    for t in ths:
        t.join(timeout=30)


def extract(path, band='r'):
    """
    Runs SourceExtractor with the band properties

    :param path: The path to the image
    :type path: str
    :param band: The name of the filter
    :type band: str
    :returns: The results from SExtractor
    :rtype: astropy.table.Table
    """
    return se[band](path)['table']


def extract_all_panstarrs(ra, dec, path, extract_id):
    """
    Extract the information with SExtractor from all PAN-STARRS images for this target.

    :param ra: The RA component of the coordinate of the target
    :type ra: float
    :param dec: The DEC component of the coordinate of the target
    :type dec: float
    :param path: Path for the image
    :type path: str
    :param extract_id: Unique id of the extract process
    :type extract_id: int
    :return: A table with all the extracted information for this source.
    :rtype: astropy.table.Table
    """
    # columns to collect
    cols = ['FWHM_IMAGE', 'FWHM_WORLD', 'KRON_RADIUS', 'PETRO_RADIUS',
            'THETA_IMAGE', 'A_IMAGE', 'B_IMAGE', 'A_WORLD', 'B_WORLD',
            'ELLIPTICITY']
    # functions to use
    funcs = [np.mean, np.median, np.std]
    stat = []

    # download all panstarrs images of this target
    download_all_bands(ra, dec, 3. / 60, path, extract_id)
    for b in se.keys():
        # extract the data from the image
        out = extract(path.format(extract_id, b), band=b)
        stat_row = [b]

        # run the functions on the data
        for c in cols:
            for f in funcs:
                stat_row.append(f(out[c]))

        # calculate the distance from the target coordinates to all detections
        # in the image
        r = np.hypot(ra - out['X_WORLD'], dec - out['Y_WORLD'])
        # if the minimal distance is larger than 6 arcsec
        if np.min(r) > 0.002:
            # ignore this source because the nearest source is to far away
            # (0.002 deg or 7.2 arcsec)
            continue

        # select the source which is the closest to the target
        p = np.where(r == np.min(r))[0][0]
        # add the data to the output
        for c in cols:
            stat_row.append(out[c][p])
        stat.append(tuple(stat_row))

    tab_cols = ['band']
    types = [('band', np.str_, 2)]
    # create the entries of the sources around the target
    for c in cols:
        for f in ['mean', 'median', 'std']:
            tab_cols.append(c + '_' + f)
            types.append((c + '_' + f, np.float64, 1))

    # create the entries of the target
    for c in cols:
        tab_cols.append(c + '_target')
        types.append((c + '_target', np.float64, 1))

    # convert the data to a astropy.table.Table
    stat = Table(data=np.array(stat, dtype=types))
    return stat


def get_times(dt):
    """
    Converts seconds into hours, minutes and seconds
    :param dt: The time in seconds
    :type dt: float
    :return: hours, minutes and seconds
    :rtype: int, int, float
    """
    mins = (dt - dt % 60) / 60
    dt = dt - 60 * mins
    hours = (mins - mins % 60) / 60
    mins = mins - 60 * hours
    return hours, mins, dt


class PanstarrsImages(Thread):
    def __init__(self, target_list, bands=None, extract_id=1, show_time=True):
        """
        PanstarrsImages provide the possibility to download images from the PAN-STARRS survey and estimate the
        properties of the source inside the image with SExtractor.
        Future features will include a binary detection from the images itself.

        :param target_list: Table with the targets
        :type target_list: astropy.table.Table
        :param bands: A list with the required bands, default are all bands of the PAN-STARRS survey
        :type bands: list
        :param extract_id:
            A unique ID for the extraction process. If you ran more than one process at the same time, you
            have to chose different ID's. Otherwise it will happen, that the same image is extracted for different
            sources.
        :param show_time: True if you want to see the needed time, else False
        :type show_time: bool
        """
        Thread.__init__(self)

        # if no bands are set, take all PAN-STARRS bands
        if bands is None:
            bands = ['g', 'r', 'i', 'z', 'y']
        self.target_list = target_list
        self.bands = bands
        self.extract_id = extract_id
        self.t0 = -1
        self.out = None
        self.show_times = show_time

    def time_estimation(self):
        """
        Estimates the time which was needed until now and the time which will be needed to finish the process.
        The results will be printed in the terminal.
        :return:
        """
        if self.show_times:
            length_out = len(self.out)

            # calculate the time since the process started
            dt = time.time() - self.t0
            hours, mins, seconds = get_times(dt)
            print('number of processed targets:', length_out)
            print('process runs since {}h {}min {} sec'.format(hours, mins, seconds))

            # estimate the time initial the process finished
            if length_out == 0:
                length_out = 1
            dt /= length_out
            dt = dt * (len(self.target_list) - length_out)
            hours, mins, seconds = get_times(dt)
            print('estimated time until the process finished {}h {}min {} sec'.format(hours, mins, seconds))

    def __do_extraction__(self, t, ra_name='ra', dec_name='dec'):
        try:
            # do the extraction of the panstarrs images
            pan = extract_all_panstarrs(t[ra_name], t[dec_name],
                                        TEMP_PANSTARRS_PATH + '{}_{}'.format(t[ra_name], t[dec_name]) +
                                        '_{:03d}_{}.fits',
                                        self.extract_id)
            # try the replace data
            try:
                # replace id, ra and dec by the input data
                pan['id'] = t['id']
                pan['ra'] = t[ra_name]
                pan['dec'] = t[dec_name]
                self.out.append(pan)
            except TypeError:
                pass
        except ValueError:
            pass

    def run(self):
        self.out = []

        # store the start time
        self.t0 = time.time()
        i = 0
        while i * 5 < len(self.target_list):
            targets = self.target_list[i * 5: (i + 1) * 5]
            th = []
            for t in targets:
                thread = Thread(target=self.__do_extraction__,
                                args=(t,))
                thread.start()
                th.append(thread)
            for t in th:
                t.join()

            self.time_estimation()

        # stack all data and store the results
        self.out = vstack(self.out)
