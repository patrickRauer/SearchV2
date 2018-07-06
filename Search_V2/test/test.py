from astropy.table import Table

from Search_V2.external.sdss import SdssSpec
from Search_V2.external.vizier import Gaia, SDSS, Wise, Apass, PanStarrs, MultiSurvey
from Search_V2.core.data import Data
from Search_V2.search.color import NSphere
from Search_V2.plot.color import ColorColor


def test_vizier():
    test_coords = Table()
    test_coords['ra'] = [50.837545, 47.233632]
    test_coords['dec'] = [-8.537334, -3.513184]

    gaia = Gaia()

    print('Gaia', gaia.query(test_coords))
    print('SDSS', SDSS().query(test_coords))
    print('Wise', Wise().query(test_coords))
    print('Apass', Apass().query(test_coords))
    print('PS', PanStarrs().query(test_coords))


def test_vizier2():
    test_coords = Table()
    test_coords['ra'] = [50.837545, 47.233632]
    test_coords['dec'] = [-8.537334, -3.513184]

    m = MultiSurvey()
    rs = m.query(test_coords)
    print(Table(rs))


def test_sdss_spec():
    sdss = SdssSpec()
    rs = sdss.download_classification()
    print(rs)
    # print(sdss.__stats__())
    sdss.show_sub_classes()


def test_data_system():
    data = Data(path='/Users/patrickr/Documents/GitHub/HCSC_Simulation/output/temp/cands_sim_-0_2.fits',
                backup=True)
    print(data)
    sphere = NSphere(['gr', 'gi', 'gz', 'gy'],
                     [1, 2, 2, 2],
                     0.75)
    sphere(data)

    cc = ColorColor(data)
    cc.color_color_grid(['gr', 'ri', 'iz', 'zy'])
    print(data.log)


if __name__ == '__main__':
    test_data_system()
