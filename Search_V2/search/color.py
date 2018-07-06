import numpy as np


class ColorSelection:
    cols = []
    center = []
    radius = -1
    comment_str = 'number of column: {}, number of left sources:{}'

    def __init__(self, cols, center, radius):
        """
        This class is the parent class of all color selection classes.
        """
        self.cols = cols
        self.center = center
        if len(self.cols) != len(self.center):
            raise ValueError('The number of cols must be the same as the' +
                             'number of centers')
        self.radius = radius

    def _process_(self, data, **kwargs):
        return None

    def __call__(self, *args, **kwargs):
        self._process_(args[0], **kwargs)


class NSphere(ColorSelection):
    def __init__(self, cols, center, radius):
        """
        N-dimensional sphere color selection class

        :param cols: Name of the cols which are included in the selection sphere
        :type cols: list
        :param center: Central coordinates which must have the same length as cols
        :type center: list
        :param radius: The selection radius
        :type radius: float
        """
        ColorSelection.__init__(self, cols, center, radius)

    def _process_(self, data, **kwargs):
        d = data.data
        cc = np.zeros(d.shape[0])
        for c, cen in zip(self.cols, self.center):
            cc += np.square(d[c] - cen)
        cc = np.sqrt(cc)
        mask = cc <= self.radius
        d = d[mask]
        data.set_data(d, 'n-sphere color selection',
                      self.comment_str.format(len(self.cols),
                                              len(d)))


class NBox(ColorSelection):
    def __init__(self, cols, center, radius):
        """
        N-dimensional box selection

        :param cols: Name of the cols which are included in the selection sphere
        :type cols: list
        :param center: Central coordinates which must have the same length as cols
        :type center: list
        :param radius: The selection radius/have of the box side
        :type radius: float
        """
        ColorSelection.__init__(self, cols, center, radius)

    def _process_(self, data, **kwargs):
        d = data.data
        for c, cen in zip(self.cols, self.center):
            cc = np.abs(d[c] - cen)
            mask = cc <= self.radius
            d = d[mask]
        data.set_data(d, 'n-box color selection',
                      self.comment_str.format(len(self.cols),
                                              len(d)))


class Grid:
    comment_str = 'number of column: {}, number of left sources:{}'
    cols = None
    grid = None
    x = None
    y = None

    def __init__(self, grid, x, y, cols):
        """
        Selection based on a grid which allows to include more complex shapes

        :param grid:
            A 2-dim array with values 0 or larger.
            Every cell with a value larger than zero will seen as
            an accepted place
        :type grid: numpy.ndarray
        :param x: The corresponding x-values of the grid
        :type x: numpy.ndarray
        :param y: The corresponding y-values of the grid
        :type y: numpy.ndarray
        :param cols: The name of the used cols
        :type cols: list
        """
        self.grid = grid
        self.x = x
        self.y = y
        self.cols = cols

    def __call__(self, data, *args, **kwargs):
        dx = __get_transformed_axis__(data, self.cols[0], self.x)
        dy = __get_transformed_axis__(data, self.cols[1], self.y)

        rs = self.grid[(dx, dy)]

        mask = rs > 0
        d = data[mask]

        data.set_data(d, 'Grid color selection',
                      self.comment_str.format(len(self.cols),
                                              len(d)))


def __get_transformed_axis__(data, c, x):
    """
    Transforms the values of data[c] to values between 0 and len(x)

    :param data: The input data set
    :type data: numpy.ndarray
    :param c: The name of the column
    :type c: str
    :param x: The values of the grid
    :type x: numpy.ndarray
    :return: The transformed column with integer values
    :rtype: numpy.ndarray
    """
    dx = data[c]
    dx -= x[0]
    dx /= x[-1] - x[0]
    dx *= x.shape[0]
    dx = np.round(x)
    return np.int32(dx)
