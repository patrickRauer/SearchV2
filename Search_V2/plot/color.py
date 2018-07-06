import pylab as pl

from Search_V2.core.data import Data


class ColorColor:
    data = None
    targets = None

    def __init__(self, data, targets=None):
        if type(data) == Data:
            self.data = data.data_original
            if len(data.data_original) != len(data.data):
                self.targets = data.data
        else:
            self.data = data
            self.targets = targets

    def color_color(self, c1, c2,
                    label_x='',
                    label_y=''):
        """
        Plots a color color histogram

        :param c1: The name of the first column
        :type c1: str
        :param c2: The name of the second column
        :type c2: str
        :param label_x: The optional x label. If it is empty, the label will be c1[0] - c1[-1]
        :type label_x: str
        :param label_y: The optional y label. If it is empty, the label will be c2[0] - c2[-1]
        :type label_y: str
        :return: 
        """
        pl.clf()
        sp = pl.subplot()
        h, x, y, scale = sp.hist2d(self.data[c1], self.data[c2],
                                   bins=200, cmap='gray_r')
        colb = pl.colorbar(scale)
        colb.set_label('counts')

        # targets are set too
        if self.targets is not None:
            sp.scatter(self.targets[c1], self.targets[c2],
                       marker='.')

        if label_x != '':
            sp.set_xlabel(label_x)
        else:
            sp.set_xlabel('{} - {}'.format(c1[0], c1[-1]))

        if label_y != '':
            sp.set_ylabel(label_y)
        else:
            sp.set_ylabel('{} - {}'.format(c2[0], c2[-1]))
        pl.show()

    def color_color_grid(self, cols, labels=None, bins=100):
        """
        Plots a color color grid and on the diagonal a histogram
        with the distribution of the colors are shown

        :param cols: The name of the color columns
        :type cols: list
        :param labels: Labels for the columns. Default is none which means that the labels will be c[0] - c[1]
        :type labels: list, None
        :param bins: The number of bins of the histograms
        :type bins: int
        :return:
        """
        if labels is None:
            labels = []
            for c in cols:
                labels.append('{} - {}'.format(c[0], c[-1]))
        length = len(cols)
        pl.clf()
        pl.subplots_adjust(wspace=0.0, hspace=0.0)
        for i, (c, l1) in enumerate(zip(cols, labels)):
            for j, (k, l2) in enumerate(zip(cols, labels)):
                if i > j:
                    continue
                # if the current plot is on the diagonal
                elif i == j:
                    sp = pl.subplot(length,
                                    length,
                                    i+length * j+1)
                    hist, x, patches = sp.hist(self.data[c],
                                               bins=bins, histtype='step', color='k')
                    if self.targets is not None:
                        sp2 = sp.twinx()
                        sp2.hist(self.targets[c], bins=bins, range=[x[0], x[-1]], color='b',
                                 histtype='step')
                        sp2.get_yaxis().set_visible(False)
                else:
                    sp = pl.subplot(length,
                                    length,
                                    i+length * j+1)
                    sp.hist2d(self.data[c], self.data[k],
                              bins=bins, cmap='gray_r')
                    if self.targets is not None:
                        sp.scatter(self.targets[c], self.targets[k],
                                   marker='.')

                if i == 0:
                    sp.set_ylabel(l2)
                else:
                    sp.get_yaxis().set_visible(False)
                if j == length-1:
                    sp.set_xlabel(l1)
                else:
                    sp.get_xaxis().set_visible(False)
        pl.subplots_adjust(wspace=0., hspace=0., top=0.9)
        pl.show()
