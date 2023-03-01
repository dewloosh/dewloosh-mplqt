# -*- coding: utf-8 -*-
from PyQt5 import QtWidgets, QtCore, QtGui
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg \
    as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as \
    NavigationToolbar
import matplotlib
matplotlib.use('Qt5Agg')


class MplFigure(QtWidgets.QMainWindow):

    def __init__(self, *args, parent=None, width=5, height=4, dpi=100, fig=None,
                 facecolor=None, name='mplWidget',  toolbar=True, **kwargs):
        super().__init__(parent)
        self.app = parent
        self.name = name
        if fig is not None:
            self.fig = fig
        else:
            if facecolor is None:
                _bcgclr = self.palette().color(QtGui.QPalette.Background)
                facecolor = [_bcgclr.red()/255, _bcgclr.green()/255,
                             _bcgclr.blue()/255]
                facecolor = [_bcgclr.red()/255, _bcgclr.green()/255,
                             _bcgclr.blue()/255]
            self.fig = Figure(figsize=(width, height), dpi=dpi,
                              edgecolor='black', facecolor=facecolor,
                              constrained_layout=False)
        self.frame = QtWidgets.QFrame()
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setFocus()
        self.canvas.setFocusPolicy(QtCore.Qt.ClickFocus)
        self.canvas.mpl_connect('pick_event', self._on_pick)
        self.canvas.mpl_connect('button_press_event', self._on_pick)
        self.canvas.mpl_connect("motion_notify_event", self._on_hover)
        if isinstance(toolbar, NavigationToolbar):
            self.toolbar = toolbar
        else:
            self.toolbar = NavigationToolbar(self.canvas, self.frame) \
                if toolbar else None
        self.fig.parentwindow = self

    def plot(self, *args, **kwargs):
        self.show()

    def show(self, *args, **kwargs):
        self._create_layout()
        self._build_plot(*args, **kwargs)
        self.canvas.draw()
        return super().show()

    @property
    def figure(self) -> Figure:
        """
        Returns the matplotlib figure.
        """
        return self.fig

    def _create_layout(self):
        self.layout = QtWidgets.QVBoxLayout()
        if self.toolbar:
            self.layout.addWidget(self.toolbar)
        self.layout.addWidget(self.canvas)
        self.canvas.setParent(self.frame)
        self.frame.setLayout(self.layout)
        self.setCentralWidget(self.frame)

    def _on_pick(self, event):
        """
        To be implemented.
        """
        pass

    def _on_hover(self, event):
        """
        To be implemented.
        """
        pass

    def _build_plot(self, *args, **kwargs):
        """
        To be implemented.
        """
        pass

    def save_to_pgf(self, filename: str = None, path: str = None):
        """
        Saves the content of the canvas to 'pgf' format.
        """
        if filename is None:
            filename = self.name
        filename = '.'.join(filename, 'pgf')
        w, h = self.fig.get_size_inches()
        w_ = 6.67
        h_ = w_ * h / w
        self.fig.set_size_inches(w_, h_)
        self.fig.savefig(filename)
        self.fig.set_size_inches(w, h)

    def _repr_latex_(self, doc=None):
        pass


if __name__ == '__main__':
    import numpy as np
    import sys

    np.random.seed(19680801)

    x = np.arange(-0.5, 10, 0.1)
    y = np.arange(4.5, 11, 0.1)
    Z = np.zeros((len(y), len(x)), dtype=np.float32)
    for i in range(10):
        Z += np.random.choice(a=np.arange(0, 1, 0.01),
                              size=(len(y), len(x)))

    app = QtWidgets.QApplication([])
    w = MplFigure()
    #w.show()

    fig = w.fig
    ax = fig.add_subplot()
    ax.pcolormesh(x, y, Z, cmap='Greys')
    #w.show()
    sys.exit(app.exec_())
