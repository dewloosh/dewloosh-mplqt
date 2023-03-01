from collections import OrderedDict, Iterable
from linkeddeepdict import LinkedDeepDict
from dewloosh.core.tools import float_to_str_sig as str_sig
from typing import TypeVar, Hashable, Callable
from dewloosh.core.typing import issequence
from typing import Iterable
from PyQt5 import QtGui, QtWidgets
from PyQt5.QtWidgets import QLineEdit, QLabel, QPushButton
from PyQt5.QtGui import QDoubleValidator
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider
from .plotwidget import MplFigure
import numpy as np
mpl.use('Qt5Agg')

__all__ = ['StressPlotter']


TScalar = TypeVar('TScalar', int, float, complex)
TReal = TypeVar('TReal', int, float)
TVector = Iterable[TScalar]
TRealVector = Iterable[TReal]
TColor = TypeVar('TColor', str, TRealVector)


scomp_to_tex = {
    'sxx': r"$\sigma_{xx}$",
    'syy': r"$\sigma_{yy}$",
    'sxy': r"$\sigma_{xy}$",
    'sxz': r"$\sigma_{xz}$",
    'syz': r"$\sigma_{yz}$",
    'sxx_e': r"$\sigma_{xx}^{\varepsilon}$",
    'syy_e': r"$\sigma_{yy}^{\varepsilon}$",
    'sxy_e': r"$\sigma_{xy}^{\varepsilon}$",
    'sxx_k': r"$\sigma_{xx}^{\kappa}$",
    'syy_k': r"$\sigma_{yy}^{\kappa}$",
    'sxy_k': r"$\sigma_{xy}^{\kappa}$",
}

    
class StressPlotter(MplFigure):

    def __init__(self, parent=None, sharelimits=False, usetex=False,
                    shell=None, separate=False, **kwargs):
        super().__init__(parent)
        mpl.rcParams['text.usetex'] = usetex
        mpl.rc('font', **{'family': "serif", 'size': 8})
        self.sharelimits = sharelimits
        self.shell = shell
        assert shell is not None, 'A shell instance must be provided!'
        if separate:
            self.scomps = ['sxx_e', 'syy_e', 'sxy_e', 'sxx_k', 'syy_k',
                            'sxy_k', 'sxz', 'syz']
        else:
            self.scomps = ['sxx', 'syy', 'sxy', 'sxz', 'syz']
        self.plotdata = LinkedDeepDict()
        self.setWindowTitle(kwargs.get('title', 'Surface Stresses'))
        self.activedatasetkey = None
        self.activeplotkey = None
        self.xticksrotation = 30

    def _create_layout(self):
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)

        _editlabel = QLabel(self)
        _editlabel.setText("z = ")
        editlabel_ = QLabel(self)
        editlabel_.setText("[m]")
        self.edit = QLineEdit(self)
        self.validator = QDoubleValidator()
        # self.edit.setValidator(self.validator)
        self.submitbutton = QPushButton(self)
        self.submitbutton.setText("Submit")
        self.submitbutton.clicked.connect(self._submit)
        self.equalbutton = QPushButton(self)
        self.equalbutton.setText("Equal Axes")
        self.equalbutton.setCheckable(True)
        self.equalbutton.clicked.connect(self._axes_range_toggle)
        self.resetbutton = QPushButton(self)
        self.resetbutton.setText("Reset")
        self.resetbutton.clicked.connect(self._reset_axes)
        hlayout = QtWidgets.QHBoxLayout()
        hlayout.addWidget(self.equalbutton)
        hlayout.addWidget(self.resetbutton)
        hlayout.addWidget(_editlabel)
        hlayout.addWidget(self.edit)
        hlayout.addWidget(editlabel_)
        hlayout.addWidget(self.submitbutton)
        layout.addLayout(hlayout)

        self.frame.setLayout(layout)
        self.setCentralWidget(self.frame)

    def _submit(self):
        try:
            z = float(self.edit.text().replace(',', '.'))
            self.slider.set_val(z)
            return True
        except Exception:
            return False

    def _axes_range_toggle(self):
        self.sharelimits = self.equalbutton.isChecked()
        self._set_axes_limits()

    def _set_axes_limits(self):
        vmin, vmax = 1e12, -1e12
        for plotdata in self.plotdata.values():
            all_lines = []
            for data in plotdata.values():
                if isinstance(data, dict) and 'lines' in data:
                    all_lines.append(data['lines'])

            visible_lines = filter(lambda l: l.get_visible(), all_lines)
            data = np.array([line.get_xdata() for line in visible_lines],
                            dtype=np.float32)
            _vmin, _vmax = np.amin(data), np.amax(data)
            if not self.sharelimits:
                self._set_xlim(plotdata['ax'], _vmin, _vmax)
            else:
                if _vmin < vmin:
                    vmin = _vmin
                if _vmax > vmax:
                    vmax = _vmax
        if self.sharelimits:
            for plotdata in self.plotdata.values():
                self._set_xlim(plotdata['ax'], vmin, vmax)
        self.canvas.draw_idle()

    def _reset_axes(self):
        self.fig.suptitle('')
        self.activedatasetkey = None
        self.slider.set_val(0.0)
        self._set_axes_limits()

    def _on_pick(self, event):
        """
        Selects datasets on clicking.
        """
        if hasattr(event, 'artist'):
            lines2D = event.artist
            plotkey, datasetkey = lines2D.plotkey, lines2D.datasetkey
            self.activedatasetkey = datasetkey
            for pkey in self.plotdata.keys():
                if 'ax' in self.plotdata[pkey]:
                    for key, data in self.plotdata[pkey].items():
                        if isinstance(data, dict) and 'lines' in data:
                            if key == datasetkey:
                                data['lines'].set_alpha(1.0)
                                data['lines'].set_linewidth(2.0)
                            else:
                                data['lines'].set_alpha(0.2)
                                data['lines'].set_linewidth(1.5)
        else:
            if event.dblclick:
                self.activedatasetkey = None
                for plotkey in self.plotdata.keys():
                    if 'ax' in self.plotdata[plotkey]:
                        for key, data in self.plotdata[plotkey].items():
                            if isinstance(data, dict) and 'lines' in data:
                                data['lines'].set_alpha(data['alpha'])
                                data['lines'].set_linewidth(
                                    data['thickness'])
        self._update_slider()

    def _on_hover(self, event):
        """
        Draws border around the axis in which the cursor is hovering.
        """
        ax = event.inaxes
        if ax is not None and hasattr(ax, 'plotkey'):
            if not ax.plotkey == self.activeplotkey:
                self.activeplotkey = ax.plotkey
                for key in self.plotdata.keys():
                    if 'ax' in self.plotdata[key]:
                        if key == self.activeplotkey:
                            self.plotdata[key]['ax'].patch.set_linewidth(
                                '2')
                            self.plotdata[key]['ax'].patch.set_edgecolor(
                                'deepskyblue')
                        else:
                            self.plotdata[key]['ax'].patch.set_linewidth(
                                '0.5')
                            self.plotdata[key]['ax'].patch.set_edgecolor(
                                'black')
        else:
            self.activeplotkey = None
            for key in self.plotdata.keys():
                if 'ax' in self.plotdata[key]:
                    self.plotdata[key]['ax'].patch.set_linewidth('0.5')
                    self.plotdata[key]['ax'].patch.set_edgecolor('black')
        self.canvas.draw_idle()

    def add_dataset(self, axisID: Hashable = None,
                    datasetID: Hashable = None, values: TRealVector = None,
                    locations: TRealVector = None, *args,
                    color: TColor = 'black', thickness: TReal = 1.0,
                    fill: bool = False, zorder: int = 0, alpha: float = 1.0,
                    **kwargs):
        assert isinstance(axisID, Hashable), "axisID must be hashable!"
        assert isinstance(
            datasetID, Hashable), "datasetID must be hashable!"
        self.plotdata[axisID][datasetID]['values'] = values
        self.plotdata[axisID][datasetID]['locations'] = locations
        self.plotdata[axisID][datasetID]['color'] = \
            self._color_to_values(color)
        self.plotdata[axisID][datasetID]['thickness'] = thickness
        self.plotdata[axisID][datasetID]['fill'] = fill
        self.plotdata[axisID][datasetID]['zorder'] = zorder
        self.plotdata[axisID][datasetID]['alpha'] = alpha

    def add_slave(self, axisID: Hashable = None, func: Callable = None,
                    *args, sourceIDs: Iterable = None,
                    locations: TRealVector = None, **kwargs):
        assert func is not None, "'func' must be specified!"
        if sourceIDs is None:
            sourceIDs = self.scomps
        for datasetID in self._common_datasetkeys(sourceIDs=sourceIDs):
            datavalues = {}
            for plotID in sourceIDs:
                datavalues[plotID] = \
                    self.plotdata[plotID][datasetID]['values']  # FIXME
            values = func(**datavalues)
            self.add_dataset(
                axisID=axisID, datasetID=datasetID, values=values)

    def _locations_of_dataset(self, datasetID: Hashable = None):
        assert datasetID is not None, "'datasetID must be provided!"
        plotkey = None
        nmax = -1
        for key, value in self.plotdata.items():
            if isinstance(value, dict) and datasetID in value:
                locations = self.plotdata[key][datasetID]['locations']
                nPoint = len(locations)
                if nPoint > nmax:
                    nmax = nPoint
                    plotkey = key
        return self.plotdata[plotkey][datasetID]['locations']

    def _common_datasetkeys(self, sourceIDs: Iterable = None) -> set:
        """
        Returns the keys of the datasets all axes have values for.
        """
        if sourceIDs is None:
            sourceIDs = self.scomps
        sets = []
        for plotkey in sourceIDs:
            if plotkey in self.plotdata:
                plotdata = self.plotdata[plotkey]
                sets.append(set())
                for key, data in plotdata.items():
                    if isinstance(data, dict) and 'values' in data:
                        sets[-1].add(key)
        return set.intersection(*sets)

    def _color_to_values(self, color: TColor = None) -> Iterable:
        """
        Returns the RGB equivalent of a color specification, or None
        if unable the recognize a valid input.
        """
        if issequence(color):
            if np.max(color) > 1.0:
                return (c/255 for c in color)
            else:
                return color
        elif isinstance(color, str):
            return mpl.colors.to_rgb(color)
        return None

    def _init_plot_data(self):
        self.layers = self.shell.layers()
        self.thickness = self.shell.tmax - self.shell.tmin
        yticks = [la.zi[0] for la in self.layers]
        yticks.append(self.layers[-1].zi[-1])
        self.yticks = np.array(yticks, dtype=np.float32)
        self.validator.setRange(self.yticks[-1], self.yticks[0], 4)

        # set limits for the datasets and store valid keys
        self.keys = OrderedDict()
        for plotkey, plotdata in self.plotdata.items():
            values, keys = [], []
            for key, data in plotdata.items():
                if isinstance(data, dict):
                    if 'values' in data:
                        keys.append(key)
                        values.append(data['values'])
                        self.keys[plotkey] = None
            a = np.vstack(values)
            minvalue = a.min()
            argmin = np.where(a == minvalue)
            plotdata['min'] = minvalue
            plotdata['minkey'] = keys[argmin[0][0]]
            maxvalue = a.max()
            argmax = np.where(a == maxvalue)
            plotdata['max'] = maxvalue
            plotdata['maxkey'] = keys[argmax[0][0]]
        self.keys = list(self.keys.keys())

        # set global min and max
        self.vmin = np.min([data['min']
                            for data in self.plotdata.values()])
        self.vmax = np.max([data['max']
                            for data in self.plotdata.values()])

    def _update_slider(self, z=None, xvalues=None):
        if z is None:
            z = self.slider.val
        for key in self.plotdata.keys():
            if 'hline' in self.plotdata[key]:
                self.plotdata[key]['hline'].set_ydata(z)
        if self.activedatasetkey is not None:
            for key in self.plotdata.keys():
                if self.activedatasetkey in self.plotdata[key]:
                    v_at_z = self._approx_at_z(
                        z, key, self.activedatasetkey)
                    self.plotdata[key]['text'].update({'visible': True,
                                                        'x': v_at_z, 'y': z,
                                                        'text':
                                                        str_sig(v_at_z)})
                else:
                    self.plotdata[key]['text'].update({'visible': False})
        else:
            for key in self.plotdata.keys():
                if 'text' in self.plotdata[key]:
                    self.plotdata[key]['text'].update({'visible': False})
        self.canvas.draw_idle()

    def _approx_at_z(self, z: TRealVector, plotkey: Hashable,
                        datasetkey: Hashable) -> TReal:
        lines2D = self.plotdata[plotkey][datasetkey]['lines']
        values = lines2D.get_xdata()
        locations = lines2D.get_ydata()
        return np.interp(z, locations, values)

    def _set_xlim(self, axs: mpl.axes, vmin: float, vmax: float,
                    offset=0.2, **kwargs):
        voffset = (vmax-vmin)*offset
        if abs(vmin-vmax) > 1e-7:
            axs.set_xlim(vmin - voffset, vmax + voffset)
        xticks = [vmin, vmax]
        axs.set_xticks(xticks)
        rotation = kwargs.get('rotation', self.xticksrotation)
        axs.set_xticklabels([str_sig(val) for val in xticks],
                            rotation=rotation)

    def _set_ylim(self, axs: mpl.axes, vmin: float, vmax: float,
                    offset=0.1, **kwargs):
        voffset = (vmax-vmin)*offset
        axs.set_ylim(vmin - voffset, vmax + voffset)

    def _build_plot(self, *args, **kwargs):
        # init
        self._init_plot_data()
        axcolor = 'lightgoldenrodyellow'
        boldfont = QtGui.QFont()
        boldfont.setBold(True)
        self.fig.clf()  # clear
        nAxes = len(self.keys) + 1  # +1 for the Slider
        width_ratios = [1 for i in range(nAxes)]
        width_ratios[-1] = 0.15
        spec = gridspec.GridSpec(ncols=nAxes, nrows=1,
                                    width_ratios=width_ratios, figure=self.fig,
                                    wspace=0.2, left=0.1)
        self.fig.suptitle('')
        # self.fig.subplots_adjust(left=0.1, bottom=0.2, right = 1.0, t
        # op = 0.85, wspace = 0.1)
        tmin, tmax = self.shell.tmin, self.shell.tmax

        # create axes
        for i, key in enumerate(self.keys):
            plotid = int("{}{}{}".format(1, nAxes, i+1))
            plotid = spec[0, i]
            ax = self.fig.add_subplot(plotid, facecolor=axcolor)
            ax.grid(False)
            ax.plotkey = key
            ax.patch.set_edgecolor('black')
            ax.patch.set_linewidth('0.5')
            if i == 0:
                ax.set_yticks(self.yticks)
                ax.set_yticklabels(["{:.3f}".format(val)
                                    for val in self.yticks])
            else:
                ax.set_yticks([])
                ax.set_yticklabels([])
            self.plotdata[key]['hline'] = \
                ax.axhline(y=0.0, color='#d62728', linewidth=1)
            bbox = dict(boxstyle="round", ec='black',
                        fc='yellow', alpha=0.8)
            self.plotdata[key]['text'] = ax.text(0.0, 0.0, "NaN", size=10,
                                                    ha="center", va="center",
                                                    visible=False, bbox=bbox)
            # layer lines
            for layer_z in self.yticks:
                ax.axhline(y=layer_z, color='black', linewidth=0.5,
                            linestyle='--')
            ax.axvline(x=0.0, color='black', linewidth=0.5, linestyle='-')
            self.plotdata[key]['ax'] = ax

        # create slider @ add_axes([left, bottom, width, height],*args,**kwargs)
        # self.sliderax = self.fig.add_axes([0.9, 0.2, 0.02, 0.65],
        # facecolor=axcolor)
        self.sliderax = self.fig.add_subplot(
            spec[0, nAxes-1], facecolor=axcolor)
        self.slider = Slider(self.sliderax, 'z [m]', valmin=tmin,
                                valmax=tmax, valinit=0.0,
                                orientation='vertical', valfmt="%.3f",
                                closedmin=True, closedmax=True)
        self.slider.on_changed(self._update_slider)

        # plot axes
        for i, plotkey in enumerate(self.keys):
            axis = self.plotdata[plotkey]['ax']
            if self.sharelimits == True:
                self._set_xlim(axis, self.vmin, self.vmax)
            else:
                self._set_xlim(axis, self.plotdata[plotkey]['min'],
                                self.plotdata[plotkey]['max'])
            self._set_ylim(axis, tmin, tmax)
            if plotkey in scomp_to_tex:
                axis.set_title(scomp_to_tex[self.scomps[i]])
            else:
                axis.set_title(str(plotkey))

            for key, data in self.plotdata[plotkey].items():
                if isinstance(data, dict):
                    if 'values' in data:
                        lines = axis.plot(data['values'], data['locations'],
                                            color=data['color'],
                                            linewidth=data['thickness'],
                                            zorder=data['zorder'],
                                            picker=5, alpha=data['alpha'])[0]
                        lines.plotkey = plotkey
                        lines.datasetkey = key
                        data['lines'] = lines
