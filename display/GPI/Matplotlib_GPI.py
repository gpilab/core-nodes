# Copyright (c) 2014, Dignity Health
# 
#     The GPI core node library is licensed under
# either the BSD 3-clause or the LGPL v. 3.
# 
#     Under either license, the following additional term applies:
# 
#         NO CLINICAL USE.  THE SOFTWARE IS NOT INTENDED FOR COMMERCIAL
# PURPOSES AND SHOULD BE USED ONLY FOR NON-COMMERCIAL RESEARCH PURPOSES.  THE
# SOFTWARE MAY NOT IN ANY EVENT BE USED FOR ANY CLINICAL OR DIAGNOSTIC
# PURPOSES.  YOU ACKNOWLEDGE AND AGREE THAT THE SOFTWARE IS NOT INTENDED FOR
# USE IN ANY HIGH RISK OR STRICT LIABILITY ACTIVITY, INCLUDING BUT NOT LIMITED
# TO LIFE SUPPORT OR EMERGENCY MEDICAL OPERATIONS OR USES.  LICENSOR MAKES NO
# WARRANTY AND HAS NOR LIABILITY ARISING FROM ANY USE OF THE SOFTWARE IN ANY
# HIGH RISK OR STRICT LIABILITY ACTIVITIES.
# 
#     If you elect to license the GPI core node library under the LGPL the
# following applies:
# 
#         This file is part of the GPI core node library.
# 
#         The GPI core node library is free software: you can redistribute it
# and/or modify it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version. GPI core node library is distributed
# in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
# the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU Lesser General Public License for more details.
# 
#         You should have received a copy of the GNU Lesser General Public
# License along with the GPI core node library. If not, see
# <http://www.gnu.org/licenses/>.


# Author: Nick Zwart
# Date: 2013 Oct 30
from __future__ import print_function
import os
import matplotlib

import gpi
from gpi import QtCore, QtGui

import numpy as np
from matplotlib.figure import Figure
#from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_qt4agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
from matplotlib.backends.backend_qt4 import SubplotToolQt

class MainWin_close(QtGui.QMainWindow):
    window_closed = gpi.Signal()
    def __init__(self):
        super(MainWin_close, self).__init__()

    def closeEvent(self, event):
        super(MainWin_close, self).closeEvent(event)
        self.window_closed.emit()


###############################################################################
# -*- coding: utf-8 -*-
#
# Copyright (C) 2009 Pierre Raybau
# Licensed under the terms of the MIT License
# see the mpl licenses directory for a copy of the license

"""Module that provides a GUI-based editor for matplotlib's figure options"""

#from __future__ import print_function
import os.path as osp

import matplotlib.backends.qt4_editor.formlayout as formlayout
#from matplotlib.backends.qt4_compat import QtGui
from matplotlib import markers

def get_icon(name):
    import matplotlib
    basedir = osp.join(matplotlib.rcParams['datapath'], 'images')
    return QtGui.QIcon(osp.join(basedir, name))

LINESTYLES = {
              '-': 'Solid',
              '--': 'Dashed',
              '-.': 'DashDot',
              ':': 'Dotted',
              'steps': 'Steps',
              'none': 'None',
              }

MARKERS = markers.MarkerStyle.markers

COLORS = {'b': '#0000ff', 'g': '#00ff00', 'r': '#ff0000', 'c': '#ff00ff',
          'm': '#ff00ff', 'y': '#ffff00', 'k': '#000000', 'w': '#ffffff'}

def col2hex(color):
    """Convert matplotlib color to hex"""
    return COLORS.get(color, color)

def figure_edit(axes, parent=None):
    """Edit matplotlib figure options"""
    sep = (None, None) # separator

    has_curve = len(axes.get_lines()) > 0

    # Get / General
    #xmin, xmax = axes.get_xlim()
    #ymin, ymax = axes.get_ylim()
    #general = [('Title', axes.get_title()),
    #           sep,
    #           (None, "<b>X-Axis</b>"),
    #           ('Min', xmin), ('Max', xmax),
    #           ('Label', axes.get_xlabel()),
    #           ('Scale', [axes.get_xscale(), 'linear', 'log']),
    #           sep,
    #           (None, "<b>Y-Axis</b>"),
    #           ('Min', ymin), ('Max', ymax),
    #           ('Label', axes.get_ylabel()),
    #           ('Scale', [axes.get_yscale(), 'linear', 'log'])
    #           ]

    if not has_curve:
        return

    if has_curve:
        # Get / Curves
        linedict = {}
        for line in axes.get_lines():
            label = line.get_label()
            if label == '_nolegend_':
                continue
            linedict[label] = line
        curves = []
        linestyles = LINESTYLES.items()
        markers = MARKERS.items()
        curvelabels = sorted(linedict.keys())
        for label in curvelabels:
            line = linedict[label]
            curvedata = [
                         ('Label', label),
                         sep,
                         (None, '<b>Line</b>'),
                         ('Style', [line.get_linestyle()] + linestyles),
                         ('Width', line.get_linewidth()),
                         ('Color', col2hex(line.get_color())),
                         sep,
                         (None, '<b>Marker</b>'),
                         ('Style', [line.get_marker()] + markers),
                         ('Size', line.get_markersize()),
                         ('Facecolor', col2hex(line.get_markerfacecolor())),
                         ('Edgecolor', col2hex(line.get_markeredgecolor())),
                         ]
            curves.append([curvedata, label, ""])

    #datalist = [(general, "Axes", "")]
    #if has_curve:
    #datalist.append((curves, "Curves", ""))
    datalist = [(curves, "Curves", "")]
        
    def apply_callback(data):
        """This function will be called to apply changes"""
        #if has_curve:
        #    general, curves = data
        #else:
        #    general, = data

        curves = data[0]
            
        # Set / General
        #title, xmin, xmax, xlabel, xscale, ymin, ymax, ylabel, yscale = general
        #axes.set_xscale(xscale)
        #axes.set_yscale(yscale)
        #axes.set_title(title)
        #axes.set_xlim(xmin, xmax)
        #axes.set_xlabel(xlabel)
        #axes.set_ylim(ymin, ymax)
        #axes.set_ylabel(ylabel)
        
        if has_curve:
            # Set / Curves
            for index, curve in enumerate(curves):
                line = linedict[curvelabels[index]]
                label, linestyle, linewidth, color, \
                    marker, markersize, markerfacecolor, markeredgecolor = curve
                line.set_label(label)
                line.set_linestyle(linestyle)
                line.set_linewidth(linewidth)
                line.set_color(color)
                if marker is not 'none':
                    line.set_marker(marker)
                    line.set_markersize(markersize)
                    line.set_markerfacecolor(markerfacecolor)
                    line.set_markeredgecolor(markeredgecolor)
            
        # Redraw
        figure = axes.get_figure()
        figure.canvas.draw()
        
    data = formlayout.fedit(datalist, title="Figure options", parent=parent,
                            icon=get_icon('qt4_editor_options.svg'), apply=apply_callback)
    if data is not None:
        apply_callback(data)
###############################################################################   


class MatplotDisplay(gpi.GenericWidgetGroup):

    """Embeds the matplotlib figure window.
    """
    valueChanged = gpi.Signal()

    def __init__(self, title, parent=None):
        super(MatplotDisplay, self).__init__(title, parent)

        # gpi interface
        self._collapsables = []
        self._subplotSettings = {}
        self._subplotPosition = {'right': 0.9, 'bottom': 0.12, 'top': 0.9, 'wspace': 0.2, 'hspace': 0.2, 'left': 0.125}
        #self._subplot_keepers = ['yscale', 'xscale']
        self._subplot_keepers = []
        self._lineSettings = []
        self._line_keepers = ['linewidth', 'linestyle', 'label', 'marker', 'markeredgecolor', 'markerfacecolor', 'markersize', 'color', 'alpha']

        # since drawing is slow, don't do it as often, use the timer as a
        # debouncer
        self._on_draw_cnt = 0
        self._updatetimer = QtCore.QTimer()
        self._updatetimer.setSingleShot(True)
        self._updatetimer.timeout.connect(self._on_draw)
        self._updatetimer.setInterval(10)

        # plot specific UI side panel
        #  -sets options for plot window so this needs to be run first
        vbox = QtGui.QVBoxLayout()
        vbox.setContentsMargins(0, 0, 0, 0)  # no spaces around this item
        vbox.setSpacing(0)

        # AUTOSCALE
        self._autoscale_btn = gpi.widgets.BasicPushButton(self)
        self._autoscale_btn.set_toggle(True)
        self._autoscale_btn.set_button_title('autoscale')
        self._autoscale_btn.valueChanged.connect(self.on_draw)
        self._collapsables.append(self._autoscale_btn)

        # GRID
        self._grid_btn = gpi.widgets.BasicPushButton(self)
        self._grid_btn.set_toggle(True)
        self._grid_btn.set_button_title('grid')
        self._grid_btn.valueChanged.connect(self.on_draw)
        self._collapsables.append(self._grid_btn)

        # X/Y LIMITS
        lims = QtGui.QGridLayout()
        self._xl = gpi.widgets.BasicDoubleSpinBox(self)
        self._xh = gpi.widgets.BasicDoubleSpinBox(self)
        self._yl = gpi.widgets.BasicDoubleSpinBox(self)
        self._yh = gpi.widgets.BasicDoubleSpinBox(self)
        self._xl.valueChanged.connect(self.on_draw)
        self._xh.valueChanged.connect(self.on_draw)
        self._yl.valueChanged.connect(self.on_draw)
        self._yh.valueChanged.connect(self.on_draw)
        self._xl.set_immediate(True)
        self._xh.set_immediate(True)
        self._yl.set_immediate(True)
        self._yh.set_immediate(True)
        self._xl.set_label('max')
        self._xh.set_label('min')
        self._xl.set_decimals(7)
        self._xh.set_decimals(7)
        self._yl.set_decimals(7)
        self._yh.set_decimals(7)
        self._xlab = QtGui.QLabel('x limits')
        self._ylab = QtGui.QLabel('y limits')
        #self._maxlab = QtGui.QLabel('max')
        #self._minlab = QtGui.QLabel('min')
        #lims.addWidget(self._maxlab,1,0,1,1)
        #lims.addWidget(self._minlab,2,0,1,1)
        lims.addWidget(self._xlab,0,1,1,1,alignment=QtCore.Qt.AlignHCenter)
        lims.addWidget(self._xh,1,1,1,1,alignment=QtCore.Qt.AlignHCenter)
        lims.addWidget(self._xl,2,1,1,1,alignment=QtCore.Qt.AlignHCenter)
        lims.addWidget(self._ylab,0,2,1,1,alignment=QtCore.Qt.AlignHCenter)
        lims.addWidget(self._yh,1,2,1,1,alignment=QtCore.Qt.AlignHCenter)
        lims.addWidget(self._yl,2,2,1,1,alignment=QtCore.Qt.AlignHCenter)
        self._collapsables.append(self._xlab)
        self._collapsables.append(self._ylab)
        self._collapsables.append(self._xl)
        self._collapsables.append(self._xh)
        self._collapsables.append(self._yl)
        self._collapsables.append(self._yh)
        #self._collapsables.append(self._minlab)
        #self._collapsables.append(self._maxlab)

        # TICK MARKS
        ticks = QtGui.QGridLayout()
        self._x_numticks = gpi.widgets.BasicSpinBox(self)
        self._x_numticks.valueChanged.connect(self.on_draw)
        self._y_numticks = gpi.widgets.BasicSpinBox(self)
        self._y_numticks.valueChanged.connect(self.on_draw)
        self._x_ticks = QtGui.QLineEdit()
        self._y_ticks = QtGui.QLineEdit()
        self._x_ticks.textChanged.connect(lambda txt: self.check_validticks(self._x_ticks))
        self._y_ticks.textChanged.connect(lambda txt: self.check_validticks(self._y_ticks))
        self._x_ticks.setPlaceholderText('comma separated list of x labels')
        self._y_ticks.setPlaceholderText('comma separated list of y labels')
        self._x_ticks.returnPressed.connect(self.on_draw)
        self._y_ticks.returnPressed.connect(self.on_draw)
        self._x_numticks.set_immediate(True)
        self._y_numticks.set_immediate(True)
        self._x_numticks.set_min(2)
        self._y_numticks.set_min(2)
        self._x_numticks.set_max(100)
        self._y_numticks.set_max(100)
        self._x_numticks.set_val(5)
        self._y_numticks.set_val(5)
        self._x_numticks.set_label('x ticks')
        self._y_numticks.set_label('y ticks')
        ticks.addWidget(self._x_numticks, 0,0,1,1)
        ticks.addWidget(self._y_numticks, 1,0,1,1)
        ticks.addWidget(self._x_ticks, 0,1,1,1)
        ticks.addWidget(self._y_ticks, 1,1,1,1)
        self._collapsables.append(self._x_numticks)
        self._collapsables.append(self._y_numticks)
        self._collapsables.append(self._x_ticks)
        self._collapsables.append(self._y_ticks)

        # TITLE, XLABEL, YLABEL
        plotlabels = QtGui.QHBoxLayout()
        self._plot_title = QtGui.QLineEdit()
        self._plot_xlab = QtGui.QLineEdit()
        self._plot_ylab = QtGui.QLineEdit()
        self._plot_title.setPlaceholderText('title')
        self._plot_xlab.setPlaceholderText('x label')
        self._plot_ylab.setPlaceholderText('y label')
        self._plot_title.returnPressed.connect(self.on_draw)
        self._plot_xlab.returnPressed.connect(self.on_draw)
        self._plot_ylab.returnPressed.connect(self.on_draw)
        plotlabels.addWidget(self._plot_title)
        plotlabels.addWidget(self._plot_xlab)
        plotlabels.addWidget(self._plot_ylab)
        self._collapsables.append(self._plot_title)
        self._collapsables.append(self._plot_xlab)
        self._collapsables.append(self._plot_ylab)

        # LEGEND
        self._legend_btn = gpi.widgets.BasicPushButton(self)
        self._legend_btn.set_toggle(True)
        self._legend_btn.set_button_title('legend')
        self._legend_btn.valueChanged.connect(self.on_draw)
        self._collapsables.append(self._legend_btn)

        # LINE OPTIONS
        self._lino_btn = gpi.widgets.BasicPushButton(self)
        self._lino_btn.set_toggle(False)
        self._lino_btn.set_button_title('line options')
        self._lino_btn.valueChanged.connect(self.lineOptionsDialog)
        self._collapsables.append(self._lino_btn)

        # SUBPLOT SPACING OPTIONS
        self._subplotso_btn = gpi.widgets.BasicPushButton(self)
        self._subplotso_btn.set_toggle(False)
        self._subplotso_btn.set_button_title('spacing options')
        self._subplotso_btn.valueChanged.connect(self.subplotSpacingOptions)
        self._collapsables.append(self._subplotso_btn)

        plot_options_layout = QtGui.QHBoxLayout()
        plot_options_layout.addWidget(self._subplotso_btn)
        plot_options_layout.addWidget(self._lino_btn)

        # panel layout
        vbox.addLayout(plotlabels)
        vbox.addWidget(self._grid_btn)
        vbox.addLayout(lims)
        vbox.addWidget(self._autoscale_btn)
        vbox.addLayout(ticks)
        vbox.addWidget(self._legend_btn)
        vbox.addLayout(plot_options_layout)
        #vbox.addWidget(self._lino_btn)
        #vbox.addWidget(self._subplotso_btn)
        vbox.insertStretch(-1,1)

        # plot window
        self._data = None
        self._plotwindow = self.create_main_frame()

        # put side panel and plot window together
        hbox = QtGui.QHBoxLayout()
        hbox.addLayout(vbox)
        hbox.addLayout(self._plotwindow)
        hbox.setStretch(0,0)
        hbox.setStretch(1,11)
        self.setLayout(hbox)

        #self._on_draw() # draw once to get initial settings
        #self.copySubplotSettings()

        # Don't hide side-panel options by default
        self.set_collapsed(False)
        self.set_grid(True)
        self.set_autoscale(True)

    # setters
    def set_val(self, data):
        '''Takes a list of npy arrays.
        '''
        if isinstance(data, list):
            self._data = data
            self.on_draw()

    def set_grid(self, val):
        self._grid_btn.set_val(val)
        self.on_draw()

    def set_autoscale(self, val):
        self._autoscale_btn.set_val(val)
        self.on_draw()

    def set_collapsed(self, val):
        """bool | Only collapse the display options, not the Plot window.
        """
        self._isCollapsed = val
        for wdg in self._collapsables:
            if hasattr(wdg, 'setVisible'):
                wdg.setVisible(not val)

    def set_xlim(self, val, quiet=False):
        '''tuple of floats (min, max)
        '''
        self._xh.set_val(val[0])
        self._xl.set_val(val[1])
        if not quiet:
            self.on_draw()

    def set_ylim(self, val, quiet=False):
        '''tuple of floats (min, max)
        '''
        self._yh.set_val(val[0])
        self._yl.set_val(val[1])
        if not quiet:
            self.on_draw()

    def set_plotOptions(self, val):
        self._subplotSettings = val

    def set_lineOptions(self, val):
        self._lineSettings = val

    def set_plotPosition(self, val):
        self._subplotPosition = val

    def set_ticks(self, s):
        self._x_numticks.set_val(s['xticknum'])
        self._y_numticks.set_val(s['yticknum'])
        self._x_ticks.setText(s['xticks'])
        self._y_ticks.setText(s['yticks'])

    def set_plotlabels(self, s):
        self._plot_title.setText(s['title'])
        self._plot_xlab.setText(s['xlabel'])
        self._plot_ylab.setText(s['ylabel'])
        self.on_draw()

    def set_legend(self, val):
        self._legend_btn.set_val(val)
        self.on_draw()

    # getters
    def get_val(self):
        return self._data

    def get_grid(self):
        return self._grid_btn.get_val()

    def get_autoscale(self):
        return self._autoscale_btn.get_val()

    def get_xlim(self):
        return (self._xh.get_val(), self._xl.get_val())

    def get_ylim(self):
        return (self._yh.get_val(), self._yl.get_val())

    def get_plotOptions(self):
        return self._subplotSettings

    def get_lineOptions(self):
        return self._lineSettings

    def get_plotPosition(self):
        return self._subplotPosition

    def get_ticks(self):
        s = {}
        s['xticknum'] = self._x_numticks.get_val()
        s['yticknum'] = self._y_numticks.get_val()
        s['xticks'] = str(self._x_ticks.text())
        s['yticks'] = str(self._y_ticks.text())
        return s

    def get_plotlabels(self):
        s = {}
        s['title'] = str(self._plot_title.text())
        s['xlabel'] = str(self._plot_xlab.text())
        s['ylabel'] = str(self._plot_ylab.text())
        return s

    def get_legend(self):
        return self._legend_btn.get_val()

    # support
    def check_validticks(self, tickwdg):
        s = tickwdg.text()
        comma_cnt = s.count(',')
        if (comma_cnt > 0) or (len(s) == 0):
            color = '#ffffff' # white
            tickwdg.setStyleSheet('QLineEdit { background-color: %s }' % color)
            return

        #color = '#fff79a' # yellow
        color = '#f6989d' # red
        tickwdg.setStyleSheet('QLineEdit { background-color: %s }' % color)
        return

    def create_main_frame(self):
        self.fig = Figure((6.0, 4.8), dpi=100, facecolor='0.98', linewidth=6.0, edgecolor='0.93')
        self.axes = None
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self)
        self.canvas.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.canvas.setFocus()

        self.canvas.setSizePolicy(QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.MinimumExpanding)

        self.mpl_toolbar = NavigationToolbar(self.canvas, self)
        self.mpl_toolbar.actionTriggered.connect(self.copySubplotSettings)

        self.canvas.mpl_connect('key_press_event', self.on_key_press)

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.canvas)  # the matplotlib canvas
        vbox.addWidget(self.mpl_toolbar)
        return vbox
        #self.setLayout(vbox)

    def lineOptionsDialog(self):
        if self.axes is None:
            return
        figure_edit(self.axes, self)
        self.copySubplotSettings()
        self.on_draw()

    def subplotSpacingOptions(self):
        if self.fig is None:
            return

        self.adj_window = MainWin_close()
        self.adj_window.window_closed.connect(self.copySubplotSettings)
        win = self.adj_window

        win.setWindowTitle("Subplot Configuration Tool")
        image = os.path.join( matplotlib.rcParams['datapath'],'images','matplotlib.png' )
        win.setWindowIcon(QtGui.QIcon( image ))

        tool = SubplotToolQt(self.fig, win)
        win.setCentralWidget(tool)
        win.setSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)

        win.show()

    def copySubplotSettings(self):
        '''Get a copy of the settings found in the 'Figure Options' editor.
        '''
        if self.axes is None:
            return

        # subplot settings
        for k in self._subplot_keepers:
            self._subplotSettings[k] = getattr(self.axes, 'get_'+k)()

        # line settings
        self._lineSettings = []
        for l in self.axes.get_lines():
            s = {}
            for k in self._line_keepers:
                s[k] = getattr(l, 'get_'+k)()
            self._lineSettings.append(s)

        # subplot position 
        self._subplotPosition = {}
        self._subplotPosition['left'] = self.fig.subplotpars.left
        self._subplotPosition['right'] = self.fig.subplotpars.right
        self._subplotPosition['top'] = self.fig.subplotpars.top
        self._subplotPosition['bottom'] = self.fig.subplotpars.bottom
        self._subplotPosition['wspace'] = self.fig.subplotpars.wspace
        self._subplotPosition['hspace'] = self.fig.subplotpars.hspace

    def applySubplotSettings(self):
        '''Everytime the plot is drawn it looses its 'Figure Options' so just
        make sure they are applied again.
        '''
        # subplot settings
        for k in self._subplot_keepers:
            if k in self._subplotSettings:
                getattr(self.axes, 'set_'+k)(self._subplotSettings[k])

        #subplot position
        self.fig.subplots_adjust(**self._subplotPosition)

    def on_draw(self):
        self._on_draw_cnt += 1
        if not self._updatetimer.isActive():
            self._updatetimer.start()

    def _on_draw(self):
        self.fig.clear()

        if self.get_autoscale():
            self.axes = self.fig.add_subplot(111, autoscale_on=self._autoscale_btn.get_val(), **self.get_plotlabels())
        else:
            self.axes = self.fig.add_subplot(111, autoscale_on=self._autoscale_btn.get_val(), xlim=self.get_xlim(), ylim=self.get_ylim(), **self.get_plotlabels())
        # self.axes.plot(self.x, self.y, 'ro')
        # self.axes.imshow(self.data, interpolation='nearest')
        # self.axes.plot([1,2,3])

        self.axes.grid(self.get_grid())
        self.axes.set_axis_bgcolor('0.97')

        if self._data is None:
            return

        self.fig.hold(True)

        # plot each set
        # print "--------------------plot the data"
        cnt = -1
        for data in self._data:
            cnt += 1

            # check for x, y data
            if cnt < len(self._lineSettings):
                s = self._lineSettings[cnt]
                if data.shape[-1] == 2:
                    self.axes.plot(data[..., 0], data[..., 1], **s)
                else:
                    self.axes.plot(data, **s)
            else:

                ln = max(data.shape)
                lw = max(5.0-np.log10(ln), 1.0)
                al = max(1.0-1.0/np.log2(ln), 0.75)

                if data.shape[-1] == 2:
                    self.axes.plot(data[..., 0], data[..., 1], alpha=al, lw=lw)
                else:
                    self.axes.plot(data, alpha=al, lw=lw)

        # LEGEND
        if self.get_legend():
            handles, labels = self.axes.get_legend_handles_labels()
            self.axes.legend(handles, labels)

        # AUTOSCALE
        if self.get_autoscale():
            self.set_xlim(self.axes.get_xlim(), quiet=True)
            self.set_ylim(self.axes.get_ylim(), quiet=True)

        # X TICKS
        xl = self._x_ticks.text().split(',')
        if len(xl) > 1:
            self.axes.set_xticks(np.linspace(*self.axes.get_xlim(), num=len(xl)))
            self.axes.set_xticklabels(xl)
        else:
            self.axes.set_xticks(np.linspace(*self.axes.get_xlim(), num=self._x_numticks.get_val()))

        # Y TICKS
        yl = self._y_ticks.text().split(',')
        if len(yl) > 1:
            self.axes.set_yticks(np.linspace(*self.axes.get_ylim(), num=len(yl)))
            self.axes.set_yticklabels(yl)
        else:
            self.axes.set_yticks(np.linspace(*self.axes.get_ylim(), num=self._y_numticks.get_val()))

        self.applySubplotSettings()

        self.canvas.draw()

        #print 'draw count: ', self._on_draw_cnt
        self._on_draw_cnt = 0

    def on_key_press(self, event):
        # print 'Matplotlib-> you pressed:' + str(event.key)
        # implement the default mpl key press events described at
        # http://matplotlib.org/users/navigation_toolbar.html#navigation-
        # keyboard-shortcuts
        try:
            from matplotlib.backend_bases import key_press_handler
            key_press_handler(event, self.canvas, self.mpl_toolbar)
        except:
            print("key_press_handler import failed. -old matplotlib version.")


class ExternalNode(gpi.NodeAPI):

    """A Qt embedded plot window using the code from: 
    http://matplotlib.org/examples/user_interfaces/embedding_in_qt4_wtoolbar.html
    keyboard shortcuts can be found here:
    http://matplotlib.org/users/navigation_toolbar.html#navigation-keyboard-shortcuts

    INPUTS
    Up to 8 data sets can be plotted simultaneously
      1D real-valued data are plotted as graph
      2D data where the 2nd dimension is 2 will be plotted as X-Y parametric plot, otherwise
      all other 2D data are plotted as series of 1D plots
    """

    def initUI(self):
        # Widgets
        self.addWidget('MatplotDisplay', 'Plot')

        # IO Ports
        self.inport_range = range(0, 8)
        for i in self.inport_range:
            self.addInPort('in' + str(i), 'NPYarray', obligation=gpi.OPTIONAL)

    def compute(self):

        # check input ports for data
        in_lst = [self.getData('in' + str(i))
                  for i in self.inport_range if self.getData('in' + str(i)) is not None]

        self.setAttr('Plot', val=in_lst)

        return 0

    def execType(self):
        '''Could be GPI_THREAD, GPI_PROCESS, GPI_APPLOOP'''
        return gpi.GPI_APPLOOP
