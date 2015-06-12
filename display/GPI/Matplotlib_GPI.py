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

import gpi
from gpi import QtCore, QtGui

import numpy as np
from matplotlib.figure import Figure
#from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_qt4agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)

class MatplotDisplay(gpi.GenericWidgetGroup):

    """Embeds the matplotlib figure window.
    """
    valueChanged = gpi.Signal()

    def __init__(self, title, parent=None):
        super(MatplotDisplay, self).__init__(title, parent)

        # gpi interface
        self._collapsables = []

        # plot specific UI side panel
        #  -sets options for plot window so this needs to be run first
        vbox = QtGui.QVBoxLayout()
        vbox.setContentsMargins(0, 0, 0, 0)  # no spaces around this item
        vbox.setSpacing(0)

        self._autoscale_btn = gpi.widgets.BasicPushButton(self)
        self._autoscale_btn.set_toggle(True)
        self._autoscale_btn.set_button_title('autoscale')
        self._autoscale_btn.valueChanged.connect(self.on_draw)
        self._collapsables.append(self._autoscale_btn)

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
        self._xlab = QtGui.QLabel('x limits')
        self._ylab = QtGui.QLabel('y limits')
        self._maxlab = QtGui.QLabel('max')
        self._minlab = QtGui.QLabel('min')
        lims.addWidget(self._maxlab,1,0,1,1)
        lims.addWidget(self._minlab,2,0,1,1)
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
        self._collapsables.append(self._minlab)
        self._collapsables.append(self._maxlab)

        vbox.addLayout(lims)
        vbox.addWidget(self._autoscale_btn)
        vbox.insertStretch(-1,1)

        # plot window
        self._data = None
        self._plotwindow = self.create_main_frame()
        self.on_draw()

        # put side panel and plot window together
        hbox = QtGui.QHBoxLayout()
        hbox.addLayout(vbox)
        hbox.addLayout(self._plotwindow)
        hbox.setStretch(0,11)
        self.setLayout(hbox)

        self.set_collapsed(True)  # hide options by default

    # setters
    def set_val(self, data):
        '''Takes a list of npy arrays.
        '''
        if isinstance(data, list):
            self._data = data
            self.on_draw()
        else:
            return

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

    # getters
    def get_val(self):
        return self._data

    def get_xlim(self):
        return (self._xh.get_val(), self._xl.get_val())

    def get_ylim(self):
        return (self._yh.get_val(), self._yl.get_val())

    def get_autoscale(self):
        return self._autoscale_btn.get_val()

    # support
    def create_main_frame(self):
        self.fig = Figure((6.0, 4.8), dpi=100, facecolor='white', linewidth=6.0, edgecolor='0.93')
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self)
        self.canvas.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.canvas.setFocus()

        self.mpl_toolbar = NavigationToolbar(self.canvas, self)

        self.canvas.mpl_connect('key_press_event', self.on_key_press)

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.canvas)  # the matplotlib canvas
        vbox.addWidget(self.mpl_toolbar)
        return vbox
        #self.setLayout(vbox)

    #def get_data2(self):
    #    return np.arange(20).reshape([4, 5]).copy()

    def on_draw(self):
        self.fig.clear()

        if self.get_autoscale():
            self.axes = self.fig.add_subplot(111, autoscale_on=self._autoscale_btn.get_val())
        else:
            self.axes = self.fig.add_subplot(111, autoscale_on=self._autoscale_btn.get_val(), xlim=self.get_xlim(), ylim=self.get_ylim())
        # self.axes.plot(self.x, self.y, 'ro')
        # self.axes.imshow(self.data, interpolation='nearest')
        # self.axes.plot([1,2,3])

        if self._data is None:
            return

        self.fig.hold(True)

        # plot each set
        # print "--------------------plot the data"
        for data in self._data:

            # check for x, y data
            if data.shape[-1] == 2:
                self.axes.plot(data[..., 0], data[..., 1], alpha=0.8, lw=2.0)
            else:
                self.axes.plot(data, alpha=0.8, lw=2.0)

        if self.get_autoscale():
            self.set_xlim(self.axes.get_xlim(), quiet=True)
            self.set_ylim(self.axes.get_ylim(), quiet=True)

        self.canvas.draw()

    def on_key_press(self, event):
        # print 'Matplotlib-> you pressed:' + str(event.key)
        # implement the default mpl key press events described at
        # http://matplotlib.org/users/navigation_toolbar.html#navigation-
        # keyboard-shortcuts
        try:
            from matplotlib.backend_bases import key_press_handler
            key_press_handler(event, self.canvas, self.mpl_toolbar)
        except:
            print "key_press_handler import failed. -old matplotlib version."


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
