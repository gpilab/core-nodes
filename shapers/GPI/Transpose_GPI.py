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


# Author: Ryan Robison
# Date: 2015aug25

import gpi
from gpi import QtGui


class OrderButtons(gpi.GenericWidgetGroup):
    """A set of reorderable tabs."""
    valueChanged = gpi.Signal()

    def __init__(self, title, parent=None):
        super(OrderButtons, self).__init__(title, parent)

        # at least one button
        wdgLayout = QtGui.QGridLayout()
        self.wdg = QtGui.QTabBar()
        self.wdg.addTab('0')
        self.wdg.setMovable(True)
        self.wdg.tabMoved.connect(self.labelsChanged)
        self._labels = ['0']
        # layout
        wdgLayout.addWidget(self.wdg)
        self.setLayout(wdgLayout)

    # setters
    def set_val(self, names):
        """list(str,str,...) | A list of labels (e.g. ['b1', 'b2',...])."""

        # clear all buttons
        self._labels = []
        while self.wdg.count() > 0:
            self.wdg.removeTab(0)

        # add buttons
        for i in xrange(len(names)):
            self.wdg.addTab(names[i])
            self._labels.append(names[i])

    def set_visible(self, val):
        self.wdg.setVisible(val)

    # getters
    def get_val(self):
        return self._labels

    # support
    def labelsChanged(self):
        self._labels = []
        for i in xrange(self.wdg.count()):
            self._labels.append(str(self.wdg.tabText(i)))
        self.valueChanged.emit()


class ExternalNode(gpi.NodeAPI):
    """Peform transpose operations on an array.
    INPUT - input array
    OUTPUT - output array

    WIDGETS:
    Info: - information on size, etc of input and output array
    Transpose - transpose any dimensions
    """

    def initUI(self):

        # Widgets
        self.shape = []
        self.trans_ind = []
        self.info_message = ""
        self.addWidget('TextBox', 'Info:')
        self.buttons = []
        # for i in range(3):
        #     self.buttons = self.buttons+[str(-i-1)]
        # self.addWidget('ExclusivePushButtons', 'Input Dim', buttons=self.buttons, val=0)
        # self.addWidget('ExclusivePushButtons', 'Output Dim', buttons=self.buttons, val=0)
        self.addWidget('PushButton', 'Transpose', toggle=True, val=0)
        self.addWidget('OrderButtons', 'Dimension Order')

        # IO Ports
        self.addInPort('in', 'NPYarray', obligation=gpi.REQUIRED)
        self.addOutPort('out', 'NPYarray')

    def validate(self):
        data = self.getData('in')
        # dim1 = self.getVal('Input Dim')
        # dim2 = self.getVal('Output Dim')

        order = self.getVal('Dimension Order')
        if len(order) != data.ndim:
            order = []
            for i in xrange(data.ndim):
                # order = order+[str(-i-1)]
                order = order+[str(i)]
        self.setAttr('Dimension Order', val=order)

        if data.ndim < 3:
            self.setAttr('Dimension Order', visible=False)
        else:
            self.setAttr('Dimension Order', visible=True)

        # if self.portEvents():
        #     if dim1 > data.ndim:
        #         dim1 = 0
        #     if dim2 > data.ndim:
        #         dim2 = 0

        # self.buttons = []
        # for i in range(data.ndim):
        #     self.buttons = self.buttons+[str(-i-1)]
        # self.setAttr('Input Dim', buttons = self.buttons, val=dim1)
        # self.setAttr('Output Dim', buttons = self.buttons, val=dim2)
        # self.trans_ind = range(-data.ndim, 0)

        # self.trans_ind[data.ndim-dim1-1] = int(self.buttons[dim2])
        # self.trans_ind[data.ndim-dim2-1] = int(self.buttons[dim1])

        return(0)

    def compute(self):

        data = self.getData('in')
        transpose = self.getVal('Transpose')
        basic_info = "Input Dimensions: "+str(data.shape)+"\n" 
        order = self.getVal('Dimension Order')
        trans_ind = data.ndim*[0]
        self.shape = list(data.shape)
        if data.ndim > 2:
            for i in xrange(len(order)):
                # trans_ind[len(order)-1-i] = int(order[i])
                trans_ind[i] = int(order[i])
                self.shape[i] = data.shape[trans_ind[i]]
        else:
            trans_ind = [1, 0]

        if transpose and data.ndim > 1:
            out = data.transpose(trans_ind)
            self.shape = out.shape
            self.setData('out', out)

        else:
            out = data

        # updata info 
        info = basic_info+"Output Dimensions: "+str(self.shape)+"\n"
        self.setAttr('Info:', val = info)

        return(0)

    def execType(self):
        '''Could be GPI_THREAD, GPI_PROCESS, GPI_APPLOOP'''
        return gpi.GPI_PROCESS
