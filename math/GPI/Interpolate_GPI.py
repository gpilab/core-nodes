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


#Author: Jim Pipe
#Date: 2014Mar

import gpi
import sys
import numpy as np
from scipy import interpolate
from math import fabs, sqrt, exp 
from numpy import linspace
from gpi import QtCore, QtGui


class Interpolate_GROUP(gpi.GenericWidgetGroup):
    """A combination of SpinBoxes, DoubleSpinBoxes, and PushButtons
    to form a unique widget suitable for interpolation options on dimensions.
    """
    valueChanged = gpi.Signal()

    def __init__(self, title, parent=None):
        super(Interpolate_GROUP, self).__init__(title, parent)

        # Widgets
        self._val = {}
        self._val['length'] = 1
        self._val['in_len'] = 1  # the original array length

        self.db = gpi.BasicDoubleSpinBox()  # factor
        self.db.set_label('factor:')
        self.db.set_min(0.001)
        self.db.set_max(gpi.GPI_FLOAT_MAX)
        self.db.set_decimals(5)
        self.db.set_singlestep(0.1)
        self.db.set_val(1)

        self.sb = gpi.BasicSpinBox()  # length
        self.sb.set_label('length:')
        self.sb.set_min(1)
        self.sb.set_val(1)
        self.sb.set_max(gpi.GPI_INT_MAX)

        self.db.valueChanged.connect(self.factChange)
        self.sb.valueChanged.connect(self.lenChange)

        vbox = QtGui.QHBoxLayout()
        vbox.addWidget(self.db)
        vbox.addWidget(self.sb)
        vbox.setStretch(0, 0)
        vbox.setStretch(1, 0)
        vbox.setStretch(2, 0)
        vbox.setContentsMargins(0, 0, 0, 0)  # we don't need margins here
        vbox.setSpacing(0)
        self.setLayout(vbox)
        
    def set_val(self, val):
        """A python-dict containing: in_len, length, and compute parms. """
        sig = False
        if 'in_len' in val:
            # otherwise this would change every time compute() was called
            if self._val['in_len'] != val['in_len']:
                self._val['in_len'] = val['in_len']
                fact = self.db.get_val()  # set len based on factor
                fact *= self._val['in_len']
                self.setLenQuietly(int(fact))
                self._val['length'] = int(fact)
        if 'length' in val:
            self._val['length'] = val['length']
            self.setLenQuietly(val['length'])
            self.setFactQuietly(float(val[
                                'length'])/float(self._val['in_len']))
            sig = True
        if 'compute' in val:
            self._val['compute'] = val['compute']
            sig = True
        if sig:
            self.valueChanged.emit()

    # getters
    def get_val(self):
        return self._val

    # support
    def factChange(self, val):
        self._val['length'] = int(self._val['in_len']*val)
        self.setLenQuietly(self._val['length'])
        self.valueChanged.emit()

    def lenChange(self, val):
        self._val['length'] = val
        self.setFactQuietly(float(val)/float(self._val['in_len']))
        self.valueChanged.emit()

    def compChange(self, val):
        self._val['compute'] = val
        self.valueChanged.emit()

    def setFactQuietly(self, val):
        self.db.blockSignals(True)
        self.db.set_val(val)
        self.db.blockSignals(False)

    def setLenQuietly(self, val):
        self.sb.blockSignals(True)
        self.sb.set_val(val)
        self.sb.blockSignals(False)

class ExternalNode(gpi.NodeAPI):
    """A node to linearly interpolate in specified directions.
    
    INPUT: input numpy array
    
    OUTPUT: interpolated data in the specified directions.
    
    WIDGETS:
    Dimension[n]:  For each dimension,
      factor - (desired output length)/(input length)
      length - desired output length
    Compute - compute
    """


    def initUI(self):

        # Widgets
        self.dim_base_name = 'Dimension['
        self.ndim = 6  # underlying c-code is only 6-dim
        for i in range(self.ndim):
            self.addWidget('Interpolate_GROUP', self.dim_base_name+str(i)+']')
        self.addWidget('PushButton', 'compute', toggle=True)

        # IO Ports
        self.addInPort('in', 'NPYarray', obligation=gpi.REQUIRED)
        self.addInPort('size', 'NPYarray', obligation=gpi.OPTIONAL)
        self.addOutPort('out1', 'NPYarray')

    def validate(self):
        if 'in' in self.portEvents() or 'size' in self.portEvents():
            data = self.getData('in')
            out_shape = self.getData('size')
           
            # visibility and bounds
            for i in xrange(self.ndim):
                if i < len(data.shape):
                    val = {'in_len': data.shape[i]}

                    # assign output lengths based on optional input array
                    if out_shape is not None:
                        offset = len(data.shape)-len(out_shape.shape)
                        if offset < 0:
                            if i < len(out_shape.shape):
                                val['length'] = out_shape.shape[i-offset]
                        else:
                            if i >= offset:
                                val['length'] = out_shape.shape[i-offset]

                    self.setAttr(self.dim_base_name+str(i)+']', visible=True, quietval=val)
                else:
                    self.setAttr(self.dim_base_name+str(i)+']',visible=False)
           
        return(0)
                       
    def compute(self):
      import sys
      import numpy as np
      import scipy as sp
      data_in = self.getData('in')
      data = data_in
      dimm = list(data_in.shape)
      self.ndim = data_in.ndim
      if self.getVal('compute'):
        for i in range(self.ndim):
          index = i
          val = self.getVal(self.dim_base_name+str(index)+']')
          interpnew = val['length']
          axisnew = dimm[i]
          dim_val = data_in.shape[i]
          new_axisnew = axisnew-1
          x = np.linspace(0, new_axisnew, axisnew)
          xnew = np.linspace(0, new_axisnew, interpnew)
          if dim_val == 1:
            m = 0
            for m in range(interpnew):
              if m>0:
                ynew = np.append(data,data,axis=i)
                data = ynew
              else:
                ynew = data
          else:
            yinterp = interpolate.interp1d(x,data,kind='linear',axis=i)
            ynew = yinterp(xnew)

          outdata = ynew
          data = ynew
          out1 = outdata

        out_dim = list(out1.shape)
        self.setData('out1', out1.astype(data_in.dtype))
      else:
        pass

      return(0)


    def execType(self):
        return gpi.GPI_PROCESS
