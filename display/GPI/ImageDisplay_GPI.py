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

# Author: Jim Pipe / Nick Zwart
# Date: 2013 Sep 01

import math
import gpi
from gpi import QtCore, QtGui
import numpy as np

# WIDGET
class WindowLevel(gpi.GenericWidgetGroup):
    """Provides an interface to the BasicCWFCSliders."""
    valueChanged = gpi.Signal()

    def __init__(self, title, parent=None):
        super(WindowLevel, self).__init__(title, parent)
        self.sl = gpi.BasicCWFCSliders()
        self.sl.valueChanged.connect(self.valueChanged)
        self.pb = gpi.BasicPushButton()
        self.pb.set_button_title('reset')
        # layout
        wdgLayout = QtGui.QVBoxLayout()
        wdgLayout.addWidget(self.sl)
        wdgLayout.addWidget(self.pb)
        self.setLayout(wdgLayout)
        # default
        self.set_min(0)
        self.set_max(100)
        self.sl.set_allvisible(True)
        self.reset_sliders()
        self.pb.valueChanged.connect(self.reset_sliders)

    # setters
    def set_val(self, val):
        """Set multiple values with a python-dict with keys:
        level, window, floor and ceiling. -Requires integer
        values for each key.
        """
        self.sl.set_center(val['level'])
        self.sl.set_width(val['window'])
        self.sl.set_floor(val['floor'])
        self.sl.set_ceiling(val['ceiling'])

    def set_min(self, val):
        """Set min for level, window, floor and ceiling (int)."""
        self.sl.set_min(val)

    def set_max(self, val):
        """Set max for level, window, floor and ceiling (int)."""
        self.sl.set_max(val)

    # getters
    def get_val(self):
        val = {}
        val['level'] = self.sl.get_center()
        val['window'] = self.sl.get_width()
        val['floor'] = self.sl.get_floor()
        val['ceiling'] = self.sl.get_ceiling()
        return val

    def get_min(self):
        return self.sl.get_min()

    def get_max(self):
        return self.sl.get_max()

    def reset_sliders(self):
        val = {}
        val['window'] = 100
        val['level'] = 50
        val['floor'] = 0
        val['ceiling'] = 100
        self.set_val(val)


class ExternalNode(gpi.NodeAPI):
    """2D image viewer for real or complex NPYarrays.

    INPUT:
    2D data, real or complex
    3D uint8 ARGB data (e.g. output of another ImageDisplay node)

    OUTPUT:
    3D data of displayed image, last dimension has length 4 for ARGB byte (uint8) data

    WIDGETS:
    Complex Display - If data are complex, allows you to show Real, Imaginary, Magnitude, Phase, or "Complex" data.
        If C (Complex) is chosen, then pixel brightness reflects value magnitude, while pixel color reflects value phase.
        If input data are real-valued, this widget is hidden

    Color Map - Chooses from a number of colormaps for real-valued data.  Not available if Scalar Display is "Sign"

    Edge Pixels - only visible for complex input data with Complex Display set to "C"
        Setting this to N creates an N-pixel color ring around the image border illustrating the phase-to-color mapping

    Black Pixels - only visible for complex input data with Complex Display set to "C"
        Setting this to N creates an N-pixel black ring around the image border (but inside the Edge Pixels ring) to
        separate the Edge pixel ring from the actual data image

    Viewport - displays the image
      Double clicking on Viewport area brings up a scaling widget to change the image size, and change graphic overlay

    L W F C - (hidden by default - double click on widget area to show sliders)
              Adjust value-to-pixel brightness mapping using Level/Window or Floor/Ceiling

    Scalar Display - visible for real data, or complex data with "Complex Display" set to R, I, M, or P
      Pass uses the real data in the data-to-pixel mapping
      Mag uses the magnitude data in the data-to-pixel mapping (i.e. affects negative values only)
      Sign will display positive values in green, and the absolutele value of negative values in magenta

    Gamma - changes gamma of display function.  Default value of 1 gives linear mapping of data to pixel value
      pixel values refect value of data^gamma

    Zero Ref - visible for real data, or complex data with "Complex Display" set to R, I, M, or P
               also invisible if "Scalar Display" set to Sign
      This is used for the data-to-pixel value mapping
      --- maps the smallest value to black, and the largest value to white
      0-> maps zero to black, and the largest value to white.  All negative numbers are black
      -0- maps zero to middle gray, with the largest magnitude set to black (if negative) or white (if positive)
      <-0 maps zero to white, and the most negative value to black.  All positive numbers are white

    Fix Range
      If Auto-Range On, the data range for pixel value mapping is rescaled whenever new data appears at input
      If Fixed-Range On, the data range is fixed (and can be changed using Range Min and Range Max)

    Range Min - shows minimum data value used for mapping to pixel values
      This value can be changed if (and only if) "Fix Range" is set to "Fixed-Ranged On"

    Range Max - shows maximum data value used for mapping to pixel values
      This value can be changed if (and only if) "Fix Range" is set to "Fixed-Ranged On"
    """

    def execType(self):
        return gpi.GPI_APPLOOP

    def initUI(self):

        # Widgets
        self.addWidget('ExclusivePushButtons','Complex Display',
                       buttons=['R','I','M','P','C'], val=4)
        self.addWidget('ExclusivePushButtons','Color Map',
                       buttons=['Gray','IceFire','Fire','Hot','HOT2','BGR'], val=0, collapsed=True)
        self.addWidget('SpinBox', 'Edge Pixels', min=0)
        self.addWidget('SpinBox', 'Black Pixels', min=0)
        self.addWidget('DisplayBox', 'Viewport:')
        self.addWidget('WindowLevel', 'L W F C:', collapsed=True)
        self.addWidget('ExclusivePushButtons','Scalar Display',
                       buttons=['Pass','Mag','Sign'], val=0)
        self.addWidget('DoubleSpinBox', 'Gamma',min=0.1,max=10,val=1,singlestep=0.05,decimals=3)
        self.addWidget('ExclusivePushButtons','Zero Ref',
                       buttons=['---','0->','-0-','<-0'], val=0)
        self.addWidget('PushButton', 'Fix Range', button_title='Auto-Range On', toggle=True)
        self.addWidget('DoubleSpinBox', 'Range Min')
        self.addWidget('DoubleSpinBox', 'Range Max')

        # IO Ports
        self.addInPort('in', 'NPYarray')
        self.addOutPort('out', 'NPYarray')

    def validate(self):

        # Complex or Scalar?
        data = self.getData('in')

        if data.ndim is not 2:
          if ((data.ndim is not 3) or
              (data.shape[-1] is not 3) and (data.shape[-1] is not 4)):
            self.log.warn("data must be 2D or 3D with the last dim=3 or 4")
            return 1

        self.setAttr('L W F C:',visible=(data.ndim==2))
        self.setAttr('Gamma',visible=(data.ndim==2))
        self.setAttr('Fix Range',visible=(data.ndim==2))
      
        if data.ndim is 3:
          self.setAttr('Complex Display',visible=False)
          self.setAttr('Color Map',visible=False)
          self.setAttr('Scalar Display',visible=False)
          self.setAttr('Edge Pixels',visible=False)
          self.setAttr('Black Pixels',visible=False)
          self.setAttr('Zero Ref',visible=False)
          self.setAttr('Range Min',visible=False)
          self.setAttr('Range Max',visible=False)

        else:

          if np.iscomplexobj(data):
            self.setAttr('Complex Display',visible=True)
            scalarvis = self.getVal('Complex Display') != 4
          else:
            self.setAttr('Complex Display',visible=False)
            scalarvis = True
          self.setAttr('Color Map',visible=scalarvis)
          self.setAttr('Scalar Display',visible=scalarvis)
          self.setAttr('Edge Pixels',visible=not scalarvis)
          self.setAttr('Black Pixels',visible=not scalarvis)

          if self.getVal('Scalar Display') == 2:
            self.setAttr('Zero Ref',visible=False)
          else:
            self.setAttr('Zero Ref',visible=scalarvis)

          self.setAttr('Range Min',visible=scalarvis)

          zval = self.getVal('Zero Ref')
          if zval == 1:
            self.setAttr('Range Min',val=0)
          elif zval == 3:
            self.setAttr('Range Max',val=0)

          if self.getVal('Fix Range'):
            self.setAttr('Fix Range',button_title="Fixed Range On")
          else:
            self.setAttr('Fix Range',button_title="Auto-Range On")

        return 0

    def compute(self):

        import math
        import numpy as np

        # make a copy for changes
        data = self.getData('in').copy()

        # Read in parameters, make a little floor:ceiling adjustment
        gamma = self.getVal('Gamma')
        lval = self.getAttr('L W F C:', 'val')
        cval = self.getVal('Complex Display')
        cmap = self.getVal('Color Map')
        sval = self.getVal('Scalar Display')
        zval = self.getVal('Zero Ref')
        fval = self.getVal('Fix Range')
        rmin = self.getVal('Range Min')
        rmax = self.getVal('Range Max')

        flor = 0.01*lval['floor']
        ceil = 0.01*lval['ceiling']
        if ceil == flor:
          if ceil == 1.:
            flor = 0.999
          else:
            ceil += 0.001

        # SHOW COMPLEX DATA
        if np.iscomplexobj(data) and cval == 4:
          mag = np.abs(data)
          phase = np.angle(data, deg=True)

          # normalize the mag
          data_min = 0.
          if fval:
            data_max = rmax
          else:
            data_max = mag.max()
            self.setAttr('Range Max',val=data_max)
          data_range = data_max-data_min
          dmask = np.ones(data.shape)
          new_min = data_range*flor + data_min
          new_max = data_range*ceil + data_min
          mag = np.minimum(np.maximum(mag,new_min*dmask),new_max*dmask)

          if new_max > new_min:
            if (gamma == 1): # Put in check for gamma=1, the common use case, just to save time
              mag = 255.*(mag - new_min)/(new_max-new_min)
            else:
              mag = 255.*pow((mag - new_min)/(new_max-new_min),gamma)
          else:
            mag = 255.*np.ones(mag.shape)

          # ADD BORDERS
          edgpix = self.getVal('Edge Pixels')
          blkpix = self.getVal('Black Pixels')
          if (edgpix + blkpix) > 0:
            # new image will be h2 x w2
            # frame defines edge pixels to paint with phase table
            h, w = mag.shape
            h2 = h + 2*(edgpix+blkpix)
            w2 = w + 2*(edgpix+blkpix)
            mag2 = np.zeros((h2,w2))
            phase2 = np.zeros((h2,w2))
            frame = np.zeros((h2,w2)) == 1
            frame[0:edgpix,:] = frame[h2-edgpix:h2,:] = True
            frame[:,0:edgpix] = frame[:,w2-edgpix:w2] = True

            mag2[edgpix+blkpix:edgpix+blkpix+h,edgpix+blkpix:edgpix+blkpix+w] = mag
            mag2[frame] = 255.

            phase2[edgpix+blkpix:edgpix+blkpix+h,edgpix+blkpix:edgpix+blkpix+w] = phase
            xloc = np.tile(np.linspace(-1.,1.,w2),(h2,1))
            yloc = np.transpose(np.tile(np.linspace(1.,-1.,h2),(w2,1)))
            phase2[frame] = np.degrees(np.arctan2(yloc[frame],xloc[frame]))

            mag = mag2
            phase = phase2

          # NOW COLORIZE
          hue = (phase+180.)/60. # hue goes from 0 to 6

          rd = np.zeros(mag.shape)
          gn = np.zeros(mag.shape)
          be = np.zeros(mag.shape)
          zmask = np.zeros(mag.shape)

          hindex0 =                         hue < 1
          hindex1 = np.logical_and(hue >= 1,hue < 2)
          hindex2 = np.logical_and(hue >= 2,hue < 3)
          hindex3 = np.logical_and(hue >= 3,hue < 4)
          hindex4 = np.logical_and(hue >= 4,hue < 5)
          hindex5 =                hue >= 5

          rd[hindex0]   = mag[hindex0]
          gn[hindex0] = (mag*hue)[hindex0]
          be[hindex0]  = zmask[hindex0]

          rd[hindex1]   = (mag*(2.-hue))[hindex1]
          gn[hindex1] = mag[hindex1]
          be[hindex1]  = zmask[hindex1]

          rd[hindex2]   = zmask[hindex2]
          gn[hindex2] = mag[hindex2]
          be[hindex2]  = (mag*(hue-2.))[hindex2]

          rd[hindex3]   = zmask[hindex3]
          gn[hindex3] = (mag*(4.-hue))[hindex3]
          be[hindex3]  = mag[hindex3]

          rd[hindex4]   = (mag*(hue-4.))[hindex4]
          gn[hindex4] = zmask[hindex4]
          be[hindex4]  = mag[hindex4]

          rd[hindex5]   = mag[hindex5]
          gn[hindex5] = zmask[hindex5]
          be[hindex5]  = (mag*(6.-hue))[hindex5]

          # Don't know what is going on!!!
          # the colors get rotated, so I am pre-swapping them to get the right color... :-/
          blue = np.uint8(rd)
          red = np.uint8(gn)
          green = np.uint8(be)
          alpha = np.uint8(255.*np.ones(blue.shape))

        # DISPLAY SCALAR DATA
        elif data.ndim is 2:

          if np.iscomplexobj(data):
            if cval == 0: # Real
              data = np.real(data)
            elif cval == 1: # Imag
              data = np.imag(data)
            elif cval == 2: # Mag
              data = np.abs(data)
            elif cval == 3: # Phase
              data = np.angle(data, deg=True)

          if sval == 1: # Mag
            data = np.abs(data)
          elif sval == 2: # Sign
            sign = np.sign(data)
            data = np.abs(data)

          # normalize the data
          if fval:
            data_min = rmin
            data_max = rmax
          else:
            data_min = data.min()
            data_max = data.max()

          if sval != 2:
            if zval == 1:
              data_min = 0.
            elif zval == 2:
              data_max = max(abs(data_min),abs(data_max))
              data_min = -data_max
            elif zval == 3:
              data_max = 0.
            data_range = data_max-data_min
            self.setAttr('Range Min',val=data_min)
            self.setAttr('Range Max',val=data_max)
          else:
            data_min = 0.
            data_max = max(abs(data_min),abs(data_max))
            data_range = data_max
            self.setAttr('Range Min',val=-data_range)
            self.setAttr('Range Max',val=data_range)

          dmask = np.ones(data.shape)
          new_min = data_range*flor + data_min
          new_max = data_range*ceil + data_min
          data = np.minimum(np.maximum(data,new_min*dmask),new_max*dmask)

          if new_max > new_min:
            if (gamma == 1): # Put in check for gamma=1, the common use case, just to save time
              data = 255.*(data - new_min)/(new_max-new_min)
            else:
              data = 255.*pow((data - new_min)/(new_max-new_min),gamma)
          else:
            data = 255.*np.ones(data.shape)

          if sval != 2: #Not Signed Data (Pass or Mag)
            # Show based on a color map
            if cmap == 0: # Grayscale
              red = green = blue = np.uint8(data)
              alpha = 255. * np.ones(blue.shape)
            else:
              rd = np.zeros(data.shape)
              gn = np.zeros(data.shape)
              be = np.zeros(data.shape)
              zmask = np.zeros(data.shape)
              fmask = np.ones(data.shape)

              if cmap == 1: # IceFire
                hue = 4.*(data/256.)
                hindex0 =                          hue < 1. 
                hindex1 = np.logical_and(hue >= 1.,hue < 2.)
                hindex2 = np.logical_and(hue >= 2.,hue < 3.)
                hindex3 = np.logical_and(hue >= 3.,hue < 4.)

                be[hindex0] = hue[hindex0]
                gn[hindex0] = zmask[hindex0]
                rd[hindex0] = zmask[hindex0]

                gn[hindex1] = (hue-1.)[hindex1]
                rd[hindex1] = (hue-1.)[hindex1]
                be[hindex1] = fmask[hindex1]

                gn[hindex2] = fmask[hindex2]
                rd[hindex2] = fmask[hindex2]
                be[hindex2] = (3.-hue)[hindex2]

                rd[hindex3] = fmask[hindex3]
                gn[hindex3] = (4.-hue)[hindex3]
                be[hindex3] = zmask[hindex3]

              elif cmap == 2: # Fire
                hue = 4.*(data/256.)
                hindex0 =                          hue < 1. 
                hindex1 = np.logical_and(hue >= 1.,hue < 2.)
                hindex2 = np.logical_and(hue >= 2.,hue < 3.)
                hindex3 = np.logical_and(hue >= 3.,hue < 4.)

                be[hindex0] = hue[hindex0]
                rd[hindex0] = zmask[hindex0]
                gn[hindex0] = zmask[hindex0]

                be[hindex1] = (2.-hue)[hindex1]
                rd[hindex1] = (hue-1.)[hindex1]
                gn[hindex1] = zmask[hindex1]

                rd[hindex2] = fmask[hindex2]
                gn[hindex2] = (hue-2.)[hindex2]
                be[hindex2] = zmask[hindex2]

                rd[hindex3] = fmask[hindex3]
                gn[hindex3] = fmask[hindex3]
                be[hindex3] = (hue-3.)[hindex3]

              elif cmap == 3: # Hot
                hue = 3.*(data/256.)
                hindex0 =                          hue < 1. 
                hindex1 = np.logical_and(hue >= 1.,hue < 2.)
                hindex2 = np.logical_and(hue >= 2.,hue < 3.)

                rd[hindex0] = hue[hindex0]
                be[hindex0] = zmask[hindex0]
                gn[hindex0] = zmask[hindex0]

                gn[hindex1] = (hue-1.)[hindex1]
                rd[hindex1] = fmask[hindex1]
                be[hindex1] = zmask[hindex1]

                rd[hindex2] = fmask[hindex2]
                gn[hindex2] = fmask[hindex2]
                be[hindex2] = (hue-2.)[hindex2]

              if cmap == 4: # Hot2, from ASIST (http://asist.umin.jp/index-e.htm)
                rindex0 = data < 20.0
                rindex1 = np.logical_and(data >=  20.0, data <= 100.0)
                rindex2 = np.logical_and(data > 100.0, data < 128.0)
                rindex3 = np.logical_and(data >= 128.0, data <= 191.0)
                rindex4 = data > 191.0
                rd[rindex0] = data[rindex0] * 4.0
                rd[rindex1] = 80.0 - (data[rindex1] - 20.0)
                rd[rindex3] = (data[rindex3] - 128.0) * 4.0
                rd[rindex4] = data[rindex4] * 0.0 + 255.0
                rd = rd/255.0;

                gindex0 = data < 45.0
                gindex1 = np.logical_and(data >= 45.0, data <= 130.0)
                gindex2 = np.logical_and(data > 130.0, data < 192.0)
                gindex3 = data >= 192.0
                gn[gindex1] = (data[gindex1] - 45.0)*3.0
                gn[gindex2] = data[gindex2] * 0.0 + 255.0
                gn[gindex3] = 252.0 - (data[gindex3] - 192.0)*4.0
                gn = gn/255.0

                bindex0 = (data < 1.0)
                bindex1 = np.logical_and(data >= 1.0, data < 86.0)
                bindex2 = np.logical_and(data >= 86.0, data <= 137.0)
                bindex3 = data > 137.0
                be[bindex1] = (data[bindex1] - 1.0)*3.0
                be[bindex2] = 255.0 - (data[bindex2] - 86.0)*5.0
                be = be/255.0

              elif cmap == 5: # BGR
                hue = 4.*(data/256.)
                hindex0 =                          hue < 1. 
                hindex1 = np.logical_and(hue >= 1.,hue < 2.)
                hindex2 = np.logical_and(hue >= 2.,hue < 3.)
                hindex3 = np.logical_and(hue >= 3.,hue < 4.)

                be[hindex0] = hue[hindex0]
                gn[hindex0] = zmask[hindex0]
                rd[hindex0] = zmask[hindex0]

                gn[hindex1] = (hue-1.)[hindex1]
                rd[hindex1] = zmask[hindex1]
                be[hindex1] = fmask[hindex1]

                gn[hindex2] = fmask[hindex2]
                rd[hindex2] = (hue-2.)[hindex2]
                be[hindex2] = (3.-hue)[hindex2]

                rd[hindex3] = fmask[hindex3]
                gn[hindex3] = (4.-hue)[hindex3]
                be[hindex3] = zmask[hindex3]

              blue = np.uint8(255.*rd)
              red = np.uint8(255.*be)
              green = np.uint8(255.*gn)
              alpha = np.uint8(255.*np.ones(blue.shape))

          else: #Signed data, positive numbers green, negative numbers magenta
            red = np.zeros(data.shape)
            green = np.zeros(data.shape)
            blue = np.zeros(data.shape)
            red[sign<=0] = data[sign<=0] 
            blue[sign<=0] = data[sign<=0] 
            green[sign>=0] = data[sign>=0] 

            red = red.astype(np.uint8)
            green = green.astype(np.uint8)
            blue = blue.astype(np.uint8)
            alpha = np.uint8(data)

        # DISPLAY RGB image
        else:
          red   = data[:,:,2].astype(np.uint8)
          green = data[:,:,1].astype(np.uint8)
          blue  = data[:,:,0].astype(np.uint8)
          if(data.ndim is 3 and data.shape[-1] is 4) :
              alpha = data[:,:,3].astype(np.uint8)
          else :
              alpha = 255.*np.ones(blue.shape)

        h, w = red.shape[:2]
        image1 = np.zeros((h, w, 4), dtype=np.uint8)
        image1[:, :, 0] = red
        image1[:, :, 1] = green
        image1[:, :, 2] = blue
        image1[:, :, 3] = alpha
        format_ = QtGui.QImage.Format_RGB32

        image = QtGui.QImage(image1.data, w, h, format_)

        #send the RGB values to the output port
        imageTru = np.zeros((h, w, 4), dtype=np.uint8)
        imageTru[:, :, 2] = red
        imageTru[:, :, 1] = green
        imageTru[:, :, 0] = blue
        imageTru[:, :, 3] = alpha
        image.ndarray = imageTru
        if image.isNull():
            self.log.warn("Image Viewer: cannot load image")

        self.setAttr('Viewport:', val=image)
        self.setData('out',imageTru)

        return 0
