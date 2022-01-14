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
# following applies:x
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

# Author: Abdulrahman Fayad
# Date: 2021 Dec
import threading 

import gpi 
from gpi import QtGui, QtWidgets, QtCore, QtOpenGL
import pyqtgraph as pg
import numpy as np

try:
    try:
        from OpenGL import gl, glu, glut
        import OpenGL.GL
        import OpenGL.GLU
        import OpenGL.GLUT
    
    # fallback for OpenGl Module import source
    except ImportError:
        # https://stackoverflow.com/questions/63475461/unable-to-import-opengl-gl-in-python-on-macos
        print('Patching OpenGL import for Big Sur (OSX)')
        from ctypes import util
        orig_util_find_library = util.find_library
        
        # retrieving the OpenGl Modules from its directory instead of cache
        def new_util_find_library (name):
            res = orig_util_find_library (name)
            if res: return res
            return '/System/Library/Frameworks/' + name + '.framework/' + name

        util.find_library = new_util_find_library
        
        from OpenGL import GLU, GLUT
        from OpenGL import GL as gl
        import OpenGL.GL
        import OpenGL.GLU
        import OpenGL.GLUT
        from OpenGL.arrays import vbo
        import gpi.GLObjects as glo

except ImportError:
    log.warn('OpenGL was not found, GL objects and windows will not be supported in this session.')
    raise


class GLWidget(QtOpenGL.QGLWidget):
    def __init__(self, parent=None):
        self.parent = parent
        QtOpenGL.QGLWidget.__init__(self, parent)

    def initializeGL(self):
        self.qglClearColor(QtGui.QColor(255, 255, 255))
        gl.glEnable(gl.GL_DEPTH_TEST)
        self.initGeometry()

    def resizeGL(self, width, height):
        gl.glViewport(0, 0, width, height)
        gl.glMatrixMode(gl.GL_PROJECTION)
        gl.glLoadIdentity()
        if height == 0: height = 1
        aspect = width / height

        GLU.gluPerspective(45.0, aspect, 1.0, 100.0)
        gl.glMatrixMode(gl.GL_MODELVIEW)

    def paintGL(self):
        gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)

        gl.glPushMatrix()    # push the current matrix to the current stack

        gl.glTranslate(0.0, 0.0, -50.0)    # third, translate cube to specified depth
        gl.glScale(20.0, 20.0, 20.0)       # second, scale cube
        gl.glTranslate(-0.5, -0.5, -0.5)   # first, translate cube center to origin

        gl.glEnableClientState(gl.GL_VERTEX_ARRAY)
        gl.glEnableClientState(gl.GL_COLOR_ARRAY)

        gl.glVertexPointer(3, gl.GL_FLOAT, 0, self.vertVBO)
        gl.glColorPointer(3, gl.GL_FLOAT, 0, self.colorVBO)

        gl.glDrawElements(gl.GL_QUADS, len(self.cubeIdxArray), gl.GL_UNSIGNED_INT, self.cubeIdxArray)

        gl.glDisableClientState(gl.GL_VERTEX_ARRAY)
        gl.glDisableClientState(gl.GL_COLOR_ARRAY)

        gl.glPopMatrix()    # restore the previous modelview matrix

    def initGeometry(self):
        self.cubeVtxArray = np.array(
            [[0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0]])
        self.vertVBO = vbo.VBO(np.reshape(self.cubeVtxArray,
                        (1, -1)).astype(np.float32))
        self.vertVBO.bind()

        self.cubeClrArray = np.array(
	    [[0.0, 0.0, 0.0],
	     [1.0, 0.0, 0.0],
	     [1.0, 1.0, 0.0],
	     [0.0, 1.0, 0.0],
	     [0.0, 0.0, 1.0],
	     [1.0, 0.0, 1.0],
	     [1.0, 1.0, 1.0],
	     [0.0, 1.0, 1.0 ]])
        self.colorVBO = vbo.VBO(np.reshape(self.cubeClrArray,
                        (1, -1)).astype(np.float32))
        self.colorVBO.bind()

        self.cubeIdxArray = np.array(
            [0, 1, 2, 3,
            3, 2, 6, 7,
            1, 0, 4, 5,
            2, 1, 5, 6,
            0, 3, 7, 4,
            7, 6, 5, 4 ])

    def mousePressEvent(self, event):
        self.init_pos = event.pos()

    def mouseMoveEvent(self, event):
        dx = event.x() - self.init_pos.x()
        dy = event.y() - self.init_pos.y()
        viewMatrix = gl.glGetFloatv(gl.GL_MODELVIEW_MATRIX)

        # init the view matrix
        gl.glLoadIdentity()

        gl.glRotatef(0.1, 0.0, 1.0, 0.0)
        gl.glMultMatrixf(viewMatrix)
        viewMatrix = gl.glGetFloatv(gl.GL_MODELVIEW_MATRIX)
        self.update()
    
class Viewer(gpi.GenericWidgetGroup):
    valueChanged = gpi.Signal(bool)

    def __init__(self, title, parent=None):
        super(Viewer, self).__init__(title, parent)

        self.glWidget = GLWidget(parent)
        # setup the image viewer widgets into a grid widget
        wdgLayout = QtWidgets.QGridLayout()
        wdgLayout.addWidget(self.glWidget, 0, 0)
        self.setLayout(wdgLayout)

    

       

class ExternalNode(gpi.NodeAPI):

    def execType(self):
        return gpi.GPI_APPLOOP

    def initUI(self):
        self.addWidget("Viewer", "Viewer")
        self.addInPort('GL Object Descriptions', 'GLOList', obligation=gpi.OPTIONAL)
        self.addOutPort('Rendered RGBA', 'NPYarray')

    def validate(self):
        return 0
        

    def compute(self):
       
        return 0
