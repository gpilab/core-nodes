import gpi

class ExternalNode(gpi.NodeAPI):
    """Specify a python boolean for use as node-data or widget-ports-parms.

    OUTPUT - boolean value

    WIDGETS:
      bool - specifies whether boolean is True or False
    """

    def initUI(self):
        # self.setTitle("Testing")
        pass

    def compute(self):
        # print("now")
        # self.setTitle("Testing")
        # self.setWindowTitle("Testing")
        # self.node._nodeIF.updateTitle()
        # print(self.node.getModuleName())
        # self.name = "Testing"

        return 0

    def execType(self):
        return gpi.GPI_THREAD
