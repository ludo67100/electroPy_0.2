# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 14:52:09 2018

@author: ludovic.spaeth
"""

from PyQt4 import QtGui,QtCore

import sys
import numpy as np
import traceback

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar

from electroPy import WinWcpIo


class Average_Popup(QtGui.QWidget):
    def __init__(self, parent = None, widget=None):    
        self.setGeometry(500,500, 800, 800)
        QtGui.QWidget.__init__(self, parent)
        layout = QtGui.QGridLayout(self)
        button = QtGui.QPushButton("Very Interesting Text Popup")
        layout.addWidget(button)
        self.move(widget.rect().bottomLeft())
        self.setWindowTitle('Signal Averaging')  