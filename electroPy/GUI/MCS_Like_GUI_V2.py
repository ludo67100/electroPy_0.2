#-*- coding: utf-8 -*-
"""
Created on Tue Dec 11 14:50:59 2018

@author: ludovic.spaeth
"""


from PyQt4 import QtGui

import sys
import numpy as np
import traceback

import pyqtgraph as pg

from electroPy.HdF5IO import HdF5IO

#---------------------------INSTANCE VARIABLES CLASS---------------------------
#This class will contain all the general variables ie the data to pass it from one
#Window to another 
#Will be usefull for filtered data etc 

class DATA():
    def __init__(self):
        print 'GetDATA'

    def get_raw_data(self,raw_data):
        
        self.raw_data = raw_data #The raw array traces, with leak, ie for CClamp exp

        
    def get_corr_data(self,corr_data):
        
        self.corr_data = corr_data #The non leaked traces, ie for VClamp exp 


#-----------------------MAIN WINDOW--------------------------------------------
class MainWindow(QtGui.QWidget):
    
    def __init__(self,parent=None):
        super(MainWindow, self).__init__()
                
        self.tag_list = np.zeros(100) #By default none traces are tagged 
        
        self.data = DATA() #To put data in 
        
        self.initUI()
        
    #------------------------GUI BASICS------------------------------------        
    def initUI(self): 
        self.setGeometry(1200,800, 1800, 800)
        self.center()
        self.setWindowTitle('This f*cking GUI')     
        
        #Grid Layout-------------
        grid = QtGui.QGridLayout()
        self.setLayout(grid)
                    
        #Canvas and Toolbar--------------------
        self.canvas = pg.GraphicsLayoutWidget()
        grid.addWidget(self.canvas, 3,0,1,6)

        
        #Import File Button----------------------------------
        btn1 = QtGui.QPushButton('Import ePhy File (HdF5)...', self)
        btn1.resize(btn1.sizeHint()) 
        btn1.clicked.connect(self.getFile)
        grid.addWidget(btn1, 0,0)
        
        #Plot Button------------------------------
        btn2 = QtGui.QPushButton('Plot (P)', self)
        btn2.resize(btn2.sizeHint())    
        btn2.clicked.connect(self.plot)
        btn2.setShortcut("p")
        grid.addWidget(btn2, 0,1)
    
        #Quit Button
        self.quit_btn = QtGui.QPushButton('Quit App')
        self.quit_btn.resize(self.quit_btn.sizeHint())
        self.quit_btn.clicked.connect(self.QuitApp)
        grid.addWidget(self.quit_btn,4,5)
            
        self.show()

    #------------------------METHODS-----------------------------------------------    
    def QuitApp(self):

        self.close()
        
        
    def getFile(self): #Penser a ajouter le channel select ici 0 par defaut
        filePath = QtGui.QFileDialog.getOpenFileName(self,'Choose ur file')       
        try : 
            path = str(filePath)
            self.file = HdF5IO(filepath=path)
            '''
            self.plot()
            self.plot() #twice, to plot directly after load w/out hitting plot btn 
            '''
        except Exception:
            traceback.print_exc()
                        
        self.setWindowTitle(str(filePath))  
        

    def plot(self):   
        try:
            
            self.canvas.clear()
            
            grid = np.array([[1,9],
                     [2,10],
                     [3,11],
                     [4,12],
                     [5,13],
                     [6,14],
                     [7,15],
                     [8,16]])
    
            def LinkAxis(a):
                '''
                Links all views from a list of view to the first 1
                Both X and Y axis
                '''    
                for i in range(len(a)):    
                    
                    if i ==0:
                        continue
                    
                    a[i].setXLink(a[0])
                    a[i].setYLink(a[0])
    
            pen = pg.mkPen(color=(255,140,0)) #Instance for color, rgb. Here should be added linethick and dotted lines
            
            Views = [] #To store the views (p) and to link them afterwards
            
            for (row,col),ch in np.ndenumerate(grid-1):
                
                time = self.file.raw_time()
                data = self.file.raw_record()[ch]
            
                p = self.canvas.addPlot(row=row,col=col)
                p.setLabel('left','Ch_%s (V)'%(ch+1))

            
                p.plot(time,data,antialias=False,downsample=50,downsampleMethod='subsample',pen=pen) #or autoDownsample=True
                
                Views.append(p)
                
            Views[0].setXLink(Views[1])
            Views[1].setXLink(Views[2])
                    
            LinkAxis(Views)
            
        except Exception : 
            traceback.print_exc()
           
        
    
    def center(self):
        qr = self.frameGeometry()
        cp = QtGui.QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())


            
    def Spike_detect(self):
        try :
            self.dialog_spike_detect.show()
            
        except Exception:
            traceback.print_exc()


if __name__ == '__main__':
    app = QtGui.QApplication.instance()
    if app is None:
        app = QtGui.QApplication(sys.argv) #For first launch
    else:
        print('QApplication instance already exists: %s' % str(app)) # In case of second launch to avoid kernel crash

    GUI = MainWindow()
    DATA = DATA()
    GUI.show() #Call GUI.object to get any self.variable from MainWindowClass
    sys.exit(app.exec_())
