# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 15:45:57 2018

GUI For electroPy

@author: ludovic.spaeth
"""

from PyQt4 import QtGui,QtCore

import sys
import numpy as np
import traceback


import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar

from electroPy import WinWcpIo, SpikeDetection

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


#---------------------------TAG LIST TABLE-------------------------------------
class Tag_List(QtGui.QWidget):
    
    def __init__(self,tag_list,parent=None):
        super(Tag_List,self).__init__()
        self.setWindowTitle('Tag List')
        self.setGeometry(50,50, 300, 800)

        #Grid Layout-------------
        grid = QtGui.QGridLayout()
        self.setLayout(grid)
        
        #Create empty table 2*100
        self.tag_table = QtGui.QTableWidget(self)

        
        self.tag_list = tag_list #Container for tag list
        
        
        grid.addWidget(self.tag_table,0,0)
        
    def array_2_table(self,array,table):
        table = self.tag_table
        self.tag_table.setRowCount(100)
        self.tag_table.setColumnCount(1)
        for row in range(100):
                table.setItem(row,0,QtGui.QTableWidgetItem(str(array[row]))) #Needs to be string to go in table
            
                
        table.resizeColumnsToContents()
        table.resizeRowsToContents()
        
        
#--------------------------SPIKE DETECTION WINDOW------------------------------
class Spike_Detect_Popup(QtGui.QWidget):
    
    def __init__(self,ephy_file,parent=None):
        super(Spike_Detect_Popup,self).__init__()
        self.setWindowTitle('Spike Detection')
        
        #Grid Layout-------------
        grid = QtGui.QGridLayout()
        self.setLayout(grid)
                    
        #Canvas and Toolbar--------------------
        self.figure = plt.figure(figsize=(1,1))  
        plt.ioff() #Avoid figure poping as extra window 
        self.canvas = FigureCanvas(self.figure)     
        self.toolbar = NavigationToolbar(self.canvas, self)
        grid.addWidget(self.canvas, 3,0,1,6)
        grid.addWidget(self.toolbar, 2,0,1,2)
        
        self.ephy_file = ephy_file #Container for trace loaded with MainWindow
               
        #Do spike detect btn 
        Do_Avg_Btn = QtGui.QPushButton('Detect these f*cking spikes', self)
        Do_Avg_Btn.resize(Do_Avg_Btn.sizeHint()) 
        #Do_Avg_Btn.clicked.connect(self.do_average)
        grid.addWidget(Do_Avg_Btn, 0,4)
        
        #Threshold textbox
        Thr_Tbox = QtGui.QLineEdit(self)
        Thr_Tbox.resize(Thr_Tbox.sizeHint())
        grid.addWidget(Thr_Tbox, 0,1)
        Thr_Tbox.setText("Threshold")      
        
        #Distance textbox
        Dist_Tbox = QtGui.QLineEdit(self)
        Dist_Tbox.resize(Dist_Tbox.sizeHint())
        grid.addWidget(Dist_Tbox, 0,2)
        Dist_Tbox.setText("Distance")     

#------------------------AVERAGING WINDOW--------------------------------------
class Average_Popup(QtGui.QWidget):
    
    def __init__(self, voltageclamp_file, parent = None):  
        super(Average_Popup,self).__init__()
        #QtGui.QWidget.__init__(self, parent)
        #self.setGeometry(500,500, 800, 800)
        self.setWindowTitle('Signal Averaging') 
        
        #Grid Layout-------------
        grid = QtGui.QGridLayout()
        self.setLayout(grid)
                    
        #Canvas and Toolbar--------------------
        self.figure = plt.figure(figsize=(1,1))  
        plt.ioff() #Avoid figure poping as extra window 
        self.canvas = FigureCanvas(self.figure)     
        self.toolbar = NavigationToolbar(self.canvas, self)
        grid.addWidget(self.canvas, 3,0,1,6)
        grid.addWidget(self.toolbar, 2,0,1,2)
        
        self.voltageclamp_file = voltageclamp_file #Container for trace loaded with MainWindow

        
        #Do average btn for VC, without leak 
        Do_Avg_Btn = QtGui.QPushButton('Do The F*cking Average (VC)', self)
        Do_Avg_Btn.resize(Do_Avg_Btn.sizeHint()) 
        Do_Avg_Btn.clicked.connect(self.do_average)
        grid.addWidget(Do_Avg_Btn, 0,0)
        
        #Do average btn for CC, WITH leak 
        Do_Avg_Btn_CC = QtGui.QPushButton('Do The F*cking Average (CC)', self)
        Do_Avg_Btn_CC.resize(Do_Avg_Btn_CC.sizeHint()) 
        Do_Avg_Btn_CC.clicked.connect(self.do_average_CC)
        grid.addWidget(Do_Avg_Btn_CC, 0,1)
        
    #Do Average method for VC
    def do_average(self,voltageclamp_file):
        try :
            voltageclamp_file = self.voltageclamp_file
                        
                    
            y = np.ravel(np.nanmean(voltageclamp_file[:,:,1],axis=0)) #The average

            x = np.ravel(voltageclamp_file[0,:,0]) #The time
            
            plt.cla()
            ax = self.figure.add_subplot(111)
            ax.clear()
            ax.plot(x,y, 'b-')
            ax.set_title('Average Trace from tagged files')
            ax.set_xlabel('Time (s)')
            ax.set_ylabel('Signal')
            self.canvas.draw()
            
            x,y = None,None #memory clear 

        except Exception:
            traceback.print_exc()
            
    #Do Average method for CC
    def do_average_CC(self,currentclamp_file):
        try :
            currentclamp_file = self.currentclamp_file
                        
                    
            y = np.ravel(np.nanmean(currentclamp_file[:,:,1],axis=0)) #The average

            x = np.ravel(currentclamp_file[0,:,0]) #The time
            
            plt.cla()
            ax = self.figure.add_subplot(111)
            ax.clear()
            ax.plot(x,y, 'b-')
            ax.set_title('Average Trace from tagged files')
            ax.set_xlabel('Time (s)')
            ax.set_ylabel('Signal')
            self.canvas.draw()
            
            x,y = None,None #memory clear 

        except Exception:
            traceback.print_exc()


#-----------------------MAIN WINDOW--------------------------------------------
class MainWindow(QtGui.QWidget):
    
    def __init__(self,parent=None):
        super(MainWindow, self).__init__()
        
        self.dialog = Average_Popup(self) #Calls average pop up window
        self.dialog_tag_list = Tag_List(self) #Calls Table Tag window
        self.dialog_spike_detect = Spike_Detect_Popup(self)
        
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
        self.figure = plt.figure(figsize=(1,1))   
        plt.ioff() #Avoid figure poping as extra window 
        self.canvas = FigureCanvas(self.figure)     
        self.toolbar = NavigationToolbar(self.canvas, self)
        grid.addWidget(self.canvas, 3,0,1,6)
        grid.addWidget(self.toolbar, 2,0,1,2)
        
        #Import File Button----------------------------------
        btn1 = QtGui.QPushButton('Import ePhy File (.wcp) ...', self)
        btn1.resize(btn1.sizeHint()) 
        btn1.clicked.connect(self.getFile)
        grid.addWidget(btn1, 0,0)
        
        #Plot Button------------------------------
        btn2 = QtGui.QPushButton('Plot (P)', self)
        btn2.resize(btn2.sizeHint())    
        btn2.clicked.connect(self.plot)
        btn2.setShortcut("p")
        grid.addWidget(btn2, 0,1)
        
        #Sweep text----------------
        sweep_text = QtGui.QLabel()
        sweep_text.setText("Sweep #")
        grid.addWidget(sweep_text,0,2)
        
        #Sweep combo box----------------
        combobox = QtGui.QComboBox(self)
        grid.addWidget(combobox,0,3)
        sweep_list = np.arange(0,100,1).astype(str).tolist()
        combobox.addItems(sweep_list)
        combobox.activated.connect(self.combo_box_update)
        current_sweep = int(combobox.currentText()) #Used in functions later very important !!
        
        #Channel text -----------------
        channel_text = QtGui.QLabel()
        channel_text.setText("Channel select")
        grid.addWidget(channel_text,0,4)
        
        #Channel combo box-----------------
        ch_combobox = QtGui.QComboBox(self)
        grid.addWidget(ch_combobox,0,5)
        ch_list = np.arange(0,2,1).astype(str).tolist()
        ch_combobox.addItems(ch_list)
        ch_combobox.activated.connect(self.ch_box_update)
        current_ch = int(ch_combobox.currentText()) #Used in functions later very important !!

        
        #Previous sweep btn---------------------------------
        Prev_btn = QtGui.QPushButton('Previous (left)', self)
        Prev_btn.resize(Prev_btn.sizeHint())    
        Prev_btn.clicked.connect(self.previous_sweep)
        Prev_btn.setShortcut("Left")
        grid.addWidget(Prev_btn, 1,0)
          
        #Next sweep btn----------------------------------
        Next_btn = QtGui.QPushButton('Next (right)', self)
        Next_btn.resize(Next_btn.sizeHint())    
        Next_btn.clicked.connect(self.next_sweep)
        grid.addWidget(Next_btn, 1,1)
        Next_btn.setShortcut("Right")
        
        #Tag button-----------------------------------
        Tag_button = QtGui.QPushButton('Tag (up)',self)
        Tag_button.resize(Tag_button.sizeHint())
        Tag_button.clicked.connect(self.tag_this)
        grid.addWidget(Tag_button,1,2)
        Tag_button.setShortcut("Up")
    
        #UnTag button----------------------------------------
        UnTag_button = QtGui.QPushButton('UnTag (down)',self)
        UnTag_button.resize(UnTag_button.sizeHint())
        UnTag_button.clicked.connect(self.untag_this)
        grid.addWidget(UnTag_button,1,3)
        UnTag_button.setShortcut("Down")
        
        #Tag All button----------------------------------
        TagAll_button = QtGui.QPushButton('Tag All',self)
        TagAll_button.resize(TagAll_button.sizeHint())
        TagAll_button.clicked.connect(self.tag_all)
        grid.addWidget(TagAll_button,1,4)
        
        #UnTag All button------------------------------------
        UnTagAll_button = QtGui.QPushButton('UnTag All',self)
        UnTagAll_button.resize(UnTagAll_button.sizeHint())
        UnTagAll_button.clicked.connect(self.untag_all)
        grid.addWidget(UnTagAll_button,1,5)

                     
        #Uniform time toggle---------------------------------
        Set_Time_Toggle = QtGui.QCheckBox('Unify Time', self)
        Set_Time_Toggle.stateChanged.connect(self.unify_time)
        grid.addWidget(Set_Time_Toggle,2,3)
        self.time = 0 #Defaut time is conserved
        #Set_Time_Toggle.toggle()  #Uncomment if u want time to be set by default              
               
        #Leak Remove toggle--------------------------------
        LeakcheckBox = QtGui.QCheckBox('Remove Leak', self)
        LeakcheckBox.stateChanged.connect(self.leak_remove)
        grid.addWidget(LeakcheckBox,2,2)
        self.leak = 0 #Defaut is leak conserved 
        #LeakcheckBox.toggle()  #Uncomment if u want leak to be removed by default

        #Tag state------------------------------------------------------

        self.tag_text = QtGui.QLabel()
        self.tag_text.setText("---")
        grid.addWidget(self.tag_text,2,4)        
        
        #Tag List btn--------------------------------------
        TagList_Btn = QtGui.QPushButton('Tag List...',self)
        TagList_Btn.resize(TagList_Btn.sizeHint())
        TagList_Btn.clicked.connect(self.get_tag_list)
        grid.addWidget(TagList_Btn,2,5)
        
        #Average button--------------------------------------------------------
        self.Av_Btn = QtGui.QPushButton('Average...')
        self.Av_Btn.resize(self.Av_Btn.sizeHint())
        self.Av_Btn.clicked.connect(self.Average)
        grid.addWidget(self.Av_Btn,4,0)
        
        #Spike detection button 
        self.Spike_detection_btn = QtGui.QPushButton('Spike Detection...')
        self.Spike_detection_btn.resize(self.Spike_detection_btn.sizeHint())
        self.Spike_detection_btn.clicked.connect(self.Spike_detect)
        grid.addWidget(self.Spike_detection_btn,4,1)
            
        
        self.file = np.ones((10,10,2)) #Default file before any load 
    
        self.show()

    #------------------------METHODS-----------------------------------------------                 
    def getFile(self): #Penser a ajouter le channel select ici 0 par defaut
        filePath = QtGui.QFileDialog.getOpenFileName(self,'Choose ur file')       
        try : 
            self.file = WinWcpIo(str(filePath)).whole_file(0)

            self.plot(0)
            self.plot(0) #twice, to plot directly after load w/out hitting plot btn 
        
        except Exception:
            traceback.print_exc()
            
            
        self.setWindowTitle(str(filePath))  
        
        
    def tag_this(self,sweep_id):
        try :
            sweep_id = self.current_sweep  #ID of current sweep 
            self.tag_list[sweep_id]=1
            self.tag_text.setText("Tagged !!")
        except Exception:
            traceback.print_exc()
            
            
    def untag_this(self,sweep_id):
        try :
            sweep_id = self.current_sweep  #ID of current sweep 
            self.tag_list[sweep_id]=0
        except Exception:
            traceback.print_exc()
            
            
    def tag_all(self):
        try:
            for i in range(len(self.tag_list)):
                self.tag_list[i]=1
        except Exception:
            traceback.print_exc()
            
            
    def untag_all(self):
        try :
            for i in range(len(self.tag_list)):
                self.tag_list[i]=0
        except Exception:
            traceback.print_exc()
            
            
    def get_tag_list(self):
        try:
            self.dialog_tag_list.array_2_table(self.tag_list,self.dialog_tag_list)
            
            self.dialog_tag_list.show()
            
        except Exception:
            traceback.print_exc()
             
            
    def leak_remove(self,state):
        try :             
            if state == QtCore.Qt.Checked:
                self.leak = 1
                self.dialog.leak = 1
                print ('Leak is removed'  )    
            else:
                self.leak = 0 
                self.dialog.leak = 0 
                print ('Leak is conserved')     
        except Exception : #For first lauch 
            self.leak = 0
        
        
    def unify_time(self,state):
        try :             
            if state == QtCore.Qt.Checked:
                self.time = 1
                print ('Time is the same for all sweeps')    
            else:
                self.time = 0 
                print ('Experimental time conserved')            
        except Exception : #For first lauch 
            self.leak = 0
      
        
    def combo_box_update(self,current_sweep):
        try:
            print ('Displaying sweep#', current_sweep)
            self.current_sweep = current_sweep
            self.plot(current_sweep)
        except Exception :
            traceback.print_exc()
      
        
    def ch_box_update(self,current_ch):
        try:
            print ('Channel #', current_ch)
            self.current_ch = current_ch
            self.plot(self.current_sweep)
        except Exception :
            traceback.print_exc()

                        
    def plot(self,sweep_id):   
        try: 
            sweep_id = self.current_sweep   
            if self.time == 0 : #Time is conserved 
                x = self.file[sweep_id,:,0] #Timepoints
                
            else : #Time is set to default for all sweeps ie from the first sweep : 
                x = self.file[0,:,0] #Timepoints


            if self.leak == 0: #LEak is conserved
               
                y = self.file[sweep_id,:,1] #Sweep, points current
                    
            else : #Implicits leak == 1 so remove leak :
                
                leak = np.mean(self.file[sweep_id,50:950,1])
                y = self.file[sweep_id,:,1]-leak #Sweep, points current

            if self.tag_list[sweep_id]==1 :
                self.tag_text.setText("Tagged !!")
                
            else : 
                self.tag_text.setText("---")

            plt.cla()
            ax = self.figure.add_subplot(111)
            ax.clear()
            ax.plot(x,y, 'b-')
            ax.set_title('AnalogSignal Sweep#%s / %s'%(self.current_sweep,self.file.shape[0]))
            ax.set_xlabel('Time (s)')
            ax.set_ylabel('Signal')
            self.canvas.draw()
            
            x = None #To clear and free memory for plot 
            y = None
                
        except Exception: 
            #Creates 0 sweep number at first launch 
            self.current_sweep = 0
                
        
    def next_sweep(self,sweep_id):
        try : 
            self.current_sweep = self.current_sweep+1
            sweep_id = self.current_sweep
            self.plot(sweep_id)
            print ('Displaying sweep#', sweep_id)
        except Exception : 
            traceback.print_exc()
      
        
    def previous_sweep(self,sweep_id):
        try : 
            self.current_sweep = self.current_sweep-1
            sweep_id = self.current_sweep
            self.plot(sweep_id)
            print ('Displaying sweep#', sweep_id)
        except Exception : 
            traceback.print_exc()

    
    def center(self):
        qr = self.frameGeometry()
        cp = QtGui.QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())


        
    def get_corrected_traces(self): #returns matrix with corrected traces, ideal for further analysis 
        self.corrected_traces = np.zeros((self.file.shape[0],self.file.shape[1],2)) #Matrix to hold non leaked traces or nans and times
        
        for trace in range(self.file.shape[0]):
            leak = np.nanmean(self.file[trace,50:950,1])
            self.file[trace,:,1] = self.file[trace,:,1]-leak


        
    def get_tagged_traces(self): #WITH LEAK !!!!!
        self.tagged_traces = np.zeros((self.file.shape[0],self.file.shape[1],2)) #Matrix to hold tagged traces or nans and times
        
        for trace in range(self.file.shape[0]):
            if self.tag_list[trace] == 1: #The trace has been tagged                
                self.tagged_traces[trace,:,1] = self.file[trace,:,1]         
                self.tagged_traces[trace,:,0] = self.file[trace,:,0]
                
            else: #Trace is rejected
                self.tagged_traces[trace,:,1] = np.nan  
                
    def get_corr_tagged_traces(self): #W/OUT LEAK !!!!!!!!!!
        self.corr_tagged_traces = np.zeros((self.file.shape[0],self.file.shape[1],2)) #Matrix to hold tagged traces or nans and times
        
        for trace in range(self.file.shape[0]):
            if self.tag_list[trace] == 1: #The trace has been tagged    
                leak = np.nanmean(self.file[trace,50:950,1])
                self.corr_tagged_traces[trace,:,1] = self.file[trace,:,1]-leak       
                self.corr_tagged_traces[trace,:,0] = self.file[trace,:,0]
                
            else: #Trace is rejected
                self.corr_tagged_traces[trace,:,1] = np.nan  
       
        
    def Average(self):
        try:
            #print ('Averaging')
            
            self.get_tagged_traces() #RAW TRACES
            self.get_corr_tagged_traces() #CORRECTED TRACES
            
            self.data.get_raw_data(self.tagged_traces)
            #self.data.get_corr_data(self.corr_tagged_traces)
            
            self.dialog.show()

                        
        except Exception:
            traceback.print_exc()
            
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
    #sys.exit(app.exec_())
