# -*- coding: utf-8 -*-
import sys, os, random
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from ukolObjekt import *

import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure


class AppForm(QMainWindow):
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.setWindowTitle('Ukol Bayes')

        self.create_main_frame()

        #self.textbox.setText('1 2 3 4')
        self.on_draw()
    

    
    def on_draw(self):

        
        if (self.jednotCov.isChecked()):
            priorcov=TypPriorCov.JEDNOTKOVA
        elif (self.kolemJCov.isChecked()):
            priorcov=TypPriorCov.KOLEM_JEDNOTKOVE
        elif (self.levoCov.isChecked()):
            priorcov=TypPriorCov.JEDNICKY_VLEVO_NAHORE
        elif (self.pravoCov.isChecked()):
            priorcov=TypPriorCov.JEDNICKY_VPRAVO_DOLE
        else:
            priorcov=TypPriorCov.RANDOM
        
        koef = 1
        if (self.kratSto.isChecked()):
            koef = koef*100
        if (self.kratDeset.isChecked()):
            koef=koef*10
        

        if (self.zeroMean.isChecked()):
            mean = TypPriorMean.VSE_NULA
        elif (self.tisicMean.isChecked()):
            mean = TypPriorMean.VSE_TIS
        elif (self.startMean.isChecked()):
            mean=TypPriorMean.ZACATEK
        else:
            mean=TypPriorMean.KONEC


        if (self.sinFce.isChecked()):
            fce = TypMnozinaFci.SINY_POSUN
        elif (self.sinPosFce.isChecked()):
            fce = TypMnozinaFci.SINY_FREKVPLUSPOSUN
        else:
            fce = TypMnozinaFci.POLYNOMY

        if (self.meanZData.isChecked()):
            typdat = TypDat.KOLEM_NULY
        elif (self.meanKData.isChecked()):
            typdat = TypDat.KOLEM_KONSTANTY
        elif (self.sinData.isChecked()):
            typdat = TypDat.KONST_KRAT_SIN_1
        else:
            typdat = TypDat.KONST_KRAT_SIN_2


        sigma=0
        alpha=0
        if (self.sigmaJRozptyl.isChecked()):
            sigma = 1.0
        elif (self.sigmaSRozptyl.isChecked()):
            sigma=0.0001
        elif (self.sigmaPRozptyl.isChecked()):
            sigma=10000.0
        elif (self.alphaJRozptyl.isChecked()):
            alpha=1.0
        elif (self.alphaSRozptyl.isChecked()):
            alpha=0.0001
        else:
            alpha=10000.0


        ukol = UkolInference(priorcov, koef, mean, fce, typdat, sigma=sigma, gammaA=alpha)


        self.axes.clear()        
       
        ukol.maluj_nahodna_data(self.axes)

        if (self.buttonPrior.isChecked()) :
            ukol.maluj_nahodne_priory(self.axes)
        elif (self.buttonBazove.isChecked()):
            ukol.maluj_base_fce(self.axes)
        else:
            ukol.maluj_nahodne_posteriory(self.axes)

        self.canvas.draw()
    
    def create_main_frame(self):
        self.main_frame = QWidget()
        
        # Create the mpl Figure and FigCanvas objects. 
        # 5x4 inches, 100 dots-per-inch
        #
        self.dpi = 100
        self.fig = Figure((5.0, 4.0), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)
        
        # Since we have only one plot, we can use add_axes 
        # instead of add_subplot, but then the subplot
        # configuration tool in the navigation toolbar wouldn't
        # work.
        #
        self.axes = self.fig.add_subplot(111)
        

        
        # Create the navigation toolbar, tied to the canvas
        #
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
        
        # Other GUI controls
        # 
        #self.textbox = QLineEdit()
        #self.textbox.setMinimumWidth(200)
       
    #def __init__(self, typPriorCov, covKoef, typPriorMean, typMnozinaFci, typDat, sigma=0, gammaA=0):

        buttony=[]
        self.vboxl = QVBoxLayout();
        self.radiogroup = QGroupBox()
        self.buttonPrior = QRadioButton(u"rnd. priory")
        self.buttonBazove = QRadioButton(u"Bázové funkce")
        self.buttonPosterior = QRadioButton(u"rnd. posteriory")
        self.buttonPosterior.setChecked(True)
        self.vboxl.addWidget(self.buttonBazove)
        self.vboxl.addWidget(self.buttonPrior)
        self.vboxl.addWidget(self.buttonPosterior)
        self.radiogroup.setLayout(self.vboxl)
        buttony+=[self.buttonPrior,self.buttonBazove,self.buttonPosterior]

        self.vboxCov = QVBoxLayout();
        self.radioboxCov = QGroupBox()
        self.jednotCov = QRadioButton(u"Cov: jednotková")
        self.jednotCov.setChecked(True)
        self.kolemJCov = QRadioButton(u"kolem jednotkové")
        self.levoCov = QRadioButton(u"vlevo nahoře")
        self.pravoCov = QRadioButton(u"vpravo dole")
        self.randomCov = QRadioButton(u"rnd")
        self.kratSto = QCheckBox(u"x100")
        self.kratDeset = QCheckBox(u"x10")
        self.vboxCov.addWidget(self.jednotCov)
        self.vboxCov.addWidget(self.kolemJCov)
        self.vboxCov.addWidget(self.levoCov)
        self.vboxCov.addWidget(self.pravoCov)
        self.vboxCov.addWidget(self.randomCov)
        self.vboxCov.addWidget(self.kratSto)
        self.vboxCov.addWidget(self.kratDeset)
        self.radioboxCov.setLayout(self.vboxCov)
        buttony+=[self.jednotCov, self.kolemJCov,self.levoCov,
                  self.pravoCov, self.randomCov, self.kratSto, self.kratDeset]

        self.vboxMean = QVBoxLayout();
        self.radioboxMean = QGroupBox()
        self.zeroMean = QRadioButton(u"mean: vše 0")
        self.zeroMean.setChecked(True)
        self.tisicMean = QRadioButton(u"vše 1000")
        self.startMean = QRadioButton(u"začátek 1000")
        self.konecMean = QRadioButton(u"konec 1000")
        self.vboxMean.addWidget(self.zeroMean)
        self.vboxMean.addWidget(self.tisicMean)
        self.vboxMean.addWidget(self.startMean)
        self.vboxMean.addWidget(self.konecMean)
        self.radioboxMean.setLayout(self.vboxMean)
        buttony += [self.zeroMean, self.tisicMean, 
            self.startMean, self.konecMean]


        self.vboxFce = QVBoxLayout();
        self.radioboxFce = QGroupBox()
        self.sinFce = QRadioButton(u"funkce:sinus posun")
        self.sinFce.setChecked(True)
        self.sinPosFce = QRadioButton(u"sinus posun+frekvence")
        self.polFce = QRadioButton(u"polynom")
        self.vboxFce.addWidget(self.sinFce)
        self.vboxFce.addWidget(self.sinPosFce)
        self.vboxFce.addWidget(self.polFce)
        self.radioboxFce.setLayout(self.vboxFce)
        buttony+=[self.sinFce, self.sinPosFce, self.polFce]


        self.vboxData = QVBoxLayout();
        self.radioboxData = QGroupBox()
        self.meanZData = QRadioButton(u"data:rnd mean 0")
        self.meanKData = QRadioButton(u"rnd mean 30")
        self.sinData = QRadioButton(u"rnd*sinus")
        self.sinData.setChecked(True)
        self.sinDData = QRadioButton(u"rnd*sinus(2x)")
        self.vboxData.addWidget(self.meanZData)
        self.vboxData.addWidget(self.meanKData)
        self.vboxData.addWidget(self.sinData)
        self.vboxData.addWidget(self.sinDData)
        self.radioboxData.setLayout(self.vboxData)
        buttony+=[self.meanZData, self.meanKData, 
            self.sinData,self.sinDData]


        self.vboxRozptyl = QVBoxLayout();
        self.radioboxRozptyl = QGroupBox()
        self.sigmaJRozptyl = QRadioButton(u"σ²:1")
        self.sigmaJRozptyl.setChecked(True)
        self.sigmaSRozptyl = QRadioButton(u"σ²:0.0001")
        self.sigmaPRozptyl = QRadioButton(u"σ²:10000")
        self.alphaJRozptyl = QRadioButton(u"α:1")
        self.alphaSRozptyl = QRadioButton(u"α:0.0001")
        self.alphaPRozptyl = QRadioButton(u"α:10000")
        self.vboxRozptyl.addWidget(self.sigmaJRozptyl)
        self.vboxRozptyl.addWidget(self.sigmaSRozptyl)
        self.vboxRozptyl.addWidget(self.sigmaPRozptyl)
        self.vboxRozptyl.addWidget(self.alphaJRozptyl)
        self.vboxRozptyl.addWidget(self.alphaSRozptyl)
        self.vboxRozptyl.addWidget(self.alphaPRozptyl)
        self.radioboxRozptyl.setLayout(self.vboxRozptyl)

        buttony+=[self.sigmaJRozptyl, self.sigmaSRozptyl, self.sigmaPRozptyl,
                  self.alphaJRozptyl, self.alphaSRozptyl, self.alphaPRozptyl]

#        self.connect(self.textbox, SIGNAL('editingFinished ()'), self.on_draw)
        
        #self.draw_button = QPushButton(u"Kresli")
        #self.connect(self.draw_button, SIGNAL('clicked()'), self.on_draw)
        
        #self.grid_cb = QCheckBox("Show &Grid")
        #self.grid_cb.setChecked(False)
        #self.connect(self.grid_cb, SIGNAL('stateChanged(int)'), self.on_draw)
        
        for button in buttony:
            #self.connect(button, SIGNAL('stateChanged(int)'), self.on_draw)
            self.connect(button, SIGNAL('clicked()'), self.on_draw)
            

        #slider_label = QLabel('Bar width (%):')
        #self.slider = QSlider(Qt.Horizontal)
        #self.slider.setRange(1, 100)
        #self.slider.setValue(20)
        #self.slider.setTracking(True)
        #self.slider.setTickPosition(QSlider.TicksBothSides)
        #self.connect(self.slider, SIGNAL('valueChanged(int)'), self.on_draw)
        
        #
        # Layout with box sizers
        # 
        hbox = QHBoxLayout()
        
        for w in [  self.radiogroup, self.radioboxCov, self.radioboxMean, 
                     self.radioboxFce, self.radioboxData,
                     self.radioboxRozptyl]:
            hbox.addWidget(w)
            hbox.setAlignment(w, Qt.AlignVCenter)
        
        vbox = QVBoxLayout()
        vbox.addWidget(self.canvas)
        vbox.addWidget(self.mpl_toolbar)
        vbox.addLayout(hbox)
        
        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)
    
   
        
    

    def add_actions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def create_action(  self, text, slot=None, shortcut=None, 
                        icon=None, tip=None, checkable=False, 
                        signal="triggered()"):
        action = QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon(":/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            self.connect(action, SIGNAL(signal), slot)
        if checkable:
            action.setCheckable(True)
        return action


def main():
    app = QApplication(sys.argv)
    form = AppForm()
    form.showMaximized()
    app.exec_()


if __name__ == "__main__":
    main()
