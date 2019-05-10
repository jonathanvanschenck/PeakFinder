import tkinter as tk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import numpy as np

def Glorz(x,mu,sig):
    return sig/((2*np.pi)*((x-mu)**2+(sig/2)**2))

def Gauss(x,mu,doubsig):
    sig = doubsig/2
    return (1/(sig*np.sqrt(2*np.pi)))*np.exp(-(x-mu)**2/(2*sig**2))

fitDic = {"Gaussian":Gauss,"Lorentzian":Glorz}

def tot(x,par,f):
    x2 = np.array(x)
    res = 0*x2+par[0]
    for i in range(len(par[1:])//3):
        res = res + par[1+i*3+0]*f(x2,par[1+i*3+2],par[1+i*3+1])
    return res

def part(x,par,f,i):
    x2 = np.array(x)
    res = 0*x2+par[0]
    res = res + par[1+i*3+0]*f(x2,par[1+i*3+2],par[1+i*3+1])
    return res

class MPL:
    def __init__(self,master,x,y,params,**kwargs):
        self.p = params
        self.frame = tk.LabelFrame(master)
        self.frame.grid(**kwargs)
        self.fig = plt.figure(figsize=(7,7))
        spec = self.fig.add_gridspec(ncols=1, nrows=5)
        self.ax = self.fig.add_subplot(spec[1:,0])
        self.axR = self.fig.add_subplot(spec[0,0])
        self.axR.get_xaxis().set_visible(False)
        plt.subplots_adjust(hspace=0,top=0.95)
        self.canvas = FigureCanvasTkAgg(self.fig, self.frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack()
        self.toolbar = NavigationToolbar2Tk(self.canvas,self.frame)
        self.toolbar.update()
        
        # Initalize Plots
        self.data = self.ax.plot(x,y)[0]
        self.axR.plot(x,len(x)*[0],color="Black")
        self.Fdatatot = self.ax.plot(x,0*np.array(x))[0]
        xL,xR = self.p.get()[0]
        mask = np.array([xi>xL and xi<xR for xi in x])
        self.Fdatatot2 = self.ax.plot(x[mask],0*np.array(x)[mask])[0]
        self.residue = self.axR.plot(x,self.data.get_ydata()-self.Fdatatot.get_ydata())[0]
        self.FdataL = []
        self.ax.set_ylim(np.min(y)*0.9,np.max(y)*1.1)
        self.updateFits()
        
        # Add click call backs
        self.canvas.mpl_connect('button_press_event',self.click)

        
    def click(self,event):
        if event.inaxes == self.ax:
            if event.button == 1:
                #print("Left click @ x=",event.xdata," y=",event.ydata)
                self.p.update(event)
            if event.button == 2:
                #print("Scroll click @ x=",event.xdata," y=",event.ydata)
                self.p.advance(event)
            if event.button == 3:
                #print("Right click @ x=",event.xdata," y=",event.ydata)
                self.p.update(event)
            self.updateFits()

        
    def updateFits(self):
        x = self.data.get_xdata()
        res = self.p.get()
        f = fitDic[res[1]]
        par = res[2]
        yL = [part(x,par,f,i) for i in range(len(par[1:])//3)]
        yL = [tot(x,par,f)]+yL
        for l in self.FdataL:
            l.remove()
        self.FdataL = [self.ax.plot(x,yL[ii+1],color=["Grey","Black"][ii==self.p.i-1],linestyle="--")[0] for ii in range(len(yL[1:]))]
        self.Fdatatot.set_ydata(yL[0])
        self.residue.set_ydata(self.data.get_ydata()-self.Fdatatot.get_ydata())
        self.axR.set_ylim(top=max(self.residue.get_ydata()),
                          bottom=min(self.residue.get_ydata()))
        xL,xR = self.p.get()[0]
        mask = np.array([xi>xL and xi<xR for xi in x])
        self.Fdatatot2.set_xdata(x[mask])
        self.Fdatatot2.set_ydata(self.Fdatatot.get_ydata()[mask])
        self.fig.canvas.draw()
        