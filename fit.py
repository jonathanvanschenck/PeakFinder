import tkinter as tk
from scipy.optimize import least_squares
import numpy as np
from tkinter import filedialog


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


class Fit:
    def __init__(self,master,params,mpl,**kwargs):
        self.p = params
        self.mpl = mpl
        self.frame = tk.LabelFrame(master,text="Fit Setup")
        self.frame.grid(**kwargs)
        self.createWidgets()
        
    def createWidgets(self):
        self.convergenceFrame = tk.LabelFrame(self.frame,text="Convergence")
        self.conv = tk.StringVar()
        self.convB = [tk.Radiobutton(self.convergenceFrame,variable=self.conv,
                                     value="Tight",text="Tight",justify=tk.LEFT),
                      tk.Radiobutton(self.convergenceFrame,variable=self.conv,
                                     value="Loose",text="Loose",justify=tk.LEFT),
                      tk.Radiobutton(self.convergenceFrame,variable=self.conv,
                                     value="Danger",text="Danger",justify=tk.LEFT)]
        self.conv.set("Tight")
        self.fitButton = tk.Button(self.frame,text="<<Fit>>",command=self.fit)
        self.saveButton = tk.Button(self.frame,text="save",command=self.save)
        
        # Pack
        self.convergenceFrame.grid(row=0,rowspan=2,column=0)
        for i in range(len(self.convB)):
            self.convB[i].grid(row=i,column=0)
        self.fitButton.grid(row=0,column=1)
        self.saveButton.grid(row=1,column=1)
            
    
    def fit(self):
        x = self.mpl.data.get_xdata()
        y = self.mpl.data.get_ydata()
        res = self.p.get()
        xL,xR = res[0]
        f = fitDic[res[1]]
        pi = np.array(res[2])
        which = np.array(res[3])
        d = np.array([max(0.01,abs(p)) for p in pi])[which]
        if self.conv.get() == "Danger":
            bound = np.transpose(sum(which)*[[-np.inf,np.inf]])
        elif self.conv.get() == "Loose":
            bound = [pi[which]-d,pi[which]+d]
        elif self.conv.get() == "Tight":
            bound = [pi[which]-0.1*d,pi[which]+0.1*d]
        mask = np.array([xi>xL and xi<xR for xi in x])
        def diff(par):
            p = np.copy(pi)
            p[which] = par
            return (y[mask]-tot(x[mask],p,f))
        self.pp = least_squares(diff,pi[which],bounds=bound)
        #print(self.pp)
        if self.pp.success:
            pi[which] = self.pp.x
            self.p.set(pi)
            self.mpl.updateFits()
        else:
            print(self.pp)
            
    def save(self):
        fl = filedialog.asksaveasfilename(title="Save As",
                                          filetypes = (("text files","*.txt"),
                                                       ("all files","*.*")))
        f = open(fl,"w")
        res = self.p.get()
        f.write("Fit Type: {}\n".format(res[1]))
        f.write("Fit Bounts: {0},{1}\n".format(*res[0]))
        f.write("Baseline={0} ({1})\n".format(res[2][0],res[3][0]))
        for i in range(len(res[2][1:])//3):
            f.write("A={0},G={1},X={2}\n".format(*res[2][1+i*3:1+3+i*3]))
        f.write("#########\n")
        f.write(",".join([str(ii) for ii in res[2]])+"\n")
        f.write(",".join([str(ii) for ii in res[3]])+"\n")
        f.write("#########\n")
        x = self.mpl.data.get_xdata()
        y = self.mpl.data.get_ydata()
        fit = self.mpl.Fdatatot.get_ydata()
        fitpiece = [l.get_ydata() for l in self.mpl.FdataL]
        comb = np.transpose([x,y,fit]+fitpiece)
        for i in range(len(x)):
            f.write(",".join([str(ii) for ii in comb[i]])+"\n")
        f.close()
            