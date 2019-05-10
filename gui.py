import tkinter as tk
import mpl
import params
import fit
from numpy import diff,mean,abs
import numpy as np

x = np.arange(0,10,0.01)
y = np.abs(np.sin(x)+np.random.normal(0,.1,len(x)))


class App(tk.Frame):
    def __init__(self,root,x,y):
        tk.Frame.__init__(self,root)
        self.root = root
        self.root.title("Bump Fit")
        self.grid(row=0,column=0)
        self.createWidgets(x,y)
        
    def createWidgets(self,x,y):
        self.p = params.FitParams(self.root,Gs=10*abs(mean(diff(x))),
                                  xL=min(x),xR=max(x),row=0,column=1)
        self.p.bl.set(mean(y))
        self.mpl = mpl.MPL(self.root,x,y,self.p,row=0,column=0,rowspan=3)
        self.fit = fit.Fit(self.root,self.p,self.mpl,row=1,column=1)
        
    def returnFitParams(self):
        return self.p.get()[2]
        
        
if __name__ == "__main__":
    root = tk.Tk()
    app = App(root,x,y)
    app.mainloop()
    print(app.returnFitParams())