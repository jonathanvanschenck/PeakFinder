import tkinter as tk
import numpy as np

def Glorz(x,mu,sig):
    return sig/((2*np.pi)*((x-mu)**2+(sig/2)**2))

def Gauss(x,mu,doubsig):
    sig = doubsig/2
    return (1/(sig*np.sqrt(2*np.pi)))*np.exp(-(x-mu)**2/(2*sig**2))

fitDic = {"Gaussian":Gauss,"Lorentzian":Glorz}

color = ["White","Yellow"]

class FitParams:
    def __init__(self,master,Gs=1.0,xL=0.0,xR=1.0,**kwargs):
        self.Gs = Gs
        self.frame = tk.LabelFrame(master,text="Parameters")
        self.frame.grid(**kwargs)
        # Setup
        self.controlFrame = tk.Frame(self.frame)
        self.addButton = tk.Button(self.controlFrame,text="Add Bump",command=self.addBump)
        self.removeButton = tk.Button(self.controlFrame,text="Remove Bump",command=self.removeBump)
        self.bounds = BoundObj(self.frame,xL=xL,xR=xR,row=1,column=0)
        self.typeFrame = tk.Frame(self.frame)
        self.typeLabel = tk.Label(self.typeFrame,text="Type:",justify=tk.RIGHT)
        self.typeVar = tk.StringVar()
        self.typeB = [tk.Radiobutton(self.typeFrame,text="Gaussian",
                                     variable=self.typeVar,value="Gaussian"),
                      tk.Radiobutton(self.typeFrame,text="Lorentzian",
                                     variable=self.typeVar,value="Lorentzian")]
        self.typeVar.set("Gaussian")
        self.bl = BaselineObj(self.frame,row=3,column=0)
        self.boList = []
        self.addBump()
        self.i = 0
        self.bl.toggleActive()
        # Pack
        self.controlFrame.grid(row=0,column=0)
        self.addButton.grid(row=0,column=0)
        self.removeButton.grid(row=0,column=1)
        self.typeFrame.grid(row=2,column=0)
        self.typeLabel.grid(row=2,column=0)
        self.typeB[0].grid(row=2,column=1)
        self.typeB[1].grid(row=2,column=2)
        self.addBump()
        tk.Button(self.frame,text="Get",command=lambda : print(self.get())).grid(row=8)
    

    def advance(self,event):
        if self.i==0:
            self.bl.toggleActive()
            self.i = (self.i+1)%(len(self.boList)+1)
            self.boList[self.i-1].toggleActive()
        else:
            self.boList[self.i-1].toggleActive()
            self.i = (self.i+1)%(len(self.boList)+1)
            if self.i==0:
                self.bl.toggleActive()
            else:
                self.boList[self.i-1].toggleActive()
#            
#        self.boList[self.i-1].toggleActive()
#        if event.dblclick:
#            self.i = (self.i-1)%(len(self.boList)+1)
#        else:
#            self.i = (self.i+1)%(len(self.boList)+1)
#        self.boList[self.i-1].toggleActive()
            
    def update(self,event):
        if self.i == 0:
            f = fitDic[self.typeVar.get()]
            bln = event.ydata
            blo = self.bl.get()[0]
            for bo in self.boList:
                res = list(bo.get())[:3]
                X = res[2]
                G = res[1]
                A = (res[0]*f(X,X,G)+blo-bln)/f(X,X,G)
                bo.set(A,G,X)
            self.bl.set(event.ydata)
        else:
            f = fitDic[self.typeVar.get()]
            res = list(self.boList[self.i-1].get())[:3]
            bl = self.bl.get()[0]
            if event.button == 1:
                X = event.xdata
                G = res[1]
                A = (event.ydata-bl)/(f(X,X,G))
            if event.button == 3:
                X = res[2]
                G = abs(X-event.xdata)*2
                A = res[0]*(f(X,X,res[1]))/(f(X,X,G))
            self.boList[self.i-1].set(A,G,X)

    
    def addBump(self):
        self.boList += [BumpObj(self.frame,label="{}:".format(1+len(self.boList)),Gs=self.Gs,row=4+len(self.boList),column=0)]
    def removeBump(self):
        if len(self.boList)>1:
            self.boList[-1].destroy()
            del self.boList[-1]
    
    def get(self):
        p = [self.bl.get()[0]]
        w = [self.bl.get()[1]]
        for bo in self.boList:
            res = bo.get()
            p += list(res[:3])
            w += list(res[3:])
        return self.bounds.get(),self.typeVar.get(),p,w
    def set(self,p):
        self.bl.set(p[0])
        for i in range(len(self.boList)):
            self.boList[i].set(*p[1+3*i:1+3+3*i])
        

class BaselineObj:
    def __init__(self,master,**kwargs):
        self.i = False
        self.frame = tk.Frame(master)
        self.frame.grid(**kwargs)
        # Setup
        self.label = tk.Label(self.frame,text="Baseline",justify=tk.RIGHT)
        self.val = tk.StringVar()
        self.val.set("0.0")
        self.entry = tk.Entry(self.frame,textvariable=self.val,
                              width=6)
        self.bool = tk.IntVar()
        self.bool.set(1)
        self.check = tk.Checkbutton(self.frame,variable=self.bool,
                                    onvalue=1,offvalue=0)
        # Pack
        self.label.grid(row=0,column=0)
        self.entry.grid(row=0,column=1)
        self.check.grid(row=0,column=2)
        
        
    def get(self):
        return float(self.val.get()), bool(self.bool.get())
    def set(self,val):#,boolean):
        self.val.set(str(val))
        #self.bool.set(int(boolean))
        
    def toggleActive(self):
        self.i = not self.i
        self.entry.config(bg=color[self.i])
        

class BoundObj:
    def __init__(self,master,xL=0.0,xR=1.0,**kwargs):
        self.frame = tk.Frame(master)
        self.frame.grid(**kwargs)
        self.label = tk.Label(self.frame,text="Fit Bounds ",justify=tk.RIGHT)
        self.labelxL = tk.Label(self.frame,text="Left:",justify=tk.RIGHT)
        self.labelxR = tk.Label(self.frame,text="Right:",justify=tk.RIGHT)
        self.xL = tk.StringVar()
        self.xR = tk.StringVar()
        self.xL.set(str(xL))
        self.xR.set(str(xR))
        self.entryxL = tk.Entry(self.frame,textvariable=self.xL,
                               width=6)
        self.entryxR = tk.Entry(self.frame,textvariable=self.xR,
                               width=6)
        # Pack
        self.label.grid(row=0,column=0)
        self.labelxL.grid(row=0,column=1)
        self.entryxL.grid(row=0,column=2)
        self.labelxR.grid(row=0,column=3)
        self.entryxR.grid(row=0,column=4)
        
    def get(self):
        return float(self.xL.get()),float(self.xR.get())
    def set(self,xL,xR):
        self.xL.set(str(xL))
        self.xR.set(str(xR))
    
class BumpObj:
    def __init__(self,master,label="",Gs=1,**kwargs):
        self.i = False
        self.frame = tk.Frame(master)
        self.frame.grid(**kwargs)
        self.label = tk.Label(self.frame,text=label,justify=tk.RIGHT)
        self.labelA = tk.Label(self.frame,text="A",justify=tk.RIGHT)
        self.labelG = tk.Label(self.frame,text="G",justify=tk.RIGHT)
        self.labelX = tk.Label(self.frame,text="X",justify=tk.RIGHT)
        self.A = tk.StringVar()
        self.G = tk.StringVar()
        self.X = tk.StringVar()
        self.A.set("0.0")
        self.G.set(str(Gs))
        self.X.set("0.0")
        self.entryA = tk.Entry(self.frame,textvariable=self.A,
                               width=6)
        self.entryG = tk.Entry(self.frame,textvariable=self.G,
                               width=6)
        self.entryX = tk.Entry(self.frame,textvariable=self.X,
                               width=6)
        self.boolA = tk.IntVar()
        self.boolG = tk.IntVar()
        self.boolX = tk.IntVar()
        self.boolA.set(1)
        self.boolG.set(1)
        self.boolX.set(1)
        self.checkA = tk.Checkbutton(self.frame,variable=self.boolA,
                                     onvalue=1,offvalue=0)
        self.checkG = tk.Checkbutton(self.frame,variable=self.boolG,
                                     onvalue=1,offvalue=0)
        self.checkX = tk.Checkbutton(self.frame,variable=self.boolX,
                                     onvalue=1,offvalue=0)
        # Pack
        self.label.grid(row=0,column=0)
        self.labelA.grid(row=0,column=1)
        self.entryA.grid(row=0,column=2)
        self.checkA.grid(row=0,column=3)
        self.labelG.grid(row=0,column=4)
        self.entryG.grid(row=0,column=5)
        self.checkG.grid(row=0,column=6)
        self.labelX.grid(row=0,column=7)
        self.entryX.grid(row=0,column=8)
        self.checkX.grid(row=0,column=9)
        
        
    def get(self):
        return float(self.A.get()),float(self.G.get()),float(self.X.get()),\
                bool(self.boolA.get()),bool(self.boolG.get()),bool(self.boolX.get())
    def set(self,A,G,X):#,boolA,boolG,boolX):
        self.A.set(str(A))
        self.G.set(str(G))
        self.X.set(str(X))
        #self.boolA.set(int(boolA))
        #self.boolG.set(int(boolG))
        #self.boolX.set(int(boolX))
        
    def destroy(self):
        listW = self.frame.grid_slaves()
        for l in listW:
            l.destroy()
        self.frame.destroy()
    
    def toggleActive(self):
        self.i = not self.i
        self.entryA.config(bg=color[self.i])
        self.entryG.config(bg=color[self.i])
        self.entryX.config(bg=color[self.i])
        