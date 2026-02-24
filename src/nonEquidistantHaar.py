# This function is taken from: https://github.com/NicoAcunaR/Non-equidistant_Haar_fluctuations/blob/main/Haar-RMS.ipynb
import pandas as pd
import numpy as np
import time
import math
from csaps import csaps

def timeSince(since):
    now = time.time()
    s = now - since
    m = math.floor(s/60)
    s -= m*60
    if m < 60:
        return '%dm %ds' % (m, s)
    else:
        h = math.floor(m/60)
        m -= h*60
        return '%dh %dm %ds' % (h, m, s)

class Haar:
    def __init__(self,t,x):
        self.t = np.array(t)
        self.x = np.array(x)

        self.deltas_t = []
        self.epsilons = []

        self.ep_min = None
        self.ep_max = None
        self.calib = None
        self.Hs = None

    def compute_deltas(self):
        for i in range(1,int(len(self.t)/2)+1):  
            self.deltas_t.append(self.t[i:]-self.t[:-i]) 

    def compute_epsilons(self):
        for i in range(1,len(self.deltas_t)+1):
            self.epsilons.append(self.deltas_t[i-1][:-i]/(self.deltas_t[i-1][:-i]+self.deltas_t[i-1][i:]))

    def fluctuations(self,ep_min=0.45,calib=2,verbose=False,prop_print=500):
        """
        Computing fluctuations

        Arguments:
            ep_min: float
                default 0.25

            calib: int
                default 2

            prop_print: int
                how often do we print the progress, default 500
        """
        self.x = self.x[:-1]
        self.t = self.t[:-1]
        self.ep_min = ep_min
        self.ep_max = 1 - self.ep_min
        self.calib = calib

        self.Hs = []
        self.delta_t = [] 
        counter = 0
        
        
        start_time = time.time()

        for H in range(2,len(self.x)-1,2):
            for start in range(len(self.x)-H-1):
                int1 = np.sum(self.x[start:start+int(H/2)]*self.deltas_t[0][start:start+int(H/2)]/self.deltas_t[int(H/2)-1][start])
                int2 = np.sum(self.x[start+int(H/2):start+H]*self.deltas_t[0][start+int(H/2):start+H]/self.deltas_t[int(H/2)-1][start+int(H/2)])
                counter += 1

                if self.epsilon_range(H,start):
                    self.Hs.append((calib*(int2 - int1))**2)  
                    self.delta_t.append(self.deltas_t[int(H/2)-1][start] + self.deltas_t[int(H/2)-1][start+int(H/2)])  
                
            if verbose:
                prop = 100*H / (len(self.x) - 1)
                if H % prop_print  == 0:
                    print("Progress: {}%, time elapsed {}".format('%.3f'%(prop),timeSince(start_time)))

        if verbose:
            print("Finished computations in {}".format(timeSince(start_time)))
            perct = (counter - len(self.Hs))/counter*100
            perct = '%.3f'%(perct)
            print("{} fluctuations removed ({}%)".format(counter - len(self.Hs),perct))
            
            
    def epsilon_range(self,H,start):
        min_condition = self.ep_min  < self.epsilons[int(H/2)-1][start]
        max_condition = self.epsilons[int(H/2)-1][start] < self.ep_max
        return min_condition and max_condition

    @property
    def data_df(self):
        if self.Hs is None:
            raise ValueError("Hs not yet defined")
        df = pd.DataFrame(data={'delta t':self.delta_t , 'Hs': self.Hs})
        return df.sort_values('delta t',axis=0).reset_index(drop=True)
    
    @property
    def csap(self):
        spiky_t=pd.DataFrame(np.log10(np.round(self.data_df['delta t'],3)), columns = ['delta t'])
        spiky_v=pd.DataFrame(self.data_df['Hs'],columns=['Hs'])
        spiky=pd.concat([spiky_t,spiky_v],axis=1)
        spiky=spiky.groupby("delta t", as_index=False).mean()
        spiky['Hs']=np.log10(np.sqrt(spiky['Hs']))
        spiky=spiky.replace([np.inf, -np.inf], np.nan).dropna()
        xi = np.linspace(min(spiky['delta t']), max(spiky['delta t']), len(spiky['delta t']))
        yi = csaps(spiky['delta t'],spiky['Hs'], xi, smooth=0.9999)
        return xi,yi

    

def Haar_RMS_fluctuations(t,x):
    H=Haar(t,x)
    H.compute_deltas()
    H.compute_epsilons()
    H.fluctuations()
    return H.csap[0],H.csap[1]