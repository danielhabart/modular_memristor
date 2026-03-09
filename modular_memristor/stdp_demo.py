import numpy as np
from matplotlib import pyplot as pt
from matplotlib.legend_handler import HandlerTuple
from multiprocessing import Pool
from pathlib import Path
import json


file_path = Path(__file__).resolve().parent / 'stdp_data.json'

with open(file_path,'r') as f:
    dtimesR, dtimesL, dataR, dataL, dtimes_longR ,dtimes_longL ,data_longR, data_longL = map(np.array, json.load(f))

dtimes = np.concatenate( (dtimesR, dtimesL) )
data = np.concatenate ( (dataR, dataL) )
dtimes_long = np.concatenate((dtimes_longR, dtimes_longL))
data_long = np.concatenate((data_longR, data_longL))


length = len(dtimesR)


# supporting functions

def ker(x):
    if x>0: return 1/np.power( eps+x,alph+1)
    else: return 0

def intker():
    x = 0.
    for t in range(int(Tmax/step)):
        x+= ker(t*step)
    return x

def pulse(t):
    if t >= 0 and t < dur*1e-3 :
        return 1
    else:
        return 0

def sgn(x):
    if x>=0: return 1
    if x<0: return -1

def pos(x):
    if x >= 0:
        return x
    else: return 0

def drel(d,rel):
    return (d-rel)/rel




#all units in basic form

alph,eps = 0.02939634033310315, 0.070

# H(s) parameters 
a = -0.58548204
eps0 = 1e-6
h11,h10= 20.78560073e-5, 184.48436321e-4
mu, g0 = 0.21597481761947285, 2.26e-03 
 
# H(w) parameters = STDP
t1,t2,a1,a2,h01,h00 = 0.176843493508, 0.619508386108, 10.962502436, 0.90766780177, 0.000101909990799, 0.00252548081238


#simulation parameters
step = 0.1e-3 
dur = 40
sep = 0.1

V0 = 0.5



#model functions, simulation evolution, dH and convolution

def conv(delt,gw,k):
    time10, time100  = (abs(delt)+10 + 2*dur)*1e-3, (abs(delt)+100 + 2*dur)*1e-3
    z = np.convolve(gw,k)[:int(Tmax/step)] 
    gw10, gw100 = z[int(time10/step)], z[int(time100/step)-1]
    return gw10 , gw100


def dH(s, w, h11,h10, h01,h00):
    return h11*np.power(abs(s)+eps0, a) + h10,   h01*np.power(abs(w)+eps0,a) +h00 #*np.exp(h0*s)

def stdp_step( x,y,w, v, tp,tn, ap,an): # sign: + is depression cycle, - is potentiation cycle

    dx = step*(-x/tp) + step * pos(v) # positive part of voltage represents activity of presynaptic
    dy = step*(-y/tn) + step * pos(-v) # negative part of volatage represents activity of postsynaptic
    dw = step*( -an*y*pos(v) + ap*x*pos(-v) ) 

    return dx,dy,dw

def stdp_ev( delt,t1,t2,a1,a2,h01,h00):

    x,y,w =  0, 0, 0 
    dgw = np.zeros(int(Tmax/step))

    if delt< 0:
        tneg, tplus = 0, (abs(delt)+dur)*1e-3
    if delt>0:
        tplus, tneg = 0, (abs(delt)+dur)*1e-3

    for t in range( int(Tmax/step)): #in s
        v = V0*pulse(t*step-tplus) - V0*pulse(t*step-tneg)

        # ds = step*mu*v
        # s += ds

        dx,dy,dw = stdp_step( x,y,w, v, t1,t2,a1,a2)
        x += dx
        y += dy
        w += dw

        dgw[t] = sum( np.array(dH(0,w,h11,h10,h01,h00)) * np.array([0/step,dw/step]) )#*G(s)
        
    gw0 = g0 
    out = dgw 
    gwGs10,gwGs100 = conv(delt,out,k)

    return  drel(gwGs10+gw0,gw0), drel(gwGs100+gw0,gw0)


def stdp_min(pa):
    # t1,t2,
    t1,t2,a1,a2,h01,h00 = pa
    dat = [(dt,t1,t2,a1,a2,h01,h00) for dt in dtimes ]

    pres10,pres100 = map(np.array, zip(*pool.starmap (stdp_ev, dat)))

    return pres10 ,pres100




Tmax = (max(np.absolute(dtimes))+100 + 2*dur)*1e-3

#normalization of kernel
Ik = intker()
k = np.array([ker((t*step))/Ik for t in range( int(Tmax/step)) ])#, (int(Tmax/step),0), mode='constant')
I_k = [ step*sum(k[:t]) for t in range(len(k)) ] 



tt = np.array([t*step for t in range( int(Tmax/step))])



if __name__ == '__main__':
    with Pool() as pool:
        
        fig, ax = pt.subplots(1,1,constrained_layout=True)
        p10,p100 = stdp_min( (t1,t2,a1,a2,h01,h00) ) 


        pm1, =ax.plot(dtimes[:length],p10[:length]*100,color='k')
        pm2, =ax.plot(dtimes[length:],p10[length:]*100 ,color='k')

        pm3, =ax.plot(dtimes[:length],p100[:length]*100,color='r')
        pm4, =ax.plot(dtimes[length:],p100[length:]*100 ,color='r')

        pd = []

        pd += [*ax.plot(dtimesR, dataR*100, linestyle='',marker ='.', color = 'b')]
        pd  += [*ax.plot(dtimes_longR, data_longR*100, linestyle='',marker ='.', color = 'r')]

        pd  += [*ax.plot(dtimesL, dataL*100, linestyle='',marker ='.', color = 'g') ]
        pd  += [*ax.plot(dtimes_longL,data_longL*100, linestyle='',marker ='.', color = 'y') ]
       

        handles = [(pd[0],pd[2]), (pd[1],pd[3]), (pm1,pm2),(pm3,pm4)]
        labels =[r'$10\,$ms',r'$100\,$ms', r'$ \frac {\Delta G(t_1)}{G_0}$', r'$ \frac {\Delta G(t_2)}{G_0}$']
        l = ax.legend(handles,labels,
               handler_map={tuple: HandlerTuple(ndivide=None)},loc='upper left')

        for text in l.get_texts(): #change font sizes
            if text.get_text() == r'$ \frac {\Delta G(t_1)}{G_0}$' or text.get_text() ==r'$ \frac {\Delta G(t_2)}{G_0}$':

                text.set_fontsize(14)  
            else:
                text.set_fontsize('11')  


        # ng.format_spines(ax)
        ax.axhline(0, color='black',linestyle=':', linewidth=1)
        ax.axvline(0, color='black',linestyle=':', linewidth=1)
        

        ax.set_xlabel('$\Delta t$ [ms]')
        ax.set_ylabel(r'$\frac{\Delta G}{G}$[%]',rotation=0,fontsize=12)

        pt.show()