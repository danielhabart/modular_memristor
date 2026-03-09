
import numpy as np
from matplotlib import pyplot as pt
import nice_graphs as ng
import json
from pathlib import Path

file_path = Path(__file__).resolve().parent / 'pulses_data.json'

with open(file_path,'r') as f:
    x,y = map(np.array, json.load(f))


def cycle(t):
    if t< 20:
        return V0
    else: return 0

def box(t):
    if t>4000 and t<5000:
        return V0
    else: return 0

def no_zeros(x, n):
    factor = 10 ** n
    x = int(x * factor) / factor             
    s = ('%f' % x).rstrip('0').rstrip('.')    
    return s




def dH(s,h1,h0):
    return h1*np.power(abs(s), a) + h0   #*np.exp(h0*s)

#due to simplicity of differential equation, evolution may be calculated simply by cumulative sum.
def sim(v,mu,h1,h0):
    dHds = np.zeros(len(v))
    
    s = np.cumsum(v)*step*mu
    dHds = dH(s,h1,h0) * mu*v

    return dHds


def conv(hs,kernel):
    z = np.convolve(hs,kernel) #[int(Tmax/step):]
    u = z[int(0.040/step)::int(datastep/step)]
    return u


#simulation parameters, in msec
datastep = 70*10**(-3)
step = 1*10**(-3)

Dt, Tr = 20,50
period = (20 + 2*Dt +Tr)*10**(-3)

V0 =    0.5
datalen = len(x)//2

Tmax = datalen*datastep
time = int(2*Tmax/step)

#model parameters
a = -0.58548204
h1,h0= 2.08560073e-6, 56.48436321e-6 #35s
mu, g0 = 0.21597481761947285, 3.0889309655883944e-05 # 35s :
gmax = 5e-5
gmin = 2*g0-gmax
slo = 4/(gmax-gmin)
alph,eps = 0.02939634033310315 , 0.07

def inputV(t,cycles):
    if t*step <= cycles* datastep:
        return cycle( (t*step)%datastep )
    else:return 0

def ker(x):
    if x>0:return 1/np.power(x+eps,alph+1)
    else: return 0

def integral(f):
    x = 0.
    for t in range(int(time/2)):
        x+= f(t*step)
    return x

def sigmoid(x):
    x = np.asarray(x)  # allows input to be list, scalar, or array
    x = np.clip(x, -500, 500)
    return 1/(1+np.exp(-x))

def nonlin(g):
    # return np.clip(g,gmin,gmax)
    return (gmax-gmin)*(sigmoid(slo*g)-1/2) +g0


Ik = integral(ker)
k = np.array([ker((t*step))/Ik for t in range(time)])
vx = np.array([inputV(t,float(500)) for t in range(time)])
tt = np.array([ t*step for t in range(time)])
tx =  np.array([t*datastep for t in range(int(2*Tmax/datastep)) ])



#graphing
pd = []  

fig, ax= pt.subplots(1,1,constrained_layout=True)
dark_gray =  '#282A2B'#'#444444'
for spine in ax.spines.values():
    spine.set_color(dark_gray)
ax.tick_params(colors= dark_gray)


ax.set_xlabel('Time [s]')
ax.set_ylabel('Conductance [S]')

ax.plot(  tx, y , marker = '.', linestyle= '',color = 'b', label='T = ' + no_zeros(float(500)*datastep,2) + r'$\,$s')

ax.plot(  tx,  conv(sim(vx,mu,h1,h0),k)[:len(y)] +g0, marker = '.',color ='tab:orange', label=r'$\ker \ast dH $' )

ax.plot(  tx, nonlin( conv(sim(vx,mu,h1,h0),k)[:len(y)] ),color ='k', label=r'sat$(\ker \ast dH) $' )

pt.legend(loc = 'upper right')
pt.show()


