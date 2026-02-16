
import numpy as np
from matplotlib import pyplot as pt
import scipy.optimize as opt

import csv

data = 'tests/data_T=35s.csv'
inputpulses= 'tests/inputpulses_T=35s.csv'

datalen =580

readV = 0.05

V0=0.5

def cycle(t):
    if t< 20:
        return V0
    else: return 0

def box(t):
    if t>4000 and t<5000:
        return V0
    else: return 0


def dH(s,h1,h0):
    return h1*np.power(abs(s), a) + h0   #*np.exp(h0*s)


def sim(v,mu,h1,h0):
    dHds = np.zeros(len(v))
    
    s = np.cumsum(v)*step*mu
    dHds = dH(s,h1,h0) * mu*v

    return dHds


def conv(hs,kernel):
    z = np.convolve(hs,kernel) #[int(Tmax/step):]
    u = z[int(0.040/step)::int(datastep/step)]
    return u


# start_vals = (1.3572192341069214e-05, 1.07103067988524e-06, 1.082, 0.41)

#in msec
datastep = 70*10**(-3)
step = 1*10**(-3)

Dt, Tr = 20,50
period = (20 + 2*Dt +Tr)*10**(-3)

Tmax = datalen*datastep


time = int(2*Tmax/step)



def inputV(t,cycles):
    if t*step <= cycles* datastep:
        return cycle( (t*step)%datastep )
    else:return 0



alph,eps = 0.02939634033310315 , 0.07
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

# def hardclip(g):
#     g = np.asarray(g)
#     if g>gmax: return gmax
#     if g<gmin: return gmin
#     else:return g


Ik = integral(ker)


k = np.array([ker((t*step))/Ik for t in range(time)])

vx = np.array([[inputV(t,float(name)) for t in range(time)] for name in names])
tt = np.array([ t*step for t in range(time)])
tx =  np.array([t*datastep for t in range(int(2*Tmax/datastep)) ])


def single_err(i,mu,g0,h1,h0):
    gres = conv(sim(vx[i],mu,h1,h0),k)[:len(y[i])] + g0
    err = (gres-y[i])**2 /np.absolute(y[i])
    return sum(err)

def single_f(params):
    mu,g0 = params

    gres = conv(sim(vx[indx],mu,h1,h0),k)[:len(y[indx])] + g0
    err = (gres-y[indx])**2 /np.absolute(y[indx])
    return sum(err)

def error_f(params):

    mu,h1,h0,g0 = params

    args = [(i,mu,g0,h1,h0) for i in range(len(names))]
    res = 0.
    for a in args:
        res += single_err(*a)
    # res = pool.starmap(single_err,args )
    # res = single_err(1,mu,Rd,R0)
    return res


indx = 2

a = -0.58548204
h1,h0= 2.08560073e-6, 56.48436321e-6 #35s


# mu =0.226
bds= [(0.12,0.3)] + [(2*10**(-5),4*10**(-5))] 

# mu, g0 = opt.differential_evolution(single_f, bounds = bds).x
# mu, g0 = 0.22651294795812438, 3.073247519189803e-05 #3.5s:
# print(mu,g0)
mu, g0 = 0.21597481761947285, 3.0889309655883944e-05 # 35s :
gmax = 5e-5
gmin = 2*g0-gmax
slo = 4/(gmax-gmin)

print(gmax, gmin)


fig, ax= pt.subplots(1,1,constrained_layout=True)
dark_gray =  '#282A2B'#'#444444'
for spine in ax.spines.values():
    spine.set_color(dark_gray)
ax.tick_params(colors= dark_gray)


# ax.set_title(f'input pulses, read after 20ms, {name}', y=1.1, pad=-14)
ax.set_xlabel('Time [s]')
ax.set_ylabel('Conductance [S]')

# ax.yaxis.set_label_coords(0.2,0.91)
# ax.xaxis.set_label_coords(0.91,0.19)

# i = 0
# for c in k:
#     i+=c

# print(i)
def no_zeros(x, n):
    factor = 10 ** n
    x = int(x * factor) / factor              # truncate
    s = ('%f' % x).rstrip('0').rstrip('.')    # remove trailing zeros
    return s

# dgds = sim(vx[indx], mu,h1,h0)

ax.plot(  tx, y[indx] , marker = '.', linestyle= '',color ='b', label='T = ' + no_zeros(float(names[indx])*datastep,2) + r'$\,$s')
ax.plot(  tx,  conv(sim(vx[indx],mu,h1,h0),k)[:len(y[indx])] +g0, marker = '.',color ='r', label=r'$\ker \ast dH $' )
ax.plot(  tx, nonlin( conv(sim(vx[indx],mu,h1,h0),k)[:len(y[indx])] ),color ='k', label=r'sat$(\ker \ast dH) $' )


pt.legend(loc = 'upper right')
pt.show()


