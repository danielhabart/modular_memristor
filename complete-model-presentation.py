
import numpy as np
from matplotlib import pyplot as pt



# ----------------------------------------------------
# 
# Modular memristor model with synaptic plasticity and volatile memory
# INPUT:
#       time            time xaxis for simulation
#       VA              external voltage applied on the device
#       IniState        Initial state of the device
#       param           device parameters
#                       [kpos apos tau beta D vth vhold Ron Roff]
#       ic              [A] compliance current                    
#                   
# OUTPUT
#       v               Applied Voltage
#       i               The current through the device
#       x               Internal State Variable
# ---------------------------------------------------- 


readV = 0.05
Pw, Dt, Tr = 10*1e-3, 20*1e-3,50*1e-3
datastep= Pw+2*Dt+Tr
cyk = 200*2
V0 = 0.4


def pos(x):
    if x >= 0:
        return x
    else: return 0

def stdp_step( x,y,w, v, tp,tn, ap,an): # sign interpretation: + is depression cycle, - is potentiation cycle

    dx = step*(-x/tp) + step * pos(v) # positive part of voltage represents presynaptic activity
    dy = step*(-y/tn) + step * pos(-v) # negative part of voltage represents postsynaptic activity
    dw = step*( -an*y*pos(v) + ap*x*pos(-v) ) 

    return dx,dy,dw

#all units in basic form

a = -0.58548204
eps0 = 1e-6
h11,h10= 10.78560073e-5, 224.48436321e-4
mu, g0 = 0.21597481761947285, 2.26e-03 

t1,t2,a1,a2,h01,h00 = 0.176843493508, 0.619508386108, 10.962502436, 0.90766780177, 0.000101909990799, 0.00252548081238

def dH(s, w, h11,h10, h01,h00):
    return h11*np.power(abs(s)+eps0, a) + h10,   h01*np.power(abs(w)+eps0,a) +h00 #*np.exp(h0*s)

def evolution(v):
    s = s0
    x,y,w =  0,0, w0
    dgw = np.zeros(int(Tmax/step))

    for n,t in enumerate(range( int(Tmax/step))): #in s

        ds = step*mu*v[n]
        s += ds

        dx,dy,dw = stdp_step( x,y,w, v[n], t1,t2,a1,a2)
        x += dx
        y += dy
        w += dw

        dgw[t] = sum( np.array(dH(s,w,h11,h10,h01,h00)) * np.array([ds/step,dw/step]) )#*G(s)
        

    return  dgw

w0 = 0#1*D/4
s0 = 0


def conv(hs,kernel):
    z = np.convolve(hs,kernel) #[int(Tmax/step):]
    u = z[int((Pw+Dt)/step)::int(datastep/step)]
    return u


# start_vals = (1.3572192341069214e-05, 1.07103067988524e-06, 1.082, 0.41)

#in msec
step = 1*10**(-3)

Tmax = datalen*datastep


time = int(2*Tmax/step)


def cycle(t):
    if t< Pw:
        return V0
    else: return 0

def inputV(t,cycles):
    if (t*step)% (2*cycles * datastep) < cycles * datastep:
        return cycle( (t*step)%datastep )

    else:return -cycle( (t*step)%datastep )



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

def nonlin(convolved,g0,gmax,a):
    return gmax*sigmoid(a*convolved) + g0



Ik = integral(ker)


k = np.array([ker((t*step))/Ik for t in range(time)])

vx = np.array([inputV(t,200) for t in range(time)])
tt = np.array([ t*step for t in range(time)])
tx =  np.array([t*datastep for t in range(int(Tmax/datastep)) ])



pd = []  


fig, ax= pt.subplots(1,1,constrained_layout=True)
dark_gray =  '#282A2B'#'#444444'
for spine in ax.spines.values():
    spine.set_color(dark_gray)
ax.tick_params(colors= dark_gray)

ng.format_spines(ax)

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

ax.plot(  x[cyk:-cyk], y[cyk:-cyk] , marker = '.', linestyle= '',color =ng.colors['b'], label='0.4V')
out = conv(evolution(vx),k)[:len(y)]
ax.plot(tx[cyk:-cyk],out[cyk:-cyk]+g0,marker = '.')
# ax.plot(  tt, dgds[:len(tt)]+g0, marker = '.',color =ng.colors['o'] )
# ax.plot(  tx, conv(evolution(vx[indx]),k)[:len(y[indx])] + g0, marker = '.',color =ng.colors['o'], label=r'$\ker \ast dH $' )
# ax.plot(  tx, conv(dgds,k)[:len(tx)]+g0*conv(gs,k)[:len(tx)], marker = '.',color =ng.colors['g'] )
# ax.plot(  tt, gs[:len(tt)] , marker = '+',color =ng.colors['b'] )
# ax.plot(  tt, np.convolve(gs,dk)[:len(tt)] , marker = '+',color =ng.colors['r'] )
# ax.plot(  tt,  s, marker = '',color =ng.colors['g'] )
# ax.plot(  tt, conv(gs))
ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 0),
                    useMathText=True, useOffset=False)
# pd += [u]

pt.legend(loc = 'upper right')
pt.show()


