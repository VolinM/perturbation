import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

hbar = 1
m = 1

a = -10
b = 10
N = 301

phis = np.load("phis.npy",allow_pickle=True)
xs = phis.item().get("xs")
ts = np.linspace(0,1,1001)

Nx = xs.shape[0]-1
Nt = ts.shape[0]-1

dx = xs[1]-xs[0]
dt = ts[1]-ts[0]

def potential(x):
    return 0
    if x>(a+2*dx) and x<(b-2*dx):
        return 0
    else:
        return 100000000
        # a = 1/10
        # return 1/(a*np.sqrt(np.pi))*np.exp(-(x/a)**2)

V = np.zeros_like(xs)
for i,row in enumerate(V):
    V[i] = potential(xs[i])

psi = pd.DataFrame(index=xs,columns=ts)
phi1 =(phis.item().get(2))
phi1 = np.append(phi1,0)
phi1 = np.insert(phi1,0,0)
phi1 = phi1/np.sqrt(np.sum(phi1**2))

# print(np.sum(phi1**2))

# phi1 = np.sqrt(1/10)*np.cos(3*np.pi*xs/20)
# phi1 = phi1/np.sqrt(np.sum(phi1**2))

# print(np.sum((phi1**2)[:],axis=0))

# plt.plot(xs,phi1)

# plt.plot(xs,phi2)
# plt.savefig("a.png")




def aniFunc(i):
    plt.clf()
    a = psi.iloc[:,i].to_numpy()
    b = (a*np.conj(a)).astype(complex).real
    plt.plot(xs,b)
    plt.ylabel(r'$|\Psi(x,t)|^2$')
    plt.xlabel(r'$x$')
    plt.title(r't = {0:.3f}'.format(np.round(ts[i],3)))
    plt.grid()


figure, ax1 = plt.subplots()

anim = animation.FuncAnimation(figure, func=aniFunc,frames=Nt,interval=1)
anim.save('anim.gif',writer='imagemagick', fps=2000)

# # # # anim = animation.FuncAnimation(figure, func=time_evol, frames=ts)
# # # # anim.save('anim.mp4',writer =animation.FFMpegWriter(fps=10))
# # # # print(V)
# # # # print(np.real(psi.iloc[:,Nt-1].to_numpy[0]))


