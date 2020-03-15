import numpy as np
import sympy as sp

from mpl_toolkits.mplot3d import axes3d 
import matplotlib.pyplot as plt 
from matplotlib import rcParams 
rcParams['savefig.transparent'] = True 
rcParams['text.usetex'] =True

from numpy import cos, sin
from numpy import pi as PI
TAU = 2 * PI
PI2 = PI / 2

# used symbols
a, b, u = sp.symbols('a b u', real=True)
ww = sp.symbols('w', positive=True)
t, t0 = sp.symbols('t t_0', positive=True)

def plot_initial_trajectory(u=[1,-1], a=1, b=0, N=1000, ax=plt.gca(projection='3d'), **kwargs):
    ucolor = ['blue', 'red', 'orange']
    for i in range(len(u)):
        w2 = a**2 + b**2 + u[i]**2
        w = np.sqrt(w2)
        Tmax = TAU / w
        t = np.linspace(0, Tmax, N)
        x = (a*u[i]*cos(t*w) - a*u[i] + b*w*sin(t*w)) / w2
        y = (a*b*cos(t*w) - a*b - u[i]*w*sin(t*w)) / w2
        z = (-a**2*cos(t*w) + a**2 + w2*cos(t*w)) / w2
        ax.plot3D(x, y, z, color=ucolor[i], label=f'$u={u[i]}$', **kwargs)


def plot_vector_field(M, R=1, a_val=1, b_val=0, u_val=[1, -1], N=20,
        ax=plt.gca(projection='3d'), plot_sphere=False, **kwargs):
    length = kwargs.pop("length", 0.1)
    ucolor = ['blue', 'red', 'orange']

    # sphere mesh
    N *= 1j
    theta, phi = np.mgrid[0:TAU:N, -PI2:PI2:N]
    X = R * np.cos(theta) * np.cos(phi)
    Y = R * np.sin(theta) * np.cos(phi)
    Z = R * np.sin(phi)
    
    # plot sphere
    if plot_sphere:
        sr = 0.9
        ax.plot_surface(sr*X, sr*Y, sr*Z, color=(0.5, 0.5, 0.5))

    # evaluate vector field
    M = M.subs({a:a_val, b:b_val})
    for ind in range(len(u_val)):
        M_val = M.subs({u:u_val[ind]})
        M_val = np.array(M_val.tolist(), dtype=float)
        U = M_val[0, 0] * X + M_val[0, 1] * Y + M_val[0, 2] * Y
        V = M_val[1, 0] * X + M_val[1, 1] * Y + M_val[1, 2] * Y
        W = M_val[2, 0] * X + M_val[2, 1] * Y + M_val[2, 2] * Y
        ax.quiver(X, Y, Z, U, V, W, normalize=True, length=length,
                  color=ucolor[ind], **kwargs)

def plot_trajectory(Xt1, Xt2, u, x0, x1, ax=plt.gca(projection='3d')):
    if x0 is not None:
        x0 = np.array(list(x0), dtype=float)
        ax.scatter([x0[0]], [x0[1]], [x0[2]], color='red', linewidths=5, label='$x(0)$', zorder=99999)
    if x1 is not None:
        x1 = np.array(list(x1), dtype=float)
        ax.scatter([x1[0]], [x1[1]], [x1[2]], color='black', linewidths=5, label='$x(T)$', zorder=99999)
    
    ax.plot3D(Xt1[0], Xt1[1], Xt1[2], color='blue', label=f'$u={u[0]}$')
    ax.plot3D(Xt2[0], Xt2[1], Xt2[2], color='red', label=f'$u={u[1]}$')
    ax.legend()