import numpy as np 
from numpy.linalg import norm
from sympy.utilities.lambdify import lambdify
import sympy as sp 
from matplotlib import pyplot as plt
from visuals import *

def numpify_x(exp_Mt, x0, a_val=1, b_val=0, u_val=1):
    # substitute parameters with values
    w_val = np.sqrt(a_val**2 + b_val**2 + u_val**2)
    exp = exp_Mt.subs({ww:w_val,a:a_val, b:b_val, u:u_val})
    if isinstance(x0, np.ndarray):
        X, Y, Z = exp * x0.reshape((3,1))
    else:
        X, Y, Z = exp * x0

    # create numpy functions
    x = lambdify(t, X, 'numpy')
    y = lambdify(t, Y, 'numpy')
    z = lambdify(t, Z, 'numpy')
    return x, y, z

def is_in_trajectory(exp_Mt, start, end, u_val=1, a_val=1, b_val=0, N=1000, epsilon=0.1):
    w_val = np.sqrt(a_val**2 + b_val**2 + u_val**2)
    Tmax = TAU / w_val
    time = np.linspace(0, Tmax, N)
    x, y, z = numpify_x(exp_Mt, start, a_val=a_val, b_val=b_val, u_val=u_val)
    X = np.array([x(time), y(time), z(time)])
    for i in range(N):
        if norm(end - X[:,i]) < epsilon:
            return time[i], X[:,:(i+1)]
    return False, None

def find_optimal(exp_Mt, X0, X1, a_val=1, b_val=0, u_first=1, N=1000, epsilon=0.01):
    w_val = np.sqrt(a_val**2 + b_val**2 + u_first**2)
    Tmax = TAU / w_val
    time = np.linspace(0, Tmax, N)
    x1 = np.array(list(X1), dtype=float)

    x, y, z = numpify_x(exp_Mt, X0, a_val=a_val, b_val=b_val, u_val=u_first)
    X = np.array([x(time), y(time), z(time)])
    for i in range(N):
        Topt, Xt2 = is_in_trajectory(exp_Mt, X[:,i], x1, u_val=-u_first,
                a_val=a_val, b_val=b_val, N=N, epsilon=epsilon)
        if Topt > 0:
            Topt += time[i]
            print("Optimal trajectory found!")
            print(f"Commutation time: {time[i]}")
            print(f"Commutation position: {X[:,i]}")
            print(f"Optimal time: {Topt}")

            Xt1 = X[:,:(i+1)] # first part of trajectory
            return Xt1, Xt2, Topt

    print("Optimal trajectory not found, try more points or consider switching fields...like chemistry or something easy")
    return False

def plot_optimal(exp_Mt, X0, X1, M, u=[1, -1], a=1, b=0, nb_points=1000, epsilon=0.01, figsize=(12, 12)):
    print(f'Finding optimal solution for u={u[0]} then u={u[0]} with a={a} b={b}')
    Xt1, Xt2, Topt = find_optimal(exp_Mt, X0, X1, a_val=a, b_val=b, u_first=u[0], epsilon=epsilon)

    fig = plt.figure(figsize=figsize)
    ax = fig.gca(projection='3d')
    plot_vector_field(M, R=float(X0.norm()), a_val=a, b_val=b, u_val=u, ax=ax, alpha=0.1)
    plot_trajectory(Xt1, Xt2, u=u, x0=X0, x1=X1, ax=ax)