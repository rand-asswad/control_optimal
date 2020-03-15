import numpy as np 
from numpy.linalg import norm
from numpy import pi as PI
import sympy as sp 
from sympy.utilities.lambdify import lambdify
from matplotlib import pyplot as plt
from visuals import plot_vector_field, plot_trajectory

TAU = 2 * PI
PI2 = PI / 2


print('Initializing system...')
# initialize matrix symbols 
a, b, u = sp.symbols('a b u', real=True)
w = sp.sqrt(a**2 + b**2 + u**2)
ww = sp.symbols('w', positive=True)

# system matrices
A = sp.Matrix([[0, a, b], [-a, 0, 0], [-b, 0, 0]])
B = sp.Matrix([[0, 0, 0], [ 0, 0,-1], [ 0, 1, 0]])
M = A + u*B
print('M')
sp.pprint(M)

# inverse of (sI - M) step by step
print("*"*180)
print('Calculating inverse matrix...')
s = sp.symbols('s') 
Ms = s * sp.eye(3) - M  # sI - M
Ms_det = sp.factor(Ms.det()).subs({w**2:ww**2})
Ms_com = Ms.cofactor_matrix()
Ms_inv = Ms_com.transpose() / Ms_det
print('(sI-M)^-1 =')

# inverse laplace transform
print("*"*180)
print('Calculating inverse laplace transform...')
t, t0 = sp.symbols('t t_0', positive=True)
exp_Mt = sp.inverse_laplace_transform(Ms_inv, s, t)
print('exp(Mt) = L-1((sI-1)^-1)(t) =')
sp.pprint(exp_Mt)

# initial conditions
X0 = sp.Matrix([0, 0, 1])
X1 = sp.Matrix([0, 1/sp.sqrt(2), 1/sp.sqrt(2)])
print("*"*180)
print('X = exp(Mt) * X0 =')
sp.pprint(exp_Mt * X0)

# define trajectory functions
def numpify_x(x0, a_val=1, b_val=0, u_val=1):
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

# check if point is in trajectory
def is_in_trajectory(start, end, u_val=1, a_val=1, b_val=0, N=1000, epsilon=0.1):
    w_val = float(w.subs({a:a_val, b:b_val, u:u_val}))
    Tmax = TAU / w_val
    time = np.linspace(0, Tmax, N)
    x, y, z = numpify_x(start, a_val=a_val, b_val=b_val, u_val=u_val)
    X = np.array([x(time), y(time), z(time)])
    for i in range(N):
        if norm(end - X[:,i]) < epsilon:
            return time[i], X[:,:(i+1)]
    return False, None

def find_optimal(a_val=1, b_val=0, u_first=1, N=1000, epsilon=0.01):
    w_val = float(w.subs({a:a_val, b:b_val, u:u_first}))
    Tmax = TAU / w_val
    time = np.linspace(0, Tmax, N)
    x1 = np.array(list(X1), dtype=float)

    x, y, z = numpify_x(X0, a_val=a_val, b_val=b_val, u_val=u_first)
    X = np.array([x(time), y(time), z(time)])
    for i in range(N):
        Topt, Xt2 = is_in_trajectory(X[:,i], x1, u_val=-u_first,
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

def plot_optimal(u=[1, -1], a=1, b=0, nb_points=1000, epsilon=0.01, figsize=(12, 12)):
    print("*"*180)
    print(f'Finding optimal solution for u={u[0]} then u={u[1]} with a={a} b={b}')
    Xt1, Xt2, Topt = find_optimal(a_val=a, b_val=b, u_first=u[0], epsilon=epsilon)

    fig = plt.figure(figsize=figsize)
    ax = fig.gca(projection='3d')
    plot_vector_field(M, R=float(X0.norm()), a_val=a, b_val=b, u_val=u, ax=ax, alpha=0.1, plot_sphere=False)
    plot_trajectory(Xt1, Xt2, u=u, x0=X0, x1=X1, ax=ax)


if __name__ == '__main__':
    plot_optimal(u=[+1, -1])
    plot_optimal(u=[-1, +1])

    plt.show()