# Contrôle Optimal

## Notions d'optimalité

Un problème de contrôle optimale consiste à trouver une commande $u^*$
qui ramène le système à un état $x^*$ en minimisant un critère donnée par
une fonction qu'on appelle **fonction coût**.

La fonction coût représente le coût du système à un état donné suivant une commande donnée.
$$ g:\R\times\R^n\times\R^m \rightarrow \R $$

Le coût total d'une commande $u:[t_0,T]\rightarrow\R^m$ se donne par
$$ C(u) = \int_{t_0}^T g(t, \x, \u) \dt $$

Afin de trouver la commande optimale, on définit le système généralisé suivant:
On pose $$x_0(t)=\int_{t_0}^t g(s, x, u) \ds$$
on remarque que $x_0(t_0)=0$ et $x_0(T)=C(u)$.
De plus, $\dot{x}_0(t) = g(t, \x, \u)$. On définit l'état généralisé du système

$$\xt(t)=\pmat{x_0(t)\\x(t)}\mat{\in\R~\\\in\R^n}
    \implies \xt(t)\in\R^{n+1}$$

On en déduit l'équation d'état généralisé
$$\xtd(t) = \pmat{\dot{x}_0(t)\\\xd}=
    \pmat{g(t, \x, \u)\\f(t, \x, \u)}=
    \ft(t, \x, \u)\in\R^{n+1}$$

Le problème revient à trouver une commande $u^*$ admissible
($\in\U$) qui ramène le système généralisé à l'état
$$\xtopt=\pmat{\min\limits_{u\in\U} C(u)\\ x^*}$$

## Temps optimalité

Afin de trouver une commande optimale en temps minimale, il suffit de définir
une fonction coût $g(t, \x, \u)=g(t)$ croissante par rapport au temps.

En pratique, il suffit de prendre $g(t)=1$ car le coût d'une commande sera donc $C(u)=T-t_0$.

## Le Hamiltonien du système

On définit le Hamiltonien du système $S=(x, u, f)$ par la fonction.

\begin{align}
H: \R^n\times\pp{\R^n\setminus\sset{0_{\R^n}}}\times\R^m &\rightarrow \R\\
    (x, p, u) &\mapsto \dotp{p(t)}{f(t,\x,\u)}
        = \dotp{p(t)}{\xd(t)}
\end{align}

où $p$ est le vecteur adjoint au système.

Pour chercher la commande optimale, on considère le Hamiltonien du système généralisé
$\tilde{S}=(\xt,u, \ft)$.

$$ \Ht(\xt, \tilde{p}, u) = \dotp{\tilde{p}}{\ft(t,x,u)}
    = p_0\cdot g(t,x,u) + \dotp{p}{f(t, x, u)}
    = \Ht(x, p_0, p, u)$$

Par abus de notation, on écriré le Hamiltonien généralisé sans le tilde
car on ne s'intéresse pas au Hamiltonien du système $S$.

On reconsidère l'équation d'état du système : $\dot x = Ax + uBx$

Il s'agit d'un système affine de la forme
$\dot x = F(x) + u\cdot G(x)$
que l'on peut identifier les fonctions $F$ et $G$, d'où

$$\begin{cases}
    F(x) = Ax &\Rightarrow \nabla F(x) = A\\
    G(x) = Bx &\Rightarrow \nabla G(x) = B
\end{cases}$$

Cette écriture nous permettra de mieux exprimer la fonction
hamiltonienne. Le Hamiltonien se donne alors par:

$$H(x, p, u) = p_0 + \dotp{p}{\xd}
    = p_0 + \dotp{p}{M\x}
    = p_0 + \dotp{p}{Ax + uBx}$$

```{python initial_conditions}
X0 = sp.Matrix([0, 0, 1])
X1 = sp.Matrix([0, 1/sp.sqrt(2), 1/sp.sqrt(2)])
```

Revenons à notre problème, on cherche à emmener le système
de l'état $x(t_0=0)$ à l'état $x(T)$ en temps minimal étant donné
$$ x(0) = `r tex(py$X0)` \qtext{et}
    x(T) = \bmat{0\\\frac{1}{\sqrt 2}\\\frac{1}{\sqrt 2}}$$

On applique le Principe du Maximum de Pontryagin (PMP),
on pose les équations du PMP:

**L'équation d'état** se donne par $\xd = \pd{H}{p}$.
En effet, $H(x,p,u) = p_0 + \dotp{p}{Mx}$ d'où
$\pd{H}{p} = Mx = \xd$.

**L'équation adjointe** se donne par $\dot{p}(t) = -\pd{H}{x}$.
\begin{align}
H(x, p, u) &= p_0 + \dotp{p}{Mx}\\
    &= p_0 + \dotp{Mx}{p}\\
    &= p_0 + \transp{\pp{Mx}}p\\
    &= p_0 + \transp{x}\transp{M}p\\
    &= p_0 + \dotp{x}{\transp{M}p}\\
\Rightarrow\pd{H}{x} &= \transp{M}p
\end{align}

Dans notre problème $M$ est anti-symétrique
(i.e. $\transp{M} = -M$), d'où $\pd{H}{x} = -Mp$.

L'équation adjointe est donc $\dot{p}(t) = M p(t)$,
on remarque que l'équation adjointe est l'équation
d'état du système dans le cas d'un système *bilinéaire
affine anti-symétrique*.

On en déduit que $p(t) = e^{Mt}p(0)$.

Il faut que $p_0 \leq 0$ également.

**La condition de maximisation :**
On cherche une commande $u^*$ qui maximise le Hamiltonien
$\forall u\in\U$.

\begin{align}
u^* &= \argmax_{u\in\U}~H(x, p, u)\\
    &= \argmax_{u\in\U}~\pp{p_0 + \dotp{p}{Ax + uBx}}\\
    &= \argmax_{u\in\U}~\pp{p_0 + \dotp{p}{Ax} + u\dotp{p}{Bx}}\\
    &= \argmax_{u\in\U}~\pp{u\dotp{p}{Bx}}\\
    &= \argmax_{u\in\U}~u\cdot\underbrace{\dotp{p}{Bx}}_{\phi(t)}\\
    &= \argmax_{u\in\U}~u\cdot\phi(t)
\end{align}

Comme le Hamiltonien est une fonction affine par rapport à $u$,
le maximum ne peut être atteint que pour $u$ bornée.
On prend $\abs{u(t)}\leq 1$, d'où
$$ u^*(t) = \sgn(\phi(t)) = \sgn(\dotp{p}{Bx})$$

On appelle $\phi(t)=\dotp{p}{Bx}$ la **fonction de commutation**.

```{python}
Xu = exp_Mt * X0

Xu_pos = Xu.subs({u:1})
Xu_neg = Xu.subs({u:-1})
```

Il s'agit donc d'une commande *bang-bang*, on cherche alors
les trajets optimales.

$$ x_{u=1}(t) = `r tex(py$Xu_pos)` \qtext{,} x_{u=-1}(t) = `r tex(py$Xu_neg)`$$

On rappelle qu'on avait trouvé que $\mathrm{Sp}(M) = \sset{0, iw, -iw}$,
on peut ainsi en déduire que :

- L'état pour une commande fixée $x_u(t)$ vit dans un hyperplan de $\R^3$
que l'on note $E_u$ tel que $\mathrm{dim}(E_u)=2$ car $M$ a une valeur propre nulle.
- Les trajectoires de $x_u(t)$ (pour une commande fixée) sont elliptiques
dans $E_u$ car les deux valeurs propres non-nulles de $M$ sont conjugées
et purement imaginaires.
- Pour $a,b$ fixés, l'état du système vit dans la surface d'un ellipsoïde,
si les trajectoires $x_{u=1}$ et $x_{u=-1}$ sont indépendants
alors tous les points sur la surface de l'ellipsoïde sont accessible à partir de $x(t_0)$.

On trace le champs de vecteur du système
à partir de $x(0)=(0, 0, 1)$ pour $a=1,b=0$
sur et les trajectoires bang-bang
à l'aide de *python*.

```{python bang_bang, result='markup'}
from visuals import plot_initial_trajectory, plot_vector_field

fig = plt.figure(figsize=(10, 10))
ax = fig.gca(projection='3d')
plot_initial_trajectory(u=[1, -1, 0], a=1, b=0, ax=ax)
plot_vector_field(M, u_val=[1, -1, 0], a_val=1, b_val=0, ax=ax, plot_sphere=False, alpha=0.1)
ax.legend()

plt.show()
```

Nous cherchons alors la commande optimale $u^*(t)$
à chaque instant $t$.

- Si $\phi(t)\geq 0$ alors $u^*(t) = 1$.
- Si $\phi(t)\leq 0$ alors $u^*(t) = -1$.
- Si $\phi(t) = 0$ alors $u^*(t) = u^*$ trajectoire singulière.

Nous identifiant les matrices $A$ et $B$ de notre système.

$$ M(u) = `r tex(py$M)`
    = \underbrace{`r tex(py$A)`}_{A} +u \underbrace{`r tex(py$B)`}_{B}$$

D'où $\phi(t) = \dotp{p(t)}{Bx(t)} = -p_2(t)x_3(t) +p_3(t)x_2(t)$.
De plus, $x(0)=\transp{(0,0,1)}$ donc $\phi(0) = -p_2(0)$.

Le temps de commutation de la trajectoire optimale
est le temps de changement de signe de $\phi$,
on cherche alors $t\geq t_0$ tel que $\phi(t) = 0$.


A l'aide du PMP, nous avons pu determiner que notre commande
est une commande *bang-bang*, nous pouvons trouver ainsi
la commande optimale en cherchant *iterativement*
le temps et lieu de commutation permettants d'aller à $x(T)$.

En effet, en prenant $u^*(0)=1$, notre commande est de la forme
$$x(t)=\begin{cases}
    e^{M(1)t}\cdot x(0) &\text{si}~t\in[0,t_c[\\
    e^{M(-1)t}\cdot x(t_c) &\text{si}~t\in[t_c,T]
\end{cases}$$
où $t_c$ est le temps de commutation et $x(t_c)$ est le lieu de commutation.

En discrétisant le temps $t_n=n\Delta t$ pour $\Delta t$
suffisament petit, nous pouvons chercher iterativement pour
chaque $t_n$ si le trajectoire $e^{M(-1)t_n}\cdot x(t_n)$
mène à $x(T)$ en vérifiant s'il existe un point $x_n$
du second trajectoire tel que $\norm{x_n-x(T)}<\eps$
pour $\eps$ suffisament petit.

```{python optimal, echo=TRUE, results='markup'}
from control import plot_vector_field, plot_trajectory, find_optimal

u = [1, -1]
a = 1
b = 0

Xt1, Xt2, Topt = find_optimal(exp_Mt, X0, X1, a_val=a, b_val=b, u_first=u[0], N=1000, epsilon=0.01)

fig = plt.figure(figsize=(10, 10))
ax = fig.gca(projection='3d')
plot_vector_field(M, R=float(X0.norm()), a_val=a, b_val=b, u_val=u, ax=ax, alpha=0.1, plot_sphere=False)
plot_trajectory(Xt1, Xt2, u=u, x0=X0, x1=X1, ax=ax)
plt.show()
```

Vérifions si le trajectoire bang-bang qui commence par $u(0)=-1$
est la trajectoire optimale.

```{python optimal2, echo=TRUE, results='markup'}
u = [-1, 1]

Xt1, Xt2, Topt = find_optimal(exp_Mt, X0, X1, a_val=a, b_val=b, u_first=u[0], N=1000, epsilon=0.01)

fig = plt.figure(figsize=(10, 10))
ax = fig.gca(projection='3d')
plot_vector_field(M, R=float(X0.norm()), a_val=a, b_val=b, u_val=u, ax=ax, alpha=0.1, plot_sphere=False)
plot_trajectory(Xt1, Xt2, u=u, x0=X0, x1=X1, ax=ax)
plt.show()
```

En effet, le trajectoire optimale est obtenu pour la commande
$u^*=1$ pour $t\in[0,t_c]$ puis $u^*=-1$ pour $t\in[t_c,T]$
où $t_c\approx 0,55$
et $x(t_c)\approx\bmat{-0,14\\-0,49\\0,87}$
en temps optimal $T\approx 1,65$.