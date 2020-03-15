Dans ce mini-projet on va étudier un système de contrôle.
On trouvera ensuite les commandes qui permettent d'aller
d'un état à l'autre du système de façon optimale
par rapport à un critère donné suivant les méthodes
étudiés en cours.

# Système de contrôle

Un système de contrôle est un système dont l'état à l'instant $t$
est décrit par $n$ variables $x_1,\ldots,x_n$. On définit le vecteur d'état
du système $x(t)=\pp{x_1(t),\ldots,x_n(t)}\in\R^n$.

On dispose de $m$ variables de contrôle $u_1,\ldots,u_m$ définissant la commande
appliquée à l'instant $t$ par le vecteur $u(t)=\pp{u_1(t),\ldots,u_m(t)}\in\R^m$.

L'évolution du système est décrit par son **équation d'état**.
$$\boxed{\dot{x}(t) = f(t, x(t), u(t))}$$

Un système contrôle est défini par une équation d'état et une condition initiale.
$$\begin{cases}
\dot{x} = f(t, x(t), u(t)) & \forall t\in [t_0,\infty[\\
x(t_0) &\text{condition initiale donnée}
\end{cases}$$

On dit qu'un état $x^*$ est accessible à partir d'un état initial $x(t_0)$
s'il existe une commande $\bar u(t)$ qui ramène le système à l'état
$x^*$ en *temps fini* $T$.
$$ \exists T<\infty, \dot x = f(t,x(t),\bar u(t)) \implies x(T)=x^*$$

# L'état du système

On considère un système bilinéaire $S$ donné par son
équation d'état.

$$ \xd = A\x + \u B\x \qou
\cas{
    x(t)\in\R^n\\
    u(t)\in\R\\
    A,B\in\Rnn
}$$

L'équation d'état se reécrit
$$ \xd = \pp{A+\u B}\x$$

Pour un contrôle constant $\u=u,\forall t\in[t_0,T]$ le système devient
$$ \xd = \underbrace{\pp{A + uB}}_{M(u)}\cdot\x \qou M(u)\in\Rnn$$

L'équation d'état devient alors une EDO linéaire.

\begin{align}
&    \xd = M\x \\
&    \xd - M\x = 0 \\
&    e^{-Mt}\xd - e^{-Mt} M\x = 0 \\
&    e^{-Mt}\cdot\frac{\x}{\dt} + \diff{e^{-Mt}}\cdot\x = 0 \\
&    \diff{e^{-Mt}\cdot\x} = 0\\
&    \int_{t_0}^t\dif{e^{-Ms}\cdot x(s)}{s}\ds = 0\\
&    \left. e^{-Ms}\cdot x(s)\right\rvert_{t_0}^t = 0\\
&    e^{-Mt}\cdot\x - e^{-Mt_0}\cdot x(t_0) = 0\\
&    \x = e^{M\pp{t-t_0}}\cdot x(t_0)
\end{align}

L'état du système est donc connu à chaque instant $t\in[t_0,T]$.
$$\boxed{\x = e^{M\pp{t-t_0}}\cdot x(t_0)}$$

Prenons un exemple pour $n=3$ la matrice $M$ donnée.

```{python initialize}
# initialize matrix symbols 
a, b, u = sp.symbols('a b u', real=True)
w = sp.sqrt(a**2 + b**2 + u**2)
ww = sp.symbols('w', positive=True)

# system matrices
A = sp.Matrix([[0, a, b], [-a, 0, 0], [-b, 0, 0]])
B = sp.Matrix([[0, 0, 0], [ 0, 0,-1], [ 0, 1, 0]])
M = A + u*B
```

$$M = `r tex(py$M)`\qtext{avec}\cas{
    a,b\in\R\\
    w = \sqrt{a^2+b^2+u^2}
}$$

Afin de connaître l'état du système, nous avons besoin
de calculer $e^{Mt}$, on utilisera la transformée
de laplace inverse.

$$e^{Mt} = \laplace{\pp{sI-M}^{-1}}(t)$$

On utilisera la formule suivante pour calculer
la matrice inverse de $sI-M$.

$$(sI-M)^{-1} = \frac{1}{\det(sI-M)}
    \transp{\mathrm{com}(sI-M)}$$

```{python inverse}
s = sp.symbols('s') 
Ms = s * sp.eye(3) - M  # sI - M
Ms_det = sp.factor(Ms.det()).subs({w**2:ww**2})
Ms_com = Ms.cofactor_matrix()
Ms_inv = Ms_com.transpose() / Ms_det
```

$$sI-M = `r tex(py$Ms)`$$

Ensuite on calcule le déterminant de la matrice résultante
$$\det(sI-M) = `r tex(py$Ms_det)` = s(s-iw)(s+iw)$$

On en déduit que $\mathrm{Sp}(M) = \sset{0, iw, -iw}$.
Comme la diagonalisation de $M$ est complexe, il vaut mieux
calculer son exponentielle par la transformée de laplace inverse.

$$(sI-M)^{-1} = \frac{1}{`r tex(py$Ms_det)`}\transp{`r tex(py$Ms_com)`}
= `r tex(py$Ms_inv)`$$

```{python laplace_inverse}
t, t0 = sp.symbols('t t_0', positive=True)
exp_Mt = sp.inverse_laplace_transform(Ms_inv, s, t)
exp_Mt0 = exp_Mt.subs({t:t-t0})
```

En appliquant la transformée de laplace inverse on obtient
$$ e^{Mt} = \laplace{\pp{sI-M}^{-1}}(t) = `r tex(py$exp_Mt)`$$

L'état du système à l'instant $t$ est donc donné par l'expression
de son vecteur d'état.

$$ \x = `r tex(py$exp_Mt0)` \cdot\bmat{x_1(t_0)\\x_2(t_0)\\x_3(t_0)}$$