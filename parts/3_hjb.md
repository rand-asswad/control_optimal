## Equation de Hamilton-Jacobi-Bellman

On définit la fonction de Bellman par

$$V(t,x(t)) = \min_{u\in[t,T[} \int_t^T g(s, x(s), u(s))\ds$$
avec $V(T, x(T)) = 0$.

Dans le cas d'un problème de temps-optimalité nous avons
$g(t, x(t), u(t)) = 1$ d'où

$$V(t, x(t)) = \min_{u\in[t,T[} (T-t) = T-t$$

L'équation de Hamilton-Jacobi-Bellman (HJB) se donne par

$$\pd{V}{t} + \min_{u\in\U}\pp{\pd{V}{x}\cdot f(t,\x,\u) + g(t,\x,\u)} = 0$$

Dans notre système,
$$\begin{cases}
    \pd{V}{t} = -1\\
    \pd{V}{x} = 0
\end{cases}\implies (-1) + \min_{u\in\U}\pp{0 + 1} = 0$$

L'équation de HJB est donc vérifiée pour tout $u\in\U$.

# Conclusion

Le système étudié represente deux cas spéciaux de système de contrôle.

- **Système affine en $u$:** $f(t, x, u) = F(x) + uG(x)$
- **Système linéaire en $x$:** $f(t, x, u) = Mx$

De plus, nous avions une matrice $M$ antisymétrique,
ce qui a présenté un cas intéressant d'un vecteur adjoint
vérifiant lui-même l'équation d'état.

Néanmoins, la recherche de la commande qui mène le système
en temps minimal de $x(0)$ à $x(T)$ représente un cas
important de problème de contrôle optimal.

Ces propriétés ont donné des cas spéciales de la forme
du Hamiltonien et ses dérivés, et de l'équation de HJB.

En conclusion, ce projet a été une opportunité d'étudier
au plus près un système de contrôle optimal en trois dimensions,
et de mieux visualiser l'état du système et son comportement,
et a ainsi permis de mieux comprendre le domaine de contrôle optimal.