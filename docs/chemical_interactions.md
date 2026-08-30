# How to include local chemical interactions in the continuum model

The chemical interactions $U$ between H are given by:

$$
U = \sum_{i,j} \phi(x_i) \varepsilon(x_i - x_j) \phi(x_j)
$$

where $\phi(x)$ is the occupation probability at site $x$ and $\varepsilon(x_i - x_j)$ is the interaction energy between site $x_i$ and $x_j$. For the sake of simplicity, I considered only the first order interaction, meaning there is no dependence on $\phi^2$. Then the local binding energy $u(x)$ (or basically the chemical potential) at site $x$ is given by:

$$
u(x) = \sum_j \varepsilon(x - x_j) \phi(x_j)
$$

We know that H atoms segregate like a cloud, meaning there is no sharp edges on the segregation. So we can take a Taylor expansion on the concentration field:

$$
\phi(x_j) \approx \phi(x) + (x_j - x) \nabla \phi(x) + \frac{1}{2} (x_j - x)^T \Delta \phi(x) (x_j - x)
$$

Then we can rewrite the binding energy equation by:

$$
u(x) = \phi(x)\sum_j \varepsilon(x - x_j) + \nabla\phi(x)\underbrace{\sum_j (x_j - x)\varepsilon(x - x_j)}_{=0} + \frac{1}{2}\Delta \phi(x)\sum_j \varepsilon(x - x_j) (x_j-x)(x_j-x)^T
$$

The second order derivative has only diagonal components, and all the components are the same due to spatial symmetry. And since $\varepsilon$ depends only on the distance between two sites, we can use the site-independent constants $A_0$ and $A_1$ to rewrite the above equation two:

$$
u(x) = A_0\phi(x) + A_1\Delta \phi(x)
$$
