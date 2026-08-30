# Kissinger equation

$$
\frac{\partial x}{\partial t}(t) = A (1 - x) e^{-\frac{E}{k(T_0 + r t)}}
$$

Via separation of variables we obtain:

$$
\int_0^x \frac{\mathrm dx'}{1 - x'} = A\int_0^t e^{-\frac{E}{k(T_0 + r t')}}\mathrm dt'
$$

We replace $E / k(T_0 + rt')$ by $\xi$, so that we obtain:

$$
\log(1 - x) = -\frac{AE}{kr}\int_{\frac{E}{kT_0}}^{\frac{E}{k(T_0 + r t)}}\frac{e^{-\xi}}{\xi^2}\mathrm d\xi
$$
