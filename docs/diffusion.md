# Diffusion equation with time-dependent temperature

I was looking at [this material](https://www.uni-muenster.de/imperia/md/content/physik_tp/lectures/ws2016-2017/num_methods_i/heat.pdf) and realized that the diffusion equation with a time-depenent diffusion coefficient can be solved analytically. The part for $X(x)$ remains the same. The tricky part is obviously for $T(t)$, because $D$ also depends on the time. I use $Z$ (= Zeit) for the time dependent part of the solution in the following, because $T$ is used for the temperature in our case.

\begin{align}
\frac{\partial Z}{\partial t} + D\left(\frac{\pi n}{L}\right)^2 Z = 0
\end{align}

Separation of variables gives:


\begin{align}
\frac{\mathrm d Z}{Z} =&- D\left(\frac{\pi n}{L}\right)^2 \mathrm dt \\
\longleftrightarrow \log \frac{Z}{Z_0} = &- D_0\left(\frac{\pi n}{L}\right)^2 \int_0^{t} \exp\left(-\frac{E}{k_{\mathrm B}(T_0+rt')}\right)\mathrm dt' \\
=& -D_0\frac{E}{rk_{\mathrm B}}\left(\frac{\pi n}{L}\right)^2 \int_{-\frac{E}{k_{\mathrm B}T_0}}^{-\frac{E}{k_{\mathrm B}(T_0+rt)}} \frac{e^{y}}{y^2}\mathrm dy \\
=& D_0\frac{E}{rk_{\mathrm B}}\left(\frac{\pi n}{L}\right)^2 \int_{-\frac{E}{k_{\mathrm B}T_0}}^{-\frac{E}{k_{\mathrm B}(T_0+rt)}} e^y\frac{\partial}{\partial y}\frac{1}{y}\mathrm dy \\
=& -D_0\frac{E}{rk_{\mathrm B}}\left(\frac{\pi n}{L}\right)^2 \left(\frac{k_{\mathrm B}(T_0+rt)}{E}e^{-\frac{E}{k_{\mathrm B}(T_0+rt)}} - \frac{k_{\mathrm B}T_0}{E}e^{-\frac{E}{k_{\mathrm B}T_0}} + \int_{-\frac{E}{k_{\mathrm B}T_0}}^{-\frac{E}{k_{\mathrm B}(T_0+rt)}} \frac{e^{y}}{y}\mathrm dy\right) \\
=& -D_0\frac{E}{rk_{\mathrm B}}\left(\frac{\pi n}{L}\right)^2 \left(\frac{k_{\mathrm B}(T_0+rt)}{E}e^{-\frac{E}{k_{\mathrm B}(T_0+rt)}} - \frac{k_{\mathrm B}T_0}{E}e^{-\frac{E}{k_{\mathrm B}T_0}} + e_i^{-\frac{E}{k_{\mathrm B}(T_0+rt)}} -  e_i^{-\frac{E}{k_{\mathrm B}T_0}}\right)\\
\longleftrightarrow Z = & Z_0\exp\left(-D_0\frac{E}{rk_{\mathrm B}}\left(\frac{\pi n}{L}\right)^2 \left(\frac{k_{\mathrm B}(T_0+rt)}{E}e^{-\frac{E}{k_{\mathrm B}(T_0+rt)}} - \frac{k_{\mathrm B}T_0}{E}e^{-\frac{E}{k_{\mathrm B}T_0}} + e_i^{-\frac{E}{k_{\mathrm B}(T_0+rt)}} -  e_i^{-\frac{E}{k_{\mathrm B}T_0}}\right)\right)\\
\end{align}
where $e_i$ is the [exponential integral Ei](https://en.wikipedia.org/wiki/Exponential_integral)

Therefore, the solution of the diffusion equation $u(x, t)$ is given by:

\begin{align}
u(x, t) = \sum_{n=1}^{\infty}\left(\frac{2}{L}\int_0^Lf(\xi)\sin\left(\frac{\pi n}{L}\xi\right)\mathrm d\xi\right)\sin\left(\frac{\pi n}{L}x\right)\exp\left(-D_0\frac{E}{rk_{\mathrm B}}\left(\frac{\pi n}{L}\right)^2 \left(\frac{k_{\mathrm B}(T_0+rt)}{E}e^{-\frac{E}{k_{\mathrm B}(T_0+rt)}} - \frac{k_{\mathrm B}T_0}{E}e^{-\frac{E}{k_{\mathrm B}T_0}} + e_i^{-\frac{E}{k_{\mathrm B}(T_0+rt)}} -  e_i^{-\frac{E}{k_{\mathrm B}T_0}}\right)\right)
\end{align}
