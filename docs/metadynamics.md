# Polynomial based free energy

$$
U(x) = w\sum_i\exp\left(-\frac{(x - x_i)^2}{2\sigma^2}\right)
$$

The idea is to replace $U$ by a polynomial:

$$
U(x) \approx U(x_0) + \left .\frac{\partial U}{\partial x}\right |_{x_0}(x-x_0) + \frac{1}{2} \left. \frac{\partial^2 U}{\partial x^2}\right |_{x_0}(x-x_0)^2 + \frac{1}{6} \left. \frac{\partial^3 U}{\partial x^3}\right |_{x_0}(x-x_0)^3 + ...
$$

Now define

$$
A(x) = \exp\left(-\frac{x^2}{2\sigma^2}\right)
$$

and

$$
\alpha(x) = \frac{x}{\sigma}
$$

Then

$$
\sigma\frac{\partial A}{\partial x} =-\alpha A
$$

$$
\sigma^2\frac{\partial^2 A}{\partial x^2} =(-1 + \alpha^2) A
$$

$$
\sigma^n\frac{\partial^n A}{\partial x^n} =\sum_k c^n_k \alpha^kA
$$
