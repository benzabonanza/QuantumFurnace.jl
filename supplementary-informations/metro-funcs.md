# Update note to this file:
- In the thesis we will mainly use a = 0, b != 0. a != 0 was a try to avoid the time singularity in the Metropolis cases, i.e. an alternative regularization to the $\eta$ regularization of Chen et al. However, this did not lead to analytically better scaling results, nor did it have a noticeable improvemenet practically in the code. (We will eventually have to change the code a bit to reflect this, since atm, we don't really have a version where only b != 0 and a = 0.)

Notation update: b -> s. But do not confuse this s-parameter here with the s-parameter of Chen et al. Here it is a shift for the lower limit of the linear combination integration to smooth out the kink and lead to some improvements for the quadrature errors. While Chen's s-parameter is for the upper limit of the integral that they send to $\infty$ in the analytical derivations.


# Gammas
$\beta := \frac{2\omega_\gamma}{\sigma^2 + \sigma_\gamma^2}$
## KMS

- Gaussian $$\gamma(\omega) = e^{-\frac{(\omega + \omega_\gamma)^2}{2\sigma_\gamma^2}}$$
- Metropolis with $g(x) = 1 / (\sqrt{2\pi} \sigma_\gamma)$ and lower limit of integration $\int_{\beta \sigma^2/ 2}^\infty dx$ $$\gamma(\omega) = e^{-\beta \max \left(\omega + \frac{\beta \sigma^2}{2}, 0\right)}$$
- Smooth with $g(x) = \frac{\kappa}{\sqrt{2\pi}\sigma_\gamma} e^{-a \beta x}$, with $\kappa = \sqrt{4a + 1} e^{a\beta^2\sigma^2 / 2} \times \kappa_{num}(b)$ and $\int_{\beta \sigma^2 (1 + b) / 2}^\infty dx$ $$\gamma_b(\omega) = \gamma_0(\omega) \frac{\kappa_{num}}{2}\left[\text{erfc}(z_-) + e^{4\sqrt{AB}}\text{erfc}(z_+) \right]$$with
	- $A = \frac{\beta}{4}(4a + 1)$
	- $B = \frac{\beta}{4}\tilde\omega^2$
	- $C = \frac{\beta}{2}(a \beta \sigma^2 + \tilde\omega)$
	- $\tilde\omega := \omega + \frac{\beta \sigma^2}{2}$
	- $z_\pm = \sqrt{A}\sqrt{b\,\beta \sigma^2 / 2} \pm \sqrt{\frac{B}{b\,\beta \sigma^2 / 2}}$
and $$\gamma_0(\omega) = e^{-\frac{\beta}{2}(\tilde\omega + \sqrt{4a + 1}|\tilde\omega|)}$$Bit more explanation: 
- $\kappa$'s first factor was chosen such that the $b=0$ case gets normalized in a way that $\|\gamma(\omega)\|_\infty = 1$. If it was less than $1$ then we would unnecessarily constrain our transitions. $\kappa_{num}(b)$ is  a numerical normalization for cases when $\beta > 0$, since there there is no simple analytical form, and it is just much easier to get this factor numerically via `gamma_norm_factor`.
- We get back the Metropolis case if we send $a,b\rightarrow 0$. We get back Glauber if $a\rightarrow 0$ and $b = 2$.
- The $b = 0$ case was derived via modified Bessel function of the 2nd kind. And for $b > 0$ it took the form of an incomplete Bessel function or with a substitution of $x = u^2$ the integral form of the heat diffusion, hence the appearance of the error functions in both cases.

## Approx. GNS
Main difference here is the fact that now $\gamma(\omega)$ has to fulfill the KMS condition in order to  provably lead to an approximate GNS DB Lindbladian. In the KMS case this was not the case. But as it turns out all this really means is:
- In the Gaussian case set $\beta := \frac{2\omega_\gamma}{\sigma_\gamma}$
- In the other cases shift the transitions functions via $\beta \sigma^2 / 2$ such that instead of having $\tilde\omega = \omega + \beta\sigma^2 / 2 \rightarrow\omega$ now we have just.


# Alphas
## KMS
$\sigma_\gamma = \sqrt{\frac{2\omega_\gamma}{\beta} - \sigma^2}$
- Gaussian $$\alpha_{\nu_1\nu_2} = \frac{\sigma_\gamma}{\sqrt{\sigma^2 + \sigma_\gamma^2}} \exp\left(-\frac{(\nu_1 + \nu_2 + 2\omega_\gamma)^2}{8(\sigma^2 + \sigma_\gamma^2)}\right)\exp\left(-\frac{(\nu_1 - \nu_2)^2}{8\sigma^2}\right)$$
- Smooth Metropolis $$\alpha_{\nu_1\nu_2} = \frac{\kappa_{num}^{b>0}}{2}e^{\frac{a \beta^2\sigma^2}{2}}e^{-\frac{(\nu_1-\nu_2)^2}{8\sigma^2}}e^{-\frac{\beta}{4}(\nu_1 + \nu_2)}e^{\frac{\beta}{4}\sqrt{4a + 1}|\nu_1 + \nu_2|}\left[\text{erfc}(z_-) + e^{4\sqrt{AB}}\text{erfc}(z_+)\right]$$with
	- $A = \frac{\beta}{4}(4a + 1)$
	- $B = \frac{\beta}{16}(\nu_1 + \nu_2)^2$
	- $C = \frac{\beta}{4}(\nu_1 + \nu_2)$
	- $z_\pm = \sqrt{A} \sqrt{\beta\sigma^2 (1 + b) / 2} \pm \sqrt{\frac{B}{\beta\sigma^2 (1 + b) / 2}}$


## GNS
Since $\gamma(\omega)$ for the approx. GNS case was just a shifted one, so will be the Kossakowski matrices, but importantly we only shift the following way:$$|\nu_1 + \nu_2| \rightarrow |\nu_1 + \nu_2 + \frac{\beta\sigma^2}{2}|$$Meaning, that only the absolute values get shifted, not the other $\nu$'s.


# Time domain functions for $B$
Only KMS here, since there is no coherent term in the approx. GNS case.
$b_1$ is the same for all: $$b_1(t) = \frac{2\sqrt{\pi}}{\beta \sigma}e^{\beta^2\sigma^2/8}\left[\frac{1}{\cosh(\frac{2\pi t}{\beta \sigma})}\star_t \sin(-\beta \sigma t) e^{-2t^2}\right]$$
- Gaussian $b_2$: $$b_2(t) = \frac{\beta \sigma_\gamma}{\pi\sqrt{\pi}}e^{-4t^2\beta \omega_\gamma - 2it\beta\omega_\gamma}$$
- Metropolis from Chen $b_2$: $$b_2(t) = \frac{1}{2\sqrt{2}\pi^2}\frac{e^{-\sigma^2\beta^2 t(2t + i)} + \mathbb{I}(|t|\leq\eta)i(2t + i)}{t(2t + i)}$$
- Smooth Metropolis $b_2$:$$b_2(t) = \kappa_{num}^{b>0}\frac{\sqrt{4a + 1}}{\sqrt{2}\pi^2}e^{-a b/2} \frac{e^{-\sigma^2\beta^2 t(2t + i)(1 + b)}}{4t^2 + a + 2it}$$


