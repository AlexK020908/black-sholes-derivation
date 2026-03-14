# From Random Walks to Black-Scholes: A Complete Derivation
---

## Part 1: The Discrete World — Random Walks

### 1.1 The Simple Symmetric Random Walk

Imagine flipping a fair coin at each time step.

Define a process $\{X_n\}$ where:

$$X_n = X_0 + \sum_{i=1}^{n} Z_i$$

where each $Z_i$ is an independent random variable:

$$Z_i = \begin{cases} +1 & \text{with probability } \frac{1}{2} \\ -1 & \text{with probability } \frac{1}{2} \end{cases}$$

**Key properties:**

- $E[Z_i] = 0$ (fair coin, no drift)
- $\text{Var}(Z_i) = E[Z_i^2] - (E[Z_i])^2 = 1 - 0 = 1$
- $E[X_n] = X_0$ (it's a martingale)
- $\text{Var}(X_n) = n$ (variance grows linearly with time, using independence of the $Z_i$)

### 1.2 Asymmetric Random Walk: Drift + Noise

Now generalize. Instead of fair $\pm 1$ steps, let:

$$Z_i = \begin{cases} +u & \text{with probability } p \\ -d & \text{with probability } 1 - p \end{cases}$$

Now $E[Z_i] = pu - (1-p)d \neq 0$ in general. The walk has **drift**.

We can always decompose any step into what we expect and the surprise:

$$Z_i = \underbrace{E[Z_i]}_{\text{expected part}} + \underbrace{(Z_i - E[Z_i])}_{\text{surprise (mean-zero noise)}}$$

**Concrete example:** Suppose $Z_i = +3$ with prob 0.6 and $Z_i = -2$ with prob 0.4.

- $E[Z_i] = 0.6(3) + 0.4(-2) = 1.0$
- The surprise is $Z_i - 1.0$, which is $+2$ or $-3$. Check: $E[\text{surprise}] = 0.6(2) + 0.4(-3) = 0$. Mean zero. ✓

**Important:** At this point, there's no time anywhere. These are just discrete steps. We'll bring in time next.

---

## Part 2: From Discrete to Continuous — Brownian Motion

### 2.1 Putting the Random Walk on a Time Grid

Take a time interval $[0, T]$ and divide it into $n$ steps of size $\Delta t = T/n$.

At each step, the process moves by $\Delta X_i = \pm\delta$ with equal probability, for some step size $\delta$ that we haven't decided yet.

After $n$ steps:

$$X(T) - X(0) = \sum_{i=1}^{n} \Delta X_i$$

Each $\Delta X_i$ has mean $0$ and variance $\delta^2$, so:

$$\text{Var}(X(T)) = n \cdot \delta^2 = \frac{T}{\Delta t} \cdot \delta^2$$

### 2.2 What Should the Step Size $\delta$ Be?

As we make the grid finer ($\Delta t \to 0$, $n \to \infty$), we need $\text{Var}(X(T))$ to stay finite and nonzero. Otherwise we get a useless limit.

The CLT tells us exactly what $\delta$ must be. Here's how.

**The CLT says:** If $\Delta X_1, \ldots, \Delta X_n$ are i.i.d. with mean $m$ and variance $v^2$, then:

$$\frac{\left(\sum_{i=1}^{n} \Delta X_i\right) - nm}{v\sqrt{n}} \xrightarrow{d} N(0, 1)$$

Equivalently: $\sum \Delta X_i \approx N(nm, \; nv^2)$.

**Our steps have:** $m = 0$ (symmetric), $v^2 = \delta^2$. So:

$$X(T) - X(0) = \sum \Delta X_i \;\approx\; N(0, \; n\delta^2)$$

**We want Brownian motion as the limit.** We haven't defined Brownian motion yet — we're *constructing* it. We want a process where the total variance over time $T$ is proportional to $T$. That is, we want:

$$n\delta^2 = (\text{some constant})^2 \times T$$

Call that constant $\sigma$ (it's just a name — it could be any positive number). Then:

$$n\delta^2 = \sigma^2 T$$

Substitute $n = T/\Delta t$:

$$\frac{T}{\Delta t} \cdot \delta^2 = \sigma^2 T$$

Cancel $T$:

$$\frac{\delta^2}{\Delta t} = \sigma^2$$

$$\delta^2 = \sigma^2 \Delta t$$

$$\boxed{\delta = \sigma\sqrt{\Delta t}}$$

**Where $\sigma$ came from:** nowhere deep. The math says "$\delta^2$ must be proportional to $\Delta t$." We call the proportionality constant $\sigma^2$. Different processes have different $\sigma$ values.

**Why the CLT derivation is better than guessing:** We didn't try different scalings and see which works. Instead, we said "the limit must be Gaussian with variance $\sigma^2 T$" (by the CLT), and this *forced* $\delta = \sigma\sqrt{\Delta t}$. No other choice is possible.

### 2.3 Standard Brownian Motion $B_t$

With $\sigma = 1$, the limiting process is called **standard Brownian motion** $B_t$:

1. $B_0 = 0$
2. **Independent increments**: $B_t - B_s$ is independent of $B_u - B_v$ for non-overlapping intervals
3. **Gaussian increments**: $B_t - B_s \sim N(0, t - s)$
4. **Continuous paths**: $t \mapsto B_t$ is continuous (but nowhere differentiable!)

Properties 1–3 come directly from the discrete walk (independence of steps + CLT). Property 4 is a feature of the continuous limit.

**Nowhere differentiable** means: you cannot write $\frac{dB_t}{dt}$ in the ordinary sense. The paths are too jagged. This is why we'll need a new calculus (Itô's).

---

## Part 3: Building $dX_t = \mu\,dt + \sigma\,dB_t$ Step by Step

This part combines the decomposition from Part 1.2 with the scaling from Part 2.2. We go slowly.

### 3.1 The Drift Part: Where $\mu\,\Delta t$ Comes From

Go back to a general (possibly asymmetric) random walk on a time grid with step size $\Delta t$. Each step $\Delta X_i$ has some expected value $E[\Delta X_i]$.

**Problem:** $E[\Delta X_i]$ depends on the grid. Example — suppose a process drifts 12 units per year total:

| Grid | $\Delta t$ | Steps per year | $E[\Delta X_i]$ |
|---|---|---|---|
| Monthly | 1/12 | 12 | 1.0 |
| Daily | 1/365 | 365 | 0.0329 |
| Hourly | 1/8760 | 8760 | 0.00137 |

Every row gives a different $E[\Delta X_i]$, but they all describe the same process.

**Fix:** Divide by $\Delta t$ to get a grid-independent **rate**:

$$\mu = \frac{E[\Delta X_i]}{\Delta t}$$

| Grid | $E[\Delta X_i]$ | $\Delta t$ | $\mu = E[\Delta X_i] / \Delta t$ |
|---|---|---|---|
| Monthly | 1.0 | 1/12 | 12 per year |
| Daily | 0.0329 | 1/365 | 12 per year |
| Hourly | 0.00137 | 1/8760 | 12 per year |

$\mu$ is the same regardless of the grid. It's a property of the *process*, not the grid.

Rearranging: $E[\Delta X_i] = \mu \cdot \Delta t$.

**Intuition:** $\mu$ is like speed. Speed × time = distance. $\mu \times \Delta t$ = expected movement per step.

### 3.2 The Noise Part: Where $\sigma\,\Delta B_i$ Comes From

After removing the expected part, the leftover noise is:

$$\text{noise}_i = \Delta X_i - E[\Delta X_i]$$

This has mean zero and some variance $\text{Var}(\Delta X_i)$.

**Same problem:** $\text{Var}(\Delta X_i)$ depends on the grid too. Example — suppose the total variance over 1 year is 5.0:

| Grid | $\Delta t$ | Steps | $\text{Var}(\Delta X_i)$ |
|---|---|---|---|
| 10 steps | 0.1 | 10 | 0.50 |
| 1000 steps | 0.001 | 1000 | 0.005 |

Different number per step, same process.

**Same fix:** Divide by $\Delta t$:

$$\sigma^2 = \frac{\text{Var}(\Delta X_i)}{\Delta t}$$

| Grid | $\text{Var}(\Delta X_i)$ | $\Delta t$ | $\sigma^2 = \text{Var}/\Delta t$ |
|---|---|---|---|
| 10 steps | 0.50 | 0.1 | 5.0 per year |
| 1000 steps | 0.005 | 0.001 | 5.0 per year |

$\sigma^2$ is grid-independent. Rearranging: $\text{Var}(\text{noise}_i) = \sigma^2 \Delta t$.

### 3.3 Factoring the Noise: Defining $\Delta B_i$

We have a noise term with variance $\sigma^2 \Delta t$. We want to separate it into:

- $\sigma$ — the grid-independent noise rate (a constant)
- something with variance $\Delta t$ — the grid-dependent random piece

**Define:**

$$\Delta B_i \equiv \pm\sqrt{\Delta t} \quad \text{with equal probability}$$

This is nothing more than a name. It's a coin flip that gives $+\sqrt{\Delta t}$ or $-\sqrt{\Delta t}$.

Check its properties:

- $E[\Delta B_i] = 0$ ✓
- $E[(\Delta B_i)^2] = (\sqrt{\Delta t})^2 = \Delta t$
- $\text{Var}(\Delta B_i) = \Delta t$ ✓

Now write:

$$\text{noise}_i = \sigma \cdot \Delta B_i$$

Check: $\text{Var}(\sigma \cdot \Delta B_i) = \sigma^2 \cdot \text{Var}(\Delta B_i) = \sigma^2 \cdot \Delta t$ ✓

The variance matches. We've factored the noise into a constant ($\sigma$) times a standard random step ($\Delta B_i$).

### 3.4 Putting It All Together

From Part 1.2, any step decomposes as:

$$\Delta X_i = \text{expected part} + \text{noise part}$$

From Section 3.1: $\text{expected part} = \mu \cdot \Delta t$

From Section 3.3: $\text{noise part} = \sigma \cdot \Delta B_i$

So:

$$\boxed{\Delta X_i = \mu \cdot \Delta t + \sigma \cdot \Delta B_i}$$

Nothing here was assumed or introduced from outside. Every piece was built:

- $\mu$ came from making the drift rate grid-independent (Section 3.1)
- $\sigma$ came from making the noise rate grid-independent (Section 3.2)
- $\Delta B_i$ is just a name for $\pm\sqrt{\Delta t}$ (Section 3.3)

### 3.5 From Discrete Sum to Integral Form

Before taking any limits, let's work purely with the discrete equation. We have:

$$\Delta X_i = \mu \, \Delta t + \sigma \, \Delta B_i$$

This holds for every step $i = 1, 2, \ldots, n$. **Sum both sides** from step 1 to step $n$:

$$\sum_{i=1}^{n} \Delta X_i = \sum_{i=1}^{n} \mu \, \Delta t + \sum_{i=1}^{n} \sigma \, \Delta B_i$$

**Left side telescopes:**

$$\sum_{i=1}^{n} \Delta X_i = (X_1 - X_0) + (X_2 - X_1) + \cdots + (X_n - X_{n-1}) = X_n - X_0$$

All the middle terms cancel.

**Right side, first sum:**

$$\sum_{i=1}^{n} \mu \, \Delta t = \mu \sum_{i=1}^{n} \Delta t = \mu \cdot T$$

(Just a constant times the total time.)

**Right side, second sum:**

$$\sum_{i=1}^{n} \sigma \, \Delta B_i = \sigma \sum_{i=1}^{n} \Delta B_i$$

(A constant $\sigma$ times the accumulated random steps.)

**Put it together (still discrete, no limits taken):**

$$X_n = X_0 + \mu \cdot T + \sigma \sum_{i=1}^{n} \Delta B_i$$

**Now take $\Delta t \to 0$:** As the grid gets infinitely fine:

- $X_n$ becomes $X_t$
- $\mu \cdot T$ becomes $\int_0^t \mu \, ds$ (an ordinary Riemann integral, just equals $\mu \cdot t$)
- $\sigma \sum \Delta B_i$ becomes $\int_0^t \sigma \, dB_s$ (this limit is called the **Itô integral** — it's literally defined as the limit of these discrete sums)

$$\boxed{X_t = X_0 + \int_0^t \mu \, ds + \int_0^t \sigma \, dB_s}$$

**Why this is called "the rigorous definition":** The shorthand $dX_t = \mu\,dt + \sigma\,dB_t$ looks like an equation involving infinitely small quantities, but $dX_t$, $dt$, $dB_t$ are not actual numbers. You can't compute with them directly. You also can't divide by $dt$ to get $dB_t/dt$, because Brownian motion is nowhere differentiable (Section 2.3). So $dX_t = \mu\,dt + \sigma\,dB_t$ is **shorthand notation** that is *defined to mean* the integral equation above. Whenever you need to prove something or compute, you go back to the integral form.

### 3.6 The Crucial Scaling: Why $(dB_t)^2 = dt$

Look at $\Delta B_i = \pm\sqrt{\Delta t}$. Square it:

$$(\Delta B_i)^2 = (\pm\sqrt{\Delta t})^2 = \Delta t$$

This is **exact**, not approximate. Whether the coin lands $+\sqrt{\Delta t}$ or $-\sqrt{\Delta t}$, the square is $\Delta t$ either way. There's no randomness left after squaring!

Sum over all steps:

$$\sum_{i=1}^{n} (\Delta B_i)^2 = \sum_{i=1}^{n} \Delta t = n \cdot \Delta t = T$$

This sum is **deterministic** — it equals $T$ with certainty, no matter what the coin flips were.

In the continuous limit, this becomes the **quadratic variation** of Brownian motion:

$$\langle B \rangle_t = t \qquad \text{or informally written as} \qquad (dB_t)^2 = dt$$

This identity is the engine of Itô calculus. In ordinary calculus, $(dx)^2 = 0$ so you ignore it. In stochastic calculus, $(dB_t)^2 = dt \neq 0$, and this changes everything.

### 3.7 The Full Multiplication Table

From the scaling $dB_t \sim \sqrt{dt}$ and the deterministic nature of $dt$:

| × | $dt$ | $dB_t$ |
|---|------|--------|
| $dt$ | $0$ | $0$ |
| $dB_t$ | $0$ | $dt$ |

**Why $(dt)^2 = 0$:** It's second order in the infinitesimal, same as ordinary calculus.

**Why $dt \cdot dB_t = 0$:** Since $dB_t$ is of order $\sqrt{dt}$, the product $dt \cdot dB_t$ is of order $dt \cdot \sqrt{dt} = (dt)^{3/2}$, which vanishes faster than $dt$.

**Why $(dB_t)^2 = dt$:** This is what we just proved in Section 3.6.

---

## Part 4: Itô's Lemma — The Chain Rule for Stochastic Calculus

### 4.1 What Is an "Itô Process"?

An Itô process is just a name for any process that has the form we built in Part 3:

$$dX_t = a(X_t, t)\,dt + b(X_t, t)\,dB_t$$

That's it. Anything with a "drift part times $dt$" plus a "noise part times $dB_t$" is an Itô process. The drift $a$ and diffusion $b$ can be constants or can depend on $X_t$ and $t$.

**You already know several Itô processes:**

- **From Part 3:** $dX_t = \mu\,dt + \sigma\,dB_t$. Here $a = \mu$ and $b = \sigma$ (both constants).
- **Brownian motion itself:** $dB_t = 0 \cdot dt + 1 \cdot dB_t$. Here $a = 0$, $b = 1$.
- **The stock model (Part 5, coming soon):** $dS_t = \mu S_t\,dt + \sigma S_t\,dB_t$. Here $a = \mu S_t$ and $b = \sigma S_t$ (both depend on $S_t$).

The term sounds fancy but it's literally just the name for the structure we spent all of Part 3 building.

### 4.2 Why the Ordinary Chain Rule Fails

Suppose $X_t$ is an Itô process and we want to know how $f(X_t)$ evolves. In ordinary calculus:

$$df = f'(X) \, dX$$

But this ignores $(dX)^2$. In ordinary calculus, $(dx)^2 = 0$ so this is fine. But now $(dX_t)^2$ contains a $(dB_t)^2 = dt$ term, which is **not** negligible. We must keep it.

### 4.3 Deriving Itô's Lemma

Let $X_t$ be a general Itô process: $dX_t = a\,dt + b\,dB_t$, where $a$ and $b$ can be anything (constants, or functions of $X_t$ and $t$).

Start with a Taylor expansion of $f(X_t, t)$ to second order:

$$df = \frac{\partial f}{\partial t} dt + \frac{\partial f}{\partial X} dX_t + \frac{1}{2} \frac{\partial^2 f}{\partial X^2} (dX_t)^2 + \text{higher order terms}$$

(The higher order terms include $dX_t \cdot dt$, $(dt)^2$, etc. — we'll show these all vanish.)

Now substitute $dX_t = a\,dt + b\,dB_t$ and compute $(dX_t)^2$:

$$(dX_t)^2 = (a\,dt + b\,dB_t)^2$$

Expand:

$$= a^2 (dt)^2 + 2ab\,dt\,dB_t + b^2 (dB_t)^2$$

Apply the multiplication table from Section 3.7:

- $a^2 (dt)^2 = 0$
- $2ab\,dt\,dB_t = 0$
- $b^2 (dB_t)^2 = b^2\,dt$ ← **this survives!**

So $(dX_t)^2 = b^2\,dt$.

The other higher-order terms vanish too:

- $dX_t \cdot dt = (a\,dt + b\,dB_t) \cdot dt = a(dt)^2 + b\,dB_t\,dt = 0$
- $(dt)^2 = 0$

Substituting back into the Taylor expansion:

$$\boxed{df = \frac{\partial f}{\partial t}\,dt + \frac{\partial f}{\partial X}\,dX_t + \frac{1}{2}\,b^2\,\frac{\partial^2 f}{\partial X^2}\,dt}$$

where $b$ is whatever multiplies $dB_t$ in the SDE for $X_t$.

Or, expanding $dX_t$:

$$df = \left(\frac{\partial f}{\partial t} + a\frac{\partial f}{\partial X} + \frac{1}{2}\,b^2\,\frac{\partial^2 f}{\partial X^2}\right) dt + b\,\frac{\partial f}{\partial X}\,dB_t$$

**This is Itô's Lemma.** The extra $\frac{1}{2}b^2 f''$ term (the "Itô correction") exists because $(dB_t)^2 = dt \neq 0$. In Part 3's simple case, $b = \sigma$ (a constant). In Part 5, $b = \sigma S_t$. The formula is the same — you just plug in whatever $b$ is.

### 4.4 Sanity Check: When $b = 0$

If there's no randomness, the Itô correction vanishes and we recover the ordinary chain rule:

$$df = \frac{\partial f}{\partial t} dt + \frac{\partial f}{\partial X} dX_t$$

---

## Part 5: Geometric Brownian Motion — Modeling Stock Prices

### 5.1 Why Not Use $dS = \mu \, dt + \sigma \, dB_t$ Directly?

If the stock price followed $dS = \mu \, dt + \sigma \, dB_t$, then $S_t$ could go negative — nonsensical for a price.

Instead, we model the **percentage** changes as having constant drift and volatility:

$$\frac{dS_t}{S_t} = \mu \, dt + \sigma \, dB_t$$

or equivalently:

$$dS_t = \mu S_t \, dt + \sigma S_t \, dB_t$$

This is an Itô process where the drift and diffusion coefficients are both proportional to $S_t$.

### 5.2 Solving the GBM SDE Using Itô's Lemma

We want to find $S_t$ explicitly. The trick: apply Itô's lemma to $f(S) = \ln S$.

**First, identify the pieces.** Our SDE is:

$$dS_t = \underbrace{\mu S_t}_{a} \, dt + \underbrace{\sigma S_t}_{b} \, dB_t$$

Comparing with the general form $dX_t = a\,dt + b\,dB_t$ from Part 4, we have $a = \mu S_t$ and $b = \sigma S_t$. The key thing: $b$ is not just $\sigma$ — it's $\sigma S_t$.

**Second, recall Itô's lemma from Part 4 (Section 4.3):**

$$df = \frac{\partial f}{\partial t}\,dt + \frac{\partial f}{\partial X}\,dX_t + \frac{1}{2}\,b^2\,\frac{\partial^2 f}{\partial X^2}\,dt$$

We're using this exact formula. We just need to plug in our specific $f$ and $b$.

**Third, compute everything we need for $f(S) = \ln S$:**

- $\frac{\partial f}{\partial t} = 0$ — because $\ln(S)$ only depends on $S$, not on $t$. If you change $t$ while keeping $S$ fixed, $\ln(S)$ doesn't change. (This would NOT be zero for a function like $S^2 e^{-t}$, which explicitly depends on $t$.)
- $\frac{\partial f}{\partial S} = \frac{1}{S}$
- $\frac{\partial^2 f}{\partial S^2} = -\frac{1}{S^2}$
- $b = \sigma S_t$ (the diffusion coefficient from our SDE)

**Fourth, plug into Itô's lemma:**

$$d(\ln S_t) = \underbrace{0}_{\partial f/\partial t}\,dt + \underbrace{\frac{1}{S_t}}_{\partial f/\partial S}\,dS_t + \frac{1}{2}\underbrace{(\sigma S_t)^2}_{b^2} \cdot \underbrace{\left(-\frac{1}{S_t^2}\right)}_{\partial^2 f/\partial S^2}\,dt$$

Simplify the last term: $\frac{1}{2}(\sigma S_t)^2 \cdot (-\frac{1}{S_t^2}) = -\frac{1}{2}\sigma^2$. So:

$$d(\ln S_t) = \frac{1}{S_t} dS_t - \frac{1}{2}\sigma^2 \, dt$$

Now substitute $dS_t = \mu S_t \, dt + \sigma S_t \, dB_t$:

$$= \frac{1}{S_t}(\mu S_t \, dt + \sigma S_t \, dB_t) - \frac{1}{2}\sigma^2 \, dt$$

$$= \mu \, dt + \sigma \, dB_t - \frac{1}{2}\sigma^2 \, dt$$

$$= \left(\mu - \frac{\sigma^2}{2}\right) dt + \sigma \, dB_t$$

**The $-\frac{\sigma^2}{2}$ is the Itô correction.** It came from the $\frac{1}{2}b^2 f''$ term — the same term we derived in Part 4 because $(dB_t)^2 = dt \neq 0$.

**Fifth, integrate both sides** (same technique as Section 3.5 — sum the discrete version, then take the limit):

$$\ln S_t - \ln S_0 = \left(\mu - \frac{\sigma^2}{2}\right)t + \sigma B_t$$

Exponentiate:

$$\boxed{S_t = S_0 \exp\left[\left(\mu - \frac{\sigma^2}{2}\right)t + \sigma B_t\right]}$$

Since $B_t \sim N(0, t)$, and looking at the equation above, $\ln S_t$ equals a **deterministic** (non-random) part plus a normal random variable. At any fixed time $t$, a deterministic number plus a normal random variable is normal. So $\ln S_t$ is normally distributed, which means $S_t$ is **log-normally** distributed — always positive (since $e^{\text{anything}} > 0$).

---

## Part 6: Deriving the Black-Scholes Equation

Parts 1–5 were pure math. Now we enter finance. This section requires some financial concepts and a new type of equation (PDEs), so we explain everything before using it.

### 6.1 Financial Background: What You Need to Know

**What is a stock?** A share of ownership in a company. Its price $S_t$ changes randomly over time. We modeled it in Part 5: $dS_t = \mu S_t\,dt + \sigma S_t\,dB_t$.

**What is a bond (risk-free asset)?** Imagine putting money in a perfectly safe bank account. If you deposit \$100 at annual interest rate $r = 5\%$, after one year you have \$105. There's no randomness — you're guaranteed the return.

Mathematically, if $M_t$ is the value of your bank account:

$$dM_t = r M_t \, dt$$

Now lets derive the M_t form with ODE:
1. let us first seperate the terms, because M_t and t are seperable. 
$$\frac{dM_t}{M_t} = rdt$$
2. take integral of both sides we have 
$$\int_{0}^{t} \frac{1}{M_t} dM_t = \int_{0}^{t}r dt$$
$$ln(M_t) - ln(M_0) = rt + C$$
$$\frac{M_t}{M_0} = e^{rt} + C$$
$$M_t = M_0e^{rt}+C$$

Notice: no $dB_t$ term. No randomness at all. This is just an ordinary differential equation. You can solve it:

$$M_t = M_0 \, e^{rt}$$

Money grows exponentially at rate $r$. The reverse is also useful: \$1 received at future time $T$ is worth $e^{-r(T-t)}$ today. This is called **discounting**.

**What is an option?** A contract that gives you the *right* (not obligation) to buy or sell a stock at a pre-agreed price.

A **European call option** says: "On date $T$ (the expiry), you can buy one share of stock at price $K$ (the strike price)."

At expiry, what's this worth?

- If $S_T = 150$ and $K = 100$: you buy at 100, sell at 150, profit = 50. The option is worth $150 - 100 = 50$.
- If $S_T = 80$ and $K = 100$: you'd be buying at 100 something worth 80. You don't exercise. The option is worth 0.

So at expiry:

$$V(S_T, T) = \max(S_T - K, \; 0)$$

**The big question:** The option's value at expiry is easy — it's just $\max(S_T - K, 0)$. But what is it worth *today*, at time $t < T$, when $S_T$ is still unknown? That's what Black-Scholes answers.

**Why $V$ depends on both $S$ and $t$:** The option's value depends on the current stock price $S$ (a higher stock price means you're more likely to profit). But it also depends on time $t$ — an option with 6 months left is worth more than one with 1 day left, because there's more time for the stock to move favorably. So $V = V(S, t)$ is a function of two variables.

**What is a portfolio?** Just a collection of financial things you hold. "I own 10 shares of stock and 3 options" is a portfolio. The total value is:

$$\Pi = 10 \cdot S + 3 \cdot V$$

**What is arbitrage?** Free money with no risk. Example: if two stores sell the same item for different prices, you buy from the cheap one and sell to the expensive one — guaranteed profit. In financial markets, we assume **no arbitrage exists** — any opportunity for free money gets exploited instantly until it disappears. This is not a mathematical theorem; it's an economic assumption.

**What does no-arbitrage imply?** If you have a portfolio with absolutely no randomness (no $dB_t$ term), it **must** earn the same return as the risk-free bank account. Why?

- If it earned more than $r$: borrow money from the bank at rate $r$, invest in the portfolio, pocket the difference. Free money. Arbitrage!
- If it earned less than $r$: sell (short) the portfolio, put the money in the bank. Free money. Arbitrage!

So: **riskless portfolio must earn rate $r$**. This simple idea drives the entire derivation.

### 6.2 Math Background: Ordinary vs Partial Differential Equations

You know ordinary differential equations (ODEs). These have one independent variable:

$$\frac{dy}{dx} = 3x + y$$

Here $y$ depends on one thing: $x$.

A **partial differential equation (PDE)** is the same idea but the unknown function depends on **multiple variables**. For example, $V(S, t)$ depends on both $S$ and $t$. So when we differentiate, we need to say *which variable* we're differentiating with respect to:

- $\frac{\partial V}{\partial t}$ = how $V$ changes when $t$ changes but $S$ is held fixed
- $\frac{\partial V}{\partial S}$ = how $V$ changes when $S$ changes but $t$ is held fixed
- $\frac{\partial^2 V}{\partial S^2}$ = the second derivative with respect to $S$

The symbol $\partial$ (instead of $d$) just means "I'm holding the other variables fixed."

**Example:** Let $V(S, t) = S^2 e^{-t}$.

- $\frac{\partial V}{\partial t} = -S^2 e^{-t}$ (differentiate with respect to $t$, treat $S$ as a constant)
- $\frac{\partial V}{\partial S} = 2S e^{-t}$ (differentiate with respect to $S$, treat $t$ as a constant)
- $\frac{\partial^2 V}{\partial S^2} = 2e^{-t}$ (differentiate $\frac{\partial V}{\partial S}$ with respect to $S$ again)

A PDE looks like an ODE but involves these partial derivatives. The Black-Scholes PDE will relate $\frac{\partial V}{\partial t}$, $\frac{\partial V}{\partial S}$, and $\frac{\partial^2 V}{\partial S^2}$ to each other.

### 6.3 Apply Itô's Lemma to $V(S_t, t)$

We use Itô's lemma from Part 4 (Section 4.3) again. Our SDE is the same as Part 5: $dS_t = \mu S_t\,dt + \sigma S_t\,dB_t$, so $b = \sigma S_t$.

But now our function is $V(S, t)$ instead of $\ln(S)$. The key difference: **$V$ depends on both $S$ and $t$**, because an option's value changes as expiry approaches even if the stock price stays still. So $\frac{\partial V}{\partial t} \neq 0$ (unlike Part 5 where $\frac{\partial(\ln S)}{\partial t} = 0$).

Applying Itô's lemma with $b = \sigma S$:

$$dV = \frac{\partial V}{\partial t} dt + \frac{\partial V}{\partial S} dS_t + \frac{1}{2} (\sigma S)^2 \frac{\partial^2 V}{\partial S^2} dt$$

Substitute $dS_t = \mu S\,dt + \sigma S\,dB_t$ into the middle term:

$$dV = \left(\frac{\partial V}{\partial t} + \mu S \frac{\partial V}{\partial S} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2}\right) dt + \sigma S \frac{\partial V}{\partial S} \, dB_t$$

Notice the structure: a deterministic drift part (everything multiplying $dt$), plus a random part (everything multiplying $dB_t$).

### 6.4 Construct a Riskless Portfolio (Delta Hedging)

Here's Black and Scholes' brilliant idea: build a portfolio where the randomness cancels out.

**The portfolio:** Hold one option and sell ("short") $\Delta$ shares of stock. The value is:

$$\Pi = V - \Delta \cdot S$$

(Why subtract? "Shorting" means you sell something you've borrowed — you owe it back later. It's like a negative holding.)

**How the portfolio changes:** Both $V$ and $S$ move, so:

$$d\Pi = dV - \Delta \, dS$$

Substitute the expressions for $dV$ (from Section 6.3) and $dS = \mu S\,dt + \sigma S\,dB_t$:

$$d\Pi = \left(\frac{\partial V}{\partial t} + \mu S \frac{\partial V}{\partial S} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2}\right) dt + \sigma S \frac{\partial V}{\partial S} \, dB_t - \Delta\left(\mu S \, dt + \sigma S \, dB_t\right)$$

Now collect the $dB_t$ terms (the random parts):

$$d\Pi = (\ldots)\,dt + \sigma S\left(\frac{\partial V}{\partial S} - \Delta\right) dB_t$$

The randomness is entirely in that $dB_t$ term. We can make it vanish by choosing:

$$\Delta = \frac{\partial V}{\partial S}$$

With this choice, the $dB_t$ term is zero:

$$d\Pi = \left(\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2}\right) dt$$

**The portfolio is now entirely deterministic.** No $dB_t$, no randomness. We know exactly how it changes. This technique is called **delta hedging** — you hold exactly $\frac{\partial V}{\partial S}$ shares of stock to cancel out the option's random fluctuations.

### 6.5 Apply No-Arbitrage

From Section 6.1: a riskless portfolio must earn the risk-free rate $r$. So:

$$d\Pi = r\Pi \, dt$$

Substitute $\Pi = V - \frac{\partial V}{\partial S} S$:

$$d\Pi = r\left(V - \frac{\partial V}{\partial S} S\right) dt$$

Now we have two expressions for $d\Pi$. Set them equal:

$$\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} = r V - rS\frac{\partial V}{\partial S}$$

Rearranging (move everything to one side):

$$\boxed{\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} + rS\frac{\partial V}{\partial S} - rV = 0}$$

**This is the Black-Scholes PDE.** It's a partial differential equation (Section 6.2) that the option price $V(S, t)$ must satisfy.

**Critical observation:** $\mu$ has completely disappeared! The stock's expected return doesn't affect the option price. Only $\sigma$ (volatility) and $r$ (risk-free rate) matter. Intuitively: the delta hedge cancels out all exposure to the stock's direction, so it doesn't matter whether the stock is expected to go up or down. This is the mathematical foundation of **risk-neutral pricing**.

---

## Part 7: Solving the PDE — The Black-Scholes Formula

### 7.1 Boundary Condition

From Section 6.1, at expiry ($t = T$) a European call option with strike $K$ is worth:

$$V(S, T) = \max(S - K, \, 0)$$

This is the **boundary condition** for our PDE — we know $V$ at the end, and we want to work backwards to find $V$ at earlier times.

### 7.2 Transformation to the Heat Equation

The Black-Scholes PDE has variable coefficients (the $S^2$ multiplying $V_{SS}$). We transform it to the standard heat equation.

**Step A: Log-coordinates and reverse time.**

Let:

$$x = \ln(S/K), \qquad \tau = \frac{1}{2}\sigma^2(T - t), \qquad V(S, t) = K \cdot v(x, \tau)$$

After computing all partial derivatives via chain rule and substituting into the PDE:

$$\frac{\partial v}{\partial \tau} = \frac{\partial^2 v}{\partial x^2} + (k - 1)\frac{\partial v}{\partial x} - kv$$

where $k = \frac{2r}{\sigma^2}$.

**Step B: Remove the lower-order terms.**

Let:

$$v(x, \tau) = e^{\alpha x + \beta \tau} \, u(x, \tau)$$

Choosing:

$$\alpha = -\frac{1}{2}(k - 1), \qquad \beta = -\frac{1}{4}(k + 1)^2$$

all first-order and zeroth-order terms cancel, leaving:

$$\frac{\partial u}{\partial \tau} = \frac{\partial^2 u}{\partial x^2}$$

**This is the heat equation** — the most classical PDE in mathematics.

### 7.3 Solving with the Gaussian Kernel

The heat equation with initial condition $u(x, 0) = u_0(x)$ has the solution:

$$u(x, \tau) = \frac{1}{2\sqrt{\pi \tau}} \int_{-\infty}^{\infty} u_0(y) \, e^{-(x - y)^2 / (4\tau)} \, dy$$

After computing this integral (it splits into two Gaussian integrals that evaluate to standard normal CDFs) and undoing all substitutions, we get:

$$\boxed{C(S, t) = S \, N(d_1) - K e^{-r(T-t)} N(d_2)}$$

where:

$$d_1 = \frac{\ln(S/K) + \left(r + \frac{\sigma^2}{2}\right)(T - t)}{\sigma\sqrt{T - t}}$$

$$d_2 = d_1 - \sigma\sqrt{T - t} = \frac{\ln(S/K) + \left(r - \frac{\sigma^2}{2}\right)(T - t)}{\sigma\sqrt{T - t}}$$

and $N(\cdot)$ is the CDF of the standard normal distribution.

---

## Part 8: Understanding the Formula

### 8.1 What Each Piece Means

| Term | Meaning |
|---|---|
| $S \cdot N(d_1)$ | Present value of receiving the stock, weighted by risk-neutral probability. Also the **delta** (hedge ratio). |
| $Ke^{-r(T-t)}$ | Present value of the strike price, discounted at the risk-free rate. |
| $N(d_2)$ | Risk-neutral probability that the option expires in the money ($S_T > K$). |
| $\sigma\sqrt{T-t}$ | Total uncertainty over the option's life — echoing the $\sqrt{\Delta t}$ scaling from Part 2. |

### 8.2 The $\sigma^2/2$ in $d_1$ vs $d_2$

From Part 5, Itô's lemma applied to $\ln S$ produced a $-\sigma^2/2$ correction. That same correction appears here:

- $d_1$ uses $(r + \sigma^2/2)$: relates to $E[S_T]$ under the risk-neutral measure
- $d_2$ uses $(r - \sigma^2/2)$: relates to the log-mean of $S_T$

The difference $d_1 - d_2 = \sigma\sqrt{T-t}$ is purely the Itô correction.

### 8.3 Put-Call Parity

For a European put:

$$P(S, t) = Ke^{-r(T-t)}N(-d_2) - S \, N(-d_1)$$

From the no-arbitrage relation:

$$C - P = S - Ke^{-r(T-t)}$$

---

## Part 9: The Complete Logical Chain

```
Part 1: Discrete random walk (coin flips, i.i.d. steps)
   │    Can decompose each step into: expected part + surprise
   │
   │    Part 2: Put walk on time grid. CLT forces step size δ = σ√(Δt).
   │    Limit is Brownian motion Bₜ.
   ▼
Part 3: Extract grid-independent rates:
   │    • μ = E[step] / Δt  (drift rate)
   │    • σ² = Var(step) / Δt  (noise rate)
   │    • Define ΔBᵢ = ±√(Δt) so noise = σ · ΔBᵢ
   │    ⟹ ΔXᵢ = μΔt + σΔBᵢ  →  dXₜ = μdt + σdBₜ
   │
   │    Key identity: (ΔBᵢ)² = Δt exactly  →  (dBₜ)² = dt
   ▼
Part 4: Itô's Lemma (corrected chain rule)
   │    Because (dBₜ)² = dt ≠ 0, Taylor expansion keeps ½σ²f'' term
   │
   │    Part 5: Model stock as GBM: dS/S = μdt + σdBₜ
   │    Solve using Itô's lemma on ln(S)
   ▼
Part 6: Apply Itô's lemma to option V(S,t)
   │    Delta hedge to kill dBₜ term → portfolio becomes riskless
   │    No-arbitrage → must earn rate r → Black-Scholes PDE
   │    (μ disappears!)
   ▼
Part 7: Transform PDE to heat equation → solve
   ▼
Black-Scholes Formula: C = S·N(d₁) - Ke^{-r(T-t)}·N(d₂)
```
