---
title: MAT1750 Project Report
author: Caroline Lin and Dmitry Paramonov
header-includes: |
    \usepackage{gensymb}
    \usepackage{float}
    \floatplacement{figure}{H}
geometry: margin=1in
---

\newcommand{\R}[0]{\mathbb{R}}
\newcommand{\N}[0]{\mathbb{N}}


## Background info.

### Laplace's equation.

Everyone knows about Laplace's equation. It's really fundamental. If you really need me to give you the equation for it, then I can give it to you:

Here is an equation.

Sometimes people shy away from explicitly defining the Laplacian, because it turns out you can define it in odd coordinate systems, and then it looks different. But if you are boring like me and coordinate-dependent then you can use this Cartesian representation of Laplace's:

Here is a second equation.

We can use it to model steady-states for heat diffusion, but the Laplacian comes up in a lot of contexts. In first-year PDEs we all learn how to solve Laplace's equation with Green's functions, but Green's functions are not practical to compute, especially when domain boundaries are not very good.

There is a full suite of numerical methods that people have come up with to deal with the problem of solving Laplace's equation with irregular boundaries. We aren't really going to talk about that here.

Instead, what we'll talk about is a slightly different way of solving Laplace's equation, based on the theory of stochastic differential equations.

### The Feynman-Kac formula.

Okay, so if you have a PhD in statistical physics then you probably know about the Feynman-Kac formula. I don't have a PhD in statistical physics, but I have heard about the Feynman-Kac formula. Here is the Feynman-Kac formula.

A formula.

If we look at it, we see that we can get evaluate the right-hand side of the formula to solve the Laplace's equation at a single point. This is pretty good!

Obviously the expectation of a stochastic process is hard to compute. However, it is pretty easy to simulate. We can use our trusty pseudo-random number generator, courtesy of the good folks at scipy, to simulate a random walk.

If we simulate enough random walks then we will get a reasonable approximation to the solution. How reasonable is the approximation? By the law of large numbers, the sampling will approach the true expectation at a rate of 1 / sqrt(N), where N is the number of samples.

It turns out that other people have written about theoretical aspects of this exact problem, so we refer you to these papers (Random Walk and the Heat
Equation, Gregory Lawler, and A Proof of the Random-Walk Method for Solving Laplace's Equation in 2-D, J. F. Reynolds) for further reading.

## Some previous work on the problem.

Generally speaking Laplace's problem has been tackled from a number of angles. For example, we can use finite elements to solve Laplace's. Here is a link to one such program (fenics).

We based a significant amount of our work on Chati et al., who did a random walk method for Laplace's and Poisson's. A quick search of related work turns up a body of "SDE" or stochastic solution methods. Here are two titles I happened to understand: "Solution of BVPs in electrodynamics by stochastic methods", and "Random Walk Method for Potential Problems".

We should talk a little more about coupling. We decided to introduce coupling in our work, based on discussions with our professor (that's you, the reader). Our coupling approach is based on a standard finite-difference method. More on this in discussion.

There's also some cost-savings work in walk-on-spheres, but we didn't look into this.

## How to implement a random walk.

### The Basic Method

For the random walk method,
we first begin by choosing the points we wish to simulate.
For each of those points, we then start simulating random walks,
by taking steps with Gaussian random noise.
Each random walk runs until it hits the boundary.
When it does, we evaluate the function on that boundary.

For each point, we perform $K$ random walks.
Then, we average the values of the function evaluations
for random walks from that point.
This averaging is the discrete equivalent of the expectation value above,
so in expectation, this gets us the value of the solution at that point.
Thus, by doing sufficiently many random walks,
we get a good unbiased estimator for the solution of the Laplacian.
And with sufficiently many random walks,
we closely can approximate the actual solution.

### The Finite Difference Modification

However, the above method has several flaws.
Notably, a random walk takes a lot of time to converge.
Furthermore, work is not reusable between points.
After solving the problem at many points,
solving the problem an another point
requires running many random walks all over again.

To avoid this, we can use finite differences
to vastly cut down on the amount of random walks we have to perform.

We begin by constructing a mesh of points within our desired boundary,
with distances between points equal to a parameter $h$.
Each point will then have 4 neighbours, the points closest to it,
one in each direction.

We classify the points within the mesh into two groups:
*interior points* and *boundary points*.
Interior points are those with four neighbours which lie within the boundary,
while boundary points are those that lie within the boundary
but which have a neighbour which lies outside the boundary.

Now, when we perform our random walks,
we only perform them on the boundary points.
We can thus obtain a good estimate
for the values of the function at the boundary points.

For the interior points, we can then just use finite differences.
We will have a sparse system,
where each interior point is entangled with at most 4 other interior points,
as well as at most 4 boundary points,
which contribute to the target values.
This lets us find all the values of the interior points
with a single linear solve,
allowing us to neatly avoid
having to perform random walks on all of those points.

The benefits of this are two-fold.
First of all, note that this reduces, asymptotically,
the amount of random walks we have to perform.
If our grid were to have spacing $h$ between the nodes,
then we would originally need to perform $O(Kh^{-2})$ random walks
(assuming a nice boundary),
since the amount of points grew as $h^{-2}$.
But by only performing the walks on the boundary points,
we cut this down to $O(Kh^{-1})$,
which can lead to huge savings.

Furthermore, note that boundary points are within $h$ of the boundary,
just by their definition.
Thus, the amount of steps necessary for a random walk to converge
should take less time for boundary points than interior points.
So many walks will converge much faster.

### The coupled set of equations.

However, the random walks can still take a while to converge.
If a random walk at a boundary point,
it can still randomly walk into the centre of the region,
and could thus take a while to actually reach the boundary.

To avoid this, we can stop our walks early if they reach the boundary.
If the random walk reaches a boundary point,
then the expectation value of that random walk
is the same as the expectation value of random walks from that boundary point.
Thus, there is no need to continue the random walks.

In practice, this can be implemented in two ways,
in an explicit or an implicit way.
We run the random walks until we hit the boundary or another boundary point.
We can then use the resulting numbers in two ways.

We can index a vector based on our boundary points,
creating a vector $\R^{B}$, where $B$ is the number of boundary points.
This lets us assign a real number to each boundary point.
Let $b_i$ be the total of the values of random walks starting at a given point $i$.

Normally, we then set $u = \frac{1}{K}b$.
This lets us solve the system directly.

Now, let $F\in\N^{B\times B}$ be a matrix of integers.
$F_{i,j}$ is the number of random walks
which start at the boundary point $i$ and end at the boundary point $j$.

The explicit way of finding $u$ is then
by taking $u=\frac{1}{K^2}Fb + \frac{1}{K}b$.
This lets us find $u$ explicitly.

We can also use an implicit form.
We can take $u=\frac{1}{K}Fu + \frac{1}{K}b$.
This requires solving a linear system,
and usually performs better than the above.

### The random walk, and boundary detection.

Overview.
The covariance matrix of the 2D gaussian is the identity. At each step we can just sample from this.

Here is some python. It's pretty simple because I abstracted away all of the details.

```python
```

#### Some details: boundary detection.

There are some subtle problems that arise from using coupled random walks.

First, we should notice that points can cross the boundary at any point between boundary points. In this case, which frequency should we update?

We decided to update the frequency of the closer of the two nodes. Something we didn't try: use fractional frequency updates based on a linear interpolant (similar to hard vs. soft k-means).

Second, random walks are discrete, so they can 'tunnel' through the boundary.

We decided to create a 'thickness' to the wall, but only on the interior-facing boundary. The thickness of the wall is such that 95% of all points starting from a boundary node will not tunnel past the wall (based on the standard deviation of the random walk).

We throw out points that go past the wall.

Something we didn't try: interpolating between successive point positions to check for a boundary collision. Based on our setup, it seems expensive.

### Other things we didn't try.

Another method for speeding up these random walks for Laplace's Problem
is to use *walk-on-spheres*.
Note that because the Laplacian is 0 on this region,
the individual steps of the random walk are identical,
except that they end up at different points.
And if we start at a given point,
if a ball of radius $R$ around that point fits within the boundary,
then a random walk has an equal probability
of ending on any point of that circle,
after every amount of steps.
Thus, we can skip a lot of steps,
and just jump to a random point on the ball.
This lets us cut down on the number of steps required,
thus vastly speeding up the algorithm.
However, we did not apply this,
as we wanted our algorithm to be more generalizable
to Poisson's Problem.

Another approach we did not use was the boundary-within-a-boundary approach.
Here, we construct another boundary on the inside of the boundary points.
Then, we can couple our boundary points with the points on this inner boundary,
which lets us cut the random walks off faster.
However, then we need to find the values on the inner boundary points.
This causes us to incur some costs to find those.

The end result of extending this is to fully couple all our points.
Now, we no longer even use our discretization of the Laplacian.
Instead, we can now fully couple both our boundary points and interior points,
which would then result in a much more dense linear system.
Solving this system would be more complicated than before,
but it would mean that we could take sufficiently fewer random walk steps.

## Performance

We ran several tests for the performance of this algorithm.

First, we considered the accuracy of this algorithm.
We ran tests on the example function $g(x,y)=x^3-3xy^2-3x^2y+y^3$.
This function satisfies Laplace's Equation,
meaning we have a target value for $u(x,y)$.

From this, we ran the above procedure to obtain values on the unit disk.
We could also calculate the values on that unit disk explictly.
From this, we could then calculate the maximum absolute error.

These runs did not use coupling,
and used constant steps for the random walks of length $0.001$.

The results are listed in the table below.
Each entry is the maximum absolute error obtained by running the algorithm
with that spacing, and with that amount of random walks.

$K \backslash h^{-1}$ | 25 | 50 | 100
--- | --- | --- | ---
25 | 0.14277 | 0.12212 | 0.11857
100 | 0.09681 | 0.07336 | 0.06178
400 | 0.03782 | 0.02709 | 0.03108
1600 | 0.01406 | 0.01657 | 0.01441

We also have the mean absolute error.

$K \backslash h^{-1}$ | 25 | 50 | 100
--- | --- | --- | ---
25 | 0.010997 | 0.006750 | 0.004682
100 | 0.011937 | 0.005103 | 0.003368
400 | 0.004551 | 0.001876 | 0.001779
1600 | 0.001479 | 0.001157 | 0.000858

In both cases, the error can be seen to decrease with $h$
and to decrease as $K$ increases.
However, the step from $25$ steps to $100$
doesn't seem to decrease the error significantly.
It appears that this number of steps just isn't significant for convergence.
The error is also linear in $h$,
with the exception of $K=400$ and $h^{-1}=100$.

We also analysed the runtime,
run on the same set of values of $K$ and $h$.
This is organized into three tables.
The first lists the runtime of the simulation steps,
the second lists the time to do the subsequent linear solve,
and the third lists the total runtime.
All times are listed in seconds.

$K \backslash h^{-1}$ | 25 | 50 | 100
--- | --- | --- | ---
25 | 1.454 | 1.963 | 4.004
100 | 6.268 | 8.477 | 25.817
400 | 23.231 | 35.812 | 68.748
1600 | 91.341 | 127.997 | 217.419

$K \backslash h^{-1}$ | 25 | 50 | 100
--- | --- | --- | ---
25 | 0.653 | 2.278 | 18.498
100 | 0.877 | 2.252 | 9.256
400 | 0.589 | 2.979 | 10.776
1600 | 0.628 | 2.415 | 10.337

$K \backslash h^{-1}$ | 25 | 50 | 100
--- | --- | --- | ---
25 | 2.107 | 4.240 | 22.503
100 | 7.146 | 10.728 | 35.0722
400 | 23.820 | 38.791 | 79.524
1600 | 91.969 | 130.412 | 227.759

The different parts scaled differently.
The linear system solve was independent of $K$, as expected,
since it happens after the simulation.
And it scaled roughly as $h^{-2}$.
This makes sense as we have $h^{-2}$ points overall,
and we are solving a sparse linear system.

The simulation steps scaled linearly in $K$, as expected.
It also scaled roughly linearly in $h^{-1}$.
This is as expected, as the number of boundary points
grows roughly linearly in $h^{-1}$.
Note that because the random walks take random amounts of time,
the resulting times are also random,
and have a pretty high variance.

We also ran some tests with a much more uneven domain.
This one was defined by a very heavily oscillatory set of sinusoidals,
resulting in a very uneven boundary.

The maximum errors of those runs are listed below.

$K \backslash h^{-1}$ | 25 | 50
--- | --- | ---
25 | 0.17144 | 0.14691
100 | 0.07071 | 0.08459
400 | 0.03723 | 0.03352

We also have the mean absolute error.

$K \backslash h^{-1}$ | 25 | 50
--- | --- | ---
25 | 0.024012 | 0.021978
100 | 0.013124 | 0.010827
400 | 0.007766 | 0.004791

As can be seen, while the mean errors did decrease slightly with $h^{-1}$,
there was not a large decrease in the error in $h^{-1}$.
Instead, the most notable decrease was with $K$,
where it appeared to decrease as $\frac{1}{\sqrt{K}}$.

We also analysed the runtime, just as before.

$K \backslash h^{-1}$ | 25 | 50
--- | --- | ---
25 | 31.307 | 51.238
100 | 144.183 | 216.124
400 | 556.048 | 910.414

$K \backslash h^{-1}$ | 25 | 50
--- | --- | ---
25 | 0.323 | 1.187
100 | 0.357 | 1.110
400 | 0.414 | 1.110

$K \backslash h^{-1}$ | 25 | 50
--- | --- | ---
25 | 31.630 | 52.425
100 | 144.539 | 217.234
400 | 556.462 | 911.524

Again, the linear solve grew as roughly $h^{-2}$,
which makes sense, just as before.
The runtime of the random walks also took time roughly proportional to $K$.

However, the random walks also took time proportional to less than $h^{-1}$.
And they took much more time than on the circular domain.
The probable reason for that is the unusual domain.
Near the centre, boundary points are very close to boundaries,
but most random walks do not quickly converge to the boundary.
Instead, because the boundary is a very thin slice,
most random walks take a lot of time.
This means that points near the centre take a long runtime,
meaning that the newer points don't add that much time,
in comparison.

## Some figures and results.

Conventions. In all of the following, $N$ is a parameter controlling the fineness of the finite-difference and boundary meshes, and $K$ is the number of random walks per point on the boundary, and $h$ is the distance between adjacent nodes in our mesh. We tested our code on the unit disk and on the level set of {#TODO: insert function name here}. For the unit disk, $h = \frac{2}{N}$ and for the TODO domain, $h = \frac{3}{N}$.

Similar to Chati et al., we test our method using a boundary value function $g(x, y) = x^3 + y^3 - 3 x y^2 - 3 y^2 x + 1$. The solution to Laplace's equation with $g(x, y)$ as the boundary condition is (trivially) $g(x, y)$ itself.

Results. Figure @fig:random-walks shows 10 random walks on the unit disk starting from a boundary point, both with and without coupling.

![10 random walks on the unit disk. Left: Without coupling. Right: With coupling. Walks were taken with a random seed of zero. The red dots indicate the final position of the walk. The small blue dots are the boundary points. The large blue dot is the initial position of the walker.](./figures/random-walks.png) {#fig:random-walks}

Table 1 shows the maximum absolute error, maximum relative error, and mean relative error for our method on the unit disk, with no coupling. Figure @fig:results-circle plots the numerical and exact solutions for the unit disk, both with and without coupling.

Table 2 shows the error for our method with coupling. Contrary to what we expected, the max error does not scale with $K$ in our coupled implementation.

\begin{table}[h!]
\begin{center}
\begin{tabular}{|c|c|c|c|c|}
	\hline
	N & K & Max abs err & Max rel err & Mean rel err \\
	\hline
	50 & 25 & 0.1167 & 21.6 & 0.033474 \\
	50 & 100 & 0.07035 & 8.1546 & 0.02415 \\
	50 & 400 & 0.03810 & 20.683& 0.0302 \\
	100 & 25 & 0.1176 & 36.625 & 0.01340 \\
	100 & 100 & 0.06578 & 7.22489 & 0.00755 \\
	100 & 400 & 0.01781 & 5.2412 & 0.002187 \\
    200 & 25 & 0.07451 & 12.5933 & 0.002244 \\
    200 & 100 & 0.04454 & 15.9934 & 0.0002541 \\
    200 & 400 & 0.01932 & 4.6903 & 0.000526\\
    50 & "$\inf$" & 1.22e-14 & 2.5057e-12 & 9.2972e-15 \\
	100 & "$\inf$" & 5.2402e-14 & 1.1273e-11 & 3.3339e-14 \\
    200 & "$\inf$" & 1.5876e-13 & 3.9257e-11 & 7.9955e-14\\
	\hline
\end{tabular}
\end{center}
\caption{ Error of the hybrid random walk method on the unit disk with no coupling for different choices of $K$ and $N$. $K = \inf$ corresponds to using the exact solution as the boundary-value data for our finite-differencing scheme.}
\end{table}

\begin{table}[h!]
\begin{center}
\begin{tabular}{|c|c|c|c|c|}
	\hline
	N & K & Max abs err & Max rel err & Mean rel err \\
	\hline
    TODO
\end{tabular}
\end{center}
\caption{ Error of the hybrid random walk method on the unit disk with coupling for different choices of $K$ and $N$. }
\end{table}

Figure @fig:squiggly plots the numerical and exact solutions for the TODO domain, as well as the domain itself. We did not use coupling to compute the numerical solution because of its poor performance on the unit disk. Table 3 shows the error for the TODO domain.

![The numerically computed solution for the TODO domain, compared against the exact solution. Left: Plot of the domain. Center: Numerically computed solution with no coupling. Right: Exact solution.](./figures/squiggly.png) {#fig:squiggly}

\begin{table}[h!]
\begin{center}
\begin{tabular}{|c|c|c|c|c|}
	\hline
	N & K & Max abs err & Max rel err & Mean rel err \\
	\hline
    TODO
\end{tabular}
\end{center}
\caption{ Error of the hybrid random walk method on the unit disk with coupling for different choices of $K$ and $N$. }
\end{table}

## Conclusions.

We have now seen that we can use random walks to solve Laplace's Problem.
While the performance is not the best,
we still get pretty good results.
The algorithm takes a lot of time to converge,
which means that in practice, it is not very efficient.

## Future Directions.

## Contributions.

Caroline did this stuff.
Dmitry did this other stuff.


