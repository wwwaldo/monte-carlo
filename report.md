---
title: Solutions to AER1316 A3
author: Caroline (Ruoyi) Lin, 1001333112
header-includes: |
    \usepackage{gensymb}
geometry: margin=1in
---

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

### Summary.

Pick a width, h.
Construct a Cartesian mesh which lies on the interior.
Perform random walks for u on the boundary of the cartesian mesh.
Solve the coupled set of equations to get u on the boundary of the mesh.
Use finite differences to compute the result.

### The finite difference operator.

Here it is.

### The coupled set of equations.

Here's the derivation.

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

Walk-on-spheres
boundary-within-a-boundary

## Performance

Here is where Dmitry puts in some data.

## Some figures and results.

Here is a test domain: A rectangle with these dimensions, with walks from a single particle.

Here is the frequency data for a coupled solution.

Some more results: Here is the error and plot for our circle data.

Here is the data for the coupled solution. Something seems to be off with the coupled solution.

Here is the error and plot for our squiggly domain data with a cusp, and no coupling.

## Conclusions.

## Future Directions.

## Contributions.

Caroline did this stuff.
Dmitry did this other stuff.


