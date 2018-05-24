---
title: 'PetIBM: toolbox and applications of the immersed-boundary method on distributed-memory architectures'
tags:
  - Computational Fluid Dynamics
  - Immersed-Boundary Method
  - PETSc
  - GPU
  - NVIDIA AmgX
authors:
 - name: Pi-Yueh Chuang
   orcid: 0000-0001-6330-2709
   affiliation: 1
 - name: Olivier Mesnard
   orcid: 0000-0001-5335-7853
   affiliation: 1
 - name: Anush Krishnan
   orcid: 0000-0001-6409-7022
   affiliation: 2
 - name: Lorena A. Barba
   orcid: 0000-0001-5812-2711
   affiliation: 1
affiliations:
 - name: Department of Mechanical and Aerospace Engineering, The George Washington University, Washington, DC, USA
   index: 1
 - name: nuTonomy Inc., Cambridge, MA, USA (previously at Boston University)
   index: 2
date: 25 January 2018
bibliography: paper.bib
---

# Summary

PetIBM is a C++ library with ready-to-use application codes to solve the two- and three-dimensional incompressible Navier-Stokes equations on fixed structured Cartesian grids with an immersed-boundary method (IBM).
PetIBM runs on distributed-memory architectures and can be used to compute the flow around multiple moving rigid immersed boundaries (with prescribed kinematics).

In the IBM framework, a collection of Lagrangian markers defines the immersed boundary (where boundary conditions are enforced) and the fluid equations are solved over the extended domain (including the body domain).
The Eulerian mesh remains unmodified when computing the flow around multiple moving immersed bodies, which removes the need for remeshing at every time step.
PetIBM discretizes the fluid equations using a second-order finite-difference scheme, various optional time-integrators, and a fully discrete projection method (@perot_1993).
It implements two immersed-boundary algorithms: the immersed-boundary projection method (@taira_colonius_2007) and its decoupled version (@li_et_al_2016).

Other open-source software packages offer immersed-boundary solvers: for example, IBAMR (@griffith_et_al_2007, @bhalla_et_al_2013) is a long-standing C++ library with MPI parallelization that also provides adaptive mesh refinement.
It can handle deforming immersed bodies and has been used in a variety of scenarios, including cardiac fluid dynamics, swimming, insect flight, and others.
PetIBM and IBAMR use different immersed-boundary schemes, however. We developed PetIBM to work with the immersed-boundary projection method, which is based on the fully discrete formulation of Perot on staggered grids and thus eliminates the need for pressure boundary conditions, which have caused many headaches for CFD practitioners (@gresho_sani_1987, @sani_gresho_1994).
PetIBM features an operator-based design, providing routines to create and manipulate discrete operators (e.g., gradient, divergence, Laplacian, convection, diffusion, etc.), so it can be used as a toolbox for researching new solution methods.
It is also capable of using graphics processing unit (GPU) architectures, a feature missing from other software, as far as we know.
A previous project implementing immersed-boundary methods on GPU architecture is cuIBM (@krishnan_et_al_2017), but it is limited to two-dimensional problems that fit on a single GPU device.

PetIBM is written in C++ and relies on the PETSc library (@petsc_1997, @petsc_user_ref_2017) for data structures and parallel routines to run on memory-distributed architectures.
PetIBM can solve one or several linear systems on multiple distributed CUDA-capable GPU devices with the NVIDIA AmgX library and AmgXWrapper (@chuang_barba_2017).
The software package includes extended documentation as well as many examples to guide users.

PetIBM has already been used to generate results published in @mesnard_barba_2017, a full replication of a study on the aerodynamics of a gliding snake species (@krishnan_et_al_2014).
PetIBM is currently used to compute the three-dimensional flow of a gliding-snake model on the cloud platform Microsoft Azure.

## Appendix: mathematical formulation

PetIBM solves the Navier-Stokes equations on an extended discretization grid that includes the interior of the immersed boundary.
To model the presence of the boundary, a forcing term is added to the momentum equation and an additional equation for the no-slip condition completes the system.
Many variants of the immersed-boundary method (IBM) depend on one models the forcing.
In PetIBM, we use regularized-delta functions to transfer data between the Eulerian grid and the Lagrangian boundary points.
The system of equation is:

\begin{equation}
\begin{cases}
\frac{\partial \mathbf{u}}{\partial t} + \mathbf{u} \cdot \nabla \mathbf{u} = -\nabla p + \frac{1}{Re} \nabla^2 \mathbf{u} + \int_{s}{\mathbf{f} \left( \mathbf{\xi} \left( \mathit{s}, \mathit{t} \right) \right) \delta \left( \mathbf{\xi} - \mathbf{x} \right)} d\mathit{s} \\
\nabla \cdot \mathbf{u} = 0 \\
\mathbf{u} \left( \mathbf{\xi} \left( \mathit{s}, t \right) \right) = \int_{\mathbf{x}}{\mathbf{u} \left( \mathbf{x} \right)} \delta \left( \mathbf{x} - \mathbf{\xi} \right) d\mathbf{x}
\end{cases}
\end{equation}

where $\mathbf{u}$ is the velocity field, $p$ is the fluid pressure, and $Re$ is the Reynolds number.

Currently, PetIBM provides two application codes implementing different versions of the IBM: (1) an immersed-boundary projection method (IBPM) based on the work of @taira_colonius_2007 and (2) a decoupled version of the IBPM proposed by @li_et_al_2016.
Those two methods fit into the framework of the projection approach of @perot_1993.
The equations are fully discretized (space and time) to form an algebraic system to be solved for the velocity $u^{n+1}$, the pressure field $\phi$, and the Lagrangian forces $\tilde{f}$.
The discretized system is:

\begin{equation}
\left[
\begin{matrix}
A & G & H \\
D & 0 & 0 \\
E & 0 & 0
\end{matrix}
\right]
\left(
\begin{matrix}
u^{n+1} \\
\phi \\
\tilde{f}
\end{matrix}
\right)
=
\left(
\begin{matrix}
r^n \\
0 \\
u_B^{n+1}
\end{matrix}
\right)
+
\left(
\begin{matrix}
bc_1 \\
bc_2 \\
0
\end{matrix}
\right)
\end{equation}

where $D$, $G$, and $A$ are the divergence, gradient, and implicit operators, respectively.
$E$ and $H$ are the interpolation and spreading operators, respectively, used to transfer the data between the Eulerian grid and the Lagrangian boundary points.
On the right-hand side, $r^n$ gathers all the explicit terms and $u_B^{n+1}$ is the known (prescribed) boundary velocity; $bc_1$ and $bc_2$ contain the boundary terms that arise from the discretization of momentum and continuity equations, respectively.

In the IBPM, we solve a modified Poisson system for the pressure field and Lagrangian forces, coupled together.
This way, the divergence-free condition and no-slip constraint are simultaneously enforced on the velocity field at the end of the time step.
The fully discretized system can be cast into the following:

\begin{equation}
\left[
\begin{matrix}
A & Q_2 \\
Q_1 & 0
\end{matrix}
\right]
\left(
\begin{matrix}
u^{n+1} \\
\lambda
\end{matrix}
\right)
=
\left(
\begin{matrix}
r_1 \\
r_2
\end{matrix}
\right)
\end{equation}

with

\begin{equation*}
Q_1 \equiv \left[ \begin{matrix} D \\ E \end{matrix} \right] ;\
Q_2 \equiv \left[ G, H \right] ;\
\lambda \equiv \left( \begin{matrix} \phi \\ \tilde{f} \end{matrix} \right) ;\
r_1 \equiv r^n + bc_1 ;\
r_2 \equiv \left( \begin{matrix} bc_2 \\ u_B^{n+1} \end{matrix} \right)
\end{equation*}

In practice, we never form the full system.
Instead, we apply a block-LU decomposition as follow:

\begin{equation}
\left[
\begin{matrix}
A & 0 \\
Q_1 & -Q_1A^{-1}Q_2
\end{matrix}
\right]
\left[
\begin{matrix}
I & A^{-1}Q_2 \\
0 & I
\end{matrix}
\right]
\left(
\begin{matrix}
u^{n+1} \\
\lambda
\end{matrix}
\right)
=
\left[
\begin{matrix}
A & 0 \\
Q_1 & -Q_1A^{-1}Q_2
\end{matrix}
\right]
\left(
\begin{matrix}
u^* \\
\lambda
\end{matrix}
\right)
=
\left(
\begin{matrix}
r_1 \\
r_2
\end{matrix}
\right)
\end{equation}

Thus, we retrieve the sequence of operations of the traditional projection method.
We solve a system for an intermediate velocity field that is corrected, after solving a modified Poisson system for the variable $\lambda$, to enforce the divergence-free condition and the no-slip constraint at the location of the immersed boundary.
The sequence is:

\begin{align}
& A u^* = r_1 \\
& Q_1A^{-1}Q_2 \lambda = Q_1 u^* - r_2 \\
& u^{n+1} = u^* - A^{-1}Q_1 \lambda
\end{align}

The IBPM implemented in PetIBM solves, at every time step, Equations (5) to (6).
(Note: the inverse of the implicit operator $A^{-1}$ is approximated by a finite Taylor series expansion.)


The IBPM requires solving, at each time step, an expensive modified Poisson system, $Q_1A^{-1}Q_2$, whose non-zero structure changes when the location of the immersed boundary is moving.
In the PetIBM implementation of the decoupled IBPM, we apply a second block-LU decomposition to decouple the pressure field from the Lagrangian forces and recover a classical Poisson system.
The fully discretized algebraic system can be cast into:

\begin{equation}
\left[
\begin{matrix}
A & H & G \\
E & 0 & 0 \\
D & 0 & 0
\end{matrix}
\right]
\left(
\begin{matrix}
u^{n+1} \\
\tilde{f} \\
\phi
\end{matrix}
\right)
=
\left(
\begin{matrix}
r^n \\
u_B^{n+1} \\
0
\end{matrix}
\right)
+
\left(
\begin{matrix}
bc_1 \\
0 \\
bc_2
\end{matrix}
\right)
\end{equation}

The velocity $u^{n+1}$ and the Lagrangian forces $\tilde{f}$ are coupled together to form a new unknown $\gamma^{n+1}$, as follows:

\begin{equation}
\left[
\begin{matrix}
\bar{A} & \bar{G} \\
\bar{D} & 0
\end{matrix}
\right]
\left(
\begin{matrix}
\gamma^{n+1} \\
\phi
\end{matrix}
\right)
=
\left(
\begin{matrix}
\bar{r_1} \\
\bar{r_2}
\end{matrix}
\right)
\end{equation}

where

\begin{equation*}
\bar{A} \equiv \left[ \begin{matrix} A & H \\ E & 0 \end{matrix} \right] ;\
\bar{G} \equiv \left[ \begin{matrix} G \\ 0 \end{matrix} \right] ;\
\bar{D} \equiv \left[ \begin{matrix} D & 0 \end{matrix} \right]
\end{equation*}

and

\begin{equation*}
\gamma^{n+1} \equiv \left( \begin{matrix} u^{n+1} \\ \tilde{f} \end{matrix}\right) ;\
\bar{r_1} \equiv \left( \begin{matrix} r_n + bc_1 \\ u_B^{n+1} \end{matrix}\right) ;\
\bar{r_2} \equiv bc_2
\end{equation*}

Two successive block-LU decompositions are applied to decouple the Lagrangian forces $\tilde{f}$ from $\gamma^{n+1}$ and to decouple the velocity from the pressure field.

The first block-LU decomposition decouples the pressure field from the new unknown $\gamma^{n+1}$:

\begin{equation}
\left[
\begin{matrix}
\bar{A} & 0 \\
\bar{D} & -\bar{D}\bar{A}^{-1}\bar{G}
\end{matrix}
\right]
\left[
\begin{matrix}
I & \bar{A}^{-1}\bar{G} \\
0 & I
\end{matrix}
\right]
\left(
\begin{matrix}
\gamma^{n+1} \\
\phi
\end{matrix}
\right)
=
\left[
\begin{matrix}
\bar{A} & 0 \\
\bar{D} & -\bar{D}\bar{A}^{-1}\bar{G}
\end{matrix}
\right]
\left(
\begin{matrix}
\gamma^* \\
\phi
\end{matrix}
\right)
=
\left(
\begin{matrix}
\bar{r_1} \\
\bar{r_2}
\end{matrix}
\right)
\end{equation}

which leads to the following sequence of operations:

\begin{align}
& \bar{A} \gamma^* = \bar{r_1} \\
& \bar{D}\bar{A}^{-1}\bar{G} \phi = \bar{D} \gamma^* - \bar{r_2} \\
& \gamma^{n+1} = \gamma^* - \bar{A}^{-1}\bar{G} \phi
\end{align}

A second block-LU decomposition is applied to the first equation above:

\begin{equation}
\left[
\begin{matrix}
A & 0 \\
E & -EA^{-1}H
\end{matrix}
\right]
\left[
\begin{matrix}
I & A^{-1}H \\
0 & I
\end{matrix}
\right]
\left(
\begin{matrix}
u^* \\
\tilde{f}
\end{matrix}
\right)
=
\left[
\begin{matrix}
A & 0 \\
E & -EA^{-1}H
\end{matrix}
\right]
\left(
\begin{matrix}
u^{* *} \\
\tilde{f}
\end{matrix}
\right)
=
\left(
\begin{matrix}
r^n + bc_1 \\
u_B^{n+1}
\end{matrix}
\right)
\end{equation}

and we end up with the following sequence:

\begin{align}
& A u^{* *} = r^n + bc_1 \\
& EA^{-1}H \tilde{f} = E u^{* *} - u_B^{n+1} \\
& u^* = u^{* *} - A^{-1}H \tilde{f}
\end{align}

The decoupled version of the IBPM implemented in PetIBM solves, at every time step, Equations (15) to (17) followed by Equations (12) and (13).

# References
