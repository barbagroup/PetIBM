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
The software package includes extended documentation as well as, many examples to guide users.

PetIBM has already been used to generate results published in @mesnard_barba_2017, a full replication of a study on the aerodynamics of a gliding snake species (@krishnan_et_al_2014).
PetIBM is currently used to compute the three-dimensional flow of a gliding-snake model on the cloud platform Microsoft Azure.

# References
