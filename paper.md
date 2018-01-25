---
title: 'PetIBM: toolbox and applications of the immersed-boundary method on distributed-memory architectures'
tags:
  - Computational Fluid Dynamics
  - Immersed-Boundary Method
  - PETSc
  - NVIDIA AmgX
authors:
 - name: Pi-Yueh Chuang
   orcid: 0000-0001-6330-2709
   affiliation: 1
 - name: Olivier Mesnard
   orcid:  0000-0001-5335-7853
   affiliation: 1
 - name: Lorena A. Barba
   orcid: 0000-0001-5812-2711
   affiliation: 1
 - name: Anush Krishnan
   affiliation: 2
affiliations:
 - name: The George Washington University
   index: 1
 - name: nuTonomy Inc. (previously at Boston University)
   index: 2
date: 25 January 2018
bibliography: paper.bib
---

# Summary

PetIBM solves the two- and three-dimensional Navier-Stokes equations with an immersed-boundary method on fixed structured Cartesian grids.
In this method, a collection of Lagrangian markers defines the immersed boundary (where boundary conditions are enforced) and the fluid equations are solved over the extended domain (including the body domain).
The Eulerian mesh remains unmodified when computing the flow around multiple moving immersed bodies, which removes the need for remeshing at every time step.
PetIBM discretizes the fluid equations using a second-order finite-difference scheme, various optional time-integrators, and a fully discrete projection method (@perot_1993).
It implements two immersed-boundary algorithms: the immersed-boundary projection method (@taira_colonius_2007) and its decoupled version (@li_et_al_2016).

PetIBM is written in C++ and relies on the PETSc library (@petsc_1997, @petsc_user_ref_2017) to run on memory-distributed architectures.
PetIBM can solve one or several linear systems on multiple distributed CUDA-capable GPU devices with the NVIDIA AmgX library and AmgXWrapper (@chuang_barba_2017).

PetIBM has already been used to generate results published in @mesnard_barba_2017, a full replication of a study on the aerodynamics of a gliding snake species (@krishnan_et_al_2014).
PetIBM is currently used to compute the three-dimensional flow of a gliding-snake model on the cloud platform Microsoft Azure.

# References
