---
title: 'AmgXWrapper: An interface between PETSc and the NVIDIA AmgX library'
tags:
  - PETSc
  - AmgX
  - multigrid
  - multi-GPU computing
authors:
 - name: Pi-Yueh Chuang
   orcid: 0000-0001-6330-2709
   affiliation: 1
 - name: Lorena A. Barba
   orcid: 0000-0001-5812-2711
   affiliation: 1
affiliations:
 - name: The George Washington University
   index: 1
date: 02 June 2017
bibliography: paper.bib
---

# Summary

This code provides access to multi-GPU computing from PETSc-based applications, using NVIDIA's AmgX library.
The PETSc (Portable, Extensible Toolkit for Scientific Computation) library,
developed by Argonne National Laboratory, is a very popular
high-performance software through for scientific applications modeled by partial differential equations [@petsc-user-ref].
It executes in distributed systems via message-passing communications with MPI.
NVIDIA's AmgX library is a multi-GPU linear-algebra package that
supports multigrid preconditioners and Krylov solvers.
Because PETSc lacks support for modern heterogeneous platforms (GPU+CPU),
AmgX can be a good option for adding GPU computing to PETSc applications.
Incorporating AmgX into PETSc applications is not straightforward, however,
due to the difference between the two libraries in design, underlying data formats, and usage.
The purpose of AmgXWrapper is to bridge PETSc and AmgX, helping PETSc-application
developers use AmgX in their software without a thorough understanding of the AmgX API.
With AmgXWrapper, developers may need only a few lines of code modification to 
add AmgX solvers in legacy PETSc applications.
The wrapper also features implicit data transfer when there are mismatched numbers
of CPU cores and GPU devices in a computing node.
This allows exploiting all possible resources on modern heterogeneous platforms.
We have presented an example of using AmgX and AmgXWrapper to accelerate an
existing PETSc-based CFD code [@chuang+barba2017].



# References
