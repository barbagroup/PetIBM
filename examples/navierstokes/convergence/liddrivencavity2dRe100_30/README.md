# Convergence test of the Navier-Stokes solver

Run the cases from the present directory:

```
petibm-navierstokes -directory 30 -options_left -log_view ascii:30/view.log
petibm-navierstokes -directory 90 -options_left -log_view ascii:90/view.log
petibm-navierstokes -directory 270 -options_left -log_view ascii:270/view.log
petibm-navierstokes -directory 810 -options_left -log_view ascii:810/view.log
```

Computes the observed orders of convergence for the velocity components and
the pressure using three consistently refined grids:

```
python scripts/getOrderConvergence.py
```
