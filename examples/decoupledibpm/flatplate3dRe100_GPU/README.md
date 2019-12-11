# 3D flow around inclined flat plate (Re=100, AR=2)

Run the example using 1 GPU:

```
export CUDA_VISIBLE_DEVICES=<idx1>
# Declare list of angles of attack to run
angles=(0 10 20 30 40 50 60 70 80 90)
for angle in ${angles[@]}; do
    echo "*** Angle of attack: $angle degrees ***"
    # Generate immersed boundary points
    python "AoA$angle/scripts/createBody.py"
    # Run PetIBM
    mpiexec -np 4 petibm-decoupledibpm \
        -directory "AoA$angle" \
        -options_left \
        -log_view ascii:"AoA$angle/view.log"
done
```

Each simulation completes in about 8 minutes when using:
- 4 CPU processes (Intel(R) Core(TM) i7-3770 CPU @ 3.40GHz),
- 1 NVIDIA K40 GPU device.

Plot the instantaneous force coefficients and compare with the experimental
results from Taira et al. (2007):

```
python scripts/plotDragCoefficient.py
```

The figure will be saved in the sub-folder `figures`.

_References:_
* Taira, K., Dickson, W. B., Colonius, T., Dickinson, M. H., & Rowley, C. W. (2007). "Unsteadiness in flow over a flat plate at angle-of-attack at low Reynolds numbers." AIAA Paper, 710, 2007.
