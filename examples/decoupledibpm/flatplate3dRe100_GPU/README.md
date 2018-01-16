# 3D flow around inclined flat plate (Re=100, AR=2)

Run the example using 1 GPU:

```
export CUDA_VISIBLE_DEVICES=0
angles = '0 10 20 30 40 50 60 70 80 90'
for angle in $angles
do
	python $angle/scripts/createBody.py
	petibm-decoupledibpm -directory $angle -options_left -log_view ascii:$angle/stdout.txt
done
```

Each simulation completes in about 10 minutes when using:
- 1 NVIDIA K40 GPU device.

Plot the instantaneous force coefficients and compare with the experimental
results from Taira et al. (2007):

```
python scripts/plotDragCoefficient.py
```

_References:_
* Taira, K., Dickson, W. B., Colonius, T., Dickinson, M. H., & Rowley, C. W. (2007). "Unsteadiness in flow over a flat plate at angle-of-attack at low Reynolds numbers." AIAA Paper, 710, 2007.
