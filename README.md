Theory:

1. Naive matrix operations with floats can sometimes be unstable, as errors propagate
2. Interval arithmetic handles errors inherently
3. HOWEVER: with regular intervals, the bounds widen rapidly, making them unsuitable for matrix operations
4. Affine arithmetic keeps track of correlations, not just the total error

Question:

Can affine arithmetic be used to provide tighter bounds in matrix operations?

Current progress:
- Currently implemented matrix operations with float and regular intervals
- Managed to observe the interval bounds growing out of control
- TODO: Research and replicate the exact circumstances where float operations prove unstable
- TODO: Research and implement affine arithmetic to see if they are at least better than interval arithmetic