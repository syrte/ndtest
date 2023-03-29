# ndtest
Multi-dimensional statistical test with python

- 2D Kolmogorov–Smirnov test (KS test)
- 1D and 2D energy distance statistics test

The code has been cited in 20 papers 
according to [Google Scholar](https://scholar.google.com/scholar?q=%22syrte%2Fndtest%22) :star_struck:

## Usage example
```
import numpy as np
import ndtest

# generate mock samples for the test
np.random.seed(42)
x1, y1 = np.random.randn(2, 100)
x2, y2 = np.random.randn(2, 100)  # same distribution as (x1, y1)
x3, y3 = np.random.randn(2, 100) * 1.5 + 0.5  # different distribution from (x1, y1)

# 2D KS
P, D = ndtest.ks2d2s(x1, y1, x2, y2, extra=True)
print(f"{P=:.3g}, {D=:.3g}")
# P=0.219, D=0.17

P, D = ndtest.ks2d2s(x1, y1, x3, y3, extra=True)
print(f"{P=:.3g}, {D=:.3g}")
# P=2.36e-05, D=0.385  # very small P, as expected.
```

See the docstring for the detailed usage and explanation.
