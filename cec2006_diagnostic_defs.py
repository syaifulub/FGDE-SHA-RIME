"""CEC2006 problems G5 and G13 (Liang et al., 2006), used for diagnostic development.

Definitions follow the CEC2006 technical report. The reference optima are those of the
pymoo implementation used in the original diagnostic experiment (G5: 5126.49811,
G13: 0.05395), so absolute errors are comparable with data/diagnostic/v2_diagnostic_fd25.csv.
"""
from __future__ import annotations
import numpy as np
from cec2006_unseen_defs import P


class Pf(P):
    def __init__(self, name, n, xl, xu, fstar):
        super().__init__(name, n, xl, xu, np.zeros(n)); self.fstar = fstar
    def pareto_front(self):
        return np.array([self.fstar])


class G5(Pf):
    def __init__(self):
        super().__init__('g5', 4, np.array([0.0, 0.0, -0.55, -0.55]),
                         np.array([1200.0, 1200.0, 0.55, 0.55]), 5126.49811)
    def _eval(self, x):
        x1, x2, x3, x4 = x.T
        f = 3*x1 + 1e-6*x1**3 + 2*x2 + (2e-6/3)*x2**3
        g = np.column_stack([-x4 + x3 - 0.55, -x3 + x4 - 0.55])
        h = np.column_stack([
            1000*np.sin(-x3 - 0.25) + 1000*np.sin(-x4 - 0.25) + 894.8 - x1,
            1000*np.sin(x3 - 0.25) + 1000*np.sin(x3 - x4 - 0.25) + 894.8 - x2,
            1000*np.sin(x4 - 0.25) + 1000*np.sin(x4 - x3 - 0.25) + 1294.8])
        return f, g, h


class G13(Pf):
    def __init__(self):
        super().__init__('g13', 5, np.array([-2.3, -2.3, -3.2, -3.2, -3.2]),
                         np.array([2.3, 2.3, 3.2, 3.2, 3.2]), 0.05395)
    def _eval(self, x):
        x1, x2, x3, x4, x5 = x.T
        f = np.exp(x1*x2*x3*x4*x5)
        h = np.column_stack([x1**2 + x2**2 + x3**2 + x4**2 + x5**2 - 10,
                             x2*x3 - 5*x4*x5, x1**3 + x2**3 + 1])
        return f, None, h


REG = {p.name: p for p in [G5(), G13()]}
def get_problem(name): return REG[name]
