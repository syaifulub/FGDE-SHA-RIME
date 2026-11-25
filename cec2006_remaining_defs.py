"""CEC2006 problems G2, G7, G10, G11, and G24 (Liang et al., 2006), added in revision.

Definitions follow the CEC2006 technical report; best-known objective values are the
reported f(x*) values and are used directly as reference optima.
"""
from __future__ import annotations
import numpy as np
from cec2006_unseen_defs import P


class Pf(P):
    """Problem with a fixed reported reference optimum."""
    def __init__(self, name, n, xl, xu, fstar):
        super().__init__(name, n, xl, xu, np.zeros(n)); self.fstar = fstar
    def pareto_front(self):
        return np.array([self.fstar])


class G2(Pf):
    def __init__(self): super().__init__('g2', 20, np.zeros(20), np.full(20, 10.0), -0.80361910412559)
    def _eval(self, x):
        c = np.cos(x); i = np.arange(1, x.shape[1] + 1)
        den = np.sqrt(np.maximum(np.sum(i * x ** 2, axis=1), 1e-300))
        f = -np.abs(np.sum(c ** 4, axis=1) - 2 * np.prod(c ** 2, axis=1)) / den
        g1 = 0.75 - np.prod(x, axis=1); g2 = np.sum(x, axis=1) - 7.5 * x.shape[1]
        return f, np.column_stack([g1, g2]), None


class G7(Pf):
    def __init__(self): super().__init__('g7', 10, np.full(10, -10.0), np.full(10, 10.0), 24.30620906818)
    def _eval(self, x):
        x1, x2, x3, x4, x5, x6, x7, x8, x9, x10 = x.T
        f = (x1**2 + x2**2 + x1*x2 - 14*x1 - 16*x2 + (x3-10)**2 + 4*(x4-5)**2 + (x5-3)**2 + 2*(x6-1)**2
             + 5*x7**2 + 7*(x8-11)**2 + 2*(x9-10)**2 + (x10-7)**2 + 45)
        g = np.column_stack([
            4*x1 + 5*x2 - 3*x7 + 9*x8 - 105,
            10*x1 - 8*x2 - 17*x7 + 2*x8,
            -8*x1 + 2*x2 + 5*x9 - 2*x10 - 12,
            3*(x1-2)**2 + 4*(x2-3)**2 + 2*x3**2 - 7*x4 - 120,
            5*x1**2 + 8*x2 + (x3-6)**2 - 2*x4 - 40,
            x1**2 + 2*(x2-2)**2 - 2*x1*x2 + 14*x5 - 6*x6,
            0.5*(x1-8)**2 + 2*(x2-4)**2 + 3*x5**2 - x6 - 30,
            -3*x1 + 6*x2 + 12*(x9-8)**2 - 7*x10])
        return f, g, None


class G10(Pf):
    def __init__(self):
        super().__init__('g10', 8, np.array([100, 1000, 1000, 10, 10, 10, 10, 10], float),
                         np.array([10000, 10000, 10000, 1000, 1000, 1000, 1000, 1000], float), 7049.24802052867)
    def _eval(self, x):
        x1, x2, x3, x4, x5, x6, x7, x8 = x.T
        f = x1 + x2 + x3
        g = np.column_stack([
            -1 + 0.0025*(x4 + x6),
            -1 + 0.0025*(x5 + x7 - x4),
            -1 + 0.01*(x8 - x5),
            -x1*x6 + 833.33252*x4 + 100*x1 - 83333.333,
            -x2*x7 + 1250*x5 + x2*x4 - 1250*x4,
            -x3*x8 + 1250000 + x3*x5 - 2500*x5])
        return f, g, None


class G11(Pf):
    def __init__(self): super().__init__('g11', 2, np.array([-1.0, -1.0]), np.array([1.0, 1.0]), 0.7499)
    def _eval(self, x):
        f = x[:, 0]**2 + (x[:, 1] - 1)**2
        h = x[:, 1] - x[:, 0]**2
        return f, None, h


class G24(Pf):
    def __init__(self): super().__init__('g24', 2, np.array([0.0, 0.0]), np.array([3.0, 4.0]), -5.50801327159536)
    def _eval(self, x):
        x1, x2 = x[:, 0], x[:, 1]
        f = -x1 - x2
        g1 = -2*x1**4 + 8*x1**3 - 8*x1**2 + x2 - 2
        g2 = -4*x1**4 + 32*x1**3 - 88*x1**2 + 96*x1 + x2 - 36
        return f, np.column_stack([g1, g2]), None


REG = {p.name: p for p in [G2(), G7(), G10(), G11(), G24()]}
def get_problem(name): return REG[name]
