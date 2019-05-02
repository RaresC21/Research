import matplotlib.pyplot as plt
import numpy as np

from Plane import Options
from Polytope import Polytope
from symmetric import *
from far_point import *
from make_cube import *

dim = 10

delta = 0.08 / np.sqrt(dim)
radius = 0.1
depth = 0

options = Options(reflection = "sphere", \
                        delta = delta, radius = radius, depth = depth)
poly = symmetric(dim, options = options)

poly.reflection_lines(100)
