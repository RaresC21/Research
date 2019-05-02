import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, draw, show, ion

import numpy as np
import time

def draw_it(marginal, c = 'blue', a = 1):
    i = 0
    for m in marginal:
        if m > 0:
            plt.scatter(i, m, color = c, alpha = a)
        i+=1

all_m_far = [[0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ,
 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ,
 0.     , 0.     , 0.     , 0.43418, 0.     , 0.     , 0.     , 0.     ,
 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ]
,[0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ,
 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ,
 0.     , 0.     , 0.     , 0.75134, 0.29063, 0.11442, 0.     , 0.     ,
 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ]
,[0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ,
 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ,
 0.     , 0.     , 0.     , 0.84159, 0.42932, 0.12877, 0.     , 0.     ,
 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ]
,[0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ,
 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ,
 0.     , 0.99116, 0.97677, 0.8937 , 0.66999, 0.46274, 0.33952, 0.26199,
 0.08577, 0.00184, 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ]
,[0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ,
 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.91966, 0.71794,
 0.60601, 0.50112, 0.46192, 0.42264, 0.31684, 0.21883, 0.16056, 0.12389,
 0.04056, 0.00087, 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ]
,[0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ,
 0.     , 0.     , 0.     , 0.     , 0.9643 , 0.91078, 0.85466, 0.74292,
 0.68154, 0.60197, 0.52229, 0.42196, 0.26803, 0.19109, 0.13236, 0.08265,
 0.02166, 0.00042, 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ]
,[0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ,
 0.     , 0.     , 0.     , 0.     , 0.97053, 0.90993, 0.84586, 0.73359,
 0.64985, 0.56462, 0.4691 , 0.32908, 0.19898, 0.12771, 0.09108, 0.04861,
 0.01054, 0.0002 , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ]
,[0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.99407, 0.98537,
 0.98114, 0.96864, 0.96377, 0.95681, 0.93359, 0.88461, 0.81589, 0.72915,
 0.65052, 0.55386, 0.46968, 0.3362 , 0.24394, 0.17178, 0.12641, 0.07985,
 0.04417, 0.02215, 0.01259, 0.00143, 0.     , 0.     , 0.     , 0.     ]
,[0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.99709, 0.99282,
 0.99075, 0.98462, 0.97966, 0.96823, 0.94257, 0.89523, 0.82232, 0.73789,
 0.64074, 0.54778, 0.44941, 0.33614, 0.24136, 0.16121, 0.10663, 0.06559,
 0.03493, 0.01708, 0.0091 , 0.0007 , 0.     , 0.     , 0.     , 0.     ]]

all_m_zero = [[0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
0. , 0. , 0. , 0.5, 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
0. , 0. , 0. ] ,[0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0.
, 0. , 0. , 0. , 0.5, 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0.
, 0. , 0. , 0. ] ,[0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
0. , 0. , 0. , 0. , 0.5, 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
0. , 0. , 0. , 0. ] ,[0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.
, 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ,
0.83003, 0.5    , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ,
0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ] ,[0.
, 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.
, 0.     , 0.     , 0.     , 0.     , 0.     , 0.75064, 0.5    , 0.     , 0.
, 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.
, 0.     , 0.     , 0.     , 0.     ] ,[0.     , 0.     , 0.     , 0.     , 0.
, 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.
, 0.     , 0.85778, 0.5    , 0.     , 0.     , 0.     , 0.     , 0.     , 0.
, 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.
] ,[0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.
, 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.92358, 0.5    , 0.
, 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.
, 0.     , 0.     , 0.     , 0.     , 0.     ] ,[0.     , 0.     , 0.     , 0.
, 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.
, 0.     , 0.     , 0.94789, 0.5    , 0.09064, 0.03754, 0.     , 0.     , 0.
, 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.
, 0.     ] ,[0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.
, 0.     , 0.     , 0.     , 0.     , 0.     , 0.98037, 0.93391, 0.87876, 0.5
, 0.36326, 0.24797, 0.1127 , 0.03231, 0.     , 0.     , 0.     , 0.     , 0.
, 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ] ,[0.     , 0.     , 0.
, 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.9929 ,
0.94964, 0.92351, 0.88083, 0.81889, 0.78154, 0.5    , 0.32875, 0.22225, 0.13626,
0.056  , 0.02142, 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ,
0.     , 0.     , 0.     ] ,[0.     , 0.     , 0.     , 0.     , 0.     , 0.
, 0.     , 0.     , 0.     , 0.     , 0.9962 , 0.97305, 0.94619, 0.87192,
0.79446, 0.66948, 0.5    , 0.387  , 0.27686, 0.1518 , 0.07034, 0.04059, 0.0033 ,
0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ]
,[0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.
, 0.     , 0.9981 , 0.9865 , 0.96832, 0.89273, 0.79864, 0.65436, 0.5    ,
0.35101, 0.25111, 0.15154, 0.06494, 0.02033, 0.00165, 0.     , 0.     , 0.     ,
0.     , 0.     , 0.     , 0.     , 0.     , 0.     ]]

plt.xlim((-1,33))
plt.ylim((-0.1,1.1))
ion()

a = 0.2

c = ['lightseagreen', 'pink', 'orange', 'red', 'green', 'mediumseagreen',
     'blue', 'midnightblue', 'mediumpurple', 'purple', 'navy', 'black']

#use all_m_zero for marginals made from reflections which started at origin
all_m = all_m_far
plt.pause(1)
for i in range(len(all_m)):
    draw_it(all_m[i], c[i], 1)
    a = a + 0.8 / len(all_m)

    plt.pause(1.0)

plt.pause(10)
