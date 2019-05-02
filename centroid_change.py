import matplotlib.pyplot as plt
import numpy as np

def centroid(pts):
    return sum(pts) / len(pts)
def cent(pts):
    return centroid(pts)

def line(p1, p2, t):
    l = []
    for tt in t:
        l.append(p1 + (p2 - p1) * tt)
    return l

def f(pts):
    p = []
    change = []
    for i in range(len(pts) - 1):
        cur = line(np.array(pts[i]), np.array(pts[i+1]), np.arange(0,1.1,1))
        for c in cur:
            p.append(c)
        change.append(len(p))
    return p, change

def draw_centroid(pts, col = 'b'):
    pts, change = f(pts)
    a = 0.2

    b = 0
    for p in range(len(pts)):
        cur = centroid(pts[:p+1])
        plt.scatter(cur[0], cur[1], color = col, alpha = a)

        '''
        if p+1 == change[b]:
            plt.scatter(cur[0], cur[1], color = 'red', alpha = 1)
            b+=1
        '''

        a += .8 / len(pts)

def draw_lines(pts):
    pts, change = f(pts)
    # print(len(pts))
    for i in range(len(pts)-1):
        #print(i)
        plt.plot([pts[i][0], pts[i+1][0]], [pts[i][1], pts[i+1][1]], color = 'green')

def test():
    plt.plot([0,2], [0.5,0.5], color = 'black')

    pts = ([0,0], [1,1], [2,0], [3,1], [4,0])
    draw_centroid(pts)
    draw_lines(pts)
    plt.show()

def symmetric_square_test():
    recur = ([0,-1], [1,0], [0,1], [-1, 0])
    pts = []
    for i in range(3):
        for r in recur:
            pts.append(r)
    draw_centroid(pts)
    draw_lines(pts)
    plt.show()

#symmetric_square_test()

def fun():

    import random

    pts = []
    for i in range(100):
        pts.append((np.random.random(), np.random.random()))
        plt.scatter(pts[-1][0], pts[-1][1], color = 'g')

    pts = sorted(pts , key=lambda k: [k[0], k[1]])
    draw_centroid(pts)
    plt.show()
