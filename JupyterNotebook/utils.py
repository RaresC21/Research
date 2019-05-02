import numpy as np

from Plane import Plane

def get_all(lines, cut):
    centroid = np.zeros(len(lines[0][0]))

    best = 0
    pt = 0

    a = 0
    b = 0
    for line in lines:
        pt1 = line[0]
        pt2 = line[1]

        cur = (pt1 + pt2) / 2
        if np.linalg.norm(cur) > best:
            best = np.linalg.norm(cur)
            pt = cur

        side1 = cut.pt_in_half_plane(pt1)
        side2 = cut.pt_in_half_plane(pt2)

        if side1 == side2:
            dist = np.linalg.norm(pt1 - pt2)
            if side1 :
                centroid += dist * (pt1 + pt2) / 2
                a += dist
            else : b += dist
        else :
            direction = pt2 - pt1
            hit = cut.line_intersection(pt1, direction)
            if side1 :
                aa = np.linalg.norm(hit - pt1)
                bb = np.linalg.norm(hit - pt2)
                centroid += aa * (hit + pt1) / 2
                a += aa
                b += bb
            else :
                aa = np.linalg.norm(hit - pt2)
                bb = np.linalg.norm(hit - pt1)
                centroid += aa * (hit + pt2) / 2
                a += aa
                b += bb

    if b == 0:
        return 0, 1, 0
    if a == 0:
        return 0, 0, 0
    return a / b, a / (a + b), centroid / a

def get_ratio(lines, cut) :
    a, b, c = get_all(lines, cut)
    return a, c

def marginal(lines, cut):
    a, b, c = get_all(lines, cut)
    return b
