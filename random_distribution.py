import numpy as np

def find_distribution(polytope, iterations, start_pos = -1):
    dimension = polytope.dimension

    standard_directions = []
    for i in range(dimension):
        cur = np.zeros(dimension)
        cur[i] = 1
        standard_directions.append(cur)
        other = np.zeros(dimension)
        other[i] = -1
        standard_directions.append(other)

    if start_pos == -1:
        p = polytope.point
    else:
        p = start_pos

    lines = []
    for i in range(iterations):
        dir_indx = np.random.randint(len(standard_directions))
        direction = standard_directions[dir_indx]

        hit, _ = polytope.get_plane_hit(p, direction)

        # choose random point from p to hit
        dist = np.linalg.norm(hit - p)
        length = np.random.uniform(0, dist)

        next_pt = p + direction * length
        lines.append([p, next_pt])
        p = next_pt

    return lines
