# %load ../Polytope.py

from Plane import Plane, Options
from scipy.optimize import linprog
import numpy as np
# from Plane import Plane, Options
from utils import marginal

import time

class Polytope:
    def __init__(self, A, b, pt, options = Options()):
        self.A = A
        self.b = b
        self.dimension = A.shape[1]
        self.num_planes = A.shape[0]
        self.options = options

        self.planes = []
        for i in range(self.num_planes):
            p = Plane(A[i], b[i], self.options, plane_id = i)
#             print(p.pt_on_plane)
            p = self.translate_plane(p, options.depth, i)
            print(p.perp_vec, p.pt_on_plane)
            self.planes.append(p)

        self.point = pt

        self.cur_lines__ = [[self.point, self.point]]
        self.pts = []
        self.averages = []
        self.centroid_distance = []
        self.centroid_change = []
        self.facet_hit = []
        self.angles = []
        #self.marginal_direction = np.zeros(self.dimension)
        #self.marginal_direction[0] = 1
        self.marginal_direction = np.random.random(self.dimension)
        self.marginal_direction /= np.linalg.norm(self.marginal_direction)
        self.last_marginal = []
        self.it = 10
        self.all_marginals = []
        self.all_marginal_diffs = []

    def translate_plane(self, plane, h, i) :
        new_pt = plane.pt_on_plane + plane.perp_vec * h
#         print(new_pt, '--')
        p = Plane(plane.perp_vec, np.dot(plane.perp_vec, new_pt), self.options, plane_id = i)
        p.pt_on_plane = new_pt
        return p

        
    def fix_marginal_direction(self, d): 
        self.marginal_direction = d / np.linalg.norm(d)
        
    def pt_in_polytope(self, pt):
        for plane in self.planes:
            if not plane.pt_in_half_plane(pt):
                return False
        return True

    def flat_plane_hit(self, pt, direction): 
        dist = 0
        indx = -1

        intersection = np.zeros(self.dimension)
        d = direction

        for i, plane in np.ndenumerate(self.planes):
            if plane.pt_on_hyperplane(pt):
                continue

            hit, v = plane.reflect(pt, direction)
            
            cur_dist = np.linalg.norm(pt - hit)
            if self.pt_in_polytope(hit) and (indx == -1 or dist > cur_dist):
                indx = i[0]
                dist = cur_dist
                intersection = hit
                d = v

#         print(direction, d)
        return intersection, d, indx

    def random_hit(self, pt, direction): 
        intersection, d, indx = self.flat_plane_hit(pt, direction)
        intersection, d = self.planes[indx].reflect_random(pt, direction)
        return intersection, d, indx
    
    def ball_plane_hit(self, pt, direction, last_plane_hit): 
        flat_intersection, _, _ = self.flat_plane_hit(pt, direction) 
        
        indx = -1
        dist = 0
        d = np.ones(self.dimension)
        intersection = np.zeros(self.dimension)
        
        for i, plane in np.ndenumerate(self.planes): 
#             if np.linalg.norm(pt - plane.projection(pt)) > self.options.radius and np.linalg.norm(flat_intersection - plane.projection(flat_intersection)) > self.options.radius:
#                 continue
            hit, v = plane.reflect_sphere(pt, direction, plane_id = last_plane_hit)
            cur_dist = np.linalg.norm(pt - hit) 
            if self.pt_in_polytope(hit) and (indx == -1 or dist > cur_dist):
                indx = i[0]
                dist = cur_dist
                intersection = hit
                d = v
        if np.linalg.norm(intersection) < 1e-5: 
            print("NOO")
        return intersection, d, indx
    
    def get_plane_hit(self, pt, direction, last_plane_hit):
        if self.options.reflection == "facet":
            return self.flat_plane_hit(pt, direction)
        elif self.options.reflection == "sphere":              
            return self.ball_plane_hit(pt, direction, last_plane_hit)
        elif self.options.reflection == "random": 
            return self.random_hit(pt, direction)
        return "Error Option" 

    def reflection_lines(self, iterations, convergence_eps = 0.01, continue_ = False):
        if not continue_:
            direction = np.random.random(self.dimension)
#             direction = np.ones(self.dimension)
#             direction = np.array([5,17])
#             print(direction / np.linalg.norm(direction))
            direction = np.zeros(self.dimension)
            direction[0] = 1
            direction += np.random.random(self.dimension) / 10
            direction = direction / np.linalg.norm(direction)
            pt = self.point
            average = np.zeros(self.dimension)
            length = 0
        else:
            direction = self.reflection_direction
            pt = self.reflection_point
            average = self.current_average
            length = self.current_length
            
        hit_average = np.zeros(self.dimension)

        indx = -1
        line = []
        last_plane_hit = -1
        for i in range(iterations):
            nxt_pt, nxt_direction, last_plane_hit = self.get_plane_hit(pt, direction, last_plane_hit)
            self.facet_hit.append(last_plane_hit)
            self.angles.append(np.pi/2 - np.arccos(np.dot(direction, self.planes[last_plane_hit].perp_vec) / np.linalg.norm(direction)))

            line.append((pt, nxt_pt))
            self.cur_lines__.append(line[-1])

            l = np.linalg.norm(nxt_pt - pt)

            mid_point = (nxt_pt + pt) / 2
            
            prev_average = average
            average = (average * length + l * mid_point) / (l + length)
            hit_average += nxt_pt
            length += l

            pt = nxt_pt

            dd = np.linalg.norm(average)
            self.pts.append(pt)
            self.centroid_distance.append(dd)
            self.centroid_change.append(np.linalg.norm(average - prev_average))
            self.averages.append(average)

            res = self.different_marginal(i, convergence_eps)
            '''
            if res == False: 
                return line
            '''

            if i % 1 == 0:
                print(i, dd)
                
            direction = nxt_direction
        self.reflection_direction = direction
        self.reflection_point = pt
        self.current_length = length
        self.current_average = average
        return line

    def compute_average(self, lines) :
        centroid = np.zeros(self.dimension)
        line_lengths = 0
        for l in lines:
            cur_len = np.linalg.norm(l[0] - l[1])
            cur_pt = (l[0] + l[1]) / 2.0
            centroid = centroid + cur_pt * cur_len
            line_lengths += cur_len
        return centroid / line_lengths

    def approximate_centroid(self, iterations = 1000, get_all = False):
        centroid = self.point
        average_centroid = centroid
        averages = [centroid]
        cents = [centroid]

        lines = []
        pt = self.point
#         print("START", pt)
        direction = np.random.uniform(-1,1,self.dimension)
#         print("initial direction", direction)

        it = 1
        average_lines = pt
        total_len = 0
        last_plane_hit = -1

        cnt = 1
        while cnt < iterations:
            cnt+=1

            next_pt, next_d, last_plane_hit = self.get_plane_hit(pt, direction, last_plane_hit)
#             if np.linalg.norm(next_pt) < 1e-5:
                # lines = lines[-10:]
                # lines.append((pt, pt + next_d))
#                 break
            lines.append([pt, next_pt])

            if it == 1:
                centroid = (lines[0][0] + lines[0][1]) / 2.0
                l = np.linalg.norm(lines[0][0] - lines[0][1])
                average_lines = centroid * l
                total_len = l
            else:
                c = (lines[-1][0] + lines[-1][1]) / 2.0
                l = np.linalg.norm(lines[-1][0] - lines[-1][1])
                centroid = (centroid * total_len + l * c) / (total_len + l)
                total_len += l

            next_centroid = (average_centroid * it + centroid) / (it + 1)
            averages.append(next_centroid)
            cents.append(centroid)

            ddd = np.sqrt(np.dot(centroid, centroid))
            self.centroid_distance.append(ddd)
            # print("pt, cnt", cnt, next_pt, next_d)
            # for plane in self.planes:
            #     ball = plane.nearest_gridpoint(next_pt)
            #     print("dist", np.linalg.norm(ball - next_pt), plane.options.radius)
            # print("direction", next_d)
#             print("centroid", cnt, ddd)

            # diff = np.linalg.norm(next_centroid - average_centroid)
            # if not self.different_marginal(lines, cnt) and cnt > 100:
                # print("iterations:", cnt)
                # break

            '''
            if cnt == iterations or (iterations > 1e10 and diff < 0.0001):
                print("iterations:", it)
                break
            '''

            average_centroid = next_centroid
            direction = next_d
            pt = next_pt
            it += 1

        if get_all:
            return averages, cents
        else:
            n = 10
            m = np.zeros(len(averages[0]))
            for a in cents[-n:]:
                m = m + a
            m = m / len(cents[-n:])
            return m, lines #self.cur_lines__

    def is_redundant(self, vec, b):
        A = np.copy(self.A)
        bb = np.copy(self.b)

        A = np.vstack([A, vec])
        bb = np.append(bb, b + 1)

        c = -1.0 * vec

        result = linprog(c, A_ub=A, b_ub=bb, bounds=(None,None))
        if result.success:
            if -result.fun <= b:
                return True
            else:
                return False

    def remove_redundant(self):
        i = 0
        while i < self.num_planes:
            tempA = np.delete(self.A, (i), axis = 0)
            tempB = np.delete(self.b, (i))

            poly = Polytope(tempA, tempB, self.point)
            poly.options = self.options

            if poly.is_redundant(self.A[i], self.b[i]):
                self = poly
            i += 1
        return self

    def different_marginal(self, it, eps = 0.01):
        if it % 100 != 0: return True
        # self.randomized_distribution(it, self.cur_lines__[-1][1])
        # m = self.cur_lines__
        # m = self.marginal_dist(self.cur_lines__)
        m = self.marginal_dist()
        if len(self.all_marginals) == 0:
            self.all_marginals.append(m)
            return True

        # np.set_printoptions(threshold=np.nan)
        # print(np.array2string(m, separator=', '))

        #print("distribution:", m)

        l = len(self.all_marginals)
#         print(self.all_marginals[l//2])
#         print(m)
#         print()
#         print()
        diff = 0
        for i in range(len(m)):
            a = m[i]
            b = self.all_marginals[l//2][i]
            diff = max(diff, abs(a - b))

#         print("cnt", it)
#         print("marginal change", diff)
        self.all_marginal_diffs.append(diff)
        self.all_marginals.append(m)
#         print("marginal diff", it, diff)
        if (diff < eps) : return False
        return True

    def marginal_dist(self):
        num_steps = 100

#         marginal_set = []
#         s = 1/(num_steps + 1)
#         for p in np.arange(s,1, s): 
#             l = -2
#             r = 2
#             while abs(l - r) > 1e-2: 
#                 mid = (l + r) / 2
#                 cut = Plane(self.marginal_direction, np.dot(self.marginal_direction * mid, self.marginal_direction))
#                 cur = marginal(self.cur_lines__, cut)
#                 if cur > p: 
#                     r = mid
#                 else: 
#                     l = mid
#             marginal_set.append(l)
#         print(marginal_set)
        dist = []
        for p in self.pts:
            dist.append(np.dot(p, self.marginal_direction))
        dist.sort()

        marginal_set = []
        for i in range(1, num_steps):
            indx = len(dist) * i // num_steps
            marginal_set.append(dist[indx])
        return np.array(marginal_set)

    def randomized_distribution(self, iterations):
        print(iterations)
        dimension = self.dimension

        standard_directions = []
        for i in range(dimension):
            cur = np.zeros(dimension)
            cur[i] = 1
            standard_directions.append(cur)
            other = np.zeros(dimension)
            other[i] = -1
            standard_directions.append(other)

        p = self.point

        lines = []
        for i in range(iterations):

            change = 0

            while change < 1e-6:

                dir_indx = np.random.randint(len(standard_directions))
                direction = standard_directions[dir_indx]

                hit, _, _ = self.get_plane_hit(p, direction, last_plane_hit = -1)

                change = np.linalg.norm(hit - p)

            dist = change
            length = np.random.uniform(0, dist)

            next_pt = p + direction * length
            lines.append([p, next_pt])
            p = next_pt
            self.pts.append(p)

            self.different_marginal(i+1)

        for l in lines:
            self.cur_lines__.append(l)
        # print(self.cur_lines__)
