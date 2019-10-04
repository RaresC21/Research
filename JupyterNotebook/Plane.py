import numpy as np
from create_basis import orthogonal_extend_basis

class Options:
    def __init__(self, reflection = "facet", delta = 0.1, radius = 0.1 * np.sqrt(3), depth = 0):
        self.reflection = reflection
        self.delta = delta
        self.radius = radius
        self.depth = depth

class Plane:
    class Ball:
        def __init__(self, center, radius):
            self.center = center
            self.radius = radius

        def pt_inside(self, pt):
            dist = np.linalg.norm(pt - self.center)

            if dist < self.radius + 1e-8:
                return True
            return False

        def does_intersect(self, o, l) :
            c = self.center
            r = self.radius

            l = l / np.linalg.norm(l)

            a = - (np.dot(l, o - c))
            b = (np.dot(l, (o - c)) ** 2) - (np.linalg.norm(o - c) ** 2) + r*r
            return b >= -1e-8

        '''
        circle centered at c, of radius r
        line: x = o + d * l (origin o, direction l)
        '''
        def line_intersection(self, o, l):
            c = self.center
            r = self.radius

            l = l / np.linalg.norm(l)

            a = - (np.dot(l, o - c))
            b = (np.dot(l, (o - c)) ** 2) - (np.linalg.norm(o - c) ** 2) + r*r
            b = np.sqrt(b)

            d1 = a + b
            d2 = a - b

            p1 = o + d1 * l
            p2 = o + d2 * l

            if np.linalg.norm(p1 - o) < np.linalg.norm(p2 - o):
                return p1
            return p2

        def reflection(self, pt, direction):
            hit_on_sphere = self.line_intersection(pt, direction)

            perp = hit_on_sphere - self.center
            b = np.dot(perp, hit_on_sphere)
            plane = Plane(perp, b)

            p, d = plane.reflect(pt, direction)
            return p, d

    def __init__(self, coef, b, options = Options(), plane_id = -1):
        norm = np.linalg.norm(coef)
        self.perp_vec = coef / norm
        self.b = b / norm
        self.dimension = coef.shape[0]
        self.options = options
        self.plane_id = plane_id
        self.basis = orthogonal_extend_basis(self.perp_vec)[1:]
        self.Q = np.transpose(np.vstack([self.basis, [self.perp_vec]]))
        self.Qinv = np.transpose(self.Q) #np.linalg.inv(self.Q)

        pt_on_plane = [0] * self.dimension
        for i, c in np.ndenumerate(coef):
            indx = i[0]
            if abs(c) > 1e-8:
                pt_on_plane[indx] = b / c
                break
        self.pt_on_plane = pt_on_plane

    def pt_in_half_plane(self, pt):
        return np.dot(pt, self.perp_vec) < self.b + 1e-8

    def pt_on_hyperplane(self, pt):
        diff = abs(np.dot(pt, self.perp_vec) - self.b)
        return diff < 1e-8

    def line_intersection(self, point, direction) :
        num = self.b - np.dot(point, self.perp_vec)
        denom = np.dot(direction, self.perp_vec)
        t = num / denom
        if t < -1e-8: t = 1e5
        return point + direction * t

    def projection(self, point):
        return self.line_intersection(point, 1 * self.perp_vec);

    def reflect(self, point, direction):
        # return self.reflect_function(point, direction)

        # move everything to origin
        intersection = self.line_intersection(point, direction)
        pp = point - self.pt_on_plane

        dd = direction

        # mirror the point across the hyperplane
        proj = self.perp_vec * np.dot(pp, self.perp_vec) \
               / np.dot(self.perp_vec, self.perp_vec)

        mirror = pp - 2 * proj + self.pt_on_plane
        direction = intersection - mirror

        direction = direction / np.linalg.norm(direction)

        return (intersection, direction)

    def reflect_random(self, point, direction): 
        intersection, _ = self.reflect(point, direction)
        direction = np.random.normal(loc = 0, scale = 1.0, size = self.dimension - 1) 
        direction /= np.linalg.norm(direction) 
        
        rad = np.random.uniform(0,1,1) ** (1.0 / (self.dimension - 1))
        direction *= rad
        
        perp = np.sqrt(1-rad*rad)
        
        direction = np.append(direction, [-perp])
        #direction = np.vstack([direction, [perp]])
    
        # change to be on hyperplane
        direction = np.matmul(self.Q, direction)
        return intersection, direction / np.linalg.norm(direction)
        
    def get_pt_in_basis(self, point) :
        projection = self.projection(point) - self.pt_on_plane;

        temp = np.vstack([self.basis, [self.perp_vec]])
        pt_plane_basis = np.linalg.solve(np.transpose(temp), projection)

        pt = np.rint(pt_plane_basis / self.options.delta)
        return pt * self.options.delta

    def reflect_function(self, point, direction): # to fix
        intersection = self.line_intersection(point, direction)
        if (np.linalg.norm(intersection) > 1e5):
            return (intersection, np.zeros(self.dimension))

        on_grid = self.get_pt_in_basis(intersection)
        pt = self.get_pt_in_basis(self.nearest_gridpoint(intersection))
        refl_direction = np.zeros(self.dimension)
        for i in range(self.dimension - 1):
            p = on_grid[i]
            x1 = pt[i]
            x2 = pt[i] + self.options.delta
            if p - x1 < x2 - p:
                d = p - x1
                if d < self.options.delta / 4:
                    refl_direction += self.basis[i] * d
                else:
                    refl_direction += self.basis[i] * ((self.options.delta / 2) - d)
            else :
                d = x2 - p
                if d < self.options.delta / 4:
                    refl_direction += self.basis[i] * (-d)
                else:
                    refl_direction += self.basis[i] * (((self.options.delta / 2) - d) * (-1))
        refl_direction -= self.perp_vec * self.options.delta;
        refl_direction /= np.linalg.norm(refl_direction)
        return (intersection, refl_direction)

    def offset(self):
        h = np.sqrt(self.options.radius**2 -(self.options.delta/2)**2)
        return h - self.options.depth

    def translate_plane(self, h) :
        new_pt = self.pt_on_plane + self.perp_vec * h
        return Plane(self.perp_vec, np.dot(self.perp_vec, new_pt), self.options)

    # q is already in coordinates of the grid
    def nearest_gridpoint(self, q):
        return np.rint(q / self.options.delta) * self.options.delta

    def standard_coords(self, point):
        return np.matmul(self.Q, point) + self.pt_on_plane

    def grid_coords(self, point):
        point = self.projection(point) - self.pt_on_plane
        return np.matmul(self.Qinv, point)

    # given point p on plane grid, find point on line (point, direction) whose projection
    # onto the plane is p itself
    def get_pt_on_line(self, p, point, direction):
        b = self.dimension
        a = b - 2
        x = p[a:b]

        d = [0,1]

        y = point[a:b]
        c = direction[a:b]

        v = abs((np.cross(x, d) - np.cross(y, d)) / np.cross(c, d))
        return point + direction * v

    # assuming line is in coordinates of the grid.
    def line_to_plane_grid(self, point, direction):
        direction[-1] = 0
        return point, direction

    # find first point along the line whose nearest ball is different from the nearest ball of the current point
    # return this ball
    def next_nearest_ball(self, ball_plane, point, direction):
        pp = point
        dd = direction
        current_ball = ball_plane.nearest_gridpoint(point)
        direction /= np.linalg.norm(direction)

        l = 0
        r = self.options.delta * self.dimension
        ans = None
        t = 0
        while (abs(r - l) > 1e-6):
            mid = (l + r) / 2
            ball = ball_plane.nearest_gridpoint(point + direction * mid)
            d = ball - current_ball
            if np.dot(d, d) > (self.options.delta / 2) ** 2:
                r = mid
                t = mid
                ans = ball
            else:
                l = mid
        return self.Ball(ans, self.options.radius), point + direction * t

    def reflect_sphere(self, point, direction, plane_id = -1):
        ball_plane = self

        onBall = False
        if plane_id == self.plane_id:
            onBall = True

        d = np.linalg.norm(point - self.projection(point))
        if d < self.options.radius + 1e-9:
            q = point
        else:
            off = self.options.radius #self.offset()
            boundary = self.translate_plane(-off)

            q = boundary.line_intersection(point, direction)
#             if np.linalg.norm(q) > 1e3:
#                 return (np.array([1e8] * self.dimension), 1e8)

        direction_ = np.matmul(self.Qinv, direction)
        point_ = np.matmul(self.Qinv, point - ball_plane.pt_on_plane)

        qq = np.matmul(self.Qinv, q - ball_plane.pt_on_plane)
        q_ = ball_plane.grid_coords(q)
        d_ = direction_.copy()
        d_[-1] = 0
        d_ /= np.linalg.norm(d_)

        ball = self.Ball(ball_plane.nearest_gridpoint(q_), self.options.radius)
        
        counter = 0
        dist = []
        while True:
            counter += 1
            dist_to_plane = ball_plane.get_pt_on_line(q_, point_, direction_)
            if abs(dist_to_plane[-1]) > 1e-4 + self.options.radius: 
                break

            if onBall:
                ball, q_ = ball_plane.next_nearest_ball(ball_plane, q_, d_)
                onBall = False
            
            if ball.does_intersect(point_, direction_):
                intersection, new_direction = ball.reflection(point_ - direction_, direction_)
                cur_intersection = intersection.copy()
                intersection[-1] = 0 # 'project'
                dif = (intersection[0] - point_[0]) / direction_[0]
                if np.linalg.norm(point_ - intersection) > 1e5:  break
                if dif > -1e-9: # checks if the intersection is in the correct direction
                    other = ball_plane.nearest_gridpoint(intersection)
                    if np.linalg.norm(other - ball.center) < 1e-5:
                        # this checked if the point we found is not inside of any other circle.
                        return self.standard_coords(cur_intersection), np.matmul(self.Q, new_direction)
#                 else: break  

            ball, q_ = self.next_nearest_ball(ball_plane, q_, d_)
            
        return (np.array([1e8] * self.dimension), 1e8)
