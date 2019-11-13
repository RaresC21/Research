import numpy as np

class Plane:
    def __init__(self, coef, b, plane_id = -1):
        norm = np.linalg.norm(coef)
        self.perp_vec = coef / norm
        self.b = b / norm
        self.dimension = coef.shape[0]
        self.plane_id = plane_id

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

class Ball:
    def __init__(self, center, radius, index):
        self.center = center
        self.radius = radius
        self.index = index

    def pt_inside(self, pt):
        dist = np.linalg.norm(pt - self.center)

        if dist < self.radius + 1e-8:
            return True
        return False

    def does_intersect(self, o, l) :
        c = self.center
        r = self.radius
        
        if np.dot(c - o, l) < 0: return False

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
    
    def angle(self, pt, direction): 
        v = pt - self.center
        v /= np.linalg.norm(v)
        d = direction / np.linalg.norm(direction)
        a = np.arccos(np.dot(v, d))
        if (np.cross(v, d)) > 0: 
            a = -a
        return a

def run(N, balls, index, pt, d, want_index = False): 
    
    indexes = []
    angles = []
    pts = []
    for _ in range(N): 
        
        dist = 1e9
        hit_ball = None
        for ball in balls: 
            if ball.index == index: 
                continue
                
            idx = -1
            if ball.does_intersect(pt, d): 
                hit = ball.line_intersection(pt, d)
                if np.linalg.norm(hit - pt) < dist: 
                    dist = np.linalg.norm(hit - pt)
                    hit_ball = ball
            
        index = hit_ball.index
        indexes.append(index)
        pt, d = hit_ball.reflection(pt, d)
                  
        pts.append(pt)
        angles.append(hit_ball.angle(pt, d))
        
    if want_index: 
        return pts, angles, indexes
    return pts, angles


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def circ(x, y, r, phi):
    return (x + r * np.cos(phi), y + r*np.sin(phi))

def plot_circle(x, y, r, h):
    phis = np.arange(0,6.28,0.01)
    plt.plot(*circ(x, y, r, phis), zs = h, c = 'black')

    # plt.xlim(-1, 1)
    # plt.ylim(-1, 1)
    # plt.gca().set_aspect('equal', adjustable='box')

def draw_balls(balls, h): 
    for ball in balls: 
        plot_circle(ball.center[0], ball.center[1], ball.radius, h)
        
def plot_orbit(pts): 
    for i in range(len(pts) - 1):
        plt.plot([pts[i][0], pts[i+1][0]], [pts[i][1], pts[i+1][1]], color = 'blue')

def unstable_curve(n, balls):
    all_pts = []
    all_angles = []
    for i in np.arange(-0.001, 0.001, 1e-5): 
        pt = np.array([-0.2, i])
        pts, angles = run(n, balls, -1, pt, np.array([1, 0]))
        all_pts.append(pts)
        all_angles.append(angles)


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
        
    for pt, angle in zip(all_pts, all_angles): 
        # if balls[0].pt_inside(pt[-1]): ## change if you want to see on a different ball
        ax.scatter(xs = [pt[-1][0]], ys = [pt[-1][1]], zs = [angle[-1]], color = 'blue')

    draw_balls(balls, -np.pi/2)
    draw_balls(balls, np.pi/2)

    # plt.gca().set_aspect('equal', adjustable='box')

    ax.set_xlim3d(-1,1)
    ax.set_ylim3d(-10,1)

    # plt.xlim(-1, 1)
    # Create cubic bounding box to simulate equal aspect ratio
    # max_range = np.array([2, 12, 3]).max()
    # Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(X.max()+X.min())
    # Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Y.max()+Y.min())
    # Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Z.max()+Z.min())
    # # Comment or uncomment following both lines to test the fake bounding box:
    # for xb, yb, zb in zip(Xb, Yb, Zb):
    #    ax.plot([xb], [yb], [zb], 'w')



    # ax.set_zlim3d(-np.pi/2,np.pi/2)

    plt.show()

radius = 0.75
# balls = [Ball(np.array([-1, 0]), radius, 0), Ball(np.array([1, 0]), radius, 1), Ball(np.array([0, -1]), radius, 2), Ball(np.array([0, 1]), radius, 3)]
balls = [Ball(np.array([-1, 0]), radius, 0), Ball(np.array([0, -1]), radius, 1), Ball(np.array([5, 5]), 7.1, 2)]
# balls = [Ball(np.array([0, -8]), 1, 0), Ball(np.array([-100, 5]), 100, 1), Ball(np.array([100, 5]), 100, 2)]


# unstable_curve(100, balls)

cols = ['red', 'green', 'blue']
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# boundary = [circ(balls[0].center[0], balls[0].center[1], balls[0].radius, np.arange(-0.4, 0.3, 1e-4))]
            # circ(balls[1].center[0], balls[1].center[1], balls[1].radius, np.arange(1.2, 2.0, 1e-1))]
def singularities(boundary, balls, N, i): 
    for pt in zip(boundary[0], boundary[1]): 
        for a in [-np.pi / 2, np.pi / 2]: 
            n = pt - balls[i].center
            d = np.array([np.cos(a) * n[0] - np.sin(a) * n[1], np.sin(a) * n[0] - np.cos(a) * n[1]])
            
            res = run(N, balls, balls[i].index, np.array(pt), d, want_index = True)

            if (res[2][-1] == 0):
                ax.scatter(res[0][-0][0], res[0][-1][1], res[1][-1], color = cols[res[2][-1]], s = 0.5)

    for a in np.arange(-np.pi/2, np.pi/2): 
        # pt = 

        n = pt - balls[i].center
        d = np.array([np.cos(a) * n[0] - np.sin(a) * n[1], np.sin(a) * n[0] - np.cos(a) * n[1]])

        try :         
            res = run(N, balls, balls[i].index, np.array(pt), d, want_index = True)

            if (res[2][-1] == 0):
                ax.scatter(res[0][-0][0], res[0][-1][1], res[1][-1], color = cols[res[2][-1]], s = 0.5)
        except: 
            continue

N = 2
for n in range(1,N+1):
    singularities(circ(balls[0].center[0], balls[0].center[1], balls[0].radius, np.arange(-0.4, 0.3, 1e-3)), balls, n, 0)
    singularities(circ(balls[1].center[0], balls[1].center[1], balls[1].radius, np.arange(1.2, 2.0, 1e-3)), balls, n, 1)
    singularities(circ(balls[2].center[0], balls[2].center[1], balls[2].radius, np.arange(3.88, 3.98, 1e-3)), balls, n, 2)

draw_balls(balls, -np.pi/2)
draw_balls(balls, np.pi/2)
ax.set_xlim3d(-1,1)
ax.set_ylim3d(-1,1)
ax.set_zlim3d(-np.pi/2,np.pi/2)

plt.show()
