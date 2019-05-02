{
  "cells": [
    {
      "metadata": {
        "trusted": true
      },
      "cell_type": "code",
      "source": "# %load ../Polytope.py\n# %run Plane.ipynb\n\nfrom Polytope import Polytope\nfrom scipy.optimize import linprog\nimport numpy as np\n# from Plane import Plane, Options\n# from Utils import marginal\n\nimport time\n\nclass Polytope:\n    def __init__(self, A, b, pt, options = Options()):\n        self.A = A\n        self.b = b\n        self.dimension = A.shape[1]\n        self.num_planes = A.shape[0]\n        self.options = options\n\n        self.planes = []\n        for i in range(self.num_planes):\n            p = Plane(A[i], b[i], self.options, plane_id = i)\n            self.planes.append(p)\n\n        self.point = pt\n\n        self.cur_lines__ = [[self.point, self.point]]\n        self.pts = []\n        self.averages = []\n        self.centroid_distance = []\n        self.centroid_change = []\n        self.facet_hit = []\n        self.angles = []\n        #self.marginal_direction = np.zeros(self.dimension)\n        #self.marginal_direction[0] = 1\n        self.marginal_direction = np.random.random(self.dimension)\n        self.last_marginal = []\n        self.it = 10\n        self.all_marginals = []\n        self.all_marginal_diffs = []\n\n    def marginal_width(self):\n        # how far can I go in the direction and still be in the polytope?\n        d = 0.01\n        while True:\n            A = np.vstack([self.A, self.marginal_direction])\n            b = np.append(self.b, -np.dot(self.marginal_direction,\n                                         self.point + d * self.marginal_direction))\n            res = linprog(self.marginal_direction, A_ub=A, b_ub=b, bounds=(None,None))\n            if res.status == 2:\n                break\n            else:\n                d = d * 2\n\n        d1 = d\n        d = 0.01\n        while True:\n            A = np.vstack([self.A, self.marginal_direction])\n            b = np.append(self.b, np.dot(self.marginal_direction,\n                                         self.point - d * self.marginal_direction))\n            res = linprog(self.marginal_direction, A_ub=A, b_ub=b, bounds=(None,None))\n            if res.status == 2:\n                break\n            else:\n                d = d * 2\n\n        return d1/2, d/2\n\n    def pt_in_polytope(self, pt):\n        for plane in self.planes:\n            if not plane.pt_in_half_plane(pt):\n                return False\n        return True\n\n    def get_plane_hit(self, pt, direction, last_plane_hit):\n        dist = 0\n        indx = -1\n\n        intersection = np.zeros(self.dimension)\n        d = direction\n\n        for i, plane in np.ndenumerate(self.planes):\n            if plane.pt_on_hyperplane(pt):\n                continue\n\n            if self.options.reflection == \"facet\":\n                hit, v = plane.reflect(pt, direction)\n            elif self.options.reflection == \"sphere\":\n                hit, v = plane.reflect_sphere(pt, direction, plane_id = last_plane_hit)\n                print(hit)\n\n            # if (np.linalg.norm(hit) > 1e5) :continue\n            cur_dist = np.linalg.norm(pt - hit)\n            # print(\"inside polytope\", self.pt_in_polytope(hit), hit)\n            if self.pt_in_polytope(hit) and (indx == -1 or dist > cur_dist):\n                indx = i[0]\n                dist = cur_dist\n                intersection = hit\n                d = v\n        return intersection, d, indx\n\n    def reflection_lines(self, iterations, continue_ = False):\n        if not continue_:\n            direction = np.random.random(self.dimension)\n            pt = self.point\n        else:\n            direction = self.reflection_direction\n            pt = self.reflection_point\n\n        average = np.zeros(self.dimension)\n        length = 0\n\n        indx = -1\n        line = []\n        last_plane_hit = -1\n        for i in range(iterations):\n            nxt_pt, nxt_direction, last_plane_hit = self.get_plane_hit(pt, direction, last_plane_hit)\n            self.facet_hit.append(last_plane_hit)\n            self.angles.append(np.pi/2 - np.arccos(np.dot(direction, self.planes[last_plane_hit].perp_vec) / np.linalg.norm(direction)))\n            if continue_:\n                line = (pt, nxt_pt)\n            else:\n                line.append((pt, nxt_pt))\n\n            l = np.linalg.norm(nxt_pt - pt)\n\n            prev_average = average\n            average = (average * length + l * nxt_pt) / (l + length)\n            length += l\n\n            pt = nxt_pt\n\n            dd = np.linalg.norm(average)\n            self.pts.append(pt)\n            self.centroid_distance.append(dd)\n            self.centroid_change.append(np.linalg.norm(average - prev_average))\n            self.averages.append(average)\n            #print(i, dd)\n            direction = nxt_direction\n        self.reflection_direction = direction\n        self.reflection_point = pt\n        return line\n\n    def compute_average(self, lines) :\n        centroid = np.zeros(self.dimension)\n        line_lengths = 0\n        for l in lines:\n            cur_len = np.linalg.norm(l[0] - l[1])\n            cur_pt = (l[0] + l[1]) / 2.0\n            centroid = centroid + cur_pt * cur_len\n            line_lengths += cur_len\n        return centroid / line_lengths\n\n    def approximate_centroid(self, iterations = 1000, get_all = False):\n        centroid = self.point\n        average_centroid = centroid\n        averages = [centroid]\n        cents = [centroid]\n\n        lines = []\n        pt = self.point\n        print(\"START\", pt)\n        #direction = np.ones(self.dimension) + np.random.randn(self.dimension)/5\n#         direction = np.random.random(self.dimension)\n        direction = np.random.uniform(-1,1,self.dimension)\n#         direction = np.array([0.05019, 0.36692])\n        # direction = np.array([-0.70499, -0.76752, 0.76134, -0.20046, -1.57117])\n        print(\"initial direction\", direction)\n\n        it = 1\n        average_lines = pt\n        total_len = 0\n        last_plane_hit = -1\n\n        cnt = 1\n        while cnt < iterations:\n            cnt+=1\n\n            next_pt, next_d, last_plane_hit = self.get_plane_hit(pt, direction, last_plane_hit)\n            if np.linalg.norm(next_pt) < 1e-5:\n                # lines = lines[-10:]\n                # lines.append((pt, pt + next_d))\n                break\n            lines.append([pt, next_pt])\n\n            if it == 1:\n                centroid = (lines[0][0] + lines[0][1]) / 2.0\n                l = np.linalg.norm(lines[0][0] - lines[0][1])\n                average_lines = centroid * l\n                total_len = l\n            else:\n                c = (lines[-1][0] + lines[-1][1]) / 2.0\n                l = np.linalg.norm(lines[-1][0] - lines[-1][1])\n                centroid = (centroid * total_len + l * c) / (total_len + l)\n                total_len += l\n\n            next_centroid = (average_centroid * it + centroid) / (it + 1)\n            averages.append(next_centroid)\n            cents.append(centroid)\n\n            ddd = np.sqrt(np.dot(centroid, centroid))\n            self.centroid_distance.append(ddd)\n            # print(\"pt, cnt\", cnt, next_pt, next_d)\n            # for plane in self.planes:\n            #     ball = plane.nearest_gridpoint(next_pt)\n            #     print(\"dist\", np.linalg.norm(ball - next_pt), plane.options.radius)\n            # print(\"direction\", next_d)\n            print(\"centroid\", cnt, ddd)\n\n            # diff = np.linalg.norm(next_centroid - average_centroid)\n            # if not self.different_marginal(lines, cnt) and cnt > 100:\n                # print(\"iterations:\", cnt)\n                # break\n\n            '''\n            if cnt == iterations or (iterations > 1e10 and diff < 0.0001):\n                print(\"iterations:\", it)\n                break\n            '''\n\n            average_centroid = next_centroid\n            direction = next_d\n            pt = next_pt\n            it += 1\n\n        if get_all:\n            return averages, cents\n        else:\n            n = 10\n            m = np.zeros(len(averages[0]))\n            for a in cents[-n:]:\n                m = m + a\n            m = m / len(cents[-n:])\n            return m, lines #self.cur_lines__\n\n    def is_redundant(self, vec, b):\n        A = np.copy(self.A)\n        bb = np.copy(self.b)\n\n        A = np.vstack([A, vec])\n        bb = np.append(bb, b + 1)\n\n        c = -1.0 * vec\n\n        result = linprog(c, A_ub=A, b_ub=bb, bounds=(None,None))\n        if result.success:\n            if -result.fun <= b:\n                return True\n            else:\n                return False\n\n    def remove_redundant(self):\n        i = 0\n        while i < self.num_planes:\n            tempA = np.delete(self.A, (i), axis = 0)\n            tempB = np.delete(self.b, (i))\n\n            poly = Polytope(tempA, tempB, self.point)\n            poly.options = self.options\n\n            if poly.is_redundant(self.A[i], self.b[i]):\n                self = poly\n            i += 1\n        return self\n\n    def different_marginal(self, lines, it, eps = 1e-1):\n        if it % 50 != 0: return True\n        # self.randomized_distribution(it, self.cur_lines__[-1][1])\n        # m = self.cur_lines__\n        # m = self.marginal_dist(self.cur_lines__)\n        m = self.marginal_dist(lines)\n        if len(self.all_marginals) == 0:\n            self.all_marginals.append(m)\n            return True\n\n        # np.set_printoptions(threshold=np.nan)\n        # print(np.array2string(m, separator=', '))\n\n        #print(\"distribution:\", m)\n\n        l = len(self.all_marginals)\n        print(self.all_marginals[l//2])\n        print(m)\n        print()\n        print()\n        diff = 0\n        for i in range(len(m)):\n            a = m[i]\n            b = self.all_marginals[l//2][i]\n            diff = max(diff, abs(a - b))\n\n        print(\"cnt\", it)\n        print(\"marginal change\", diff)\n        self.all_marginal_diffs.append(diff)\n        self.all_marginals.append(m)\n        if (diff < eps) : return False\n        return True\n\n    def marginal_dist(self, lines):\n        num_steps = 10\n\n        marginal_set = []\n        for p in range(1,num_steps+1):\n            amnt = p / num_steps\n\n            l = 0\n            r = 1e2\n            start = - self.marginal_direction * 1e1\n            while (abs(l - r) > 1e-4):\n                mid = (l + r) / 2\n                pos = start + self.marginal_direction * mid\n                cut = Plane(self.marginal_direction, -np.dot(self.marginal_direction, pos))\n                m = marginal(lines, cut)\n\n                #print(m)\n                if m < amnt:\n                    r = mid\n                else:\n                    l = mid\n            #print(m, amnt, l)\n            finish = (l * self.marginal_direction + start)\n            dist = np.linalg.norm(finish)\n            if np.linalg.norm(start) < np.linalg.norm(start - finish):\n                dist *= -1\n            marginal_set.append(dist)\n\n        # step = abs(self.marginal_ends[0] + self.marginal_ends[1]) / num_steps\n        #\n        # pos = self.point - self.marginal_direction * self.marginal_ends[1]\n        # for i in range(num_steps):\n        #     cut = Plane(self.marginal_direction, -np.dot(self.marginal_direction, pos))\n        #     m = marginal(lines, cut)\n        #     marginal_set.append(m)\n        #\n        #     pos += self.marginal_direction * step\n        return np.array(marginal_set)\n\n\n    def randomized_distribution(self, iterations, start_pos = [None, None]):\n        print(iterations)\n        dimension = self.dimension\n\n        standard_directions = []\n        for i in range(dimension):\n            cur = np.zeros(dimension)\n            cur[i] = 1\n            standard_directions.append(cur)\n            other = np.zeros(dimension)\n            other[i] = -1\n            standard_directions.append(other)\n\n        if start_pos[0] != None:\n            p = self.point\n        else:\n            p = start_pos\n\n        lines = []\n        for i in range(iterations):\n            dir_indx = np.random.randint(len(standard_directions))\n            direction = standard_directions[dir_indx]\n\n            hit, _ = self.get_plane_hit(p, direction)\n\n            # choose random point from p to hit\n            dist = np.linalg.norm(hit - p)\n            length = np.random.uniform(0, dist)\n\n            next_pt = p + direction * length\n            lines.append([p, next_pt])\n            p = next_pt\n\n        for l in lines:\n            self.cur_lines__.append(l)\n        # print(self.cur_lines__)\nprint(\"GOOD\")",
      "execution_count": null,
      "outputs": []
    },
    {
      "metadata": {
        "trusted": true
      },
      "cell_type": "code",
      "source": "",
      "execution_count": null,
      "outputs": []
    },
    {
      "metadata": {
        "trusted": true
      },
      "cell_type": "code",
      "source": "",
      "execution_count": null,
      "outputs": []
    },
    {
      "metadata": {
        "trusted": true
      },
      "cell_type": "code",
      "source": "",
      "execution_count": null,
      "outputs": []
    },
    {
      "metadata": {
        "trusted": true
      },
      "cell_type": "code",
      "source": "",
      "execution_count": null,
      "outputs": []
    },
    {
      "metadata": {
        "trusted": true
      },
      "cell_type": "code",
      "source": "",
      "execution_count": null,
      "outputs": []
    },
    {
      "metadata": {
        "trusted": true
      },
      "cell_type": "code",
      "source": "",
      "execution_count": null,
      "outputs": []
    }
  ],
  "metadata": {
    "kernelspec": {
      "name": "python3",
      "display_name": "Python 3",
      "language": "python"
    },
    "language_info": {
      "mimetype": "text/x-python",
      "nbconvert_exporter": "python",
      "name": "python",
      "pygments_lexer": "ipython3",
      "version": "3.5.4",
      "file_extension": ".py",
      "codemirror_mode": {
        "version": 3,
        "name": "ipython"
      }
    }
  },
  "nbformat": 4,
  "nbformat_minor": 2
}