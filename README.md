
# Volume Approx. Implementation Details

### Hyperplanes

The convex hull will be represented as the intersection of half-spaces.
So, we have a handy class, `Plane`, which will give us some basic utility
functions for dealing with hyperplanes such as

1. Determine if a given point is on the hyperplane
2. Determine if a point is in the half-space
3. Given a line (defined as a point along with a direction vector),
   determine the intersection of the line with the plane.
4. Reflect a point given the direction it travels in.

Implementation details can be found in the file `Plane.py`


### The Polytope

We describe the polytope as

$$Ax \leq b$$

The polytope class provides the following functions we can make use of.

1. Given an initial point and a direction, bounce this point around the
   polytope, reflecting it wen it hits a wall.

   a. First, we must determine which hyperplane the point will intersect first.
      We determine this by finding the intersection with each,
      and we take the one which lies within the polytope.

   b. Currently, the initial direction is arbitrary.
2. We approximate the centroid by simply taking the average of all the
   midpoints of each line segment we created through the reflection process.

##### Removing Redundant Inequalities

We can remove linear inequalities if the resulting polytope represents
the same convex set as before. We do so by iterating through each inequality
and determining if its remove has any effect as follows:

Say we wish to check the $i^{th}$ inequality in $A$. let $A'$ be $A$ with
the $i^{th}$ row removed, and $b'$ be $b$ with the $i^{th}$ element removed.
Then, we solve the following linear program:

\newpage

maixmize

$$A_i x$$

subject to

$$ A'x \leq b' \\  A_i x \leq b_i + 1 $$

If the maximum objective value if at most $b_i$, then $A'x \leq b'$ inherently
implies that $A_i x \leq b_i$ and so, we can safely remove the $i^{th}$
inequality.

Implementation-wise, we simply use a built-in linear program solver from the
scipy library.

### Computing Volume

We first perform the reflection procedure to approximate the centroid and we
draw an axis-aligned hyperplane through this point (axis is chosen round-robin).
This divides the polytope in two. The ratio of the volumes of these is
approximated by the ratio of the lengths of the lines which lie in each half.
We recursively perform the above steps on one of the resulting polytopes.

##### When to Stop

We can reliably compute the exact volume of any n-dimensional cube, given the
side lengths. Eventually, after some number of iterations, the resulting
polytope will be defined exclusively by hyperplanes that we inserted
during the course of the algorithm. By construction, these define a cube. Having
this volume, and the ratios, we may approximate the volume of the original
polytope.

However, to do the above, we must know that our polytope indeed consists only
of the axis-aligned hyperplanes. While simply adding our half-space at each
iteration is sufficient to describe the divided polytope, we must have a
minimal description of it to be able to decide that it is a cube. So here is
where we use the previously defined method to remove redundant inequalities. (We
might think of a more efficient way to do this in the future).

Additionally, note that we are required to know of a point inside the polytope
in order to be able to perform the reflections. Initially, we'll assume that
we are given some point (eg. the origin). At successive iterations, we will
use the approximated centroid we computed in the step prior.

## Testing

##### Cubes

We can compute the exact volume of a cube ourselves, so we have an easy way
to determine the accuracy of the algorithm. We force the algorithm to run for
at least 2 * (# of dimensions) times.

The accuracy depends greatly on the number of reflections we perform at each
step, and needs to increase with the number of dimensions to keep the same level
of accuracy as would be needed for lower dimensions.

You can run the `test(d, i)` function in `cube_test.py` to run
the algorithm on a cube of side length 2 in `d` dimensions, using `i` iterations
of reflection to approximate the centroid.

##### Visualizing Reflections

We can visualize the reflections for a square by running the file
`draw_square.py`.

##### Some Arbitrary 2D Shapes

We can test this on some polygons - we can either compute the areas exactly using pen and paper, or we can get the computer to approximate the area by essentially creating a lattice with squares of some small side length and add up the number of squares within the polygon. See `test2d.py`

## Sphere Reflections

Essentially, it does not work well. We'll first discuss the implementation, and then the pitfalls that arose, and compare it to the normal reflection method.

Each plane is essentially covered with spheres. We set up a grid on the plane with each cube of side length $\delta$.

![image](/home/rares/Research/pics/grid.png =300x300)

Here, the 'plane' is the red line, with $\delta = \sqrt2$. The black is the line we are reflecting. We do not create all of the spheres (as this may be an exponential number), but rather create them as needed. We create spheres as follows:

  1. Find grid point $c$ closest to the intersection between plane and line.
  2. Set center of circle along perpendicular to the plane at $c$, offset by $r$
  3. Set radius $r \cdot \alpha$ for some chosen constant $\alpha$.

![image](/home/rares/Research/pics/circle.png =300x270)

The green line is the perpendicular through $c$, and the distance between the center and $c$ is $r$. Having the sphere (circle in this case), it only remains to reflect the line about this sphere. This is simply reflecting across the tangent at the point of intersection:

![image](/home/rares/Research/pics/sphere_reflection.png =300x200)

There is an issue here, however. If we had indeed created every circle at each grid point, the line may have intersected a different circle before the intersection we have found with the method above. In fact, we may even find an intersection that lies within another circle. This brings the problem: what happens on the next reflection? There is no sphere to intersect, since we are already inside of a sphere.

![image](/home/rares/Research/pics/problem.png =300x200)

In the above example, the plane the circle lie on is simply the x-axis. We would have reflected against the purple circle. However, the actual first intersection occurs with the red circle. This is not entirely problematic, however, this becomes an issue once the next intersection is with a sphere which contains the previous point of intersection. Moreover, it is impossible to cover the plane with disjoint circles.

This approach does work to some extent, however, the issue mentioned above makes consistently produces errors when computing volumes. Nonetheless, we can compare the results of both reflection methods.

Here is the result of reflecting lines 100 times within a square using facet reflections:

![image](/home/rares/Research/pics/normal_reflect.png =300x300)

and here it is using reflections on spheres:

![image](/home/rares/Research/pics/reflection_with_sphere.png =300x300)

These would clearly change depending on the initial direction chosen, however, the facet reflection always has a great amount of symmetry. With spherical reflections, we have a more hectic arrangement, however, interestingly, I have found that there is often a concentration around the corners of the square. That is, lines go from corner to corner.
