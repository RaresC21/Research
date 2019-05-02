# Useful Aspects of Using a Notebook 

In short, this is nice because it's interactive. For example, we were looking at the distance of the average from the origin (the centroid in the case we were considering) and found some unexpected peaks. We wanted to see what was happening during one of those peaks. To do so, we would need to perform the reflections, plot it to see where a peak was, and then save the points in some separated file, then read the same points again from another script to analyze what has happening during that specific interval. 

This notebook is very helfpul since it remembers everything for us - the polytope, the reflection points etc. Any variable that we may want saved. Then, in a separated cell, we may use these variables after the fact. 

Also, this is all online, so I don't have to send you a lot of zipped files each time I change something. 

### Warning 

If you want to run 10^5 iterations or more, I would recommend running this locally. This runs everything on a virtual machine somwhere which isn't very performant, and would take significantly longer than if run locally. Also, if you try storing that many points in the notebook, everything becomes very slow, so I advise against it. 


# Automatically Updating Changes

If you make a change in a script file (like Plane.py), the changes do not get automatically imported into a running ipython kernel. By default, it will use the old version unless you restart the kernel. To get around this, include the following two lines at the top of any file: 

` %reload_ext autoreload `
` %autoreload 2 `

# Inline Plots 

I've found using 

` %matplotlib inline `

works well to show plots (this may be a default so you may not need to write this) 

Additionally, to enlarge the plots, since by default they are fairly small, use the following: 

` plt.figure(figsize=(width, height)) `

