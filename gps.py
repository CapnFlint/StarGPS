import numpy as np
import math

def distance(a, b):
    return np.linalg.norm(a - b)

class Marker:
    def __init__(self, name="unknown", x=0, y=0, z=0, r=0):
        self.name = name
        self.x = x
        self.y = y
        self.z = z
        self.range = r
        self.vec = np.array([x,y,z])

    def __str__(self):
        return "\t{0}: {1},{2},{3}".format(self.name, self.x,self.y,self.z)

###[ Config Values ]################################

bounds = 929 #distance between adjacent OMs

# Marker X/Y locations, relative to area size

OM1 = [0.5,0.5]
OM2 = [0.5,0.5]
OM3 = [0,0.5]
OM4 = [1,0.5]
OM5 = [0.5,0]
OM6 = [0.5,1]

ranges = [387.0,674.5,652.5,336.0]

####################################################

def init_markers():
    om = [
        Marker("OM3", OM3[0] * bounds, OM3[1] * bounds, 0, 387.0),
        Marker("OM4", OM4[0] * bounds, OM4[1] * bounds, 0, 674.5),
        Marker("OM5", OM5[0] * bounds, OM5[1] * bounds, 0, 652.5)
        ]

    print "Locations:"
    for m in om:
        print m
    return om

def compute():
    om = init_markers()
    dist = [
        distance(om[1].vec, om[0].vec),
        distance(om[2].vec, om[1].vec),
        distance(om[0].vec, om[2].vec)
        ]

    print "Distances:"
    print dist



    xa = (np.square(om[0].range) + np.square(dist[0]) - np.square(om[1].range)) / (2 * dist[0])
    xb = (np.square(om[1].range) + np.square(dist[1]) - np.square(om[2].range)) / (2 * dist[1])
    xc = (np.square(om[2].range) + np.square(dist[2]) - np.square(om[0].range)) / (2 * dist[2])

print "x values:"
print xa
print xb
print xc

ya = math.sqrt(np.square(om[0].range) - np.square(xa))
yb = math.sqrt(np.square(om[1].range) - np.square(xb))
yc = math.sqrt(np.square(om[2].range) - np.square(xc))

print "y values:"
print ya
print yb
print yc

exa = (om[1].x - om[0].x) / dist[0]
eya = (om[1].y - om[0].y) / dist[0]
print om[1]
print om[0]
print dist[0]
print eya
exb = (om[2].x - om[1].x) / dist[1]
eyb = (om[2].y - om[1].y) / dist[1]

exc = (om[0].x - om[2].x) / dist[2]
eyc = (om[0].y - om[2].y) / dist[2]

nxa = -eya
nxb = -eyb
nxc = -eyc

nya = exa
nyb = exb
nyc = exc

Q1xa = om[0].x + xa * exa
Q1xb = om[1].x + xb * exb
Q1xc = om[2].x + xc * exc

print "---"
print Q1xa
print ya
print nxa
print "---"
Q1ya = om[0].y + xa * eya
Q1yb = om[1].y + xb * eyb
Q1yc = om[2].y + xb * eyc

if (ya == 0 or yb == 0 or yc == 0):
    res = [
            [Q1xa, Q1ya],
            [Q1xa, Q1ya],
            [Q1xb, Q1yb],
            [Q1xb, Q1yb],
            [Q1xc, Q1yc],
            [Q1xc, Q1yc]
            ]
else:
            
    Q2xa = Q1xa - ya * nxa
    Q2xb = Q1xb - yb * nxb
    Q2xc = Q1xc - yc * nxc
    print Q2xa
    print "---"
    Q2ya = Q1ya - ya * nya
    Q2yb = Q1yb - yb * nyb
    Q2yc = Q1yc - yc * nyc

    Q1xa += ya * nxa
    Q1xb += yb * nxb
    Q1xc += yc * nxc
    print Q1xa
    print "---"
    Q1ya += ya * nya
    Q1yb += yb * nyb
    Q1yc += yc * nyc

    print Q2ya
    print Q1ya
    print Q2xa
    print Q1xa
    ka = (Q2ya - Q1ya) / (Q2xa - Q1xa)
    da = Q2ya - (ka * Q2xa)
    yamax = ka * bounds + da
    yamin = ka * 0 + da

    kb = (Q2yb - Q1yb) / (Q2xb - Q1xb)
    db = Q2yb - (kb * Q2xb)
    ybmax = kb * bounds + db
    ybmin = kb * 0 + db

    kc = (Q2yc - Q1yc) / (Q2xc - Q1xc)
    dc = Q2yc - (kc * Q2xc)
    ycmax = kc * bounds + dc
    ycmin = kc * 0 + dc

    res = [
            [Q1xa, Q1ya],
            [Q2xa, Q2ya],
            [Q1xb, Q1yb],
            [Q2xb, Q2yb],
            [Q1xc, Q1yc],
            [Q2xc, Q2yc],
            [yamax, yamin],[ybmax,ybmin],[ycmax,ycmin]]

print res
