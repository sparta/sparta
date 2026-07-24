import math

# 10-degree half-angle sphere-cone from issue #598
Rn    = 0.01905          # nose radius, 1.905 cm
theta = math.radians(10) # half angle
L     = 0.03             # total axial length, 3 cm (nose tip -> base plane)

st, ct, tt = math.sin(theta), math.cos(theta), math.tan(theta)

# sphere center on axis at (Rn, 0); nose tip at (0,0)
# tangency point where sphere meets cone
x_t = Rn*(1.0 - st)
r_t = Rn*ct
r_base = r_t + (L - x_t)*tt

pts = []

# 1) sphere arc: beta from 180deg (tip) down to (90+theta) (tangent)
n_arc = 60
b0, b1 = math.pi, math.pi/2 + theta
for i in range(n_arc+1):
    b = b0 + (b1-b0)*i/n_arc
    x = Rn + Rn*math.cos(b)
    r = Rn*math.sin(b)
    pts.append((x, r))

# 2) cone: tangent point -> base corner (skip first, dup of arc end)
n_cone = 30
for i in range(1, n_cone+1):
    f = i/n_cone
    pts.append((x_t + (L-x_t)*f, r_t + (r_base-r_t)*f))

# 3) base face: base corner -> axis (skip first, dup)
n_base = 30
for i in range(1, n_base+1):
    f = i/n_base
    pts.append((L, r_base*(1.0-f)))

# snap axis endpoints exactly to r=0
pts[0]  = (0.0, 0.0)
pts[-1] = (L, 0.0)

# de-dup consecutive identical
clean = [pts[0]]
for p in pts[1:]:
    if abs(p[0]-clean[-1][0])>1e-15 or abs(p[1]-clean[-1][1])>1e-15:
        clean.append(p)
pts = clean

n = len(pts)
with open("data.spherecone", "w") as f:
    f.write("surf file: 10-deg half-angle sphere-cone (issue #598), axisymmetric\n\n")
    f.write("%d points\n" % n)
    f.write("%d lines\n\n" % (n-1))
    f.write("Points\n\n")
    for i,(x,r) in enumerate(pts):
        f.write("%d %.10g %.10g\n" % (i+1, x, r))
    f.write("\nLines\n\n")
    for i in range(n-1):
        f.write("%d %d %d\n" % (i+1, i+1, i+2))

print("points=%d lines=%d" % (n, n-1))
print("tangent (x_t,r_t) = (%.5g, %.5g)" % (x_t, r_t))
print("base radius r_base = %.5g" % r_base)
print("axis endpoints:", pts[0], pts[-1])
