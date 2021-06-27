import ROOT as RT
from argparse import ArgumentParser as ap
import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt

parser = ap()
parser.add_argument("-i", type=str, help='Input')
parser.add_argument("-o", type=str, help='Output string')
args = parser.parse_args()

def get_disp_points(hx, hy, hz):
  results = []
  for i in range(1, hx.GetNbinsX()+1):
    for j in range(1, hx.GetNbinsY()+1):
      for k in range(1, hx.GetNbinsZ()+1):
        results.append(
            [hx.GetBinContent(i, j, k),
             hy.GetBinContent(i, j, k),
             hz.GetBinContent(i, j, k)])
  return results

def get_bin_points(h):
  results = []
  for i in range(1, h.GetNbinsX()+1):
    x = h.GetXaxis().GetBinCenter(i)
    for j in range(1, h.GetNbinsY()+1):
      y = h.GetYaxis().GetBinCenter(j)
      for k in range(1, h.GetNbinsZ()+1):
        z = h.GetZaxis().GetBinCenter(k)
        results.append([x, y, z])
  return results
'''
fIn = RT.TFile(args.i, "OPEN")
bkwd_z_pos = fIn.Get("RecoBkwd_Displacement_Z_Pos")
bkwd_z_neg = fIn.Get("RecoBkwd_Displacement_Z_Neg")
bkwd_x_pos = fIn.Get("RecoBkwd_Displacement_X_Pos")
bkwd_x_neg = fIn.Get("RecoBkwd_Displacement_X_Neg")
bkwd_y_pos = fIn.Get("RecoBkwd_Displacement_Y_Pos")
bkwd_y_neg = fIn.Get("RecoBkwd_Displacement_Y_Neg")

bkwd_pos_points = get_bin_points(bkwd_z_pos)
bkwd_neg_points = get_bin_points(bkwd_z_neg)

bkwd_pos_displaced = get_disp_points(bkwd_x_pos, bkwd_y_pos, bkwd_z_pos)
bkwd_neg_displaced = get_disp_points(bkwd_x_neg, bkwd_y_neg, bkwd_z_neg)

fIn.Close()
'''


points = np.array([[0, 0], [0, 1.0], [1, 0], [1.0, 1.0]])
corr_points = np.array([[.15, .15], [.2, 1.3], [1.3, .3], [1.3, 1.4]])

tri = Delaunay(points)
tri_corr = Delaunay(corr_points)
plt.triplot(points[:,0], points[:,1], tri_corr.simplices)
plt.plot(points[:,0], points[:,1], 'o')
plt.grid()

#plt.triplot(corr_points[:,0], corr_points[:,1], tri_disp.simplices)
plt.triplot(corr_points[:,0], corr_points[:,1], tri_corr.simplices)
plt.plot(corr_points[:,0], corr_points[:,1], 'o')

#true_pts = []
#for i in range(0, 12):
#  for j in range(0, 12):
#    true_pts.append([i*.1, j*.1])
#true_grid = np.array(true_pts)
true_grid = np.array([[.25, .5], [.75, .5]])
plt.plot(true_grid[:,0], true_grid[:,1], 'o')

print(points)
print(tri_corr.simplices)

distorted_grid = []

for pt in true_grid:
  i = tri_corr.find_simplex(pt)
  print(i)
  s = tri_corr.simplices[i]
  print(pt, i)
  print(s)
  r = tri_corr.transform[i,2]
  b = tri_corr.transform[i,:2].dot(np.transpose(pt - r))
  b = np.append(b, 1. - b.sum())
  print(b)

  new_point_x = 0.
  new_point_y = 0.
  for j, bi in zip(s,b):
    new_point_x += points[j][0]*bi
    new_point_y += points[j][1]*bi
  print([new_point_x, new_point_y])
  distorted_grid.append([new_point_x, new_point_y])

distorted_array = np.array(distorted_grid)

deltas = [ [d[0] - t[0], d[1] - t[1]] for d, t in zip(distorted_grid, true_grid)]
#plt.plot(distorted_array[:,0], distorted_array[:,1], 'o')
delta_array = np.array(deltas)
for pt, d in zip(true_grid, deltas):
  plt.arrow(pt[0], pt[1], d[0], d[1], head_width=.01)
plt.show()
