import ROOT as RT
from argparse import ArgumentParser as ap
from array import array

def line_to_coords(line):
  return [float(i.strip()) for i in line.split()]

parser = ap()

parser.add_argument('-o', type=str, required=True)
parser.add_argument('-i', type=str, required=True)

args = parser.parse_args()

with open(args.i, 'r') as f:
  coords = [line_to_coords(l) for l in f.readlines() if '#' not in l]
print(coords)

tree = RT.TTree('tree', '')
x = array('d', [0])
y = array('d', [0])
z = array('d', [0])
theta = array('d', [0])
phi = array('d', [0])

tree.Branch('x', x, 'x/D')
tree.Branch('y', y, 'y/D')
tree.Branch('z', z, 'z/D')
tree.Branch('theta', theta, 'theta/D')
tree.Branch('phi', phi, 'phi/D')

for c in coords:
  (x[0], y[0], z[0], theta[0], phi[0]) = c
  tree.Fill() 
fout = RT.TFile(args.o, 'recreate')
tree.Write()
fout.Close()
