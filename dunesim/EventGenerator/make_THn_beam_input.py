import ROOT as RT
from array import array 
from argparse import ArgumentParser as ap
parser = ap()

parser.add_argument("-i", type=str, help='Input file', default="")
parser.add_argument("-o", type=str, help='Output file', default="beam_pdf.root")
parser.add_argument("-n", type=int, help='Max number of entries', default=-1)
parser.add_argument("-p", type=int, help='Which momentum', default=1)
args = parser.parse_args()


if not args.i:
  exit()

inputFile = RT.TFile(args.i, "OPEN")
tree = inputFile.Get("beamreco/tree")
outputFile = RT.TFile(args.o, "RECREATE")
                   #p,  fvu, fhu, fvd, fhd
nBins = array("i", [100, 32, 32, 32, 32])

if args.p == 1:
  min_p = 0.
  max_p = 2.
elif args.p == 2:
  min_p = 1.
  max_p = 3.
elif args.p == 3:
  min_p = 2.
  max_p = 4.
elif args.p == 6:
  min_p = 4.
  max_p = 8.
elif args.p == 7:
  min_p = 6.
  max_p = 8.

mins = array("d", [min_p, 0., 0., 0., 0.]) 
maxes = array("d", [max_p, 192., 192., 192., 192.]) 

Pions = RT.THnSparseD("Pions", "", 5, nBins, mins, maxes)
Protons = RT.THnSparseD("Protons", "", 5, nBins, mins, maxes)
Electrons = RT.THnSparseD("Electrons", "", 5, nBins, mins, maxes)
Kaons = RT.THnSparseD("Kaons", "", 5, nBins, mins, maxes)

if args.p in [1, 2, 3]:
  pdfs = [Pions, Protons, Electrons]
elif args.p in [6, 7]:
  pdfs = [Pions, Protons, Kaons]

counter = 0

if args.n < 0:
  max_entries = tree.GetEntries()
else:
  max_entries = args.n

for e in tree:
  if counter >= max_entries: break
  f_v_up = [i for i in e.fibers_v_upstream]
  f_h_up = [i for i in e.fibers_h_upstream]
  f_v_down = [i for i in e.fibers_v_downstream]
  f_h_down = [i for i in e.fibers_h_downstream]
  perfectP = e.perfectP
  Momentum = e.Momentum
  pdgs = [i for i in e.possible_pdg]
 

  if not (perfectP and len(f_v_up) == 1 and len(f_h_up) == 1 and
          len(f_v_down) == 1 and len(f_h_down) == 1 and len(pdgs) > 0):
    continue
  
  if Momentum < min_p or Momentum > max_p: continue 
  
  data = array("d", [Momentum, f_v_up[0], f_h_up[0], f_v_down[0], f_h_down[0]])

  if args.p in [1, 2]:
    if pdgs[0] == 2212:
      Protons.Fill(data)
    elif pdgs[0] == 13:
      Pions.Fill(data)   
    elif pdgs[0] == 11:
      Electrons.Fill(data)

  elif args.p == 3:
    if 13 in pdgs:
      Pions.Fill(data)
    elif 2212 in pdgs:
      Protons.Fill(data)
    elif 11 in pdgs:
      Electrons.Fill(data)

  elif args.p in [6, 7]:
    if 13 in pdgs:
      Pions.Fill(data)
    elif 2212 in pdgs:
      Protons.Fill(data)
    elif 321 in pdgs:
      Kaons.Fill(data)

  counter += 1 

outputFile.cd()

for pdf in pdfs:
  if pdf.Projection(0).Integral() > 0:
    pdf.Scale(1./pdf.Projection(0).Integral())
  pdf.Write()

outputFile.Close()
