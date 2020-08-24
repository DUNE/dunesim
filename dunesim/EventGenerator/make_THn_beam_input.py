import ROOT as RT
from array import array 
from argparse import ArgumentParser as ap
parser = ap()

parser.add_argument( "-i", type=str, help='Input file', default="")
parser.add_argument( "-o", type=str, help='Output file', default="beam_pdf.root")
parser.add_argument( "-n", type=int, help='Max number of entries', default=-1)
args = parser.parse_args()


if not args.i:
  exit()

inputFile = RT.TFile(args.i, "OPEN")
tree = inputFile.Get("beamreco/tree")

''' fibers_h_upstream = NULL
 fibers_v_upstream = NULL
  fibers_h_downstream = NULL
   fibers_v_downstream = NULL
'''

outputFile = RT.TFile(args.o, "RECREATE")
                   #p,  fvu, fhu, fvd, fhd
nBins = array("i", [100, 32, 32, 32, 32])
mins = array("d", [0., 0., 0., 0., 0.]) 
maxes = array("d", [2., 192., 192., 192., 192.]) 

Pions = RT.THnSparseD("Pions", "", 5, nBins, mins, maxes)
Protons = RT.THnSparseD("Protons", "", 5, nBins, mins, maxes)
Electrons = RT.THnSparseD("Electrons", "", 5, nBins, mins, maxes)
#kaon?

counter = 0
nProtons = 0
nPions = 0
nElectrons = 0
#kaon?

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
  
  if Momentum < 0. or Momentum > 2.: continue 
  
  data = array("d", [Momentum, f_v_up[0], f_h_up[0], f_v_down[0], f_h_down[0]])

  if pdgs[0] == 2212:
    Protons.Fill(data)
    nProtons += 1
  elif pdgs[0] == 13:
    Pions.Fill(data)   
    nPions += 1
  elif pdgs[0] == 11:
    Electrons.Fill(data)
    nElectrons += 1
  #else kaon?

  counter += 1 

outputFile.cd()

Pions.Scale(1. / nPions)
Pions.Write()

Protons.Scale(1. / nProtons)
Protons.Write()

if nElectrons > 0: Electrons.Scale(1. / nElectrons)
Electrons.Write()

##Saving the number of particles
nElec = RT.TVectorD(1)
nElec[0] += nElectrons
nElec.Write("nElectrons")

nProt = RT.TVectorD(1)
nProt[0] += nProtons
nProt.Write("nProtons")

nPi = RT.TVectorD(1)
nPi[0] += nPions
nPi.Write("nPions")
###############################

outputFile.Close()
