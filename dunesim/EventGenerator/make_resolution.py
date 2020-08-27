import ROOT as RT
from array import array 
from argparse import ArgumentParser as ap
parser = ap()

parser.add_argument("-i", type=str, help='Input file', default="")
parser.add_argument("-o", type=str, help='Output file', default="beam_res.root")
parser.add_argument("-b", type=str, help='Binning', default="75, -.15, .15")
parser.add_argument("--b2", type=str, help='Binning', default="75, -.15, .15")
parser.add_argument("--eb", type=str, help='Binning for electrons',
                    default="75, -.15, .5")
parser.add_argument("--eb2", type=str, help='Binning for electrons',
                    default="75, -.15, .5")
args = parser.parse_args()

if not args.i:
  exit()

RT.gROOT.SetBatch(1)

inputFile = RT.TFile(args.i, "OPEN")
tree = inputFile.Get("tree")

outputFile = RT.TFile(args.o, "RECREATE")

tree.Draw("(reco_p*1.e3 - NP04front_p)/NP04front_p>>hPionsRes(" + args.b + ")", "PDG == 211 || PDG == -13")
tree.Draw("(reco_p*1.e3 - NP04front_p)/NP04front_p>>hProtonsRes(" + args.b + ")", "PDG == 2212")
tree.Draw("(reco_p*1.e3 - NP04front_p)/NP04front_p>>hElectronsRes(" + args.eb + ")", "PDG == -11")
tree.Draw("(reco_p*1.e3 - NP04front_p)/NP04front_p>>hKaonsRes(" + args.b + ")", "PDG == 321")

tree.Draw("NP04front_p*1.e-3:reco_p>>hPionsRes2D(" + args.b2 + ")", "PDG == 211 || PDG == -13")
tree.Draw("NP04front_p*1.e-3:reco_p>>hProtonsRes2D(" + args.b2 + ")", "PDG == 2212")
tree.Draw("NP04front_p*1.e-3:reco_p>>hElectronsRes2D(" + args.eb2 + ")", "PDG == -11")
tree.Draw("NP04front_p*1.e-3:reco_p>>hKaonsRes2D(" + args.b2 + ")", "PDG == 321")

hPionsRes = RT.gDirectory.Get("hPionsRes")
hProtonsRes = RT.gDirectory.Get("hProtonsRes")
hElectronsRes = RT.gDirectory.Get("hElectronsRes")
hKaonsRes = RT.gDirectory.Get("hKaonsRes")

hPionsRes2D = RT.gDirectory.Get("hPionsRes2D")
hProtonsRes2D = RT.gDirectory.Get("hProtonsRes2D")
hElectronsRes2D = RT.gDirectory.Get("hElectronsRes2D")
hKaonsRes2D = RT.gDirectory.Get("hKaonsRes2D")

outputFile.cd()
hPionsRes.Write()
hProtonsRes.Write()
hElectronsRes.Write()
hKaonsRes.Write()
hPionsRes2D.Write()
hProtonsRes2D.Write()
hElectronsRes2D.Write()
hKaonsRes2D.Write()
outputFile.Close()
