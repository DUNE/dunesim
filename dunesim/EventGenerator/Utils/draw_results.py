import ROOT as RT
from argparse import ArgumentParser as ap
RT.gStyle.SetOptStat(0)
RT.gROOT.SetBatch()

def save_c(c):
  c.SaveAs(c.GetName() + '.png')

def draw_trigger_pdgs(fOut, t, args):
  h = RT.TH1D('hPDG', 'Triggered Particle Yields;;Entries', 5, 0, 5)
  parts = [-11, -13, 211, 2212, -999]
  #parts = [-11, -13, 211, 321, 2212, -999]
  #labels = ['e^{+}', '#mu^{+}', '#pi^{+}', 'K^{+}', 'p', 'Other']
  labels = ['e^{+}', '#mu^{+}', '#pi^{+}', 'p', 'Other']
  for i in range(1, len(parts)+1):
    if parts[i-1] == -999:
      #cut = ('trigger_pdg != -11 && trigger_pdg != -13 && '
      #       'trigger_pdg != 211 && trigger_pdg != 321 && '
      #       'trigger_pdg != 2212')
      cut = ('trigger_pdg != -11 && trigger_pdg != -13 && '
             'trigger_pdg != 211 && '
             'trigger_pdg != 2212')
      h.SetBinContent(i, t.GetEntries(cut))
    else:
      h.SetBinContent(i, t.GetEntries(f'trigger_pdg == {parts[i-1]}'))
    h.GetXaxis().SetBinLabel(i, labels[i-1])

  fOut.cd()

  c = RT.TCanvas('cPDG', '')
  h.SetMinimum(0)
  h.SetLineWidth(2)
  c.SetTicks()
  h.Draw()
  c.Write()
  if args.save: save_c(c)
  h.Write()

  c = RT.TCanvas('cPDG_rel', '')
  h2 = h.Clone('hPDG_rel')
  h2.SetMinimum(0)
  h2.SetMaximum(1.1)
  h2.SetLineWidth(2)

  for i in range(len(parts)):
    h2.SetBinContent(i+1, h2.GetBinContent(i+1)/t.GetEntries())
  c.SetTicks()
  h2.Draw()
  c.Write()
  if args.save: save_c(c)

  if args.print:

    for i in range(h.GetNbinsX()):
      print(h.GetXaxis().GetBinLabel(i+1), ' ')
    for i in range(h.GetNbinsX()):
      print(h.GetBinContent(i+1), ' ')
    print(h.Integral(), t.GetEntries())

    for i in range(h2.GetNbinsX()):
      print(h2.GetXaxis().GetBinLabel(i+1), ' ')
    for i in range(h2.GetNbinsX()):
      print(h2.GetBinContent(i+1), ' ')
    print()

  h2.Write()

  return h

def draw_background_yields(fOut, t, args):
  h = RT.TH1D('hBG_PDG', 'Triggered Particle Yields;;Entries', 7, 0, 7)
  parts = [11, 13, 211, 321, 2212, 22, -999]
  labels = ['e', '#mu', '#pi', 'K', 'p', '#gamma', 'Other']
  for i in range(1, len(parts)+1):
    if parts[i-1] == -999:
      cut = ('abs(background_pdgs) != 11 && abs(background_pdgs) != 13 && '
             'abs(background_pdgs) != 211 && abs(background_pdgs) != 321 && '
             'abs(background_pdgs) != 2212 && abs(background_pdgs) != 22')
      h.SetBinContent(i, t.GetEntries(cut))
    else:
      h.SetBinContent(i, t.GetEntries(f'abs(background_pdgs) == {parts[i-1]}'))
    h.GetXaxis().SetBinLabel(i, labels[i-1])

  fOut.cd()
  c = RT.TCanvas('cBG_PDG', '')
  h.SetMinimum(0)
  h.SetLineWidth(2)
  c.SetTicks()
  h.Draw()
  c.Write()
  if args.save: save_c(c)
  h.Write()

  return h


def draw_trigger_profiles(fOut, t, args):
  parts = [-11, -13, 211, 321, 2212]
  names = {
    -11:'positron',
    -13:'muplus',
    211:'piplus',
    321:'kplus',
    2212:'proton',
  }
  labels = ['e^{+}', '#mu^{+}', '#pi^{+}', 'K^{+}', 'p']
  binnings = {
    -11: '(50, -25, 10, 50, 430, 470)',
    -13: '(50, -25, 15, 50, 430, 470)',
    211: '(50, -20, 10, 50, 430, 470)',
    321: '(50, -15, 5, 50, 430, 470)',
   2212: '(50, -20, 10, 50, 430, 470)',
  }
  for i in range(len(parts)):
    #h = RT.TH1D(f'hXY_{names[parts[i]}')
    #binning = '(50, -15, 5, 50, 430, 470)'
    binning = binnings[parts[i]]
    t.Draw(f'trigger_front_y:trigger_front_x>>hXY_{names[parts[i]]}{binning}',
           f'trigger_pdg == {parts[i]}')
    h = RT.gDirectory.Get(f'hXY_{names[parts[i]]}')
    h.SetTitle(f'Front Profile -- Triggered {labels[i]};Front X (cm);Front Y (cm)')

    fOut.cd()
    c = RT.TCanvas(f'cXY_{names[parts[i]]}', '')
    c.SetLogz()
    c.SetTicks()
    h.Draw('colz')
    c.Write()
    if args.save: save_c(c)
    h.Write()

  return h

'''background_front_p", "background_pdgs != 22 && abs(background_pdgs) != 2112 background_front_y < 600 && abs(background_front_x) < 350'''

class Values:
  def __init__(self, e):
    self.bg_p = [i for i in e.background_front_p]
    self.bg_pdgs = [i for i in e.background_pdgs]
    self.bg_y = [i for i in e.background_front_y]
    self.bg_x = [i for i in e.background_front_x]

    self.trigger_pdg = e.trigger_pdg

    self.get_valid_bg()
    self.get_plug_bg()
  
  def get_plug_bg(self):
    self.plug_bg = []
    for i in self.valid_bg:
      if (self.bg_y[i] < 465. and self.bg_y[i] > 435. and
          self.bg_x[i] > -20. and self.bg_x[i] < 10.):
        self.plug_bg.append(i)
      
  def get_valid_bg(self, limit_neutron=False):
    self.valid_bg = []
    for i in range(len(self.bg_pdgs)):
      if self.bg_pdgs[i] == 22: continue
      if limit_neutron and self.bg_pdgs[i] == 2112: continue
      if self.bg_y[i] > 600 or self.bg_y[i] < 0.: continue
      if abs(self.bg_x[i]) > 350: continue

      self.valid_bg.append(i)
      
def tree(t, args):
  hN_BG = RT.TH1D('hN_BG', 'All Good Background Particles;N Background Particles', 10,0,10)
  hN_BG_plug = RT.TH1D('hN_BG_plug', 'Background Particles Near Beam Plug;N Background Particles', 10,0,10)
  hP_BG = RT.TH1D('hP_BG', 'All Good Background Particles;BG Particle Momentum [GeV/c]', 100,0,1.5)
  hP_BG_plug = RT.TH1D('hP_BG_plug', 'Background Particles Near Beam Plug;BG Particle Momentum [GeV/c]', 100,0,1.5)
  for e in t:
    #if e.interacted: continue

    vals = Values(e) 
    hN_BG.Fill(len(vals.valid_bg))
    hN_BG_plug.Fill(len(vals.plug_bg))
    #print(vals.valid_bg) 
    for i in vals.valid_bg:
      hP_BG.Fill(vals.bg_p[i])
      #print(vals.bg_pdgs[i])

    for i in vals.plug_bg:
      hP_BG_plug.Fill(vals.bg_p[i])

  hN_BG.Write()
  hN_BG_plug.Write()
  hP_BG.Write()
  hP_BG_plug.Write()

  if args.save:
    cN_BG = RT.TCanvas('cN_BG')
    hN_BG.Draw()
    save_c(cN_BG)

    cN_BG_plug = RT.TCanvas('cN_BG_plug')
    hN_BG_plug.Draw()
    save_c(cN_BG_plug)

    cP_BG = RT.TCanvas('cP_BG')
    hP_BG.Draw()
    save_c(cP_BG)

    cP_BG_plug = RT.TCanvas('cP_BG_plug')
    hP_BG_plug.Draw()
    save_c(cP_BG_plug)

def process(args):
  fOut = RT.TFile(args.o, 'recreate')

  fIn = RT.TFile.Open(args.i)
  t = fIn.Get('tree')

  fOut.cd()
  tree(t, args)
  hPDG = draw_trigger_pdgs(fOut, t, args)
  hBG_PDG = draw_background_yields(fOut, t, args)

  hProfile = draw_trigger_profiles(fOut, t, args)

  fOut.Close()
  fIn.Close()

if __name__ == '__main__':
  parser = ap()
  parser.add_argument('-i', required=True, help='Input file')
  parser.add_argument('-o', default='yield_results.root', help='Output file')
  parser.add_argument('--save', action='store_true')
  parser.add_argument('--print', action='store_true')
  args = parser.parse_args()

  process(args)


