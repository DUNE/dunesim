import subprocess
from math import ceil
import sys
import argparse
import os 
from glob import glob as ls

'''Method to add fcl file to the submit cmd'''
def add_fcl(cmd, args):
  #Check if the fcl is valid
  if '.fcl' not in args.fcl:
    print('Warning. Invalid fcl', args.fcl)
    sys.exit(1)

  #check it exists
  thefile = ls(args.fcl)
  print(thefile)
  if len(thefile) == 0:
    print('Could not find fcl')
    sys.exit(1)

  #Isolate the name
  splitted = args.fcl.split('/')
  fcl_name = splitted[-1]
  print(fcl_name)

  #Check if a directory was supplied -- if so bring it along in the job
  if len(splitted) > 1:
    fcl_dir = '/'.join(splitted[:-1])
    print('At', fcl_dir)
    cmd += [f'-Osubmit.f_1=dropbox://{args.fcl}']
    cmd += [f'-Oglobal.fcl_name=\\\\\\${{CONDOR_DIR_INPUT}}/{fcl_name}']
  else:
    cmd += [f'-Oglobal.fcl_name={fcl_name}']



def add_output_dir(cmd, args):
  if args.output_dir is None: return
  cmd += [f'-Oglobal.output_dir={args.output_dir}']

def add_tar(cmd, args):
  if args.tar is None: return
  cmd += [f'-Osubmit.tar_file_name=dropbox://{args.tar}']

def submit_loop(cmd, args):

  ##Check start file is >= 0
  if args.start_file < 0:
    print('Error need to provide start_file >= 0')
    sys.exit(1)

  ##Get the number of jobs to submit
  njobs = ceil((args.n_files)/args.n_per_job)
  print('N jobs:', njobs)
  
  start_file = args.start_file
  base_cmd = cmd
  for i in range(njobs):
    cmd = base_cmd
    cmd += [f'-Oglobal.start_file={start_file}', f'-Oglobal.nfiles={args.n_per_job}']
    start_file += args.n_per_job
    subprocess.run(cmd)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description = 'Submission script for beam studies')
  parser.add_argument('--config', type=str, help='Which config',
                      default='%s/cfg_files/h4_vle_studies.cfg'%(os.environ['DUNESIM_DIR']))
  parser.add_argument('--output_dir', type=str, help='Output top dir', default=None)
  parser.add_argument('--dry_run', action='store_true', help='Tell fife_launch to do a dry_run')
  parser.add_argument('--lifetime', type=str, default='1h')
  parser.add_argument('--fcl', type=str, required=True)
  parser.add_argument('--n_files', type=int, required=True)
  parser.add_argument('--n_per_job', type=int, default=1)
  parser.add_argument('--start_file', type=int, default=0)
  parser.add_argument('--tar', type=str, default=None)

  args = parser.parse_args()

  cmd = ['fife_launch', '-c', args.config]

  add_fcl(cmd, args)
  add_output_dir(cmd, args)
  add_tar(cmd, args)
  if args.dry_run: cmd += ['--dry_run']

  print(cmd)
  submit_loop(cmd, args)
