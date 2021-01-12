#!/usr/bin/python

import argparse
parser = argparse.ArgumentParser(\
  description = 'Submits MD simulations to Bridges.' + \
    'Loads a .PDB file from from ../1-model-water/[system_name].pdb.' + \
    'Stores simulation results in [repetition]/.')
parser.add_argument('--type', \
  choices=['shared','small','test'], \
  default='test', \
  help="The type of job." + \
    "'shared' is 48 hrs on GPU-shared, " + \
    "'small' is 8 hrs on GPU-small, " + \
    "and 'test' is 0.5 hrs on GPU-small.")
parser.add_argument('--system_name', \
  default='1ubq', \
  help='The name of the system.')
parser.add_argument('--repetition', \
  type=int, \
  default=0, \
  help='The repetition of the simulation.')
parser.add_argument('--simulation_time',
  type=int, \
  default=100, \
  help='The amount of time to simulate, in nanoseconds.')
args = parser.parse_args()

if args.type=='shared':
  partition = 'GPU-shared'
  wall_time = '48:00:00'
elif args.type=='small':
  partition = 'GPU-small'
  wall_time = '08:00:00'
elif args.type=='test':
  partition = 'GPU-small'
  wall_time = '00:30:00'

import os.path
if not os.path.isdir(repr(args.repetition)):
  os.mkdir(repr(args.repetition))
script_path = os.path.join(os.getcwd(), repr(args.repetition))
os.chdir(script_path)

slurm_script = '''#!/bin/bash
#SBATCH -N 1
#SBATCH -p {0}
#SBATCH --ntasks-per-node 1
#SBATCH --gres=gpu:p100:1
#SBATCH -t {1}

#echo commands to stdout
set -x

# Set up the environment
source ~/.bashrc
module load anaconda3
conda activate openmm

# Run MD
cd {2}
python ../run_simulation.py --system_name {3} --simulation_time {4}
'''.format(partition, wall_time, script_path, \
  args.system_name, args.simulation_time)

print(slurm_script)
F = open('MD_simulation.job','w')
F.write(slurm_script)
F.close()

os.system('sbatch MD_simulation.job')
