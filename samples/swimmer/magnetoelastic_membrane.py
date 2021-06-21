# Copyright (C) 2010-2019 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
Simulation of magnetoelastic membranes
"""

import espressomd
import espressomd.magnetostatics

print(espressomd.features())
required_features = ["EXTERNAL_FORCES", "DIPOLES", "WCA"]
espressomd.assert_features(required_features)

import numpy as np
import os
import object_in_fluid as oif

import ipdb

# system setup
boxX = 100.0
boxY = 100.0
boxZ = 100.0

periodicX = True
periodicY = True
periodicZ = True

time_step = 0.001

system = espressomd.System(box_l=(boxX, boxY, boxZ))
system.periodicity = (periodicX, periodicY, periodicZ)
system.time_step = time_step
system.cell_system.skin = 0.2

# cuda devices
#print(system.cuda_init_handle.list_devices())

# input files
input_folder = "./input/"
nodes_file = input_folder + "1082.6.6.nodes"
triangles_file = input_folder + "1082.6.6.triangles"

# output files
output_folder = "./output/"
if not os.path.isdir(output_folder):
    os.makedirs(output_folder)
else:
    print("output files already existed and will be overwritten.")

# particle types
type_A = 0
type_B = 1

# stretching constants
ks_A = 4000.0
ks_B = 4000.0

# bending constants
kb_A = 50.0
kb_B = 50.0

# define the elastic object
rescale=(5.0, 5.0, 5.0)
elastic_object_type = oif.ElasticObjectType(nodes_file=nodes_file,
        triangles_file=triangles_file, rescale=rescale, ks_A=ks_A, ks_B=ks_B, kb_A=kb_A, kb_B=kb_B)
elastic_object_type.initialize()

# create the elastic object
translate = (0.0, 0.0, 0.0)
rotate = (0.0, 0.0, 0.0)
membrane = oif.ElasticObject(system, elastic_object_type, object_id=0, particle_type_A=type_A, particle_type_B=type_B, translate=translate, rotate=rotate)
membrane.initialize()

# define magnetic interactions
mu_A = 1.0
mu_B = 1.0
for p in system.part:
    if p.type == type_A:
        p.dip = (0.0, 0.0, 1.0)
        p.dipm = mu_A
    if p.type == type_B:
        p.dip = (0.0, 0.0, 1.0)
        p.dipm = mu_B

H_field = [1,0,0]
H_constraint = espressomd.constraints.HomogeneousMagneticField(H=H_field)
system.constraints.add(H_constraint)

# magnetic dipole-dipole interactions
dipolar_direct_sum = espressomd.magnetostatics.DipolarDirectSumCpu(prefactor=1)
system.actors.add(dipolar_direct_sum)

# dipolar p3m
#p3m = espressomd.magnetostatics.DipolarP3M(prefactor=1, accuracy=1e-3, r_cut=2.0)
#system.actors.add(p3m)

# define lj interactions
wca_epsilon = 1.0
wca_sigma = 0.4

system.non_bonded_inter[0, 0].wca.set_params(epsilon=wca_epsilon, sigma=wca_sigma)
system.non_bonded_inter[1, 1].wca.set_params(epsilon=wca_epsilon, sigma=wca_sigma)
system.non_bonded_inter[0, 1].wca.set_params(epsilon=wca_epsilon, sigma=wca_sigma)

# thermostat
#system.thermostat.turn_off()
system.thermostat.set_langevin(kT=0.1, gamma=1.0, seed=42)

# intergrator
#system.integrator.set_steepest_descent(f_max=0, gamma=0.1, max_displacement=0.01)
system.integrator.set_vv()

output_attributes = ['type', 'mean_curvature', 'gaussian_curvature', 'normal', 'area', 'dipole']
membrane.output_vtk_data(output_folder+"swimmer_0.vtk", output_attributes)

steps = 500
cycles = 5

i=0
#system.force_cap = 1
while i < cycles:
    print(system.analysis.min_dist())
    energy = system.analysis.energy()
    #ipdb.set_trace()
    print('kinetic = {} bonded = {} non_bonded = {} dipolar = {}'.format(energy['kinetic'], energy['bonded'], energy['non_bonded'], energy['dipolar']))
    system.integrator.run(steps)
    ipdb.set_trace()
    i += 1
    membrane.output_vtk_data(output_folder+"swimmer_{}.vtk".format(i), output_attributes)
    print('current steps = {}'.format(i*steps))







