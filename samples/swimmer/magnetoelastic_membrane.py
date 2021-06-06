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

required_features = ["EXTERNAL_FORCES"]
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

time_step = 0.01

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

# stretching constants
ks_A = 1.0
ks_B = 1.0

# bending constants
kb_A = 1.0
kb_B = 1.0

# define the elastic object
rescale=(1.0, 1.0, 1.0)
elastic_object_type = oif.ElasticObjectType(nodes_file=nodes_file,
        triangles_file=triangles_file, rescale=rescale, ks_A=ks_A, ks_B=ks_B, kb_A=kb_A, kb_B=kb_B)
elastic_object_type.initialize()

# create the elastic object
translate = (0.0, 0.0, 0.0)
rotate = (0.0, 0.0, 0.0)
membrane = oif.ElasticObject(system, elastic_object_type, object_id=0, particle_type_A=0, particle_type_B=1, translate=translate, rotate=rotate)
membrane.initialize()
output_attributes = ['type', 'mean_curvature', 'gaussian_curvature', 'principal_curvatures', 'normal', 'area']
membrane.output_vtk_data(output_folder+"swimmer.vtk", output_attributes)
ipdb.set_trace()







