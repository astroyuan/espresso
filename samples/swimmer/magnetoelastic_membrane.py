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

###########################################################
# System Parameters
###########################################################

###########
# general #
###########
boxX = 100.0
boxY = 100.0
boxZ = 100.0

periodicX = True
periodicY = True
periodicZ = True

ncycles = 50
steps_per_cycle = 1000
time_step = 0.001

equi_steps = 10000
total_steps = ncycles*steps_per_cycle
total_time = time_step*total_steps

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

#######################
# membrane properties #
#######################
radius          = 0.25         # particle size
diameter        = 2.0 * radius
membrane_radius = 5.0          # membrane radius

# particle types
type_A = 0
type_B = 1

# stretching constants
ks_A = 400.0
ks_B = 4000.0
ks_AA = ks_A / 2.0
ks_BB = ks_B / 2.0
ks_AB = ks_A*ks_B / (ks_A + ks_B)

# bending constants
kb_A = 50.0
kb_B = 50.0
kb_AA = kb_A
kb_BB = kb_B
kb_AB = (kb_A + kb_B) / 2.0

# particle dipole moments
mu_A = 0.0
mu_B = 1.0

####################
# fluid properties #
####################
dx = 1.0 # lattice spacing
dt = time_step # LB time step
dens = 1.0 # fluid density
visc = 1.0 # fluid kinematic viscosity
ext_force_density = (0.000, 0.0, 0.0)
kT = 1.0 # temperature
seed = 42 # random seed

###################
# external fields #
###################
field_on = True
precession_angle = 90/180 * np.pi
#precession_angle = np.arccos((1./3.)**(1./2.))
precession_period = 1000*time_step
precession_freq = 1 / precession_period
precession_phi0 = 0.0
precession_axis = 0

###########################################################
# System Configuration
###########################################################

##########
# basics #
##########
system = espressomd.System(box_l=(boxX, boxY, boxZ))
system.periodicity = (periodicX, periodicY, periodicZ)
system.time_step = time_step
system.cell_system.skin = 0.2

# MPI
print("MPI configuration: ", system.cell_system.node_grid)

# CUDA devices
print("Available GPU devices: ", system.cuda_init_handle.device_list)

########################
# elastic interactions #
########################

# define the elastic object
rescale = np.ones(3) * membrane_radius
elastic_object_type = oif.ElasticObjectType(nodes_file=nodes_file,
        triangles_file=triangles_file, rescale=rescale, ks_A=ks_A, ks_B=ks_B, kb_A=kb_A, kb_B=kb_B, refState='flat')
elastic_object_type.initialize()

# create the elastic object
translate = (0.0, 0.0, 0.0)
rotate = (0.0, 0.0, 0.0)
membrane = oif.ElasticObject(system, elastic_object_type, object_id=0, particle_type_A=type_A, particle_type_B=type_B, translate=translate, rotate=rotate)
membrane.initialize()

#########################
# magnetic interactions #
#########################
# dip = (0.0, 0.0, 0.0) leads to undefined quat
for p in system.part:
    if p.type == type_A:
        p.dipm = mu_A
        p.quat = (1.0, 0.0, 0.0, 0.0)
    if p.type == type_B:
        p.dipm = mu_B
        p.quat = (1.0, 0.0, 0.0, 0.0)

# magnetic dipole-dipole interactions
#dipolar_direct_sum = espressomd.magnetostatics.DipolarDirectSumCpu(prefactor=1)
#system.actors.add(dipolar_direct_sum)

#dipolar_direct_sum = espressomd.magnetostatics.DipolarDirectSumGpu(prefactor=1)
#system.actors.add(dipolar_direct_sum)

# dipolar p3m
p3m = espressomd.magnetostatics.DipolarP3M(prefactor=1, accuracy=1e-3, mesh=32)
system.actors.add(p3m)

#######################
# steric interactions #
#######################
wca_epsilon = 1.0
wca_sigma = diameter

system.non_bonded_inter[0, 0].wca.set_params(epsilon=wca_epsilon, sigma=wca_sigma)
system.non_bonded_inter[1, 1].wca.set_params(epsilon=wca_epsilon, sigma=wca_sigma)
system.non_bonded_inter[0, 1].wca.set_params(epsilon=wca_epsilon, sigma=wca_sigma)

##########################
# impose external fields #
##########################
system.dipoleseter.enabled = True
# field strength
system.dipoleseter.amp = 1.0
# precession angle
system.dipoleseter.angle = precession_angle
# precession frequency
system.dipoleseter.freq = precession_freq
# initial phase
system.dipoleseter.phase0 = precession_phi0
# precession axis
system.dipoleseter.axis = precession_axis
# particle types to apply
system.dipoleseter.types = [type_A, type_B]
# particle dipoles per type
system.dipoleseter.mus = [mu_A, mu_B]

#############################
# hydrodynamic interactions #
#############################
lbf = espressomd.lb.LBFluidGPU(agrid=dx, dens=dens, visc=visc, tau=dt, kT=kT, seed=seed, ext_force_density=ext_force_density)
system.actors.add(lbf)
system.thermostat.set_lb(LB_fluid=lbf, gamma=10.0, seed=123)

#system.thermostat.set_langevin(kT=kT, gamma=10.0, seed=41)

###########################################################
# System Configuration Summary
###########################################################
vertices_num = elastic_object_type.mesh.vertices_num
edges_num = elastic_object_type.mesh.edges_num
faces_num = elastic_object_type.mesh.faces_num
dihedrals_num = elastic_object_type.mesh.dihedrals_num
l0 = elastic_object_type.mesh.avg_edge_length

# dimensionless parameters

### von karman number
# microscopic
para_e_A = ks_AA*l0**2/kb_AA
para_e_B = ks_BB*l0**2/kb_BB
# continuum
para_e_A_c = ks_AA*membrane_radius**2/kb_AA
para_e_B_c = ks_BB*membrane_radius**2/kb_BB

### magnetoelastic parameter
angle_factor = np.cos(precession_angle)**2 - 1.0/3.0
mesh_factor = 2*np.sqrt(3)
# microscopic
para_m_A = (3*mu_A)**2/l0**3 * angle_factor / kb_AA / 4.0
para_m_B = (3*mu_B)**2/l0**3 * angle_factor / kb_BB / 4.0
# continuum
para_m_A_c = (3*mu_A)**2/l0**5 * mesh_factor * angle_factor * membrane_radius**2 / kb_AA
para_m_B_c = (3*mu_B)**2/l0**5 * mesh_factor * angle_factor * membrane_radius**2 / kb_BB

print('\n')
print('{:#^50}'.format('Configuration Summary'))

print('{:-^50}'.format('Simulation Setup'))
print('steps: {:<8} time setp: {:<8} total time: {:<8}'.format(total_steps, time_step, total_time))

print('{:-^50}'.format('Membrane parameters'))
print('particle radius: {:<8} membrane radius: {:<8}'.format(radius,membrane_radius))
print('vertices number: {:<8} edges number: {:<8} faces number: {:<8}'.format(vertices_num, edges_num, faces_num))
print('surface nodes density: {:<8}'.format(vertices_num/(4*np.pi*membrane_radius**2)))

print('{:-^50}'.format('Fluid parameters'))

print('{:-^50}'.format('External Field Setup'))
print('precession angle: {:<8} period: {:<8} freq: {:<8} phi0: {:<8}'.format(precession_angle/np.pi*180, precession_period, precession_freq, precession_phi0))

print('{:-^50}'.format('Dimensionless Parameters'))
print('microscopic:')
print('Foppl-von Karman parameter AA: {:.8} BB:{:.8}'.format(para_e_A, para_e_B))
print('magnetoelastic parameter AA: {:.8} BB:{:.8}'.format(para_m_A, para_m_B))
print('continuum:')
print('Foppl-von Karman parameter AA: {:.8} BB:{:.8}'.format(para_e_A_c, para_e_B_c))
print('magnetoelastic parameter AA: {:.8} BB:{:.8}'.format(para_m_A_c, para_m_B_c))

###########################################################
# Equilibration
###########################################################
print("Start equilibration")
#system.timestep = time_step
#system.integrator.run(equi_steps)
print("equilibration done.")

###########################################################
# RUN
###########################################################

output_attributes = ['type', 'mean_curvature', 'gaussian_curvature', 'normal', 'area', 'dipole']

# initial frame
print('min distance: ', system.analysis.min_dist())
energy = system.analysis.energy()
print('time = {} kinetic = {} bonded = {} non_bonded = {} dipolar = {}'.format(system.time, energy['kinetic'], energy['bonded'], energy['non_bonded'], energy['dipolar']))
membrane.output_vtk_data(output_folder+"swimmer_0.vtk", output_attributes)

ipdb.set_trace()

i=0
system.time = 0
#system.force_cap = 1
while i < ncycles:
    # execute the simulation
    system.integrator.run(steps_per_cycle)
    i += 1

    # control external fields
    #if i % 20 < 10:
    #    system.dipoleseter.angle = 90/180 * np.pi
    #else:
    #    system.dipoleseter.angle = 0.0

    # output
    print('min distance: ', system.analysis.min_dist())

    energy = system.analysis.energy()
    print('time = {} kinetic = {} bonded = {} non_bonded = {} dipolar = {}'.format(system.time, energy['kinetic'], energy['bonded'], energy['non_bonded'], energy['dipolar']))

    membrane.output_vtk_data(output_folder+"swimmer_{}.vtk".format(i), output_attributes)







