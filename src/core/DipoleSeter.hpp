/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef DIPOLESETER_HPP
#define DIPOLESETER_HPP

#include "ParticleRange.hpp"
#include "rotation.hpp"

class DipoleSeter {
public:
  double precession_angle;
  double precession_freq;

  Utils::Vector3d lsdip;

private:
  Utils::Vector4d quat;
  double dipm;

  void apply(ParticleRange &particles) {
      // get current dipole
      dip[0] = 0;
      dip[1] = 0;
      dip[2] = 1;

      // convert dipole into dipm and quat
      std::tie(quat, dipm) = convert_dip_to_quat(dip);

      // set dipole strength and quaternion
      for (auto &p : particles) {
          p.p.dipm = dipm;
          for (int i = 0; i < 4; i++) {
              p.r.quat[i] = quat[i];
          }
      }
  }
};

#endif
