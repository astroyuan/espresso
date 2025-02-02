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
/** \file
 *  force cap calculation.
 *
 *  For more information see \ref forcecap.hpp "forcecap.hpp".
 */

#include "forcecap.hpp"

#include "Particle.hpp"
#include "communication.hpp"
#include "event.hpp"

#include <utils/Vector.hpp>

#include <cmath>

double force_cap = 0.0;

double forcecap_get() { return force_cap; }

void forcecap_cap(ParticleRange particles) {
  if (force_cap <= 0) {
    return;
  }

  auto const fc2 = force_cap * force_cap;

  for (auto &p : particles) {
    auto const f2 = p.f.f.norm2();
    if (f2 > fc2) {
      auto const scale = force_cap / std::sqrt(f2);

      for (int i = 0; i < 3; i++) {
        p.f.f[i] *= scale;
      }
    }
  }
}

void mpi_set_forcecap_local(double force_cap) {
  ::force_cap = force_cap;
  on_forcecap_change();
}

REGISTER_CALLBACK(mpi_set_forcecap_local)

void mpi_set_forcecap(double force_cap) {
  mpi_call_all(mpi_set_forcecap_local, force_cap);
}