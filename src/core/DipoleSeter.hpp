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
#include "integrate.hpp"
#include <stdexcept>
#include <vector>
#include <algorithm>

extern double sim_time;

class DipoleSeter {
public:
  bool &enabled() { return m_enabled; }

  double &amp() { return m_amp; }
  double &angle() { return m_angle; }
  double &freq() { return m_freq; }
  double &phase0() { return m_phase0; }
  int &axis() { return m_axis; }

  double &time() { return m_time; }
  Utils::Vector3d &H() { return m_H; }

  void set_types(std::vector<int> const &types) {
      m_types.clear();
      for (auto const &type : types) {
          m_types.push_back(type);
      }
  }
  std::vector<int> get_types() const { return m_types; }

  void set_mus(std::vector<double> const &mus) {
      m_mus.clear();
      for (auto const &mu : mus) {
          m_mus.push_back(mu);
      }
  }
  std::vector<double> get_mus() const { return m_mus; }

  void set_field_strength(double t) {
    // current field strength
    m_Hm = m_amp;
  }

  void set_field_direction(double t) {
    // current phase angle
    m_phi = t * m_freq * M_PI + m_phase0;

    switch (m_axis) {
        case 0:
            // x-axis
            m_Hdir[0] = std::cos(m_angle);
            m_Hdir[1] = std::sin(m_angle) * std::cos(m_phi);
            m_Hdir[2] = std::sin(m_angle) * std::sin(m_phi);
            break;
        case 1:
            // y-axis
            m_Hdir[0] = std::sin(m_angle) * std::sin(m_phi);
            m_Hdir[1] = std::cos(m_angle);
            m_Hdir[2] = std::sin(m_angle) * std::cos(m_phi);
            break;
        case 2:
            // z-axis
            m_Hdir[0] = std::sin(m_angle) * std::cos(m_phi);
            m_Hdir[1] = std::sin(m_angle) * std::sin(m_phi);
            m_Hdir[2] = std::cos(m_angle);
            break;
        default:
            throw std::runtime_error("invalid precession axis: " + std::to_string(m_axis));
    }
  }

  void apply(ParticleRange &particles) {
      if (m_enabled == false) return;

      // setup fields
      //printf("current simulation time: %f\n", sim_time);
      m_time = sim_time;
      set_field_strength(sim_time);
      set_field_direction(sim_time);

      // current fields
      m_H = m_Hdir * m_Hm;
      //printf("current fields: %f %f %f\n", m_H[0], m_H[1], m_H[2]);

      // get dipole direction - align with external fields
      convert_director_to_quat(m_Hdir, m_quat);

      // set dipole strength and quaternion
      for (auto &p : particles) {
          auto it_type = std::find(m_types.begin(), m_types.end(), p.p.type);
          // only apply to interested particle types
          if ( it_type != m_types.end() ) {
              auto mu = m_mus.at(it_type - m_types.begin());
              //printf("id: %d type: %d mu: %f\n", p.p.identity, *it_type, mu);

              p.p.dipm = mu * m_Hm;
              for (int i = 0; i < 4; i++) {
                p.r.quat[i] = m_quat[i];
              }
          }
      }
  }

private:
  bool m_enabled = false;

  // external field setup
  double m_amp = 0.0;
  double m_angle = 0.0;
  double m_freq = 0.0;
  double m_phi = 0.0;
  double m_phase0 = 0.0;
  int m_axis = 0;

  double m_time = 0.0;
  // field strength and direction
  Utils::Vector3d m_H;
  Utils::Vector3d m_Hdir;
  double m_Hm;

  // particle types
  std::vector<int> m_types;
  // particle dipole moments
  std::vector<double> m_mus;

  // dipole moments
  //Utils::Vector3d m_dip;
  Utils::Vector4d m_quat;
  //double m_dipm;

};

#endif
