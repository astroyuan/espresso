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

#ifndef SCRIPT_INTERFACE_DIPOLESETER_HPP
#define SCRIPT_INTERFACE_DIPOLESETER_HPP

#include "ScriptInterfaceBase.hpp"
#include "auto_parameters/AutoParameters.hpp"

#include "core/dipoleseter_global.hpp"

namespace ScriptInterface {

class DipoleSeter : public AutoParameters<DipoleSeter> {
public:
  DipoleSeter() {
    add_parameters({{"enabled", dipoleseter.enabled()},
                    {"amp", dipoleseter.amp()},
                    {"angle", dipoleseter.angle()},
                    {"freq", dipoleseter.freq()},
                    {"phase0", dipoleseter.phase0()},
                    {"axis", dipoleseter.axis()},
                    {"H", dipoleseter.H()},
                    {"types", [](Variant const &v){ dipoleseter.set_types(get_value<std::vector<int>>(v));}, [](){ return dipoleseter.get_types(); }},
                    {"mus", [](Variant const &v){ dipoleseter.set_mus(get_value<std::vector<double>>(v));}, [](){ return dipoleseter.get_mus(); }}
                    });
  }
};
} // namespace ScriptInterface
#endif
