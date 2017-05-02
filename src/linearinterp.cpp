/*!
 *
 * \brief Definitions of some routines from the linear interp library
 * \details TODO This is a more detailed description of the linear interp routines.
 * \author Modified by Andres Payo and David Favis-Mortlock
 * \date 2017
 * \copyright GNU Lesser General Public License
 *
 * \file linearinterp.cpp
 * \brief Contains definitions of routines from the linear interp library
 *
 */

/*==============================================================================================================================

 This file is part of CoastalME, the Coastal Modelling Environment.

 CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

==============================================================================================================================*/
#include "linearinterp.h"

//template<typename Real>
int nearestNeighbourIndex(std::vector<double> &x, double &value)
{
   double dist = std::numeric_limits<double>::max();
   double newDist = dist;
   size_t idx = 0;

   for (size_t i = 0; i < x.size(); ++i) 
   {
      newDist = std::abs(value - x[i]);
      if (newDist <= dist) 
      {
         dist = newDist;
         idx = i;
      }
   }

   return idx;
}


//template<typename Real>
std::vector<double> interp1(std::vector<double> &x, std::vector<double> &y, std::vector<double> &x_new)
{
   std::vector<double> y_new;
   double dx, dy;
   size_t x_max_idx = x.size() - 1;
   size_t x_new_size = x_new.size();

   y_new.reserve(x_new_size);

   for (size_t i = 0; i < x_new_size; ++i)
   {
      size_t idx = nearestNeighbourIndex(x, x_new[i]);

      if (x[idx] > x_new[i])
      {
         dx = idx > 0 ? (x[idx] - x[idx - 1]) : (x[idx + 1] - x[idx]);
         dy = idx > 0 ? (y[idx] - y[idx - 1]) : (y[idx + 1] - y[idx]);
      }
      else
      {
         dx = idx < x_max_idx ? (x[idx + 1] - x[idx]) : (x[idx] - x[idx - 1]);
         dy = idx < x_max_idx ? (y[idx + 1] - y[idx]) : (y[idx] - y[idx - 1]);
      }

      double m = dy / dx;
      double b = y[idx] - x[idx] * m;

      y_new.push_back(x_new[i] * m + b);
   }

   return y_new;
}
