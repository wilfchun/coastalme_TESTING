/*!
 *
 * \file calc_external_forcing.cpp
 * \brief Calculates external forcings
 * \details TODO A more detailed description of these routines.
 * \author David Favis-Mortlock
 * \author Andres Payo

 * \date 2017
 * \copyright GNU General Public License
 *
 */

/*==============================================================================================================================

 This file is part of CoastalME, the Coastal Modelling Environment.

 CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

==============================================================================================================================*/
#include <iostream>
using std::cerr;
using std::endl;

#include "cme.h"
#include "simulation.h"


/*===============================================================================================================================

 Calculate external forcings: change in still water level, TODO tide

===============================================================================================================================*/
int CSimulation::nCalcExternalForcing(void)
{
   // Increment SWL (note that increment may be zero)
   m_dThisTimestepSWL += m_dDeltaSWLPerTimestep;


//    int nSize = m_VdTideData.size();
//
//    if (nSize != 0)
//    {
//       // We have tide data
//       static int nTideDataCount = 0;
//
//       // Wrap the tide data, i.e. start again with the first record if we do not have enough
//       if (nTideDataCount > nSize-1)
//       {
//          // TODO If tide file is not long enough, then make this error message more informative and also give a warning message at start of simulation
//          string const str1 = "reached end of tide data, starting again from first line";
//          cerr << WARN << str1 << endl;
//          LogStream << WARN << str1 << endl;
//
//          nTideDataCount = 0;
//       }
//
//       m_dThisTimestepSWL = m_dOrigSWL + m_VdTideData[nTideDataCount];
//       nTideDataCount++;
//    }

   // Update min and max still water levels
   m_dMaxSWL = tMax(m_dThisTimestepSWL, m_dMaxSWL);
   m_dMinSWL = tMin(m_dThisTimestepSWL, m_dMinSWL);

   return RTN_OK;
}
