/*!
 *
 * \file calc_external_forcing.cpp
 * \brief Calculates external forcings
 * \details TODO A more detailed description of these routines.
 * \author David Favis-Mortlock
 * \author Andres Payo

 * \date 2018
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
using std::cout;

#include "cme.h"
#include "simulation.h"


/*===============================================================================================================================

 Calculate external forcings: change in still water level, tide level and deep water waves height, orientation and period

===============================================================================================================================*/
int CSimulation::nCalcExternalForcing(void)
{
   // Increment SWL (note that increment may be zero)
   m_dAccumulatedSeaLevelChange += m_dDeltaSWLPerTimestep;

   unsigned int nSize = static_cast<unsigned int>(m_VdTideData.size());
   if (nSize != 0)
   {
      // We have tide data
      static unsigned int nTideDataCount = 0;

      // Wrap the tide data, i.e. start again with the first record if we do not have enough
      if (nTideDataCount > nSize-1)
         nTideDataCount = 0;

      m_dThisTimestepSWL = m_dOrigSWL + m_VdTideData[nTideDataCount] + m_dAccumulatedSeaLevelChange;
   
      // cout << m_dThisTimestepSWL << endl;
      nTideDataCount++;
   }

   // Update min and max still water levels
   m_dMaxSWL = tMax(m_dThisTimestepSWL, m_dMaxSWL);
   m_dMinSWL = tMin(m_dThisTimestepSWL, m_dMinSWL);
   
   // Update the wave height, orientation and period for this time step and start again with the first record if we do not have enough
   if (! m_bSingleDeepWaterWaveValues)
   {
      // We have wave time series data: the number of time steps is total size divided by the number of points
      unsigned int nWaveTimeSteps = static_cast<unsigned int>(m_VdDeepWaterWavePointHeightTS.size()) / static_cast<unsigned int>(m_VnDeepWaterWavePointID.size());     
      static unsigned int nWaveDataCount = 0;

       if (nWaveDataCount > nWaveTimeSteps-1)
       {
          // Wrap the tide data, i.e. start again with the first record if we do not have enough
          nWaveDataCount = 0;
       }
      
      // Update this time step deep water wave values: the order on the vector is determined by the points ID i.e. to ensure that stations match with time series
      unsigned int 
      nNumberDeepWaterWaveStations = static_cast<unsigned int>(m_VnDeepWaterWavePointID.size()),
         nTot = nNumberDeepWaterWaveStations * nWaveDataCount;
         
      for (unsigned int j = 0; j < nNumberDeepWaterWaveStations; j++)
      {
         m_VdDeepWaterWavePointHeight[j] = m_VdDeepWaterWavePointHeightTS[(m_VnDeepWaterWavePointID[j]-1) + nTot];
         m_VdDeepWaterWavePointAngle[j]  = m_VdDeepWaterWavePointAngleTS[(m_VnDeepWaterWavePointID[j]-1) + nTot];
         m_VdDeepWaterWavePointPeriod[j] = m_VdDeepWaterWavePointPeriodTS[(m_VnDeepWaterWavePointID[j]-1) + nTot];
      }
      
      nWaveDataCount++;
   }
   
   return RTN_OK;
}
