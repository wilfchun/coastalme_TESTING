/*!
 *
 * \file do_sediment_input_event.cpp
 * \brief Deposits sediment onto the grid
 * \details TODO A more detailed description of these routines.
 * \author David Favis-Mortlock
 * \author Andres Payo
 *
 * \date 2021
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
using std::endl;

#include "cme.h"
#include "cell.h"
#include "sediment_input_event.h"


/*===============================================================================================================================

 Check to see if we have any sediment input events this timestep, if so then do the event(s)

===============================================================================================================================*/
int CSimulation::nCheckForSedimentInputEvent(void)
{
   // Go through all sediment input events, check for any this timestep
   int nEvents = static_cast<int>(m_pVSedInputEvent.size());

   for (int n = 0; n < nEvents; n++)
   {
      if (m_pVSedInputEvent[n]->ulGetEventTimeStep() == m_ulIter)
      {
         int nRet = nDoSedimentInputEvent(n);
         if (nRet != RTN_OK)
            return nRet;
      }
   }

   return RTN_OK;
}


/*===============================================================================================================================

 Do a sediment input event

===============================================================================================================================*/
int CSimulation::nDoSedimentInputEvent(int const nEvent)
{
   // Get values for the sediment input event
   int nLocID = m_pVSedInputEvent[nEvent]->nGetLocationID();
   double
      dFineSedVol = m_pVSedInputEvent[nEvent]->dGetFineSedVol(),
      dSandSedVol = m_pVSedInputEvent[nEvent]->dGetSandSedVol(),
      dCoarseSedVol = m_pVSedInputEvent[nEvent]->dGetCoarseSedVol(),
      dLen = m_pVSedInputEvent[nEvent]->dGetLen(),
      dWidth = m_pVSedInputEvent[nEvent]->dGetWidth();

   // Now get the location from values read from the shapefile
   int
      nEvents = static_cast<int>(m_VnSedimentInputLocationID.size()),
      nEventGridX = -1,
      nEventGridY = -1;
   for (int n = 0; n < nEvents; n++)
   {
      if (m_VnSedimentInputLocationID[n] == nLocID)
      {
         nEventGridX = nRound(m_VdSedimentInputLocationX[n]);
         nEventGridY = nRound(m_VdSedimentInputLocationY[n]);
      }
   }

   // Should never get here
   if (nEventGridX == -1)
      return RTN_ERR_SEDIMENT_INPUT_EVENT;

   // Is this sediment input event at a pre-specified location, or at a point on a coast?
   if (m_bSedimentInputLocationIsExact)
   {
      // Sediment input is at a pre-specified point
      LogStream << "Sediment input event " << nEvent << " at pre-specified point [" << nEventGridX << "][" << nEventGridY << "] = {" << dGridXToExtCRSX(nEventGridX) << ", " << dGridYToExtCRSY(nEventGridY) << "] with Location ID " << nLocID;

      int nTopLayer = m_pRasterGrid->m_Cell[nEventGridX][nEventGridY].nGetTopLayerAboveBasement();

      double dFineDepth = dFineSedVol / m_dCellArea;
      m_pRasterGrid->m_Cell[nEventGridX][nEventGridY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddFine(dFineDepth);

      double dSandDepth = dSandSedVol / m_dCellArea;
      m_pRasterGrid->m_Cell[nEventGridX][nEventGridY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddSand(dSandDepth);

      double dCoarseDepth = dCoarseSedVol / m_dCellArea;
      m_pRasterGrid->m_Cell[nEventGridX][nEventGridY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddCoarse(dCoarseDepth);

      LogStream << ", depth of fine sediment = " << dFineDepth << " m, depth of sand sediment = " << dSandDepth << " m, depth of coarse sediment = " << dCoarseDepth << " m" << endl;
   }
   else
   {
      // Is at a point on a coast

   }





   return RTN_OK;
}





