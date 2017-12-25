/*!
 *
 * \file do_beach_actual_sediment_movement.cpp
 * \brief Does between-polygon actual (supply-limited) redistribution of transported beach sediment
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
// #include <assert.h>
#include <cmath>
#include <cfloat>
#include <iostream>
using std::cout;
using std::endl;

#include <algorithm>
using std::stable_sort;

#include "cme.h"
#include "simulation.h"
#include "coast.h"


/*===============================================================================================================================

 Function used to sort polygons before doing the polygon-to-polygon sediment budget

 For both LH and RH arguments, the first value is the polygon coast ID, the second value is the down- or up-coast direction, and subsequent numbers are adjacent polygon coastIDs in that direction. If the first argument must be ordered before the second, return true

===============================================================================================================================*/
bool bPolygonAndAdjCompare(const vector<int>& nVLeft, const vector<int>& nVRight)
{
   // For safety, check that the LHS polygon has at least one adjacent polygon (it should have, apart from the bad situation where just one big polygon is created)
   if ((nVLeft.size() >= 3) && (nVRight.size() >= 3))
   {
      // Polygons at the grid edge are processed last, so put LHS grid-edge polygons on the RHS
      if (nVLeft[2] == INT_NODATA)
         return false;

      // Polygons at the grid edge are processed last, so keep RHS grid-edge polygons where they are
      if (nVRight[2] == INT_NODATA)
         return true;

      // Now sort out polygon-to-polygon dependencies. We need to put 'target' polygons after 'source' polygons, so that the source id processed before the target. So does the LHS polygon have the RHS polygon as one of its adjacent polygons?
      for (unsigned int n = 2; n < nVLeft.size(); n++)
      {
         if (nVRight[0] == nVLeft[n])
            // It does, so keep the existing sequence
            return true;
      }

      // Does the RHS polygon have the LHS polygon as one of its adjacent polygons?
      for (unsigned int n = 2; n < nVRight.size(); n++)
      {
         if (nVLeft[0] == nVRight[n])
            // It does, so swap them
            return false;
      }
   }

   bool bDownCoast = nVLeft[1];
   if (bDownCoast)
      // Sediment going down-coast
      return nVLeft < nVRight;
   else
      // Sediment going up-coast
      return nVLeft > nVRight;

   // Default return value, should never get here
   return true;
}


/*===============================================================================================================================

 Does between-polygon and within-polygon actual (supply-limited) redistribution of transported beach sediment

===============================================================================================================================*/
int CSimulation::nDoAllActualBeachErosionAndDeposition(void)
{
   int nRet = RTN_OK;

   LogStream << "ADJACENT POLYGONS, POLYGON-TO-POLYGON SHARES, AND POLYGON D50 VALUES" << endl;
   LogStream << "Num \tGlobal\tCoast\t\tUncons\t(Dirn Adj Share)..." << endl;
   LogStream << "    \tID    \tID   \t\td50" << endl;
   for (unsigned int n = 0; n < m_pVCoastPolygon.size(); n++)
   {
      LogStream << n << "\t\t" << m_pVCoastPolygon[n]->nGetGlobalID() << "\t\t\t" << m_pVCoastPolygon[n]->nGetCoastID() << "\t\t\t" << m_pVCoastPolygon[n]->dGetAvgUnconsD50() << "\t\t\t";
      for (int m = 0; m < m_pVCoastPolygon[n]->nGetNumUpCoastAdjacentPolygons(); m++)
      {
         if (! m_pVCoastPolygon[n]->bDownCoastThisTimestep())
            LogStream << "(UP  \t" << m_pVCoastPolygon[n]->nGetUpCoastAdjacentPolygon(m) << "\t" << m_pVCoastPolygon[n]->dGetUpCoastAdjacentPolygonBoundaryShare(m) << ")\t";
      }
      for (int m = 0; m < m_pVCoastPolygon[n]->nGetNumDownCoastAdjacentPolygons(); m++)
      {
         if (m_pVCoastPolygon[n]->bDownCoastThisTimestep())
            LogStream << "(DOWN\t" << m_pVCoastPolygon[n]->nGetDownCoastAdjacentPolygon(m) << "\t" << m_pVCoastPolygon[n]->dGetDownCoastAdjacentPolygonBoundaryShare(m) << ")\t";
      }
      LogStream << endl;
   }
   LogStream << endl;

   LogStream << "PER-POLYGON POTENTIAL BEACH EROSION, AND ACTUAL DEPOSITION FROM SHORE PLATFORM EROSION" << endl;
   LogStream << "Num \tGlobal\tCoast\t\tPotential\tActual\t\tFine\t\t\tSand\t\t\tCoarse" << endl;
   LogStream << "    \tID    \tID   \t\tErosion\t\tDeposition\tDeposition\tDeposition\tDeposition" << endl;
   for (unsigned int n = 0; n < m_pVCoastPolygon.size(); n++)
      LogStream << n << "\t\t" << m_pVCoastPolygon[n]->nGetGlobalID() << "\t\t\t" << m_pVCoastPolygon[n]->nGetCoastID() << "\t\t\t" << -m_pVCoastPolygon[n]->dGetDeltaPotentialErosion() << "\t\t\t\t" << m_pVCoastPolygon[n]->dGetDeltaActualTotalSediment() << "\t\t\t" << m_pVCoastPolygon[n]->dGetDeltaActualUnconsFine() <<  "\t\t\t" << m_pVCoastPolygon[n]->dGetDeltaActualUnconsSand() << "\t\t\t" << m_pVCoastPolygon[n]->dGetDeltaActualUnconsCoarse() << endl;
   LogStream << endl;

   // OK, we know the potential depth of erosion on each polygon, but we do not yet know the actual (supply-limited) depth. Nor do we know how much of the sediment which is to be removed is fine, sand or coarse. So estimate actual (supply-limited) sediment removal, but don't actually erode the polygon
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      for (int nPoly = 0; nPoly < m_VCoast[nCoast].nGetNumPolygons(); nPoly++)
      {
         // Get the potential erosion on this polygon (is a depth in m)
         double dSedChange = m_VCoast[nCoast].pGetPolygon(nPoly)->dGetDeltaPotentialErosion();
         m_dThisTimestepPotentialBeachErosion -= dSedChange;

         // And determine actual erosion on this polygon, in sediment size categories
         nRet = nEstimateActualBeachErosionOnPolygon(nCoast, nPoly, -dSedChange);
         if (nRet != RTN_OK)
            return nRet;
      }
   }

   LogStream << "BEFORE BETWEEN-POLYGON BEACH SEDIMENT ROUTING:" << endl;
   LogStream << "PER-POLYGON ESTIMATED ACTUAL BEACH EROSION (-ve) AND DEPOSITION (+ve). This does not include any unconsolidated sediment from platform erosion" << endl;
   LogStream << "Num \t\tGlobal\tCoast\t\tEstimated\t\tEstimated\tEstimated\tEstimated" << endl;
   LogStream << "    \t\tID    \tID   \t\tTotal    \t\tFine     \tSand     \tCoarse" << endl;
   double
      dCheckTotErosion = 0,
      dCheckFineErosion = 0,
      dCheckSandErosion = 0,
      dCheckCoarseErosion = 0,
      dCheckTotDeposition = 0,
      dCheckFineDeposition = 0,
      dCheckSandDeposition = 0,
      dCheckCoarseDeposition = 0;
   for (unsigned int n = 0; n < m_pVCoastPolygon.size(); n++)
   {
      LogStream << n << "\t\t\t" << m_pVCoastPolygon[n]->nGetGlobalID() << "\t\t\t" << m_pVCoastPolygon[n]->nGetCoastID() << "\t\t\t" << m_pVCoastPolygon[n]->dGetDeltaEstimatedUnconsFine() + m_pVCoastPolygon[n]->dGetDeltaEstimatedUnconsSand() + m_pVCoastPolygon[n]->dGetDeltaEstimatedUnconsCoarse() << "\t\t\t" << m_pVCoastPolygon[n]->dGetDeltaEstimatedUnconsFine() <<  "\t\t\t" << m_pVCoastPolygon[n]->dGetDeltaEstimatedUnconsSand() <<  "\t\t\t" << m_pVCoastPolygon[n]->dGetDeltaEstimatedUnconsCoarse() << endl;

      if ((m_pVCoastPolygon[n]->dGetDeltaEstimatedUnconsFine() + m_pVCoastPolygon[n]->dGetDeltaEstimatedUnconsSand() + m_pVCoastPolygon[n]->dGetDeltaEstimatedUnconsCoarse()) < 0)
         dCheckTotErosion += (m_pVCoastPolygon[n]->dGetDeltaEstimatedUnconsFine() + m_pVCoastPolygon[n]->dGetDeltaEstimatedUnconsSand() + m_pVCoastPolygon[n]->dGetDeltaEstimatedUnconsCoarse());
      else
         dCheckTotDeposition += (m_pVCoastPolygon[n]->dGetDeltaEstimatedUnconsFine() + m_pVCoastPolygon[n]->dGetDeltaEstimatedUnconsSand() + m_pVCoastPolygon[n]->dGetDeltaEstimatedUnconsCoarse());

      if (m_pVCoastPolygon[n]->dGetDeltaEstimatedUnconsFine() < 0)
         dCheckFineErosion += m_pVCoastPolygon[n]->dGetDeltaEstimatedUnconsFine();
      else
         dCheckFineDeposition += m_pVCoastPolygon[n]->dGetDeltaEstimatedUnconsFine();

      if (m_pVCoastPolygon[n]->dGetDeltaEstimatedUnconsSand() < 0)
         dCheckSandErosion += m_pVCoastPolygon[n]->dGetDeltaEstimatedUnconsSand();
      else
         dCheckSandDeposition += m_pVCoastPolygon[n]->dGetDeltaEstimatedUnconsSand();

      if (m_pVCoastPolygon[n]->dGetDeltaEstimatedUnconsCoarse() < 0)
         dCheckCoarseErosion += m_pVCoastPolygon[n]->dGetDeltaEstimatedUnconsCoarse();
      else
         dCheckCoarseDeposition += m_pVCoastPolygon[n]->dGetDeltaEstimatedUnconsCoarse();
   }
   double dActualSedimentDeliveryRatio = 0;
   if (dCheckTotErosion != 0)
      dActualSedimentDeliveryRatio = (-dCheckTotErosion - dCheckTotDeposition) / -dCheckTotErosion;

   LogStream << endl << m_ulIteration<< ": total estimated erosion = " << -dCheckTotErosion << " total estimated deposition = " << dCheckTotDeposition << " (sediment delivery ratio = " << dActualSedimentDeliveryRatio << ")" << endl;

   LogStream << m_ulIteration<< ": estimated fine erosion = " << -dCheckFineErosion << " estimated fine deposition = " << dCheckFineDeposition << endl;
   LogStream << m_ulIteration<< ": estimated sand erosion = " << -dCheckSandErosion << " estimated sand deposition = " << dCheckSandDeposition << endl;
   LogStream << m_ulIteration<< ": estimated coarse erosion = " << -dCheckCoarseErosion << " estimated coarse deposition = " << dCheckCoarseDeposition << endl << endl;;

   // Now route actually-eroded sand/coarse sediment to adjacent polygons (or off-grid)
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      // Sort first
      vector<vector<int> > nVVPolyAndAdjacent;
//      for (int nPoly = 0; nPoly < m_VCoast[nCoast].nGetNumPolygons(); nPoly++)
      for (int nPoly = m_VCoast[nCoast].nGetNumPolygons()-1; nPoly >= 0; nPoly--)
      {
         CGeomCoastPolygon* pPoly = m_VCoast[nCoast].pGetPolygon(nPoly);

         vector<int> nVPolyAndAdj;

         // The first array item is the polygon coast ID
         nVPolyAndAdj.push_back(pPoly->nGetCoastID());

         if (pPoly->bDownCoastThisTimestep())
         {
            // Sediment is leaving this polygon in a down-coast direction: so set this as the second array item
            nVPolyAndAdj.push_back(true);

            // Set subsequent array items to be the IDs of adjacent polygons
            for (int nAdj = 0; nAdj < pPoly->nGetNumDownCoastAdjacentPolygons(); nAdj++)
            {
               int nAdjPolyID = pPoly->nGetDownCoastAdjacentPolygon(nAdj);
               nVPolyAndAdj.push_back(nAdjPolyID);
            }
         }
         else
         {
            // Sediment is leaving this polygon in an up-coast direction: so set this as the second array item
            nVPolyAndAdj.push_back(false);

            // Set subsequent array items to be the IDs of adjacent polygons
            for (int nAdj = 0; nAdj < pPoly->nGetNumUpCoastAdjacentPolygons(); nAdj++)
            {
               int nAdjPolyID = pPoly->nGetUpCoastAdjacentPolygon(nAdj);
               nVPolyAndAdj.push_back(nAdjPolyID);
            }
         }

         nVVPolyAndAdjacent.push_back(nVPolyAndAdj);
      }

//       LogStream << "UNSORTED SEQUENCE OF POLYGON PROCESSING" << endl;
//       for (unsigned int n = 0; n < nVVPolyAndAdjacent.size(); n++)
//       {
//          for (unsigned int m = 0; m < nVVPolyAndAdjacent[n].size(); m++)
//             LogStream << nVVPolyAndAdjacent[n][m] << " ";
//          LogStream << endl;
//       }
//       LogStream << endl;

      // Sort the array using bPolygonAndAdjCompare(), so that 'target' polygons are processed after 'source' polygons
      stable_sort(nVVPolyAndAdjacent.begin(), nVVPolyAndAdjacent.end(), bPolygonAndAdjCompare);

//       LogStream << "SORTED SEQUENCE OF POLYGON PROCESSING" << endl;
//       for (unsigned int n = 0; n < nVVPolyAndAdjacent.size(); n++)
//       {
//          for (unsigned int m = 0; m < nVVPolyAndAdjacent[n].size(); m++)
//             LogStream << nVVPolyAndAdjacent[n][m] << " ";
//          LogStream << endl;
//       }
//       LogStream << endl;

      // Now go through the polygons in the sorted sequence
      for (int n = 0; n < m_VCoast[nCoast].nGetNumPolygons(); n++)
      {
         // And route this to adjacent polygons
         int nPoly = nVVPolyAndAdjacent[n][0];
         nRet = nRouteActualBeachErosionToAdjacentPolygons(nCoast, nPoly);
         if (nRet != RTN_OK)
            return nRet;
      }
   }

   // TEST for mass balance
   LogStream << endl << "AFTER BETWEEN-POLYGON SEDIMENT BUDGET BUT BEFORE BEACH EROSION AND DEPOSITION ON CELLS: PER-POLYGON ACTUAL EROSION/DEPOSITION" << endl << "(including any deposition of uncons sediment from shore platform erosion)" << endl;
   LogStream << "Num \tGlobal\tCoast\t\tActual\t\tFine  \t\tSand  \t\tCoarse" << endl;
   LogStream << "    \tID    \tID   \t\tTotal \t\tTotal \t\tTotal \t\tTotal" << endl;
   double
      dCheckActualErosion = 0,
      dCheckActualDeposition = 0,
      dCheckActualFineErosion = 0,
      dCheckActualFineDeposition = 0,
      dCheckActualSandErosion = 0,
      dCheckActualSandDeposition = 0,
      dCheckActualCoarseErosion = 0,
      dCheckActualCoarseDeposition = 0;

   for (unsigned int n = 0; n < m_pVCoastPolygon.size(); n++)
   {
      LogStream << n << "\t\t" << m_pVCoastPolygon[n]->nGetGlobalID() << "\t\t\t" << m_pVCoastPolygon[n]->nGetCoastID() << "\t\t\t" << m_pVCoastPolygon[n]->dGetDeltaActualTotalSediment() << "\t\t\t" << m_pVCoastPolygon[n]->dGetDeltaActualUnconsFine() <<  "\t\t\t" << m_pVCoastPolygon[n]->dGetDeltaActualUnconsSand() << "\t\t\t" << m_pVCoastPolygon[n]->dGetDeltaActualUnconsCoarse() << endl;

      if (m_pVCoastPolygon[n]->dGetDeltaActualTotalSediment() > 0)
         dCheckActualDeposition += m_pVCoastPolygon[n]->dGetDeltaActualTotalSediment();

      if (m_pVCoastPolygon[n]->dGetDeltaActualTotalSediment() < 0)
         dCheckActualErosion += m_pVCoastPolygon[n]->dGetDeltaActualTotalSediment();

      if (m_pVCoastPolygon[n]->dGetDeltaActualUnconsFine() > 0)
         dCheckActualFineDeposition += m_pVCoastPolygon[n]->dGetDeltaActualUnconsFine();

      if (m_pVCoastPolygon[n]->dGetDeltaActualUnconsFine() < 0)
         dCheckActualFineErosion += m_pVCoastPolygon[n]->dGetDeltaActualUnconsFine();

      if (m_pVCoastPolygon[n]->dGetDeltaActualUnconsSand() > 0)
         dCheckActualSandDeposition += m_pVCoastPolygon[n]->dGetDeltaActualUnconsSand();

      if (m_pVCoastPolygon[n]->dGetDeltaActualUnconsSand() < 0)
         dCheckActualSandErosion += m_pVCoastPolygon[n]->dGetDeltaActualUnconsSand();

      if (m_pVCoastPolygon[n]->dGetDeltaActualUnconsCoarse() > 0)
         dCheckActualCoarseDeposition += m_pVCoastPolygon[n]->dGetDeltaActualUnconsCoarse();

      if (m_pVCoastPolygon[n]->dGetDeltaActualUnconsCoarse() < 0)
         dCheckActualCoarseErosion += m_pVCoastPolygon[n]->dGetDeltaActualUnconsCoarse();
   }
   LogStream << endl;

   LogStream << m_ulIteration << ": all-polygon actual beach erosion = " << dCheckActualErosion << " all-polygon actual beach deposition = " << dCheckActualDeposition << " all-polygon actual beach loss from grid = " << m_dThisTimestepActualFineSedLostBeachErosion + m_dThisTimestepActualSandSedLostBeachErosion + m_dThisTimestepActualCoarseSedLostBeachErosion << endl;

   LogStream << m_ulIteration << ": all-polygon fine beach erosion = " << dCheckActualFineErosion << " all-polygon actual fine beach deposition = " << dCheckActualFineDeposition << " all-polygon actual fine beach loss from grid = " << m_dThisTimestepActualFineSedLostBeachErosion << endl;

   LogStream << m_ulIteration << ": all-polygon sand beach erosion = " << dCheckActualSandErosion << " all-polygon actual sand beach deposition = " << dCheckActualSandDeposition << " all-polygon actual sand beach loss from grid = " << m_dThisTimestepActualSandSedLostBeachErosion << endl;

   LogStream << m_ulIteration << ": all-polygon coarse beach erosion = " << dCheckActualCoarseErosion << " all-polygon actual coarse beach deposition = " << dCheckActualCoarseDeposition << " all-polygon actual coarse beach loss from grid = " << m_dThisTimestepActualCoarseSedLostBeachErosion << endl << endl;

   // We have an actual sediment budget, in sediment size categories, for all polygons: so process all polygons and do either erosion or deposition on cells within each polygon
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      for (int nPoly = 0; nPoly < m_VCoast[nCoast].nGetNumPolygons(); nPoly++)
      {
         nRet = nDoWithinPolygonBeachRedistribution(nCoast, nPoly);
         if (nRet != RTN_OK)
            return nRet;
      }
   }

   return RTN_OK;
}


/*===============================================================================================================================

 Distribute the change in actual (supply-limited) unconsolidated beach sediment from this polygon to adjacent polygons

===============================================================================================================================*/
int CSimulation::nRouteActualBeachErosionToAdjacentPolygons(int const nCoast, int const nPoly)
{
   int nNumPolygons = m_VCoast[nCoast].nGetNumPolygons() ;

   CGeomCoastPolygon* pPolygon = m_VCoast[nCoast].pGetPolygon(nPoly);

   double
      dTotFineEroded = 0,
      dTotSandEroded = 0,
      dTotCoarseEroded = 0,
      dTotFineToPoly = 0,
      dTotSandToPoly = 0,
      dTotCoarseToPoly = 0;

   // 'Actual' includes sediment from platform erosion and cliff collapse, 'estimated' was estimated in nEstimateActualBeachErosionOnPolygon()
   double
      dTotFineChange = pPolygon->dGetDeltaEstimatedUnconsFine() + pPolygon->dGetDeltaActualUnconsFine(),
      dTotSandChange = pPolygon->dGetDeltaEstimatedUnconsSand() + pPolygon->dGetDeltaActualUnconsSand(),
      dTotCoarseChange = pPolygon->dGetDeltaEstimatedUnconsCoarse() + pPolygon->dGetDeltaActualUnconsCoarse();

//    LogStream << m_ulIteration << ": polygon " << nPoly << " has actual delta unconsolidated sediment (from platform erosion and cliff collapse): fine = " << pPolygon->dGetDeltaActualUnconsFine() << " sand = " << pPolygon->dGetDeltaActualUnconsSand() << " coarse = " << pPolygon->dGetDeltaActualUnconsCoarse() << " TOTAL = " << pPolygon->dGetDeltaActualUnconsFine() + pPolygon->dGetDeltaActualUnconsSand() + pPolygon->dGetDeltaActualUnconsCoarse() << endl;
//    LogStream << m_ulIteration << ": polygon " << nPoly << " has estimated delta unconsolidated sediment (from beach erosion): fine = " << pPolygon->dGetDeltaEstimatedUnconsFine() << " sand = " << pPolygon->dGetDeltaEstimatedUnconsSand() << " coarse = " << pPolygon->dGetDeltaEstimatedUnconsCoarse() << " TOTAL = " << pPolygon->dGetDeltaEstimatedUnconsFine() + pPolygon->dGetDeltaEstimatedUnconsSand() + pPolygon->dGetDeltaEstimatedUnconsCoarse() << " (potential erosion = " << -pPolygon->dGetDeltaPotentialErosion() << ")" << endl;
//    LogStream << m_ulIteration << ": polygon " << nPoly << " has total delta unconsolidated sediment (actual plus estimated): fine = " << dTotFineChange << " sand = " << dTotSandChange << " coarse = " << dTotCoarseChange << " TOTAL = " << dTotFineChange + dTotSandChange + dTotCoarseChange << endl;

   // If any of these total change values are +ve, it means deposition on this polygon. So ignore these values in terms of routing sediment off this polygon
   if (dTotFineChange < 0)
      dTotFineEroded = -dTotFineChange;

   if (dTotSandChange < 0)
      dTotSandEroded = -dTotSandChange;

   if (dTotCoarseChange < 0)
      dTotCoarseEroded = -dTotCoarseChange;

//    LogStream << m_ulIteration << ": dTotFineEroded = " << dTotFineEroded << " dTotSandEroded = " << dTotSandEroded << " dTotCoarseEroded = " << dTotCoarseEroded << endl;

   // Now calculate actual unconsolidated sediment movement to adjacent polygons
   if (pPolygon->bDownCoastThisTimestep())
   {
      // Moving eroded sediment down-coast
      int nNumAdjPoly = pPolygon->nGetNumDownCoastAdjacentPolygons();
      for (int n = 0; n < nNumAdjPoly; n++)
      {
         int nAdjPoly = pPolygon->nGetDownCoastAdjacentPolygon(n);
         double
            dBoundaryShare = pPolygon->dGetDownCoastAdjacentPolygonBoundaryShare(n),
            dFineToPoly = dTotFineEroded * dBoundaryShare,
            dSandToPoly = dTotSandEroded * dBoundaryShare,
            dCoarseToPoly = dTotCoarseEroded * dBoundaryShare;

         dTotFineToPoly += dFineToPoly;
         dTotSandToPoly += dSandToPoly;
         dTotCoarseToPoly += dCoarseToPoly;

         if (nAdjPoly == INT_NODATA)
         {
            if (nPoly == 0)
            {
               // This is the polygon at the up-coast end of the coastline: uncons sediment movement is down-coast but there is no adjacent polygon!
               LogStream << m_ulIteration<< ": " << ERR << "polygon " << nPoly << " is at the up-coast end of the coastline, actual sediment movement is DOWN-COAST. But there is no adjacent coast-end polygon!" << endl;

//                return RTN_ERR_NO_ADJACENT_POLYGON;
            }
            else if (nPoly == nNumPolygons-1)
            {
               // TEST DFM 
               pPolygon->AddDeltaActualUnconsSand(1);               
               
               // This is the polygon at the down-coast end of the coastline, and uncons sediment movement is down-coast. Decide what to do based on the user setting m_nUnconsSedimentHandlingAtGridEdges
               if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_CLOSED)
               {
                  // Closed grid edges: no uncons sediment moves off-grid, nothing is removed from this polygon
                  LogStream << m_ulIteration << ": polygon " << nPoly << " is at the down-coast end of the coastline, actual sediment movement is DOWN-COAST. With closed grid edges, no sand or coarse unconsolidated sediment goes off-grid" << endl;
               }
               else if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_OPEN)
               {
                  // Open grid edges, so this sediment goes off-grid
                  m_dThisTimestepActualFineSedLostBeachErosion += dFineToPoly;
                  m_dThisTimestepActualSandSedLostBeachErosion += dSandToPoly;
                  m_dThisTimestepActualCoarseSedLostBeachErosion += dCoarseToPoly;

                  // Remove the sediment from this polygon
                  pPolygon->AddDeltaActualUnconsFine(-dTotFineToPoly);
                  pPolygon->AddDeltaActualUnconsSand(-dTotSandToPoly);
                  pPolygon->AddDeltaActualUnconsCoarse(-dTotCoarseToPoly);

                  LogStream << m_ulIteration << ": polygon " << nPoly << " is at the down-coast end of the coastline, actual sediment movement is DOWN-COAST. With open grid edges, sediment goes off-grid: fine = " << dFineToPoly << " sand = " << dSandToPoly << " coarse = " << dCoarseToPoly << endl;
               }
               else if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_RECIRCULATE)
               {
                  // Re-circulating grid edges, so send the sediment to the polygon at the up-coast end of this coastline TODO Check whether this causes mass balance problems, depending on the sequence of polygon processing
                  int nOtherEndPoly = 0;
                  CGeomCoastPolygon* pOtherEndPoly = m_VCoast[nCoast].pGetPolygon(nOtherEndPoly);
                  pOtherEndPoly->AddDeltaActualUnconsFine(dFineToPoly);
                  pOtherEndPoly->AddDeltaActualUnconsSand(dSandToPoly);
                  pOtherEndPoly->AddDeltaActualUnconsCoarse(dCoarseToPoly);

                  // Add to the off-grid totals even tho' it is not really moving off-grid
                  m_dThisTimestepActualFineSedLostBeachErosion += dFineToPoly;
                  m_dThisTimestepActualSandSedLostBeachErosion += dSandToPoly;
                  m_dThisTimestepActualCoarseSedLostBeachErosion += dCoarseToPoly;

                  // And remove the sediment from this polygon
                  pPolygon->AddDeltaActualUnconsFine(-dTotFineToPoly);
                  pPolygon->AddDeltaActualUnconsSand(-dTotSandToPoly);
                  pPolygon->AddDeltaActualUnconsCoarse(-dTotCoarseToPoly);

                  LogStream << m_ulIteration << ": polygon " << nPoly << " is at the down-coast end of the coastline, actual sediment movement is DOWN-COAST. With re-circulating option this goes to up-coast coast-end polygon " << nOtherEndPoly << ": fine = " << dFineToPoly << " sand = " << dSandToPoly << " coarse = " << dCoarseToPoly << endl;
               }
            }
         }
         else
         {
            // This isn't a coast-end polygon, so just move the sediment to the adjacent polygon
            CGeomCoastPolygon* pAdjPoly = m_VCoast[nCoast].pGetPolygon(nAdjPoly);
            pAdjPoly->AddDeltaActualUnconsFine(dFineToPoly);
            pAdjPoly->AddDeltaActualUnconsSand(dSandToPoly);
            pAdjPoly->AddDeltaActualUnconsCoarse(dCoarseToPoly);

            // And remove the sediment from this polygon
            pPolygon->AddDeltaActualUnconsFine(-dTotFineToPoly);
            pPolygon->AddDeltaActualUnconsSand(-dTotSandToPoly);
            pPolygon->AddDeltaActualUnconsCoarse(-dTotCoarseToPoly);

            LogStream << m_ulIteration << ": polygon " << nPoly << " has actual sediment movement DOWN-COAST to polygon " << nAdjPoly << " fine = " << dFineToPoly << " sand = " << dSandToPoly << " coarse = " << dCoarseToPoly << endl;
         }
      }
   }
   else
   {
      // Moving eroded sediment up-coast
      int nNumAdjPoly = pPolygon->nGetNumUpCoastAdjacentPolygons();
      for (int n = 0; n < nNumAdjPoly; n++)
      {
         int nAdjPoly = pPolygon->nGetUpCoastAdjacentPolygon(n);
         double
            dBoundaryShare = pPolygon->dGetUpCoastAdjacentPolygonBoundaryShare(n),
            dFineToPoly = dTotFineEroded * dBoundaryShare,
            dSandToPoly = dTotSandEroded * dBoundaryShare,
            dCoarseToPoly = dTotCoarseEroded * dBoundaryShare;

         dTotFineToPoly += dFineToPoly;
         dTotSandToPoly += dSandToPoly;
         dTotCoarseToPoly += dCoarseToPoly;

         if (nAdjPoly == INT_NODATA)
         {
            if (nPoly == nNumPolygons-1)
            {
               // This is the polygon at the down-coast end of the coastline: uncons sediment movement is up-coast but there is no adjacent polygon!
               LogStream << m_ulIteration << ": " << ERR << "polygon " << nPoly << " is at the down-coast end of the coastline, actual sediment movement is UP-COAST. But there is no adjacent coast-end polygon!" << endl;

//                return RTN_ERR_NO_ADJACENT_POLYGON;
            }
            else if (nPoly == 0)
            {
               // TEST DFM 
               pPolygon->AddDeltaActualUnconsSand(1);
               
               // This is the polygon at the up-coast end of the coastline, and uncons sediment movement is up-coast. Decide what to do based on the user setting m_nUnconsSedimentHandlingAtGridEdges
               if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_CLOSED)
               {
                  // Closed grid edges: no sand or coarse uncons sediment moves off-grid, nothing is removed from this polygon
                  LogStream << m_ulIteration << ": pPolygon " << nPoly << " is at the up-coast end of the coastline, actual sediment movement is UP-COAST. With closed grid edges, no sand or coarse unconsolidated sediment goes off-grid" << endl;
               }
               else if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_OPEN)
               {
                  // Open grid edges, so this sediment goes off-grid
                  m_dThisTimestepActualFineSedLostBeachErosion += dFineToPoly;
                  m_dThisTimestepActualSandSedLostBeachErosion += dSandToPoly;
                  m_dThisTimestepActualCoarseSedLostBeachErosion += dCoarseToPoly;

                  // Remove the sediment from this polygon
                  pPolygon->AddDeltaActualUnconsFine(-dTotFineToPoly);
                  pPolygon->AddDeltaActualUnconsSand(-dTotSandToPoly);
                  pPolygon->AddDeltaActualUnconsCoarse(-dTotCoarseToPoly);

                  LogStream << m_ulIteration << ": polygon " << nPoly << " is at the up-coast end of the coastline, actual sediment movement is UP-COAST. With open grid edges, sediment goes off-grid: fine = " << dFineToPoly << " sand = " << dSandToPoly << " coarse = " << dCoarseToPoly << endl;
               }
               else if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_RECIRCULATE)
               {
                  // Re-circulating grid edges, so send the sediment to the polygon at the down-coast end of this coastline TODO Check whether this causes mass balance problems, depending on the sequence of polygon processing
                  int nOtherEndPoly = nNumPolygons-1;
                  CGeomCoastPolygon* pOtherEndPoly = m_VCoast[nCoast].pGetPolygon(nOtherEndPoly);
                  pOtherEndPoly->AddDeltaActualUnconsFine(dFineToPoly);
                  pOtherEndPoly->AddDeltaActualUnconsSand(dSandToPoly);
                  pOtherEndPoly->AddDeltaActualUnconsCoarse(dCoarseToPoly);

                  // Add to the off-grid totals even tho' it is not really moving off-grid
                  m_dThisTimestepActualFineSedLostBeachErosion += dFineToPoly;
                  m_dThisTimestepActualSandSedLostBeachErosion += dSandToPoly;
                  m_dThisTimestepActualCoarseSedLostBeachErosion += dCoarseToPoly;

                  // And remove the sediment from this polygon
                  pPolygon->AddDeltaActualUnconsFine(-dTotFineToPoly);
                  pPolygon->AddDeltaActualUnconsSand(-dTotSandToPoly);
                  pPolygon->AddDeltaActualUnconsCoarse(-dTotCoarseToPoly);

                  LogStream << m_ulIteration << ": polygon " << nPoly << " is at the up-coast end of the coastline, actual sediment movement is UP-COAST. With re-circulating option this goes to down-coast coast-end polygon " << nOtherEndPoly << ": fine = " << dFineToPoly << " sand = " << dSandToPoly << " coarse = " << dCoarseToPoly << endl;
               }
            }
         }
         else
         {
            // This isn't a coast-end polygon, so just move the sediment to the adjacent polygon
            CGeomCoastPolygon* pAdjPoly = m_VCoast[nCoast].pGetPolygon(nAdjPoly);
            pAdjPoly->AddDeltaActualUnconsFine(dFineToPoly);
            pAdjPoly->AddDeltaActualUnconsSand(dSandToPoly);
            pAdjPoly->AddDeltaActualUnconsCoarse(dCoarseToPoly);

            // And remove the sediment from this polygon
            pPolygon->AddDeltaActualUnconsFine(-dTotFineToPoly);
            pPolygon->AddDeltaActualUnconsSand(-dTotSandToPoly);
            pPolygon->AddDeltaActualUnconsCoarse(-dTotCoarseToPoly);

            LogStream << m_ulIteration << ": polygon " << nPoly << " has actual sediment movement UP-COAST to polygon " << nAdjPoly << " fine = " << dFineToPoly << " sand = " << dSandToPoly << " coarse = " << dCoarseToPoly << endl;
         }
      }
   }

//    LogStream << "\tPolygon " << nPoly << " now has delta unconsolidated sediment: fine = " << pPolygon->dGetDeltaActualUnconsFine() << " sand = " << pPolygon->dGetDeltaActualUnconsSand() << " coarse = " << pPolygon->dGetDeltaActualUnconsCoarse() << endl << endl;

   return RTN_OK;
}
