/*!
 *
 * \file do_beach_within_polygon.cpp
 * \brief Does within-polygon actual erosion and distribution of transported beach sediment
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
#include <assert.h>

#include <cmath>
#include <cfloat>
#include <iostream>
using std::cout;
using std::endl;

#include "cme.h"
#include "simulation.h"
#include "coast.h"


/*===============================================================================================================================

 Does within-polygon redistribution (erosion or deposition) of unconsolidated beach sediment

===============================================================================================================================*/
int CSimulation::nDoWithinPolygonBeachRedistribution(int const nCoast, int const nPoly)
{
   // Get the depth of sediment to be redistributed on this polygon. Is a depth in m: -ve for erosion, +ve for deposition
   double dSedChange = m_VCoast[nCoast].pGetPolygon(nPoly)->dGetDeltaActualTotalSediment();

   int nRet = RTN_OK;
   if (tAbs(dSedChange) < SEDIMENT_ELEV_TOLERANCE)
      // Do nothing for tiny amounts
      return nRet;

   if (dSedChange < 0)
   {
      // Net erosion on this polygon, so calculate a net decrease in depth of unconsolidated sediment on the cells within the polygon. Note that some cells may gain in elevation (i.e. have some unconsolidated sediment deposition) however
      nRet = nDoBeachErosionOnCells(nCoast, nPoly, -dSedChange);
   }
   else
   {
      // Net deposition on this polygon, however we will only be depositing sand and coarse sediment here
      dSedChange = m_VCoast[nCoast].pGetPolygon(nPoly)->dGetDeltaActualUnconsSand() + m_VCoast[nCoast].pGetPolygon(nPoly)->dGetDeltaActualUnconsCoarse();

      // Calculate a net increase in depth of unconsolidated sediment on the cells within the polygon by depositing fine and sand. Note that some cells may decrease in elevation (i.e. have some unconsolidated sediment erosion) however
      nRet = nDoBeachDepositionOnCells(nCoast, nPoly, dSedChange);
   }

   return nRet;
}


/*===============================================================================================================================

 Erodes unconsolidated beach sediment on the cells within a polygon

 As with the estimation of actual erosion, this is done by working down the coastline and constructing profiles which are parallel to the up-coast polygon boundary; then reversing direction and going up-coast, constructing profiles parallel to the down-coast boundary. Then iteratively fit a Dean equilibrium profile until the normal's share of the change in total depth of unconsolidated sediment is accommodated under the revised profile. For erosion, this reduces the beach volume

===============================================================================================================================*/
int CSimulation::nDoBeachErosionOnCells(int const nCoast, int const nPoly, double const dTotTargetToErode)
{
   // Don't bother with tiny amounts
   if (dTotTargetToErode < SEDIMENT_ELEV_TOLERANCE)
      return RTN_OK;

   CGeomCoastPolygon* pPolygon = m_VCoast[nCoast].pGetPolygon(nPoly);

   double dSedToErodeOnThisPolygon = dTotTargetToErode;

   // Get the grid cell co-ords of this polygon's up-coast and down-coast profiles
   int
      nUpCoastProfile = pPolygon->nGetUpCoastProfile(),
      nDownCoastProfile = pPolygon->nGetDownCoastProfile();

   CGeomProfile* pUpCoastProfile = m_VCoast[nCoast].pGetProfile(nUpCoastProfile);
   CGeomProfile* pDownCoastProfile = m_VCoast[nCoast].pGetProfile(nDownCoastProfile);

   // We are using only part of each profile, seaward as far as the depth of closure. First find the seaward end point of the up-coast part-profile, this does not change as the landwards offset changes
   int nIndex = pUpCoastProfile->nGetCellGivenDepth(m_pRasterGrid, m_dDepthOfClosure);
   if (nIndex == INT_NODATA)
   {
      LogStream << m_ulIteration << ": " << ERR << "in nDoBeachErosionOnCells() for polygon " << nPoly << ", could not find the seaward end point of the up-coast profile (" << nUpCoastProfile << ") for depth of closure = " << m_dDepthOfClosure << ". Lengthen the coastline normals." << endl;

      return RTN_ERR_BAD_BEACH_EROSION_PROFILE;
   }

   // The part-profile length is one greater than nIndex, since pPtiGetCellGivenDepth() returns the index of the cell at depth of closure
   int nUpCoastPartProfileLen = nIndex + 1;

//    assert(bIsWithinValidGrid(&PtiUpCoastPartProfileSeawardEnd));

   // Store the cell co-ordinates of the boundary part-profile in reverse (sea to coast) order so we can append to the coastward end as we move inland (i.e. as nInlandOffset increases)
   vector<CGeom2DIPoint> PtiVUpCoastPartProfileCell;
   for (int n = 0; n < nUpCoastPartProfileLen; n++)
      PtiVUpCoastPartProfileCell.push_back(*pUpCoastProfile->pPtiGetCellInProfile(nUpCoastPartProfileLen - n - 1));

   int
      nUpCoastProfileCoastPoint = pUpCoastProfile->nGetNumCoastPoint(),
      nDownCoastProfileCoastPoint = pDownCoastProfile->nGetNumCoastPoint(),
      nXUpCoastProfileExistingCoastPoint = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nUpCoastProfileCoastPoint)->nGetX(),
      nYUpCoastProfileExistingCoastPoint = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nUpCoastProfileCoastPoint)->nGetY(),
      nCoastSegLen;

   // Store the coast point numbers for this polygon so that we can shuffle them
   vector<int> nVCoastPoint;
   if (nDownCoastProfileCoastPoint == m_VCoast[nCoast].nGetCoastlineSize()-1)
   {
      // This is the final down-coast polygon, so also include the down-coast polygon boundary
      nCoastSegLen = nDownCoastProfileCoastPoint - nUpCoastProfileCoastPoint + 1;
      for (int nCoastPoint = nUpCoastProfileCoastPoint; nCoastPoint <= nDownCoastProfileCoastPoint; nCoastPoint++)
         nVCoastPoint.push_back(nCoastPoint);
   }
   else
   {
      // This is not the final down-coast polygon, so do not include the polygon's down-coast boundary
      nCoastSegLen = nDownCoastProfileCoastPoint - nUpCoastProfileCoastPoint;
      for (int nCoastPoint = nUpCoastProfileCoastPoint; nCoastPoint < nDownCoastProfileCoastPoint; nCoastPoint++)
         nVCoastPoint.push_back(nCoastPoint);
   }

   // Shuffle the coast points, this is necessary so that leaving the loop does not create sequence-related artefacts
   Rand1Shuffle(&(nVCoastPoint.at(0)), nCoastSegLen);

   // Estimate the volume of sediment which is to be eroded from each parallel profile
   double
      dAllSedimentTargetPerProfile = dTotTargetToErode / nCoastSegLen,
      dTotFineEroded = 0,             // Grand totals for the whole polygon
      dTotSandEroded = 0,
      dTotCoarseEroded = 0;

   // Now traverse the polygon's existing coastline in a random (but broadly DOWN-COAST i.e. increasing coastpoint indices) sequence, at each coast point fitting a Dean profile which is parallel to the up-coast polygon boundary
   // DOWN-COAST (increasing coastpoint indices) ================================================================================
   for (int n = 0; n < nCoastSegLen; n++)
   {
      // Pick a random coast point
      int nCoastPoint = nVCoastPoint[n];
      CGeom2DIPoint PtiCoastPoint = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint);
      int
         nCoastX = PtiCoastPoint.nGetX(),
         nCoastY = PtiCoastPoint.nGetY();

      // Now calculate the x-y offset between this coast point, and the coast point of the up-coast normal
      int
         nXOffset = nCoastX - PtiVUpCoastPartProfileCell.back().nGetX(),
         nYOffset = nCoastY - PtiVUpCoastPartProfileCell.back().nGetY();

      // Get the x-y coords of a profile starting from this coast point and parallel to the up-coast polygon boundary profile (these are in reverse sequence, like the boundary part-profile)
      vector<CGeom2DIPoint> PtiVParProfile;
      for (int n = 0; n < nUpCoastPartProfileLen; n++)
      {
         CGeom2DIPoint PtiTmp(PtiVUpCoastPartProfileCell[n].nGetX() + nXOffset, PtiVUpCoastPartProfileCell[n].nGetY() + nYOffset);
         PtiVParProfile.push_back(PtiTmp);
      }

      // Get the elevations of the start and end points of the parallel profiles (as we extend the profile inland, the elevation of the new coast point of the Dean profile is set to the elevation of the original coast point)
      int
         nParProfEndX = PtiVParProfile[0].nGetX(),
         nParProfEndY = PtiVParProfile[0].nGetY();

      // Safety check
      if (! bIsWithinValidGrid(nParProfEndX, nParProfEndY))
      {
//          LogStream << WARN << "01 @@@@ while eroding polygon " << nPoly << " in DOWN-COAST direction, hit edge of grid at [" << nParProfEndX << "][" << nParProfEndY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

         KeepWithinValidGrid(nCoastX, nCoastY, nParProfEndX, nParProfEndY);
         PtiVParProfile[0].SetX(nParProfEndX);
         PtiVParProfile[0].SetY(nParProfEndY);
      }

      double
         dParProfCoastElev = m_pRasterGrid->m_Cell[nCoastX][nCoastY].dGetSedimentTopElev(),
         dParProfEndElev = m_pRasterGrid->m_Cell[nParProfEndX][nParProfEndY].dGetSedimentTopElev();

      int
         nParProfLen,
         nInlandOffset = -1;
      vector<double> VdParProfileDeanElev;

      // OK, loop until we can erode sufficient unconsolidated sediment
      bool
         bHitEdge = false,
         bEndProfile = false;
      while (true)
      {
         // Move inland by one cell
         nInlandOffset++;

         if (nInlandOffset > 0)
         {
            if (nInlandOffset > (pUpCoastProfile->nGetNumCellsInProfile()-1))
            {
//                LogStream << m_ulIteration << ": reached end of up-coast profile " << nUpCoastProfile << " during down-coast ACTUAL erosion of unconsolidated sediment for coast " << nCoast << " polygon " << nPoly << " (nCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << ")" << endl;

               bEndProfile = true;
               break;
            }

            // Extend the parallel profile by one cell in the coastward direction. First get the coords of the cell that is nInlandOffset cells seaward from the existing up-coast part-profile start point
            CGeom2DIPoint PtiUpCoastTmp = *pUpCoastProfile->pPtiGetCellInProfile(nInlandOffset);

            // Then get the offset between this PtiUpCoastTmp cell and the existing up-coast part-profile start point, and use the reverse of this offset to get the co-ords of the cell that extends the existing up-coast part-profile landwards
            int
               nXUpCoastStartOffset = PtiUpCoastTmp.nGetX() - nXUpCoastProfileExistingCoastPoint,
               nYUpCoastStartOffset = PtiUpCoastTmp.nGetY() - nYUpCoastProfileExistingCoastPoint,
//                nXUpCoastThisStart = nXUpCoastProfileExistingCoastPoint - nXUpCoastStartOffset,
//                nYUpCoastThisStart = nYUpCoastProfileExistingCoastPoint - nYUpCoastStartOffset;
               nXUpCoastThisStart = nCoastX - nXUpCoastStartOffset,
               nYUpCoastThisStart = nCoastY - nYUpCoastStartOffset;
               
            // Is the new landwards point within the raster grid?
            if (! bIsWithinValidGrid(nXUpCoastThisStart, nYUpCoastThisStart))
            {
               // It isn't
//                LogStream << WARN << "reached edge of grid at [" << nXUpCoastThisStart << "][" << nYUpCoastThisStart << "] during DOWN-COAST erosion of unconsolidated sediment for coast " << nCoast << " polygon " << nPoly << " (nCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << ")" << endl;

               // TODO Need to improve this: at present we just abandon erosion on this coast point and move to another coast point
               bHitEdge = true;
               break;

//               return RTN_ERR_EDGEOFGRID;
            }

            CGeom2DIPoint PtiThisUpCoastStart(nXUpCoastThisStart, nYUpCoastThisStart);

            // Append this new landward cell to the up-coast part-profile
//             PtiVUpCoastPartProfileCell.push_back(PtiThisUpCoastStart);

            // And calculate the co-ords of a new landwards cell for the parallel profile
            int
               nXParNew = nXUpCoastThisStart + nXOffset,
               nYParNew = nYUpCoastThisStart + nYOffset;

            // Safety check
            if (! bIsWithinValidGrid(nXParNew, nYParNew))
            {
//                LogStream << WARN << "02 @@@@ while eroding polygon " << nPoly << " in DOWN-COAST direction (nInlandOffset = " << nInlandOffset << "), hit edge of grid at [" << nParProfEndX << "][" << nParProfEndY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its landward end" << endl;

               KeepWithinValidGrid(nCoastX, nCoastY, nXParNew, nYParNew);
            }

            // Append this to the parallel profile
            CGeom2DIPoint PtiTmp(nXParNew, nYParNew);
            PtiVParProfile.push_back(PtiTmp);
         }

         // Get the distance between the start and end of the parallel profile, in external CRS units. Note that the parallel profile co-ords are in reverse sequence
         CGeom2DPoint
            PtStart = PtGridCentroidToExt(&PtiVParProfile.back()),
            PtEnd = PtGridCentroidToExt(&PtiVParProfile[0]);

         double dParProfileLen = dGetDistanceBetween(&PtStart, &PtEnd);

         // Solve for dA so that the existing elevations at the end of the parallel profile, and at the end of a Dean equilibrium profile on that part-normal, are the same
         double dParProfA = (dParProfCoastElev - dParProfEndElev) /  pow(dParProfileLen, DEAN_POWER);

         nParProfLen = PtiVParProfile.size();
         VdParProfileDeanElev.resize(nParProfLen, 0);

         double dInc = dParProfileLen / (nParProfLen-1);

         // For this eroding parallel profile, calculate the Dean equilibrium profile of the unconsolidated sediment h(y) = A * y^(2/3) where h(y) is the distance below the highest point in the profile at a distance y from the landward start of the profile
         CalcDeanProfile(&VdParProfileDeanElev, dInc, dParProfCoastElev, dParProfA, false, 0, 0);

         double dParProfTotDiff = 0;

         for (int n = 0; n < nParProfLen; n++)
         {
            CGeom2DIPoint PtiTmp = PtiVParProfile[nParProfLen - n - 1];
            int
               nX = PtiTmp.nGetX(),
               nY = PtiTmp.nGetY();

            // Safety check
            if (! bIsWithinValidGrid(nX, nY))
            {
//                LogStream << WARN << "03 @@@@ while constructing parallel profile to assess DOWN-COAST erosion on polygon " << nPoly << ", hit edge of grid at [" << nParProfEndX << "][" << nParProfEndY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

               KeepWithinValidGrid(nCoastX, nCoastY, nX, nY);
               PtiTmp.SetX(nX);
               PtiTmp.SetY(nY);
            }

            // Don't do anything to intervention cells
            if (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION)
               continue;

            // Don't do cells twice
            if (! m_pRasterGrid->m_Cell[nX][nY].bBeachErosionOrDepositionThisTimestep())
            {
               double
                  dTmpElev = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev(),
                  dDiff = dTmpElev - VdParProfileDeanElev[n];

               dParProfTotDiff += dDiff;
            }
         }

//       // DEBUG STUFF -----------------------------------------------------
//          if (m_ulIteration == 16)
//          {
//             LogStream << "\tFor polygon " << nPoly << " doing DOWN-COAST erosion, parallel profile from [" << PtiVParProfile.back().nGetX() << "][" << PtiVParProfile.back().nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiVParProfile.back().nGetX()) << ", " <<  dGridCentroidYToExtCRSY(PtiVParProfile.back().nGetY()) << "} to [" << PtiVParProfile[0].nGetX() << "][" << PtiVParProfile[0].nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiVParProfile[0].nGetX()) << ", " <<  dGridCentroidYToExtCRSY(PtiVParProfile[0].nGetY()) << "}, nParProfLen = " << nParProfLen << " dParProfileLen = " << dParProfileLen << " dParProfCoastElev = " << dParProfCoastElev << " dParProfEndElev = " << dParProfEndElev << " dParProfA = " << dParProfA << endl;
//
//             LogStream << "\tPresent profile = ";
//             for (int n = 0; n < nParProfLen; n++)
//             {
//                CGeom2DIPoint PtiTmp = PtiVParProfile[nParProfLen - n - 1];
//                int
//                   nX = PtiTmp.nGetX(),
//                   nY = PtiTmp.nGetY();
//
//                // Safety check
//                if (! bIsWithinValidGrid(nX, nY))
//                {
//                   KeepWithinValidGrid(nX, nY);
//                }
//
//                // Don't do anything to intervention cells
//                if (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION)
//                   continue;
//
//                // Don't do cells twice
//                if (! m_pRasterGrid->m_Cell[nX][nY].bBeachErosionOrDepositionThisTimestep())
//                {
//                   LogStream << m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev() << " ";
//                }
//             }
//             LogStream << endl;
//
//             LogStream << "\tParallel Dean profile for erosion = ";
//             for (int n = 0; n < nParProfLen; n++)
//             {
//                LogStream << VdParProfileDeanElev[n] << " ";
//             }
//             LogStream << endl;
//
//             LogStream << "\tnCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << endl << "Difference = ";
//             for (int n = 0; n < nParProfLen; n++)
//             {
//                CGeom2DIPoint PtiTmp = PtiVParProfile[nParProfLen - n - 1];
//                int
//                   nX = PtiTmp.nGetX(),
//                   nY = PtiTmp.nGetY();
//
//                // Safety check
//                if (! bIsWithinValidGrid(nX, nY))
//                {
//                   KeepWithinValidGrid(nX, nY);
//                }
//
//                // Don't do anything to intervention cells
//                if (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION)
//                   continue;
//
//                // Don't do cells twice
//                if (! m_pRasterGrid->m_Cell[nX][nY].bBeachErosionOrDepositionThisTimestep())
//                {
//                   double
//                      dTmpElev = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev(),
//                      dDiff = dTmpElev - VdParProfileDeanElev[n];
//                   LogStream << dDiff << " ";
//                }
//             }
//             LogStream << endl << endl;
//          }
//       // DEBUG STUFF -----------------------------------------------------

         // So will we be able to erode as much as is needed?
         if (dParProfTotDiff > dAllSedimentTargetPerProfile)
         {
//             LogStream << "\tDOWN-COAST nCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << " dParProfTotDiff = " << dParProfTotDiff << " dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << " ENOUGH TO BE ERODED" << endl;

            break;
         }
      }

      // TODO Improve this, see above
      if (bHitEdge || bEndProfile)
         continue;

//       assert(dParProfTotDiff > 0);

      // OK, this value of nInlandOffset gives us enough erosion. So start eroding the parallel profile from the landward end (which is inland from the existing coast, if nInlandOffset > 0)
      double
         dFineEroded = 0,                 // Totals for this parallel profile
         dSandEroded = 0,
         dCoarseEroded = 0,
         dSedToErodeOnThisProfile = tMin(dAllSedimentTargetPerProfile, dSedToErodeOnThisPolygon);

//       LogStream << "dRatio = " << dRatio << endl;

      for (int nDistSeawardFromNewCoast = 0; nDistSeawardFromNewCoast < nParProfLen; nDistSeawardFromNewCoast++)
      {
         // Don't bother with tiny amounts
         if (dSedToErodeOnThisPolygon < SEDIMENT_ELEV_TOLERANCE)
            dSedToErodeOnThisPolygon = 0;
         if (dSedToErodeOnThisProfile < SEDIMENT_ELEV_TOLERANCE)
            dSedToErodeOnThisProfile = 0;

         // Leave the loop if we have eroded enough for this polygon
         if (dSedToErodeOnThisPolygon <= 0)
         {
//             LogStream << "\tIn nDoBeachErosionOnCells() going DOWN-COAST, nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << " leaving loop at start of timestep (" << nDistSeawardFromNewCoast << " / " << nParProfLen << ") because enough erosion for polygon, dSedToErodeOnThisPolygon = " << dSedToErodeOnThisPolygon << endl;
            break;
         }

         // Leave the loop if we have done enough erosion for this profile
         if (dSedToErodeOnThisProfile <= 0)
         {
//             LogStream << "\tIn nDoBeachErosionOnCells() going DOWN-COAST, nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << " leaving loop at start of timestep (" << nDistSeawardFromNewCoast << " / " << nParProfLen << ") because enough erosion for profile, dSedToErodeOnThisProfile = " << dSedToErodeOnThisProfile << " dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << endl;

            break;
         }

         CGeom2DIPoint PtiTmp = PtiVParProfile[nParProfLen - nDistSeawardFromNewCoast - 1];
         int
            nX = PtiTmp.nGetX(),
            nY = PtiTmp.nGetY();

         // Safety check
         if (! bIsWithinValidGrid(nX, nY))
         {
//                LogStream << WARN << "04 @@@@ while eroding polygon " << nPoly << " in DOWN-COAST direction, hit edge of grid at [" << nParProfEndX << "][" << nParProfEndY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

            KeepWithinValidGrid(nCoastX, nCoastY, nX, nY);
            PtiTmp.SetX(nX);
            PtiTmp.SetY(nY);
         }

         // Don't do anything to intervention cells
         if (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION)
            continue;

         // Don't do cells twice
         if (! m_pRasterGrid->m_Cell[nX][nY].bBeachErosionOrDepositionThisTimestep())
         {
            // Get this cell's current elevation
            double dThisElevNow = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev();

//             LogStream << "\tnPoly = " << nPoly << " going DOWN-COAST, [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "}  nCoastPoint = " << nCoastPoint << " nDistSeawardFromNewCoast = " << nDistSeawardFromNewCoast << " dThisElevNow = " << dThisElevNow << " Dean Elev = " << VdParProfileDeanElev[nDistSeawardFromNewCoast] << endl;

            // Subtract the two elevations
            double dElevDiff = dThisElevNow - VdParProfileDeanElev[nDistSeawardFromNewCoast];
            if (dElevDiff > SEDIMENT_ELEV_TOLERANCE)
            {
               // The current elevation is higher than the Dean elevation, so we have potential beach erosion (i.e. not constrained by availability of unconsolidated sediment) here
               m_ulThisTimestepNumPotentialBeachErosionCells++;

               m_pRasterGrid->m_Cell[nX][nY].SetPotentialBeachErosion(dElevDiff);

//                LogStream << "\tnPoly = " << nPoly << " doing DOWN-COAST, potential beach erosion = " << dElevDiff << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "} nCoastPoint = " << nCoastPoint << " nDistSeawardFromNewCoast = " << nDistSeawardFromNewCoast << " dThisElevNow = " << dThisElevNow << " Dean Elev = " << VdParProfileDeanElev[nDistSeawardFromNewCoast] << endl;

               // Now get the number of the highest layer with non-zero thickness
               int nThisLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();

               // Safety check
               if (nThisLayer == INT_NODATA)
                  return RTN_ERR_NO_TOP_LAYER;

               if (nThisLayer != NO_NONZERO_THICKNESS_LAYERS)
               {
                  // We still have at least one layer left with non-zero thickness (i.e. we are not down to basement), and the cell's current elevation is higher than the Dean equilibrium profile elevation. So do some beach erosion
                  double
                     dToErode = tMin(dElevDiff, dSedToErodeOnThisProfile, dSedToErodeOnThisPolygon),
                     dFine = 0,
                     dSand = 0,
                     dCoarse = 0;
   //                assert(dToErode > 0);
                  ErodeBeachConstrained(nX, nY, nThisLayer, dToErode, dFine, dSand, dCoarse);

                  // Update totals for this parallel profile
                  dFineEroded += dFine;
                  dSandEroded += dSand;
                  dCoarseEroded += dCoarse;

                  // Update totals for the polygon
                  dTotFineEroded += dFine;
                  dTotSandEroded += dSand;
                  dTotCoarseEroded += dCoarse;

                  double dTmpTot = dFine + dSand + dCoarse;
                  dSedToErodeOnThisProfile -= dTmpTot;
                  dSedToErodeOnThisPolygon -= dTmpTot;

                  // Update this-timestep totals
                  m_ulThisTimestepNumActualBeachErosionCells++;
                  m_dThisTimestepActualBeachErosionFine += dFine;
                  m_dThisTimestepActualBeachErosionSand += dSand;
                  m_dThisTimestepActualBeachErosionCoarse += dCoarse;

//                   LogStream << "\tIn nDoBeachErosionOnCells(), nPoly = " << nPoly << " going DOWN-COAST, actual beach erosion = " << dTmpTot << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "}  nCoastPoint = " << nCoastPoint << " nDistSeawardFromNewCoast = " << nDistSeawardFromNewCoast << endl;
               }
            }
            else if ((dElevDiff < -SEDIMENT_ELEV_TOLERANCE) && ((m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea()) || (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_DRIFT)))
            {
               // The current elevation is below the Dean elevation, so we have can have beach deposition here provided that we have some previously-eroded unconsolidated sediment to deposit
               double dTotSandAndCoarseEroded = dTotSandEroded + dTotCoarseEroded;
               if (dTotSandAndCoarseEroded > SEDIMENT_ELEV_TOLERANCE)
               {
                  // Assume that the sediment size fractions which are deposited are in the same ratios as the previously-eroded sediment
                  double
                     dTotToDeposit = tMin(-dElevDiff, dTotSandAndCoarseEroded),
                     dSandToDeposit = dTotToDeposit * dTotSandEroded / dTotSandAndCoarseEroded,
                     dCoarseToDeposit = dTotToDeposit * dTotCoarseEroded / dTotSandAndCoarseEroded;

                  int nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopLayerAboveBasement();

                  // Safety check
                  if (nTopLayer == INT_NODATA)
                     return RTN_ERR_NO_TOP_LAYER;

                  if (dSandToDeposit > SEDIMENT_ELEV_TOLERANCE)
                  {
                     dSandToDeposit = tMin(dSandToDeposit, dTotSandEroded);

                     double dSandNow = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetSand();
                     m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetSand(dSandNow + dSandToDeposit);

                     // Set the changed-this-timestep switch
                     m_bUnconsChangedThisTimestep[nTopLayer] = true;

                     dTotSandEroded -= dSandToDeposit;

                     dSedToErodeOnThisProfile += dSandToDeposit;
                     dSedToErodeOnThisPolygon += dSandToDeposit;
                  }

                  if (dCoarseToDeposit > SEDIMENT_ELEV_TOLERANCE)
                  {
                     dCoarseToDeposit = tMin(dCoarseToDeposit, dTotCoarseEroded);

                     double dCoarseNow = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetCoarse();
                     m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetCoarse(dCoarseNow + dCoarseToDeposit);

                     // Set the changed-this-timestep switch
                     m_bUnconsChangedThisTimestep[nTopLayer] = true;

                     dTotCoarseEroded -= dCoarseToDeposit;

                     dSedToErodeOnThisProfile += dCoarseToDeposit;
                     dSedToErodeOnThisPolygon += dCoarseToDeposit;
                  }

                  // Now update the cell's layer elevations
                  m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();

                  // Update the cell's sea depth
                  m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();

                  // Update the cell's beach deposition, and total beach deposition, values
                  m_pRasterGrid->m_Cell[nX][nY].IncrBeachDeposition(dSandToDeposit + dCoarseToDeposit);

                  // And set the landform category
                  m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetLFSubCategory(LF_SUBCAT_DRIFT_BEACH);

                  // Update per-timestep totals
                  m_ulThisTimestepNumBeachDepositionCells++;
                  m_dThisTimestepBeachDepositionSand += dSandToDeposit;
                  m_dThisTimestepBeachDepositionCoarse += dCoarseToDeposit;

//                   LogStream << "\tIn nDoBeachErosionOnCells(), nPoly = " << nPoly << " going DOWN-COAST, beach deposition = " << dTotToDeposit << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "}  nCoastPoint = " << nCoastPoint << " nDistSeawardFromNewCoast = " << nDistSeawardFromNewCoast << endl;
               }
            }
         }
      }

//       LogStream << "\tIn nDoBeachErosionOnCells() going DOWN-COAST, nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << " dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << " dSedToErodeOnThisProfile = " <<  dAllSedimentTargetPerProfile << " dSedToErodeOnThisPolygon = " << dSedToErodeOnThisPolygon << " dTotTargetToErode = " << dTotTargetToErode << endl;
   }

   // Have we eroded the full potential-erosion amount?
   if (dSedToErodeOnThisPolygon > SEDIMENT_ELEV_TOLERANCE)
   {
      // No, so do the same in a broadly UP-COAST direction, but only if each cell has not previously been eroded this timestep
      // UP-COAST (decreasing coastpoint indices) ==================================================================================

      // Start by finding the seaward end point of the up-coast part-profile, this does not change as the landwards offset changes
      int nIndex = pDownCoastProfile->nGetCellGivenDepth(m_pRasterGrid, m_dDepthOfClosure);
      if (nIndex == INT_NODATA)
      {
         LogStream << m_ulIteration << ": " << ERR << "in nDoBeachErosionOnCells() for polygon " << nPoly << ", could not find the seaward end point of the down-coast profile (" << nUpCoastProfile << ") for depth of closure = " << m_dDepthOfClosure << ". Lengthen the coastline normals." << endl;

         return RTN_ERR_BAD_BEACH_EROSION_PROFILE;
      }

      // The part-profile length is one greater than nIndex, since pPtiGetCellGivenDepth() returns the index of the cell at depth of closure
      int nDownCoastPartProfileLen = nIndex + 1;

//       assert(bIsWithinValidGrid(&PtiDownCoastPartProfileSeawardEnd));

      // Store the cell co-ordinates of the boundary part-profile in reverse (sea to coast) order so we can append to the coastward end as we move inland (i.e. as nInlandOffset increases)
      vector<CGeom2DIPoint> PtiVDownCoastPartProfileCell;
      for (int n = 0; n < nDownCoastPartProfileLen; n++)
         PtiVDownCoastPartProfileCell.push_back(*pDownCoastProfile->pPtiGetCellInProfile(nDownCoastPartProfileLen - n - 1));

      int
         nXDownCoastProfileExistingCoastPoint = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nDownCoastProfileCoastPoint)->nGetX(),
         nYDownCoastProfileExistingCoastPoint = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nDownCoastProfileCoastPoint)->nGetY();

      // Store the coast point numbers for this polygon so that we can shuffle them
      nVCoastPoint.resize(0);
      if (nUpCoastProfileCoastPoint == 0)
      {
         // This is the first up-coast polygon, so also include the up-coast polygon boundary
         nCoastSegLen = nDownCoastProfileCoastPoint - nUpCoastProfileCoastPoint + 1;
         for (int nCoastPoint = nUpCoastProfileCoastPoint; nCoastPoint <= nDownCoastProfileCoastPoint; nCoastPoint++)
            nVCoastPoint.push_back(nCoastPoint);
      }
      else
      {
         // This is not the final down-coast polygon, so do not include the polygon's down-coast boundary
         nCoastSegLen = nDownCoastProfileCoastPoint - nUpCoastProfileCoastPoint;
         for (int nCoastPoint = nUpCoastProfileCoastPoint; nCoastPoint < nDownCoastProfileCoastPoint; nCoastPoint++)
            nVCoastPoint.push_back(nCoastPoint);
      }

      // Shuffle the coast points, this is necessary so that leaving the loop does not create sequence-related artefacts
      Rand1Shuffle(&(nVCoastPoint.at(0)), nCoastSegLen);

      // Now traverse the polygon's existing coastline in a random sequence, fitting a Dean profile at each coast point
      for (int n = 0; n < nCoastSegLen; n++)
      {
         // Pick a random coast point
         int nCoastPoint = nVCoastPoint[n];
         CGeom2DIPoint PtiCoastPoint = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint);
         int
            nCoastX = PtiCoastPoint.nGetX(),
            nCoastY = PtiCoastPoint.nGetY();

         // Now calculate the x-y offset between this coast point, and the coast point of the down-coast normal
         int
            nXOffset = nCoastX - PtiVDownCoastPartProfileCell.back().nGetX(),
            nYOffset = nCoastY - PtiVDownCoastPartProfileCell.back().nGetY();

         // Get the x-y coords of a profile starting from this coast point and parallel to the down-coast polygon boundary profile (these are in reverse sequence, like the boundary part-profile)
         vector<CGeom2DIPoint> PtiVParProfile;
         for (int n = 0; n < nDownCoastPartProfileLen; n++)
         {
            CGeom2DIPoint PtiTmp(PtiVDownCoastPartProfileCell[n].nGetX() + nXOffset, PtiVDownCoastPartProfileCell[n].nGetY() + nYOffset);
            PtiVParProfile.push_back(PtiTmp);
         }

         // Get the elevations of the start and end points of the parallel profiles (as we extend the profile inland, the elevation of the new coast point of the Dean profile is set to the elevation of the original coast point)
         int
            nParProfEndX = PtiVParProfile[0].nGetX(),
            nParProfEndY = PtiVParProfile[0].nGetY();

         // Safety check
         if (! bIsWithinValidGrid(nParProfEndX, nParProfEndY))
         {
//             LogStream << WARN << "05 @@@@ while eroding polygon " << nPoly << " in UP-COAST direction, hit edge of grid at [" << nParProfEndX << "][" << nParProfEndY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

            KeepWithinValidGrid(nCoastX, nCoastY, nParProfEndX, nParProfEndY);
            PtiVParProfile[0].SetX(nParProfEndX);
            PtiVParProfile[0].SetY(nParProfEndY);
         }

         double
            dParProfCoastElev = m_pRasterGrid->m_Cell[nCoastX][nCoastY].dGetSedimentTopElev(),
            dParProfEndElev = m_pRasterGrid->m_Cell[nParProfEndX][nParProfEndY].dGetSedimentTopElev();

         int
            nParProfLen,
            nInlandOffset = -1;
         vector<double> VdParProfileDeanElev;

         // OK, loop until we can erode sufficient unconsolidated sediment
         bool
            bEndProfile = false,
            bHitEdge = false;
         while (true)
         {
            // Move inland by one cell
            nInlandOffset++;

            if (nInlandOffset > 0)
            {
               if (nInlandOffset > (pDownCoastProfile->nGetNumCellsInProfile()-1))
               {
//                   LogStream << m_ulIteration << ": reached end of down-coast profile " << nDownCoastProfile << " during up-coast ACTUAL erosion of unconsolidated sediment for coast " << nCoast << " polygon " << nPoly << " (nCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << ")" << endl;

                  bEndProfile = true;
                  break;
               }

               // Extend the parallel profile by one cell in the coastward direction. First get the coords of the cell that is nInlandOffset cells seaward from the existing down-coast part-profile start point
               CGeom2DIPoint PtiDownCoastTmp = *pDownCoastProfile->pPtiGetCellInProfile(nInlandOffset);

               // Then get the offset between this PtiDownCoastTmp cell and the existing down-coast part-profile start point, and use the reverse of this offset to get the co-ords of the cell that extends the existing up-coast part-profile landwards
               int
                  nXDownCoastStartOffset = PtiDownCoastTmp.nGetX() - nXDownCoastProfileExistingCoastPoint,
                  nYDownCoastStartOffset = PtiDownCoastTmp.nGetY() - nYDownCoastProfileExistingCoastPoint,
//                   nXDownCoastThisStart = nXDownCoastProfileExistingCoastPoint - nXDownCoastStartOffset,
//                   nYDownCoastThisStart = nYDownCoastProfileExistingCoastPoint - nYDownCoastStartOffset;
                  nXDownCoastThisStart = nCoastX - nXDownCoastStartOffset,
                  nYDownCoastThisStart = nCoastY - nYDownCoastStartOffset;
                  
               // Is the new landwards point within the raster grid?
               if (! bIsWithinValidGrid(nXDownCoastThisStart, nYDownCoastThisStart))
               {
                  // It isn't
//                   LogStream << m_ulIteration << ": " << WARN << "reached edge of grid at [" << nXDownCoastThisStart << "][" << nYDownCoastThisStart << "] during UP-COAST erosion of unconsolidated sediment for coast " << nCoast << " polygon " << nPoly << " (nCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << ")" << endl;

                  // TODO Need to improve this: at present just abandons erosion on this coast point and moves to another coast point
                  bHitEdge = true;
                  break;

//                  return RTN_ERR_EDGEOFGRID;
               }

               CGeom2DIPoint PtiThisDownCoastStart(nXDownCoastThisStart, nYDownCoastThisStart);

               // Append this new landward cell to the down-coast part-profile
//                PtiVDownCoastPartProfileCell.push_back(PtiThisDownCoastStart);

               // And calculate the co-ords of a new landwards cell for the parallel profile
               int
                  nXParNew = nXDownCoastThisStart + nXOffset,
                  nYParNew = nYDownCoastThisStart + nYOffset;

               // Safety check
               if (! bIsWithinValidGrid(nXParNew, nYParNew))
               {
//                   LogStream << WARN << "06 @@@@ while eroding polygon " << nPoly << " in UP-COAST direction, hit edge of grid at [" << nParProfEndX << "][" << nParProfEndY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its landward end" << endl;

                  KeepWithinValidGrid(nCoastX, nCoastY, nXParNew, nYParNew);
               }

               // Append this to the parallel profile
               CGeom2DIPoint PtiTmp(nXParNew, nYParNew);
               PtiVParProfile.push_back(PtiTmp);
            }

            // Get the distance between the start and end of the parallel profile, in external CRS units. Note that the parallel profile co-ords are in reverse sequence
            CGeom2DPoint
               PtStart = PtGridCentroidToExt(&PtiVParProfile.back()),
               PtEnd = PtGridCentroidToExt(&PtiVParProfile[0]);

            double dParProfileLen = dGetDistanceBetween(&PtStart, &PtEnd);

            // Solve for dA so that the existing elevations at the end of the parallel profile, and at the end of a Dean equilibrium profile on that part-normal, are the same
            double dParProfA = (dParProfCoastElev - dParProfEndElev) /  pow(dParProfileLen, DEAN_POWER);

            nParProfLen = PtiVParProfile.size();
            VdParProfileDeanElev.resize(nParProfLen, 0);

            double dInc = dParProfileLen / (nParProfLen-1);

            // For this eroding parallel profile, calculate the Dean equilibrium profile of the unconsolidated sediment h(y) = A * y^(2/3) where h(y) is the distance below the highest point in the profile at a distance y from the landward start of the profile
            CalcDeanProfile(&VdParProfileDeanElev, dInc, dParProfCoastElev, dParProfA, false, 0, 0);

            double dParProfTotDiff = 0;
            for (int n = 0; n < nParProfLen; n++)
            {
               CGeom2DIPoint PtiTmp = PtiVParProfile[nParProfLen - n - 1];
               int
                  nX = PtiTmp.nGetX(),
                  nY = PtiTmp.nGetY();

               // Safety check
               if (! bIsWithinValidGrid(nX, nY))
               {
//                   LogStream << WARN << "07 @@@@ while constructing parallel profile to assess UP-COAST erosion on polygon " << nPoly << ", hit edge of grid at [" << nParProfEndX << "][" << nParProfEndY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

                  KeepWithinValidGrid(nCoastX, nCoastY, nX, nY);
                  PtiTmp.SetX(nX);
                  PtiTmp.SetY(nY);
               }

               // Don't do anything to intervention cells
               if (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION)
                  continue;

               // Don't do cells twice
               if (! m_pRasterGrid->m_Cell[nX][nY].bBeachErosionOrDepositionThisTimestep())
               {
                  double
                     dTmpElev = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev(),
                     dDiff = dTmpElev - VdParProfileDeanElev[n];

                  dParProfTotDiff += dDiff;
               }
            }

//          // DEBUG STUFF -----------------------------------------------------
//             if (m_ulIteration == 16)
//             {
//                LogStream << "\tFor polygon " << nPoly << " doing UP-COAST erosion, parallel profile from [" << PtiVParProfile.back().nGetX() << "][" << PtiVParProfile.back().nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiVParProfile.back().nGetX()) << ", " <<  dGridCentroidYToExtCRSY(PtiVParProfile.back().nGetY()) << "} to [" << PtiVParProfile[0].nGetX() << "][" << PtiVParProfile[0].nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiVParProfile[0].nGetX()) << ", " <<  dGridCentroidYToExtCRSY(PtiVParProfile[0].nGetY()) << "}, nParProfLen = " << nParProfLen << " dParProfileLen = " << dParProfileLen << " dParProfCoastElev = " << dParProfCoastElev << " dParProfEndElev = " << dParProfEndElev << " dParProfA = " << dParProfA << endl;
//
//                LogStream << "\tPresent profile = ";
//                for (int n = 0; n < nParProfLen; n++)
//                {
//                   CGeom2DIPoint PtiTmp = PtiVParProfile[nParProfLen - n - 1];
//                   int
//                      nX = PtiTmp.nGetX(),
//                      nY = PtiTmp.nGetY();
//
//                   // Safety check
//                   if (! bIsWithinValidGrid(nX, nY))
//                   {
//                      KeepWithinValidGrid(nX, nY);
//                   }
//
//                      // Don't do anything to intervention cells
//                      if (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION)
//                         continue;
//
//                   // Don't do cells twice
//                   if (! m_pRasterGrid->m_Cell[nX][nY].bBeachErosionOrDepositionThisTimestep())
//                   {
//                      LogStream << m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev() << " ";
//                   }
//                }
//                LogStream << endl;
//
//                LogStream << "\tParallel Dean profile for erosion = ";
//                for (int n = 0; n < nParProfLen; n++)
//                {
//                   LogStream << VdParProfileDeanElev[n] << " ";
//                }
//                LogStream << endl;
//
//                LogStream << "\tnCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << endl << "Difference = ";
//                for (int n = 0; n < nParProfLen; n++)
//                {
//                   CGeom2DIPoint PtiTmp = PtiVParProfile[nParProfLen - n - 1];
//                   int
//                      nX = PtiTmp.nGetX(),
//                      nY = PtiTmp.nGetY();
//
//                   // Safety check
//                   if (! bIsWithinValidGrid(nX, nY))
//                   {
//                      KeepWithinValidGrid(nX, nY);
//                   }
//
//                      // Don't do anything to intervention cells
//                      if (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION)
//                         continue;
//
//                   // Don't do cells twice
//                   if (! m_pRasterGrid->m_Cell[nX][nY].bBeachErosionOrDepositionThisTimestep())
//                   {
//                      double
//                         dTmpElev = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev(),
//                         dDiff = dTmpElev - VdParProfileDeanElev[n];
//                      LogStream << dDiff << " ";
//                   }
//                }
//                LogStream << endl << endl;
//             }
//             // DEBUG STUFF -----------------------------------------------------

            // So will we be able to erode as much as is needed?
            if (dParProfTotDiff >= dAllSedimentTargetPerProfile)
            {
//                LogStream << "\tUP-COAST nCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << " dParProfTotDiff = " << dParProfTotDiff << " dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << " ENOUGH TO BE ERODED" << endl;

               break;
            }
         }

         // TODO Improve this, see above
         if (bHitEdge || bEndProfile)
            continue;

//          assert(dParProfTotDiff > 0);

         // OK, this value of nInlandOffset gives us enough erosion. So start eroding the parallel profile from the landward end (which is inland from the existing coast, if nInlandOffset > 0)
         double
            dFineEroded = 0,                 // Totals for this parallel profile
            dSandEroded = 0,
            dCoarseEroded = 0,
            dSedToErodeOnThisProfile = tMin(dAllSedimentTargetPerProfile, dSedToErodeOnThisPolygon);

         for (int nDistSeawardFromNewCoast = 0; nDistSeawardFromNewCoast < nParProfLen; nDistSeawardFromNewCoast++)
         {
            // Don't bother with tiny amounts
            if (dSedToErodeOnThisPolygon < SEDIMENT_ELEV_TOLERANCE)
               dSedToErodeOnThisPolygon = 0;
            if (dSedToErodeOnThisProfile < SEDIMENT_ELEV_TOLERANCE)
               dSedToErodeOnThisProfile = 0;

            // Leave the loop if we have eroded enough for this polygon
            if (dSedToErodeOnThisPolygon <= 0)
            {
//                LogStream << "\tIn nDoBeachErosionOnCells() going UP-COAST, nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << " leaving loop at start of timestep (" << nDistSeawardFromNewCoast << " / " << nParProfLen << ") because enough erosion for polygon, dSedToErodeOnThisPolygon = " << dSedToErodeOnThisPolygon << endl;

               break;
            }

            // Leave the loop if we have done enough erosion for this profile
            if (dSedToErodeOnThisProfile <= 0)
            {
//                LogStream << "\tIn nDoBeachErosionOnCells() going UP-COAST, nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << " leaving loop at start of timestep (" << nDistSeawardFromNewCoast << " / " << nParProfLen << ") because enough erosion for profile, dSedToErodeOnThisProfile = " << dSedToErodeOnThisProfile << " dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << endl;

               break;
            }

            CGeom2DIPoint PtiTmp = PtiVParProfile[nParProfLen - nDistSeawardFromNewCoast - 1];
            int
               nX = PtiTmp.nGetX(),
               nY = PtiTmp.nGetY();

            // Safety check
            if (! bIsWithinValidGrid(nX, nY))
            {
//                LogStream << WARN << "08 @@@@ while constructing parallel profile for UP-COAST erosion on polygon " << nPoly << ", hit edge of grid at [" << nParProfEndX << "][" << nParProfEndY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;
               KeepWithinValidGrid(nCoastX, nCoastY, nX, nY);
               PtiTmp.SetX(nX);
               PtiTmp.SetY(nY);
            }

            // Don't do anything to intervention cells
            if (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION)
               continue;

            // Don't do cells twice
            if (! m_pRasterGrid->m_Cell[nX][nY].bBeachErosionOrDepositionThisTimestep())
            {
               // Get this cell's current elevation
               double dThisElevNow = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev();

//                LogStream << "\tnPoly = " << nPoly << " going UP-COAST, [" << nX << "][" << nY << "]  = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "} nCoastPoint = " << nCoastPoint << " nDistSeawardFromNewCoast = " << nDistSeawardFromNewCoast << " dThisElevNow = " << dThisElevNow << " Dean Elev = " << VdParProfileDeanElev[nDistSeawardFromNewCoast] << endl;

               // Subtract the two elevations
               double dElevDiff = dThisElevNow - VdParProfileDeanElev[nDistSeawardFromNewCoast];
               if (dElevDiff > SEDIMENT_ELEV_TOLERANCE)
               {
                  // The current elevation is higher than the Dean elevation, so we have potential beach erosion (i.e. not constrained by availability of unconsolidated sediment) here
                  m_ulThisTimestepNumPotentialBeachErosionCells++;

                  m_pRasterGrid->m_Cell[nX][nY].SetPotentialBeachErosion(dElevDiff);

//                   if (m_ulIteration == 16)
//                      LogStream << "\tnPoly = " << nPoly << " going UP-COAST, potential beach erosion = " << dElevDiff << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "} nCoastPoint = " << nCoastPoint << " nDistSeawardFromNewCoast = " << nDistSeawardFromNewCoast << " dThisElevNow = " << dThisElevNow << " Dean Elev = " << VdParProfileDeanElev[nDistSeawardFromNewCoast] << endl;

                  // Now get the number of the highest layer with non-zero thickness
                  int nThisLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();

                  // Safety check
                  if (nThisLayer == INT_NODATA)
                     return RTN_ERR_NO_TOP_LAYER;

                  if (nThisLayer != NO_NONZERO_THICKNESS_LAYERS)
                  {
                     // We still have at least one layer left with non-zero thickness (i.e. we are not down to basement), and the cell's current elevation is higher than the Dean equilibrium profile elevation. So do some beach erosion
                     double
                        dToErode = tMin(dElevDiff, dSedToErodeOnThisProfile, dSedToErodeOnThisPolygon),
                        dFine = 0,
                        dSand = 0,
                        dCoarse = 0;
   //                   assert(dToErode > 0);
                     ErodeBeachConstrained(nX, nY, nThisLayer, dToErode, dFine, dSand, dCoarse);

                     // Update totals for this parallel profile
                     dFineEroded += dFine;
                     dSandEroded += dSand;
                     dCoarseEroded += dCoarse;

                     // Update totals for the polygon
                     dTotFineEroded += dFine;
                     dTotSandEroded += dSand;
                     dTotCoarseEroded += dCoarse;

                     double dTmpTot = dFine + dSand + dCoarse;
                     dSedToErodeOnThisProfile -= dTmpTot;
                     dSedToErodeOnThisPolygon -= dTmpTot;

                     // Update this-timestep totals
                     m_ulThisTimestepNumActualBeachErosionCells++;
                     m_dThisTimestepActualBeachErosionFine += dFine;
                     m_dThisTimestepActualBeachErosionSand += dSand;
                     m_dThisTimestepActualBeachErosionCoarse += dCoarse;

//                      LogStream << "\tIn nDoBeachErosionOnCells(), nPoly = " << nPoly << " going UP-COAST, actual beach erosion = " << dTmpTot << " at [" << nX << "][" << nY << "]  = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "} nCoastPoint = " << nCoastPoint << " nDistSeawardFromNewCoast = " << nDistSeawardFromNewCoast << endl;
                  }
               }
               else if ((dElevDiff < -SEDIMENT_ELEV_TOLERANCE) && ((m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea()) || (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_DRIFT)))
               {
                  // The current elevation is below the Dean elevation, so we can have beach deposition here provided that we have some previously-eroded unconsolidated sediment to deposit
                  double dTotSandAndCoarseEroded = dTotSandEroded + dTotCoarseEroded;
                  if (dTotSandAndCoarseEroded > SEDIMENT_ELEV_TOLERANCE)
                  {
                     // Assume that the sediment size fractions which are deposited are in the same ratios as the previously-eroded sediment
                     double
                        dTotToDeposit = tMin(-dElevDiff, dTotSandAndCoarseEroded),
                        dSandToDeposit = dTotToDeposit * dTotSandEroded / dTotSandAndCoarseEroded,
                        dCoarseToDeposit = dTotToDeposit * dTotCoarseEroded / dTotSandAndCoarseEroded;

                     int nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopLayerAboveBasement();

                     // Safety check
                     if (nTopLayer == INT_NODATA)
                        return RTN_ERR_NO_TOP_LAYER;

                     if (dSandToDeposit > SEDIMENT_ELEV_TOLERANCE)
                     {
                        dSandToDeposit = tMin(dSandToDeposit, dTotSandEroded);

                        double dSandNow = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetSand();
                        m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetSand(dSandNow + dSandToDeposit);

                        // Set the changed-this-timestep switch
                        m_bUnconsChangedThisTimestep[nTopLayer] = true;

                        dTotSandEroded -= dSandToDeposit;

                        dSedToErodeOnThisProfile += dSandToDeposit;
                        dSedToErodeOnThisPolygon += dSandToDeposit;
                     }

                     if (dCoarseToDeposit > SEDIMENT_ELEV_TOLERANCE)
                     {
                        dCoarseToDeposit = tMin(dCoarseToDeposit, dTotCoarseEroded);

                        double dCoarseNow = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetCoarse();
                        m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetCoarse(dCoarseNow + dCoarseToDeposit);

                        // Set the changed-this-timestep switch
                        m_bUnconsChangedThisTimestep[nTopLayer] = true;

                        dTotCoarseEroded -= dCoarseToDeposit;

                        dSedToErodeOnThisProfile += dCoarseToDeposit;
                        dSedToErodeOnThisPolygon += dCoarseToDeposit;
                     }

                     // Now update the cell's layer elevations
                     m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();

                     // Update the cell's sea depth
                     m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();

                     // Update the cell's beach deposition, and total beach deposition, values
                     m_pRasterGrid->m_Cell[nX][nY].IncrBeachDeposition(dSandToDeposit + dCoarseToDeposit);

                     // And set the landform category
                     m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetLFSubCategory(LF_SUBCAT_DRIFT_BEACH);

                     // Update this-timestep totals
                     m_ulThisTimestepNumBeachDepositionCells++;
                     m_dThisTimestepBeachDepositionSand += dSandToDeposit;
                     m_dThisTimestepBeachDepositionCoarse += dCoarseToDeposit;

//                      LogStream << "\tIn nDoBeachErosionOnCells(), nPoly = " << nPoly << " going UP-COAST, beach deposition = " << dTotToDeposit << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "} nCoastPoint = " << nCoastPoint << " nDistSeawardFromNewCoast = " << nDistSeawardFromNewCoast << endl;
                  }
               }
            }
         }

//          LogStream << "\tIn nDoBeachErosionOnCells() going UP-COAST, nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << " dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << " dSedToErodeOnThisProfile = " <<  dAllSedimentTargetPerProfile << " dSedToErodeOnThisPolygon = " << dSedToErodeOnThisPolygon << " dTotTargetToErode = " << dTotTargetToErode << endl;
      }
   }

   // How much have we been able to erode?
//    LogStream << "\tIn nDoBeachErosionOnCells() going UP-COAST, nPoly = " << nPoly << " dTotFineEroded = " << dTotFineEroded << " dTotSandEroded = " << dTotSandEroded << " dTotCoarseEroded = " << dTotCoarseEroded << endl;

   double dActualErosion = dTotFineEroded + dTotSandEroded + dTotCoarseEroded;

   // Is the depth eroded within TOLERANCE of the target depth-equivalent?
   if (bFPIsEqual(dActualErosion, dTotTargetToErode, TOLERANCE))
      LogStream << "\tPolygon " << nPoly << " actual beach erosion approx equal to target beach erosion: actual = " << dActualErosion << " target = " << dTotTargetToErode << endl;
   else
   {
      LogStream << ERR << "on polygon " << nPoly << " actual beach erosion LESS THAN target beach erosion: actual = " << dActualErosion << " target = " << dTotTargetToErode << endl;

      m_dThisTimestepMassBalanceErosionError += (dTotTargetToErode - dActualErosion);
   }

   return RTN_OK;
}


/*===============================================================================================================================

 Erodes the unconsolidated beach sediment on a single cell, returns the depth-equivalents of fine, sand and coarse sediment removed

===============================================================================================================================*/
void CSimulation::ErodeBeachConstrained(int const nX, int const nY, int const nThisLayer, double const dMaxToErode, double& dFine, double& dSand, double& dCoarse)
{
   // Find out how much unconsolidated sediment we have available on this cell
   double
      dExistingAvailableFine = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->dGetFine(),
      dExistingAvailableSand = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->dGetSand(),
      dExistingAvailableCoarse = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->dGetCoarse();

   // Is there any unconsolidated sediment on this cell?
   if ((dExistingAvailableFine + dExistingAvailableSand + dExistingAvailableCoarse) < SEDIMENT_ELEV_TOLERANCE)
      return;

   // We have some unconsolidated sediment, so partition the total lowering for this cell between the three size fractions: do this by relative erodibility
   int
      nFineWeight = (dExistingAvailableFine > 0 ? 1 : 0),
      nSandWeight = (dExistingAvailableSand > 0 ? 1 : 0),
      nCoarseWeight = (dExistingAvailableCoarse > 0 ? 1 : 0);

   double
      dTotErodibility = (nFineWeight * m_dFineErodibilityNormalized) + (nSandWeight * m_dSandErodibilityNormalized) + (nCoarseWeight * m_dCoarseErodibilityNormalized),
      dTotActualErosion = 0;

   if (nFineWeight)
   {
      // Erode some fine-sized consolidated sediment
      double dFineLowering = (m_dFineErodibilityNormalized * dMaxToErode) / dTotErodibility;

      // Make sure we don't get -ve amounts left on the cell
      dFine = tMin(dExistingAvailableFine, dFineLowering);
      double dRemaining = dExistingAvailableFine - dFine;

      dTotActualErosion += dFine;

      // Set the value for this layer
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->SetFine(dRemaining);

      // And set the changed-this-timestep switch
      m_bUnconsChangedThisTimestep[nThisLayer] = true;
   }

   if (nSandWeight)
   {
      // Erode some sand-sized consolidated sediment
      double dSandLowering = (m_dSandErodibilityNormalized * dMaxToErode) / dTotErodibility;

      // Make sure we don't get -ve amounts left on the source cell
      dSand = tMin(dExistingAvailableSand, dSandLowering);
      double dRemaining = dExistingAvailableSand - dSand;

      dTotActualErosion += dSand;

      // Set the value for this layer
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->SetSand(dRemaining);

      // Set the changed-this-timestep switch
      m_bUnconsChangedThisTimestep[nThisLayer] = true;
   }

   if (nCoarseWeight)
   {
      // Erode some coarse-sized consolidated sediment
      double dCoarseLowering = (m_dCoarseErodibilityNormalized * dMaxToErode) / dTotErodibility;

      // Make sure we don't get -ve amounts left on the source cell
      dCoarse = tMin(dExistingAvailableCoarse, dCoarseLowering);
      double dRemaining = dExistingAvailableCoarse - dCoarse;

      dTotActualErosion += dCoarse;

      // Set the value for this layer
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->SetCoarse(dRemaining);

      // Set the changed-this-timestep switch
      m_bUnconsChangedThisTimestep[nThisLayer] = true;
   }

   // Set the actual erosion value for this cell
   m_pRasterGrid->m_Cell[nX][nY].SetActualBeachErosion(dTotActualErosion);

   if (dTotActualErosion > 0)
   {
      // Recalculate the elevation of every layer
      m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();

      // And update the cell's sea depth
      m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();
   }

//    LogStream << "\t[" << nX << "][" << nY << "] dPotentialErosion = " << dPotentialErosion << " dTotActualErosion = " << dTotActualErosion << " dFine = " << dFine << " dSand = " << dSand << " dCoarse = " << dCoarse << endl;
}


/*===============================================================================================================================

 Deposits unconsolidated beach sediment on the cells within a polygon

 As with the estimation of actual erosion, this is done by working down the coastline and constructing profiles which are parallel to the up-coast polygon boundary; then reversing direction and going up-coast, constructing profiles parallel to the down-coast boundary. Then iteratively fit a Dean equilibrium profile until the normal's share of the change in total depth of unconsolidated sediment is accommodated under the revised profile. For deposition, this adds to the beach volume

===============================================================================================================================*/
int CSimulation::nDoBeachDepositionOnCells(int const nCoast, int const nPoly, double const dTotTargetToDeposit)
{
   // Don't bother with tiny amounts
   if (dTotTargetToDeposit < SEDIMENT_ELEV_TOLERANCE)
      return RTN_OK;

   CGeomCoastPolygon* pPolygon = m_VCoast[nCoast].pGetPolygon(nPoly);

   // Get the grid cell co-ords of this polygon's up-coast and down-coast profiles
   int
      nUpCoastProfile = pPolygon->nGetUpCoastProfile(),
      nDownCoastProfile = pPolygon->nGetDownCoastProfile();

   CGeomProfile* pUpCoastProfile = m_VCoast[nCoast].pGetProfile(nUpCoastProfile);
   CGeomProfile* pDownCoastProfile = m_VCoast[nCoast].pGetProfile(nDownCoastProfile);

   // We are using only part of each profile, seaward as far as the depth of closure. First find the seaward end point of the up-coast part-profile
//    CGeom2DIPoint PtiUpCoastPartProfileSeawardEnd;
//    int nIndex =  pUpCoastProfile->nGetCellGivenDepth(m_pRasterGrid, m_dDepthOfClosure, &PtiUpCoastPartProfileSeawardEnd);
   int nIndex =  pUpCoastProfile->nGetCellGivenDepth(m_pRasterGrid, m_dDepthOfClosure);
   if (nIndex == INT_NODATA)
   {
      LogStream << m_ulIteration << ": " << ERR << "in nDoBeachDepositionOnCells() for polygon " << nPoly << ", could not find the seaward end point of the up-coast profile (" << nUpCoastProfile << ") for depth of closure = " << m_dDepthOfClosure << ". Lengthen the coastline normals." << endl;

      return RTN_ERR_BAD_BEACH_EROSION_PROFILE;
   }

   // The part-profile length is one greater than nIndex, since pPtiGetCellGivenDepth() returns the index of the cell at depth of closure. This will be the number of cells in the Dean profile portion of every parallel profile
   int nUpCoastDeanLen = nIndex + 1;

//    assert(bIsWithinValidGrid(&PtiUpCoastPartProfileSeawardEnd));

   // Get the distance between the start and end of the part-profile (the Dean length), in external CRS units
//    CGeom2DPoint
//       PtUpCoastProfileStart = *pUpCoastProfile->pPtGetPointInProfile(0),
//       PtUpCoastProfileEnd = PtGridCentroidToExt(&PtiUpCoastPartProfileSeawardEnd);

//   double dUpCoastDeanLen = dGetDistanceBetween(&PtUpCoastProfileStart, &PtUpCoastProfileEnd);

   int
      nUpCoastProfileCoastPoint = pUpCoastProfile->nGetNumCoastPoint(),
      nDownCoastProfileCoastPoint = pDownCoastProfile->nGetNumCoastPoint(),
      nXUpCoastProfileExistingCoastPoint = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nUpCoastProfileCoastPoint)->nGetX(),
      nYUpCoastProfileExistingCoastPoint = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nUpCoastProfileCoastPoint)->nGetY(),
      nCoastSegLen;

   // Store the coast point numbers for this polygon so that we can shuffle them
   vector<int> nVCoastPoint;
   if (nDownCoastProfileCoastPoint == m_VCoast[nCoast].nGetCoastlineSize()-1)
   {
      // This is the final down-coast polygon, so also include the down-coast polygon boundary
      nCoastSegLen = nDownCoastProfileCoastPoint - nUpCoastProfileCoastPoint + 1;
      for (int nCoastPoint = nUpCoastProfileCoastPoint; nCoastPoint <= nDownCoastProfileCoastPoint; nCoastPoint++)
         nVCoastPoint.push_back(nCoastPoint);
   }
   else
   {
      // This is not the final down-coast polygon, so do not include the polygon's down-coast boundary
      nCoastSegLen = nDownCoastProfileCoastPoint - nUpCoastProfileCoastPoint;
      for (int nCoastPoint = nUpCoastProfileCoastPoint; nCoastPoint < nDownCoastProfileCoastPoint; nCoastPoint++)
         nVCoastPoint.push_back(nCoastPoint);
   }

   // Shuffle the coast points, this is necessary so that leaving the loop does not create sequence-related artefacts
   Rand1Shuffle(&(nVCoastPoint.at(0)), nCoastSegLen);

   // Get the volume of sediment which is to be deposited on the polygon and on each parallel profile
   double
      dSandToDepositOnPoly = pPolygon->dGetDeltaActualUnconsSand(),
      dCoarseToDepositOnPoly = pPolygon->dGetDeltaActualUnconsCoarse(),
      dAllSedimentTargetPerProfile = dTotTargetToDeposit / nCoastSegLen,
      dSandTargetPerProfile = dSandToDepositOnPoly / nCoastSegLen,
      dCoarseTargetPerProfile = dCoarseToDepositOnPoly / nCoastSegLen,
      dTotSandDeposited = 0,              // Grand totals for the whole polygon
      dTotCoarseDeposited = 0;

   // Now traverse the polygon's existing coastline in a random (but broadly DOWN-COAST i.e. increasing coast point indices) sequence, fitting a Dean profile at each coast point
   // DOWN-COAST (increasing coastpoint indices) ================================================================================
   for (int n = 0; n < nCoastSegLen; n++)
   {
      // Pick a random coast point
      int nCoastPoint = nVCoastPoint[n];
      CGeom2DIPoint PtiCoastPoint = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint);
      int
         nCoastX = PtiCoastPoint.nGetX(),
         nCoastY = PtiCoastPoint.nGetY();

      // Calculate the x-y offset between this coast point, and the coast point of the up-coast normal
      int
         nXOffset = nCoastX - nXUpCoastProfileExistingCoastPoint,
         nYOffset = nCoastY - nYUpCoastProfileExistingCoastPoint,
         nSeawardOffset = -1;
//          nParProfLen;
      vector<CGeom2DIPoint> PtiVParProfile;
      vector<double> VdParProfileDeanElev;

      // OK, loop until we can deposit sufficient unconsolidated sediment
      while (true)
      {
         // Move seaward by one cell
         nSeawardOffset++;

         // And lengthen the parallel profile
         int nParProfLen = nUpCoastDeanLen + nSeawardOffset;

         if (nParProfLen > (pUpCoastProfile->nGetNumCellsInProfile()))
         {
            // We've reached the seaward end of the up-coast profile, and still cannot deposit sufficient sediment. Need to quit, since mass balance will not be preserved (TODO find a way round this)
            LogStream << m_ulIteration << ": " << WARN << "reached seaward end of up-coast profile during DOWN-COAST deposition of unconsolidated sediment for coast " << nCoast << " polygon " << nPoly << " (nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << ")" << endl;

            break;
         }

         // Get the x-y coords of a profile starting from this coast point and parallel to the up-coast polygon boundary profile (these are in natural sequence, like the boundary part-profile)
         PtiVParProfile.resize(0);
         for (int n = 0; n < nParProfLen; n++)
         {
            CGeom2DIPoint PtiProf = *pUpCoastProfile->pPtiGetCellInProfile(n);
            CGeom2DIPoint PtiTmp(PtiProf.nGetX() + nXOffset, PtiProf.nGetY() + nYOffset);
            PtiVParProfile.push_back(PtiTmp);
         }

         // Get the existing elevation of the seaward end of the parallel profile
         int
            nSeaEndX = PtiVParProfile.back().nGetX(),
            nSeaEndY = PtiVParProfile.back().nGetY();

         // Safety check
         if (! bIsWithinValidGrid(nSeaEndX, nSeaEndY))
         {
            LogStream << WARN << "09 @@@@ while doing DOWN-COAST deposition on polygon " << nPoly << ", hit edge of grid at [" << nSeaEndX << "][" << nSeaEndY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

            KeepWithinValidGrid(nCoastX, nCoastY, nSeaEndX, nSeaEndY);
            PtiVParProfile.back().SetX(nSeaEndX);
            PtiVParProfile.back().SetY(nSeaEndY);
         }

         double dParProfEndElev = m_pRasterGrid->m_Cell[nSeaEndX][nSeaEndY].dGetSedimentTopElev();

         // Set the start elevation for the Dean profile just a bit above SWL for this timestep, so that it is a Bruun profile
         double dParProfStartElev = m_dThisTimestepSWL + m_dDeanProfileStartAboveSWL;

         // Calculate the total length of the parallel profile, including any seaward offset
         double dParProfLen = dGetDistanceBetween(&PtiVParProfile.front(), &PtiVParProfile.back());

         // Now calculate the length of the Dean profile-only part i.e. without any seaward offset. The approach used here is approximate but probably OK
         double dParProfDeanLen = dParProfLen - (nSeawardOffset * m_dCellSide);

         // Solve for dA so that the existing elevations at the end of the parallel profile, and at the end of a Dean equilibrium profile on that part-normal, are the same
         double dParProfA = (dParProfStartElev - dParProfEndElev) /  pow(dParProfDeanLen, DEAN_POWER);

         nParProfLen = PtiVParProfile.size();
         VdParProfileDeanElev.resize(nParProfLen, 0);

//          for (int m = 0; m < static_cast<int>(PtiVParProfile.size()); m++)
//             LogStream << "[" << PtiVParProfile[m].nGetX() << "][" << PtiVParProfile[m].nGetY() << "] ";
//          LogStream << endl;

         double dInc = dParProfDeanLen / (nParProfLen - nSeawardOffset - 2);

         // The elevation of the coast point in the Dean profile is the same as the elevation of the current coast point TODO is this correct? Should it be dParProfStartElev?
         double dCoastElev = m_pRasterGrid->m_Cell[nCoastX][nCoastY].dGetSedimentTopElev();

         // For this depositing parallel profile, calculate the Dean equilibrium profile of the unconsolidated sediment h(y) = A * y^(2/3) where h(y) is the distance below the highest point in the profile at a distance y from the landward start of the profile
         CalcDeanProfile(&VdParProfileDeanElev, dInc, dParProfStartElev, dParProfA, true, nSeawardOffset, dCoastElev);

         double dParProfTotDiff = 0;
         for (int n = 0; n < nParProfLen; n++)
         {
            CGeom2DIPoint PtiTmp = PtiVParProfile[n];
            int
               nX = PtiTmp.nGetX(),
               nY = PtiTmp.nGetY();

            // Safety check
            if (! bIsWithinValidGrid(nX, nY))
            {
//                LogStream << WARN << "10 @@@@ while doing DOWN-COAST deposition on polygon " << nPoly << ", hit edge of grid at [" << nX << "][" << nY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

               KeepWithinValidGrid(nCoastX, nCoastY, nX, nY);
               PtiTmp.SetX(nX);
               PtiTmp.SetY(nY);
            }

            // Don't do anything to intervention cells
            if (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION)
               continue;

            // Don't do cells twice
            if (! m_pRasterGrid->m_Cell[nX][nY].bBeachErosionOrDepositionThisTimestep())
            {
               double
                  dTmpElev = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev(),
                  dDiff = VdParProfileDeanElev[n] - dTmpElev;

               dParProfTotDiff += dDiff;
            }
         }

//          // DEBUG STUFF -----------------------------------------------------
//          LogStream << endl << "\tFor polygon " << nPoly << " doing DOWN-COAST deposition, nSeawardOffset = " << nSeawardOffset << ", parallel profile from [" << PtiVParProfile[0].nGetX() << "][" << PtiVParProfile[0].nGetY() << "] to [" << PtiVParProfile.back().nGetX() << "][" << PtiVParProfile.back().nGetY() << "], nUpCoastDeanLen = " << nUpCoastDeanLen << " dUpCoastDeanLen = " << dUpCoastDeanLen << " nParProfLen = " << nParProfLen << " dParProfDeanLen = " << dParProfDeanLen << " dInc = " << dInc << " dParProfStartElev = " << dParProfStartElev << " dParProfEndElev = " << dParProfEndElev << " dParProfA = " << dParProfA << endl;
//
//          LogStream << "\tExisting profile for deposition = ";
//          for (int n = 0; n < nParProfLen; n++)
//          {
//             CGeom2DIPoint PtiTmp = PtiVParProfile[n];
//             int
//                nX = PtiTmp.nGetX(),
//                nY = PtiTmp.nGetY();
//
//             // Safety check
//             if (! bIsWithinValidGrid(nX, nY))
//                KeepWithinValidGrid(nX, nY);
//
//             LogStream << m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev() << " ";
//          }
//          LogStream << endl;
//          LogStream << "\tParallel Dean equilibrium profile for deposition = ";
//          for (int n = 0; n < nParProfLen; n++)
//          {
//             LogStream << VdParProfileDeanElev[n] << " ";
//          }
//          LogStream << endl;
//
//          LogStream << "\tnCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " difference = ";
//          for (int n = 0; n < nParProfLen; n++)
//          {
//             CGeom2DIPoint PtiTmp = PtiVParProfile[n];
//             int
//                nX = PtiTmp.nGetX(),
//                nY = PtiTmp.nGetY();
//
//             // Safety check
//             if (! bIsWithinValidGrid(nX, nY))
//                KeepWithinValidGrid(nX, nY);
//
//             double
//                dTmpElev = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev(),
//                dDiff = VdParProfileDeanElev[n] - dTmpElev;
//
//             LogStream << dDiff << " ";
//          }
//          LogStream << endl;
//          // END DEBUG STUFF -----------------------------------------------------


         // So will we be able to deposit as much as is needed?
         if (dParProfTotDiff >= dAllSedimentTargetPerProfile)
         {
//             LogStream << "\tDOWN-COAST nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " dParProfTotDiff = " << dParProfTotDiff << " dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << " ENOUGH TO BE DEPOSITED" << endl;

            break;
         }
      }

//      assert(dParProfTotDiff > 0);

      // OK, this value of nSeawardOffset gives us enough deposition. So start depositing on the parallel profile from the landward end
      double
         dSandDeposited = 0,        // Totals for this parallel profile
         dCoarseDeposited = 0,
         dSandToDepositOnProf = dSandTargetPerProfile,
         dCoarseToDepositOnProf = dCoarseTargetPerProfile,
         dSandRatio = dSandToDepositOnProf / (dSandToDepositOnProf + dCoarseToDepositOnProf),
         dCoarseRatio = 1 - dSandRatio;
//       assert(dRatio >= 0);

      for (unsigned int nSeawardFromCoast = 0; nSeawardFromCoast < PtiVParProfile.size(); nSeawardFromCoast++)
      {
         // Don't bother with tiny amounts
         if (dSandToDepositOnPoly < SEDIMENT_ELEV_TOLERANCE)
            dSandToDepositOnPoly = 0;
         if (dSandToDepositOnProf < SEDIMENT_ELEV_TOLERANCE)
            dSandToDepositOnProf = 0;
         if (dCoarseToDepositOnPoly < SEDIMENT_ELEV_TOLERANCE)
            dCoarseToDepositOnPoly = 0;
         if (dCoarseToDepositOnProf < SEDIMENT_ELEV_TOLERANCE)
            dCoarseToDepositOnProf = 0;

         // Leave the loop if we have deposited enough for this polygon
         if ((dSandToDepositOnPoly <= 0) && (dCoarseToDepositOnPoly <= 0))
         {
//             LogStream << "\tIn nDoWithinPolygonBeachRedistribution() going DOWN-COAST, nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " leaving loop at start of timestep (" << nSeawardFromCoast << " / " << nParProfLen << ")  because enough deposition for polygon, dSandToDepositOnPoly = " << dSandToDepositOnPoly << " dCoarseToDepositOnPoly = " << dCoarseToDepositOnPoly << endl;

            break;
         }

         // Leave the loop if we have done enough deposition for this profile
         if ((dSandToDepositOnProf <= 0) && (dCoarseToDepositOnProf <= 0))
         {
//             LogStream << "\tIn nDoWithinPolygonBeachRedistribution() going DOWN-COAST, nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " leaving loop at start of timestep (" << nSeawardFromCoast << " / " << nParProfLen << ") because enough deposition for profile, dSandToDepositOnProf = " << dSandToDepositOnProf << " dCoarseToDepositOnProf = " << dCoarseToDepositOnProf << " dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << endl;
            break;
         }

//          assert(nSeawardFromCoast < PtiVParProfile.size());
         CGeom2DIPoint PtiTmp = PtiVParProfile[nSeawardFromCoast];
         int
            nX = PtiTmp.nGetX(),
            nY = PtiTmp.nGetY();

         // Safety check
         if (! bIsWithinValidGrid(nX, nY))
         {
//             LogStream << WARN << "11 @@@@ while doing DOWN-COAST deposition on polygon " << nPoly << ", hit edge of grid at [" << nX << "][" << nY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

            KeepWithinValidGrid(nCoastX, nCoastY, nX, nY);
            PtiTmp.SetX(nX);
            PtiTmp.SetY(nY);
         }

         // Don't do anything to intervention cells
         if (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION)
            continue;

         // Don't do cells twice
         if (! m_pRasterGrid->m_Cell[nX][nY].bBeachErosionOrDepositionThisTimestep())
         {
            // Get this cell's current elevation
            double dThisElevNow = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev();

//             LogStream << "\tnPoly = " << nPoly << " going DOWN-COAST, [" << nX << "][" << nY << "] nCoastPoint = " << nCoastPoint << " nSeawardFromCoast = " << nSeawardFromCoast << " dThisElevNow = " << dThisElevNow << " Dean Elev = " << VdParProfileDeanElev[nSeawardFromCoast] << endl;

            // Subtract the two elevations
//             assert(nSeawardFromCoast < VdParProfileDeanElev.size());
            double dElevDiff = VdParProfileDeanElev[nSeawardFromCoast] - dThisElevNow;
            if (dElevDiff > SEDIMENT_ELEV_TOLERANCE)
            {
               bool
                  bSandDeposited = false,
                  bCoarseDeposited = false;
               double
                  dSandToDeposit = 0,
                  dCoarseToDeposit = 0;

               // The current elevation is below the Dean elevation, so we have can have beach deposition here
               if (dSandToDepositOnProf > SEDIMENT_ELEV_TOLERANCE)
               {
                  dSandToDeposit = tMin(dElevDiff * dSandRatio, dSandToDepositOnProf, dSandToDepositOnPoly);

                  int nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopLayerAboveBasement();

                  // Safety check
                  if (nTopLayer == INT_NODATA)
                     return RTN_ERR_NO_TOP_LAYER;

                  if (dSandToDeposit > SEDIMENT_ELEV_TOLERANCE)
                  {
                     bSandDeposited = true;

                     double dSandNow = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetSand();
                     m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetSand(dSandNow + dSandToDeposit);

                     // Set the changed-this-timestep switch
                     m_bUnconsChangedThisTimestep[nTopLayer] = true;

                     dSandDeposited += dSandToDeposit;
                     dTotSandDeposited += dSandToDeposit;

                     dSandToDepositOnProf -= dSandToDeposit;
                     dSandToDepositOnPoly -= dSandToDeposit;

                     // Update the cell's beach deposition, and total beach deposition, values
                     m_pRasterGrid->m_Cell[nX][nY].IncrBeachDeposition(dSandToDeposit);

//                      LogStream << "\tIn nDoWithinPolygonBeachRedistribution(), nPoly = " << nPoly << " going DOWN-COAST, sand deposition = " << dSandToDeposit << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "}  nCoastPoint = " << nCoastPoint << " nSeawardFromCoast = " << nSeawardFromCoast << endl;
                  }
               }

               if (dCoarseToDepositOnProf > SEDIMENT_ELEV_TOLERANCE)
               {
                  dCoarseToDeposit = tMin(dElevDiff * dCoarseRatio, dCoarseToDepositOnProf, dCoarseToDepositOnPoly);

                  int nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopLayerAboveBasement();

                  // Safety check
                  if (nTopLayer == INT_NODATA)
                     return RTN_ERR_NO_TOP_LAYER;

                  if (dCoarseToDeposit > SEDIMENT_ELEV_TOLERANCE)
                  {
                     bCoarseDeposited = true;

                     double dCoarseNow = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetCoarse();
                     m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetCoarse(dCoarseNow + dCoarseToDeposit);

                     // Set the changed-this-timestep switch
                     m_bUnconsChangedThisTimestep[nTopLayer] = true;

                     dCoarseDeposited += dCoarseToDeposit;
                     dTotCoarseDeposited += dCoarseToDeposit;

                     dCoarseToDepositOnProf -= dCoarseToDeposit;
                     dCoarseToDepositOnPoly -= dCoarseToDeposit;

                     // Update the cell's beach deposition, and total beach deposition, values
                     m_pRasterGrid->m_Cell[nX][nY].IncrBeachDeposition(dCoarseToDeposit);

   //                   LogStream << "\tIn nDoWithinPolygonBeachRedistribution(), nPoly = " << nPoly << " going DOWN-COAST, coarse deposition = " << dCoarseToDeposit << " at [" << nX << "][" << nY << "] nCoastPoint = " << nCoastPoint << " nSeawardFromCoast = " << nSeawardFromCoast << endl;
                  }
               }

               if (bSandDeposited || bCoarseDeposited)
               {
                  // Update the cell's layer elevations
                  m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();

                  // Update the cell's sea depth
                  m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();

                  // And set the landform category
                  m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetLFSubCategory(LF_SUBCAT_DRIFT_BEACH);

                  // Update this-timestep totals
                  m_ulThisTimestepNumBeachDepositionCells++;
                  m_dThisTimestepBeachDepositionSand += dSandToDeposit;
                  m_dThisTimestepBeachDepositionCoarse += dCoarseToDeposit;

               }
            }
            else if ((dElevDiff < -SEDIMENT_ELEV_TOLERANCE) && ((m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea()) || (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_DRIFT)))
            {
               // The current elevation is higher than the Dean elevation, so we have potential beach erosion (i.e. not constrained by availability of unconsolidated sediment) here
               m_ulThisTimestepNumPotentialBeachErosionCells++;

               m_pRasterGrid->m_Cell[nX][nY].SetPotentialBeachErosion(-dElevDiff);

//                   LogStream << "\tIn nDoWithinPolygonBeachRedistribution(), nPoly = " << nPoly << " going UP-COAST, potential beach erosion = " << -dElevDiff << " at [" << nX << "][" << nY << "] nCoastPoint = " << nCoastPoint << " nSeawardFromCoast = " << nSeawardFromCoast << " dThisElevNow = " << dThisElevNow << " Dean Elev = " << VdParProfileDeanElev[nSeawardFromCoast] << endl;

               // Now get the number of the highest layer with non-zero thickness
               int nThisLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();

               // Safety check
               if (nThisLayer == INT_NODATA)
                  return RTN_ERR_NO_TOP_LAYER;

               if (nThisLayer != NO_NONZERO_THICKNESS_LAYERS)
               {
                  // We still have at least one layer left with non-zero thickness (i.e. we are not down to basement), and the cell's current elevation is higher than the Dean equilibrium profile elevation. So do some beach erosion
                  double
                     dFine = 0,
                     dSand = 0,
                     dCoarse = 0;

                  ErodeBeachConstrained(nX, nY, nThisLayer, -dElevDiff, dFine, dSand, dCoarse);

                  // Totals for this parallel profile
                  dSandDeposited -= dSand;
                  dCoarseDeposited -= dCoarse;

//                   assert(dSandDeposited >= 0);
//                   assert(dCoarseDeposited >= 0);

                  // Update totals for the polygon (note that these may become slightly -ve if erosion occurs early during this routine, but this is not a problem since the values will eventually become +ve again
                  dTotSandDeposited -= dSand;
                  dTotCoarseDeposited -= dCoarse;

                  dSandToDepositOnProf += dSand;
                  dSandToDepositOnPoly += dSand;
                  dCoarseToDepositOnProf += dCoarse;
                  dCoarseToDepositOnPoly += dCoarse;

                  // Update this-timestep totals
                  m_ulThisTimestepNumActualBeachErosionCells++;
                  m_dThisTimestepActualBeachErosionFine += dFine;
                  m_dThisTimestepActualBeachErosionSand += dSand;
                  m_dThisTimestepActualBeachErosionCoarse += dCoarse;

//                   LogStream << "\tPolygon net deposition " << nPoly << " going UP-COAST, actual beach erosion = " << dFine + dSand + dCoarse << " at [" << nX << "][" << nY << "] nCoastPoint = " << nCoastPoint << " nSeawardFromCoast = " << nSeawardFromCoast << endl;
               }
            }
         }
      }

//       LogStream << "\tIn nDoWithinPolygonBeachRedistribution() going DOWN-COAST, nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << " dSandToDepositOnProf = " << dSandToDepositOnProf << " dCoarseToDepositOnProf = " << dCoarseToDepositOnProf << " dSandToDepositOnPoly = " << dSandToDepositOnPoly << " dCoarseToDepositOnPoly = " << dCoarseToDepositOnPoly << endl;
   }

   // Have we deposited the full amount?
   double dTotDeposited = dTotSandDeposited + dTotCoarseDeposited;
   if ((! bFPIsEqual(dTotDeposited, dTotTargetToDeposit, TOLERANCE)) && (dTotDeposited < dTotTargetToDeposit))
   {
      // No, so do the same in an UP-COAST (i.e. decreasing coastpoint indices) direction
      // UP-COAST (decreasing coastpoint indices) ==================================================================================

      // Start by finding the seaward end point of the down-coast part-profile, as before we are using only part of each profile, seaward as far as the depth of closure
//       CGeom2DIPoint PtiDownCoastPartProfileSeawardEnd;
//       int nIndex = pDownCoastProfile->nGetCellGivenDepth(m_pRasterGrid, m_dDepthOfClosure, &PtiDownCoastPartProfileSeawardEnd);
      int nIndex = pDownCoastProfile->nGetCellGivenDepth(m_pRasterGrid, m_dDepthOfClosure);
      if (nIndex == INT_NODATA)
      {
         LogStream << m_ulIteration << ": " << ERR << "in nDoBeachDepositionOnCells() for polygon " << nPoly << ", could not find the seaward end point of the down-coast profile (" << nUpCoastProfile << ") for depth of closure = " << m_dDepthOfClosure << ". Lengthen the coastline normals." << endl;
         return RTN_ERR_BAD_BEACH_EROSION_PROFILE;
      }

      // The part-profile length is one greater than nIndex, since pPtiGetCellGivenDepth() returns the index of the cell at depth of closure. This will be the number of cells in the Dean profile portion of every parallel profile
      int nDownCoastDeanLen = nIndex + 1;

//       assert(bIsWithinValidGrid(&PtiDownCoastPartProfileSeawardEnd));

      // Get the distance between the start and end of the part-profile (the Dean length), in external CRS units
//       CGeom2DPoint
//          PtDownCoastProfileStart = *pDownCoastProfile->pPtGetPointInProfile(0),
//          PtDownCoastProfileEnd = PtGridCentroidToExt(&PtiDownCoastPartProfileSeawardEnd);

//       double dDownCoastDeanLen = dGetDistanceBetween(&PtDownCoastProfileStart, &PtDownCoastProfileEnd);

      int
         nXDownCoastProfileExistingCoastPoint = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nDownCoastProfileCoastPoint)->nGetX(),
         nYDownCoastProfileExistingCoastPoint = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nDownCoastProfileCoastPoint)->nGetY();

      int nCoastSegLen;

      // Store the coast point numbers for this polygon so that we can shuffle them
      nVCoastPoint.resize(0);
      if (nUpCoastProfileCoastPoint == 0)
      {
         // This is the final up-coast polygon, so also include the up-coast polygon boundary
         nCoastSegLen = nDownCoastProfileCoastPoint - nUpCoastProfileCoastPoint + 1;
         for (int nCoastPoint = nUpCoastProfileCoastPoint; nCoastPoint <= nDownCoastProfileCoastPoint; nCoastPoint++)
            nVCoastPoint.push_back(nCoastPoint);
      }
      else
      {
         // This is not the final down-coast polygon, so do not include the polygon's down-coast boundary
         nCoastSegLen = nDownCoastProfileCoastPoint - nUpCoastProfileCoastPoint;
         for (int nCoastPoint = nUpCoastProfileCoastPoint; nCoastPoint < nDownCoastProfileCoastPoint; nCoastPoint++)
            nVCoastPoint.push_back(nCoastPoint);
      }

      // Shuffle the coast points, this is necessary so that leaving the loop does not create sequence-related artefacts
      Rand1Shuffle(&(nVCoastPoint.at(0)), nCoastSegLen);

      // Recalc the targets for deposition per profile
      dSandTargetPerProfile = dSandToDepositOnPoly / nCoastSegLen;
      dCoarseTargetPerProfile = dCoarseToDepositOnPoly / nCoastSegLen;

      // Now traverse the polygon's existing coastline in a random (but broadly UP-COAST, i.e. decreasing coastpoint indices) sequence, fitting a Dean profile at each coast point
      for (int n = 0; n < nCoastSegLen; n++)
      {
         // Pick a random coast point
         int nCoastPoint = nVCoastPoint[n];
         CGeom2DIPoint PtiCoastPoint = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint);
         int
            nCoastX = PtiCoastPoint.nGetX(),
            nCoastY = PtiCoastPoint.nGetY();

         // Calculate the x-y offset between this coast point, and the coast point of the down-coast normal
         int
            nXOffset = nCoastX - nXDownCoastProfileExistingCoastPoint,
            nYOffset = nCoastY - nYDownCoastProfileExistingCoastPoint,
            nSeawardOffset = -1;
//             nParProfLen;
         vector<CGeom2DIPoint> PtiVParProfile;
         vector<double> VdParProfileDeanElev;

         // OK, loop until we can deposit sufficient unconsolidated sediment
         while (true)
         {
            // Move seaward by one cell
            nSeawardOffset++;

            // And lengthen the parallel profile
            int nParProfLen = nDownCoastDeanLen + nSeawardOffset;

            if (nParProfLen > (pDownCoastProfile->nGetNumCellsInProfile()))
            {
               // We've reached the seaward end of the down-coast profile, and still cannot deposit sufficient sediment. Need to quit, since mass balance will not be preserved (TODO find a way round this)
               LogStream << m_ulIteration << ": " << WARN << "reached seaward end of down-coast profile during UP-COAST deposition of unconsolidated sediment for coast " << nCoast << " polygon " << nPoly << " (nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << ")" << endl;

               break;
            }

            // Get the x-y coords of a profile starting from this coast point and parallel to the down-coast polygon boundary profile (these are in natural sequence, like the boundary part-profile)
            PtiVParProfile.resize(0);
            for (int n = 0; n < nParProfLen; n++)
            {
               CGeom2DIPoint PtiProf = *pDownCoastProfile->pPtiGetCellInProfile(n);
               CGeom2DIPoint PtiTmp(PtiProf.nGetX() + nXOffset, PtiProf.nGetY() + nYOffset);
               PtiVParProfile.push_back(PtiTmp);
            }

            // Get the existing elevation of the seaward end of the parallel profile
            int
               nSeaEndX = PtiVParProfile.back().nGetX(),
               nSeaEndY = PtiVParProfile.back().nGetY();

            // Safety check
            if (! bIsWithinValidGrid(nSeaEndX, nSeaEndY))
            {
//                LogStream << WARN << "12 @@@@ while doing UP-COAST deposition on polygon " << nPoly << ", hit edge of grid at [" << nSeaEndX << "][" << nSeaEndY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

               KeepWithinValidGrid(nCoastX, nCoastY, nSeaEndX, nSeaEndY);
               PtiVParProfile.back().SetX(nSeaEndX);
               PtiVParProfile.back().SetY(nSeaEndY);
            }

            double dParProfEndElev = m_pRasterGrid->m_Cell[nSeaEndX][nSeaEndY].dGetSedimentTopElev();

            // Set the start elevation for the Dean profile just a bit above SWL for this timestep, so that it is a Bruun profile
            double dParProfStartElev = m_dThisTimestepSWL + m_dDeanProfileStartAboveSWL;

            // Calculate the total length of the parallel profile, including any seaward offset
            double dParProfLen = dGetDistanceBetween(&PtiVParProfile.front(), &PtiVParProfile.back());

            // Now calculate the length of the Dean profile-only part i.e. without any seaward offset. The approach used here is approximate but probably OK
            double dParProfDeanLen = dParProfLen - (nSeawardOffset * m_dCellSide);

            // Solve for dA so that the existing elevations at the end of the parallel profile, and at the end of a Dean equilibrium profile on that part-normal, are the same
            double dParProfA = (dParProfStartElev - dParProfEndElev) /  pow(dParProfDeanLen, DEAN_POWER);

            nParProfLen = PtiVParProfile.size();
            VdParProfileDeanElev.resize(nParProfLen, 0);

            double dInc = dParProfDeanLen / (nParProfLen - nSeawardOffset - 2);

            // The elevation of the coast point in the Dean profile is the same as the elevation of the current coast point TODO is this correct? Should it be dParProfStartElev?
            double dCoastElev = m_pRasterGrid->m_Cell[nCoastX][nCoastY].dGetSedimentTopElev();

            // For this depositing parallel profile, calculate the Dean equilibrium profile of the unconsolidated sediment h(y) = A * y^(2/3) where h(y) is the distance below the highest point in the profile at a distance y from the landward start of the profile
            CalcDeanProfile(&VdParProfileDeanElev, dInc, dParProfStartElev, dParProfA, true, nSeawardOffset, dCoastElev);

            double dParProfTotDiff = 0;
            for (int n = 0; n < nParProfLen; n++)
            {
               CGeom2DIPoint PtiTmp = PtiVParProfile[n];
               int
                  nX = PtiTmp.nGetX(),
                  nY = PtiTmp.nGetY();

               // Safety check
               if (! bIsWithinValidGrid(nX, nY))
               {
//                   LogStream << WARN << "13 @@@@ while constructing parallel profile for UP-COAST deposition on polygon " << nPoly << ", hit edge of grid at [" << nX << "][" << nY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

                  KeepWithinValidGrid(nCoastX, nCoastY, nX, nY);
                  PtiTmp.SetX(nX);
                  PtiTmp.SetY(nY);
               }

               // Don't do anything to intervention cells
               if (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION)
                  continue;

               // Don't do cells twice
               if (! m_pRasterGrid->m_Cell[nX][nY].bBeachErosionOrDepositionThisTimestep())
               {
                  double
                     dTmpElev = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev(),
                     dDiff = VdParProfileDeanElev[n] - dTmpElev;

                  dParProfTotDiff += dDiff;
               }
            }

//             // DEBUG STUFF -----------------------------------------------------
//             LogStream << endl << "\tFor polygon " << nPoly << " doing UP-COAST deposition, nSeawardOffset = " << nSeawardOffset << ", parallel profile from [" << PtiVParProfile[0].nGetX() << "][" << PtiVParProfile[0].nGetY() << "] to [" << PtiVParProfile.back().nGetX() << "][" << PtiVParProfile.back().nGetY() << "], nDownCoastDeanLen = " << nDownCoastDeanLen << " dDownCoastDeanLen = " << dDownCoastDeanLen << " nParProfLen = " << nParProfLen << " dParProfDeanLen = " << dParProfDeanLen << " dInc = " << dInc << " dParProfStartElev = " << dParProfStartElev << " dParProfEndElev = " << dParProfEndElev << " dParProfA = " << dParProfA << endl;
//
//             LogStream << "\tExisting profile for deposition = ";
//             for (int n = 0; n < nParProfLen; n++)
//             {
//                CGeom2DIPoint PtiTmp = PtiVParProfile[n];
//                int
//                   nX = PtiTmp.nGetX(),
//                   nY = PtiTmp.nGetY();
//
//                // Safety check
//                if (! bIsWithinValidGrid(nX, nY))
//                   KeepWithinValidGrid(nX, nY);
//
//                LogStream << m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev() << " ";
//             }
//             LogStream << endl;
//             LogStream << "\tParallel Dean equilibrium profile for deposition = ";
//             for (int n = 0; n < nParProfLen; n++)
//             {
//                LogStream << VdParProfileDeanElev[n] << " ";
//             }
//             LogStream << endl;
//
//             LogStream << "\tnCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " difference = ";
//             for (int n = 0; n < nParProfLen; n++)
//             {
//                CGeom2DIPoint PtiTmp = PtiVParProfile[n];
//                int
//                   nX = PtiTmp.nGetX(),
//                   nY = PtiTmp.nGetY();
//
//                // Safety check
//                if (! bIsWithinValidGrid(nX, nY))
//                   KeepWithinValidGrid(nX, nY);
//
//                double
//                   dTmpElev = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev(),
//                   dDiff = VdParProfileDeanElev[n] - dTmpElev;
//
//                LogStream << dDiff << " ";
//             }
//             LogStream << endl;
//             // END DEBUG STUFF -----------------------------------------------------

            // So will we be able to deposit as much as is needed?
            if (dParProfTotDiff >= dAllSedimentTargetPerProfile)
            {
//                LogStream << "\tUP-COAST nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " dParProfTotDiff = " << dParProfTotDiff << " dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << " ENOUGH TO BE DEPOSITED" << endl;

               break;
            }
         }

   //      assert(dParProfTotDiff > 0);

         // OK, this value of nSeawardOffset gives us enough deposition. So start depositing on the parallel profile from the landward end
         double
//            dRatio = tMin(1.0, dAllSedimentTargetPerProfile / dParProfTotDiff),
            dSandDeposited = 0,        // Totals for this parallel profile
            dCoarseDeposited = 0,
            dSandToDepositOnProf = dSandTargetPerProfile,
            dCoarseToDepositOnProf = dCoarseTargetPerProfile,
            dSandRatio = dSandToDepositOnProf / (dSandToDepositOnProf + dCoarseToDepositOnProf),
            dCoarseRatio = 1 - dSandRatio;


         for (unsigned int nSeawardFromCoast = 0; nSeawardFromCoast < PtiVParProfile.size(); nSeawardFromCoast++)
         {
            // Don't bother with tiny amounts
            if (dSandToDepositOnPoly < SEDIMENT_ELEV_TOLERANCE)
               dSandToDepositOnPoly = 0;
            if (dSandToDepositOnProf < SEDIMENT_ELEV_TOLERANCE)
               dSandToDepositOnProf = 0;
            if (dCoarseToDepositOnPoly < SEDIMENT_ELEV_TOLERANCE)
               dCoarseToDepositOnPoly = 0;
            if (dCoarseToDepositOnProf < SEDIMENT_ELEV_TOLERANCE)
               dCoarseToDepositOnProf = 0;

            // Leave the loop if we have deposited enough for this polygon
            if ((dSandToDepositOnPoly <= 0) && (dCoarseToDepositOnPoly <= 0))
            {
//                LogStream << "\tIn nDoWithinPolygonBeachRedistribution() going UP-COAST, nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " leaving loop at start of timestep (" << nSeawardFromCoast << " / " << nParProfLen << ") because enough deposition for polygon, dSandToDepositOnPoly = " << dSandToDepositOnPoly << " dCoarseToDepositOnPoly = " << dCoarseToDepositOnPoly << endl;
               break;
            }

            // Leave the loop if we have done enough deposition for this profile
            if ((dSandToDepositOnProf <= 0) && (dCoarseToDepositOnProf <= 0))
            {
//                LogStream << "\tIn nDoWithinPolygonBeachRedistribution() going UP-COAST, nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " leaving loop at start of timestep (" << nSeawardFromCoast << " / " << nParProfLen << ") because enough deposition for profile, dSandToDepositOnProf = " << dSandToDepositOnProf << " dCoarseToDepositOnProf = " << dCoarseToDepositOnProf << " dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << endl;
               break;
            }

            CGeom2DIPoint PtiTmp = PtiVParProfile[nSeawardFromCoast];
            int
               nX = PtiTmp.nGetX(),
               nY = PtiTmp.nGetY();

            // Safety check
            if (! bIsWithinValidGrid(nX, nY))
            {
//                LogStream << WARN << "14 @@@@ while constructing parallel profile for UP-COAST deposition on polygon " << nPoly << ", hit edge of grid at [" << nX << "][" << nY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

               KeepWithinValidGrid(nCoastX, nCoastY, nX, nY);
               PtiTmp.SetX(nX);
               PtiTmp.SetY(nY);
            }

            // Don't do anything to intervention cells
            if (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION)
               continue;

            // Don't do cells twice
            if (! m_pRasterGrid->m_Cell[nX][nY].bBeachErosionOrDepositionThisTimestep())
            {
               // Get this cell's current elevation
               double dThisElevNow = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev();

      //          LogStream << "\tnPoly = " << nPoly << " going UP-COAST, [" << nX << "][" << nY << "] nCoastPoint = " << nCoastPoint << " nSeawardFromCoast = " << nSeawardFromCoast << " dThisElevNow = " << dThisElevNow << " Dean Elev = " << VdParProfileDeanElev[nSeawardFromCoast] << endl;

               // Subtract the two elevations
//                assert(nSeawardFromCoast < VdParProfileDeanElev.size());
               double dElevDiff = VdParProfileDeanElev[nSeawardFromCoast] - dThisElevNow;
               if (dElevDiff > SEDIMENT_ELEV_TOLERANCE)
               {
                  bool
                     bSandDeposited = false,
                     bCoarseDeposited = false;
                  double
                     dSandToDeposit = 0,
                     dCoarseToDeposit = 0;

                  // The current elevation is below the Dean elevation, so we have can have beach deposition here
                  if (dSandToDepositOnProf > SEDIMENT_ELEV_TOLERANCE)
                  {
                     dSandToDeposit = tMin(dElevDiff * dSandRatio, dSandToDepositOnProf, dSandToDepositOnPoly);

                     int nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopLayerAboveBasement();

                     // Safety check
                     if (nTopLayer == INT_NODATA)
                        return RTN_ERR_NO_TOP_LAYER;

                     if (dSandToDeposit > SEDIMENT_ELEV_TOLERANCE)
                     {
                        bSandDeposited = true;

                        double dSandNow = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetSand();
                        m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetSand(dSandNow + dSandToDeposit);

                        // Set the changed-this-timestep switch
                       m_bUnconsChangedThisTimestep[nTopLayer] = true;

                        dSandDeposited += dSandToDeposit;
                        dTotSandDeposited += dSandToDeposit;

                        dSandToDepositOnProf -= dSandToDeposit;
                        dSandToDepositOnPoly -= dSandToDeposit;

                        // Update the cell's beach deposition, and total beach deposition, values
                        m_pRasterGrid->m_Cell[nX][nY].IncrBeachDeposition(dSandToDeposit);

//                         LogStream << "\tIn nDoWithinPolygonBeachRedistribution() going UP-COAST, nPoly = " << nPoly << " going up-coast, sand deposition = " << dSandToDeposit << " at [" << nX << "][" << nY << "] nCoastPoint = " << nCoastPoint << " nSeawardFromCoast = " << nSeawardFromCoast << endl;
                     }
                  }

                  if (dCoarseToDepositOnProf > SEDIMENT_ELEV_TOLERANCE)
                  {
                     dCoarseToDeposit = tMin(dElevDiff * dCoarseRatio, dCoarseToDepositOnProf, dCoarseToDepositOnPoly);

                     int nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopLayerAboveBasement();

                     // Safety check
                     if (nTopLayer == INT_NODATA)
                        return RTN_ERR_NO_TOP_LAYER;

                     if (dCoarseToDeposit > SEDIMENT_ELEV_TOLERANCE)
                     {
                        bCoarseDeposited = true;

                        double dCoarseNow = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetCoarse();
                        m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetCoarse(dCoarseNow + dCoarseToDeposit);

                        // Set the changed-this-timestep switch
                        m_bUnconsChangedThisTimestep[nTopLayer] = true;

                        dCoarseDeposited += dCoarseToDeposit;
                        dTotCoarseDeposited += dCoarseToDeposit;

                        dCoarseToDepositOnProf -= dCoarseToDeposit;
                        dCoarseToDepositOnPoly -= dCoarseToDeposit;

                        // Update the cell's beach deposition, and total beach deposition, values
                        m_pRasterGrid->m_Cell[nX][nY].IncrBeachDeposition(dCoarseToDeposit);

//                         LogStream << "\tIn nDoWithinPolygonBeachRedistribution(), nPoly = " << nPoly << " going UP-COAST, coarse deposition = " << dCoarseToDeposit << " at [" << nX << "][" << nY << "] nCoastPoint = " << nCoastPoint << " nSeawardFromCoast = " << nSeawardFromCoast << endl;
                     }
                  }

                  if (bSandDeposited || bCoarseDeposited)
                  {
                     // Now update the cell's layer elevations
                     m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();

                     // Update the cell's sea depth
                     m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();

                     // And set the landform category
                     m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetLFSubCategory(LF_SUBCAT_DRIFT_BEACH);

                     // Update this-timestep totals
                     m_ulThisTimestepNumBeachDepositionCells++;
                     m_dThisTimestepBeachDepositionSand += dSandToDeposit;
                     m_dThisTimestepBeachDepositionCoarse += dCoarseToDeposit;

                  }
               }
               else if ((dElevDiff < -SEDIMENT_ELEV_TOLERANCE) && ((m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea()) || (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_DRIFT)))
               {
                  // The current elevation is higher than the Dean elevation, so we have potential beach erosion (i.e. not constrained by availability of unconsolidated sediment) here
                  m_ulThisTimestepNumPotentialBeachErosionCells++;

                  m_pRasterGrid->m_Cell[nX][nY].SetPotentialBeachErosion(-dElevDiff);

//                   LogStream << "\tIn nDoWithinPolygonBeachRedistribution(), nPoly = " << nPoly << " going UP-COAST, potential beach erosion = " << -dElevDiff << " at [" << nX << "][" << nY << "] nCoastPoint = " << nCoastPoint << " nSeawardFromCoast = " << nSeawardFromCoast << " dThisElevNow = " << dThisElevNow << " Dean Elev = " << VdParProfileDeanElev[nSeawardFromCoast] << endl;

                  // Now get the number of the highest layer with non-zero thickness
                  int nThisLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();

                  // Safety check
                  if (nThisLayer == INT_NODATA)
                     return RTN_ERR_NO_TOP_LAYER;

                  if (nThisLayer != NO_NONZERO_THICKNESS_LAYERS)
                  {
                     // We still have at least one layer left with non-zero thickness (i.e. we are not down to basement), and the cell's current elevation is higher than the Dean equilibrium profile elevation. So do some beach erosion
                     double
                        dFine = 0,
                        dSand = 0,
                        dCoarse = 0;

                     ErodeBeachConstrained(nX, nY, nThisLayer, -dElevDiff, dFine, dSand, dCoarse);

                     // Totals for this parallel profile
                     dSandDeposited -= dSand;
                     dCoarseDeposited -= dCoarse;

//                      assert(dSandDeposited >= 0);
//                      assert(dCoarseDeposited >= 0);

                     // Totals for the polygon (note that these may become slightly -ve if erosion occurs early during this routine, but this is not a problem since the values will eventually become +ve again
                     dTotSandDeposited -= dSand;
                     dTotCoarseDeposited -= dCoarse;

                     dSandToDepositOnProf += dSand;
                     dSandToDepositOnPoly += dSand;
                     dCoarseToDepositOnProf += dCoarse;
                     dCoarseToDepositOnPoly += dCoarse;

                     // Update this-timestep totals
                     m_ulThisTimestepNumActualBeachErosionCells++;
                     m_dThisTimestepActualBeachErosionFine += dFine;
                     m_dThisTimestepActualBeachErosionSand += dSand;
                     m_dThisTimestepActualBeachErosionCoarse += dCoarse;

   //                   LogStream << "\tPolygon net deposition " << nPoly << " going UP-COAST, actual beach erosion = " << dFine + dSand + dCoarse << " at [" << nX << "][" << nY << "] nCoastPoint = " << nCoastPoint << " nSeawardFromCoast = " << nSeawardFromCoast << endl;
                  }
               }
            }
         }

//          LogStream << "\tIn nDoWithinPolygonBeachRedistribution() going UP-COAST, nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nSeawardOffset = " << nSeawardOffset << " dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << " dSandToDepositOnProf = " << dSandToDepositOnProf << " dCoarseToDepositOnProf = " << dCoarseToDepositOnProf << " dSandToDepositOnPoly = " << dSandToDepositOnPoly << " dCoarseToDepositOnPoly = " << dCoarseToDepositOnPoly << endl;
      }
   }

   // How much have we been able to deposit on this polygon?
   double dTmpTot = dTotSandDeposited + dTotCoarseDeposited;
   string strMsg;
   if (bFPIsEqual(dTmpTot, dTotTargetToDeposit, SEDIMENT_ELEV_TOLERANCE))
      strMsg = "deposition OK";
   else
   {
      if (dTmpTot > dTotTargetToDeposit)
      {
         strMsg = "TOO MUCH deposition";

         m_dThisTimestepMassBalanceDepositionError += (dTmpTot - dTotTargetToDeposit);
      }
      else
      {
         strMsg = "NOT ENOUGH deposition";

         m_dThisTimestepMassBalanceDepositionError += (dTotTargetToDeposit - dTmpTot);
      }
   }

   LogStream << m_ulIteration << ": polygon " << nPoly << " " << strMsg << " deposited = " << (dTotSandDeposited + dTotCoarseDeposited) << " target = " << dTotTargetToDeposit << " (sand = " << dTotSandDeposited << " sand remaining = " << dSandToDepositOnPoly << " coarse = " << dTotCoarseDeposited << " coarse remaining = " << dCoarseToDepositOnPoly << ")" << endl;

   return RTN_OK;
}
