/*!
 *
 * \file do_beach_estimate_actual_erosion.cpp
 * \brief Estimates actual (supply-limited) beach erosion on polygons, and constructs a between-polygon budget of actual beach sediment movement
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
//#include <assert.h>

#include <cmath>
#include <cfloat>
#include <iostream>
using std::cout;
using std::endl;

#include "cme.h"
#include "simulation.h"
#include "coast.h"


/*===============================================================================================================================

 We know the potential erosion for this polygon, so determine actual (supply-limited) erosion in sediment size categories. But do not do any actual erosion yet

===============================================================================================================================*/
int CSimulation::nEstimateActualBeachErosionOnPolygon(int const nCoast, int const nPoly, double const dPotentialErosion)
{
   CGeomCoastPolygon* pPolygon = m_VCoast[nCoast].pGetPolygon(nPoly);

   double
      dSedToErodeOnThisPolygon = dPotentialErosion,
      dTotFineEroded = 0,              // Totals for this polygon
      dTotSandEroded = 0,
      dTotCoarseEroded = 0;

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
      LogStream << m_ulTimestep << ": " << ERR << "in nEstimateActualBeachErosionOnPolygon() for polygon " << nPoly << ", could not find the seaward end point of the up-coast profile (" << nUpCoastProfile << ") for depth of closure = " << m_dDepthOfClosure << ". Lengthen the coastline normals." << endl;

      return RTN_ERR_BAD_BEACH_EROSION_PROFILE;
   }

   // The part-profile length is one greater than nIndex, since pPtiGetCellGivenDepth() returns the index of the cell at depth of closure
   int nUpCoastPartProfileLen = nIndex + 1;

//    assert(bIsWithinGrid(&PtiUpCoastPartProfileSeawardEnd));

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

   // Estimate the volume of sediment which is to be eroded from each parallel profile
   double
      dAllSedimentTargetPerProfile = dPotentialErosion / nCoastSegLen;

   // Traverse the polygon's existing coastline in a DOWN-COAST (i.e. increasing coastpoint indices) sequence, at each coast point fitting a Dean profile which is parallel to the up-coast polygon boundary
   for (int n = 0; n < nCoastSegLen; n++)
   {
      // Get the coast point
      int nCoastPoint = nVCoastPoint[n];
      
      CGeom2DIPoint PtiCoastPoint = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint);
      int
         nCoastX = PtiCoastPoint.nGetX(),
         nCoastY = PtiCoastPoint.nGetY();
         
      // Is the coast cell an intervention structure?
      if (m_pRasterGrid->m_Cell[nCoastX][nCoastY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION)
         // No erosion possible here, so move on
         continue;

      // Not an intervention structure, so calculate the x-y offset between this coast point, and the coast point of the up-coast normal
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
      if (! bIsWithinGrid(nParProfEndX, nParProfEndY))
      {
//          LogStream << WARN << "01 @@@@ while estimating actual beach erosion for coast " << nCoast << " polygon " << nPoly << " in DOWN-COAST direction, hit edge of grid at [" << nParProfEndX << "][" << nParProfEndY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

         KeepWithinGrid(nCoastX, nCoastY, nParProfEndX, nParProfEndY);
         PtiVParProfile[0].SetX(nParProfEndX);
         PtiVParProfile[0].SetY(nParProfEndY);
      }

      double
         dParProfCoastElev = m_pRasterGrid->m_Cell[nCoastX][nCoastY].dGetSedimentTopElev(),
         dParProfEndElev = m_pRasterGrid->m_Cell[nParProfEndX][nParProfEndY].dGetSedimentTopElev();

      int
         nParProfLen,
         nInlandOffset = -1;
      vector<double> dVParProfileDeanElev;

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
//                LogStream << m_ulTimestep << ": reached end of up-coast profile " << nUpCoastProfile << " during down-coast ESTIMATION of actual beach erosion for coast " << nCoast << " polygon " << nPoly << " (nCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << ")" << endl;

               bEndProfile = true;
               break;
            }

            // Extend the parallel profile by one cell in the coastward direction. First get the coords of the cell that is nInlandOffset cells seaward from the existing up-coast part-profile start point
            CGeom2DIPoint PtiUpCoastTmp = *pUpCoastProfile->pPtiGetCellInProfile(nInlandOffset);

            // Then get the offset between this PtiUpCoastTmp cell and the existing up-coast part-profile start point, and use the reverse of this offset to get the co-ords of the cell that extends the existing up-coast part-profile landwards
            int
               nXUpCoastStartOffset = PtiUpCoastTmp.nGetX() - nXUpCoastProfileExistingCoastPoint,
               nYUpCoastStartOffset = PtiUpCoastTmp.nGetY() - nYUpCoastProfileExistingCoastPoint,
               nXUpCoastThisStart = nXUpCoastProfileExistingCoastPoint - nXUpCoastStartOffset,
               nYUpCoastThisStart = nYUpCoastProfileExistingCoastPoint - nYUpCoastStartOffset;

            // Is the new landwards point within the raster grid?
            if (! bIsWithinGrid(nXUpCoastThisStart, nYUpCoastThisStart))
            {
               // It isn't
//                LogStream << WARN << "reached edge of grid at [" << nXUpCoastThisStart << "][" << nYUpCoastThisStart << "] during DOWN-COAST estimation of beach erosion for coast " << nCoast << " polygon " << nPoly << " (nCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << ")" << endl;

               // TODO Need to improve this: at present we just abandon erosion on this coast point and move to another coast point
               bHitEdge = true;
               break;
            }

            CGeom2DIPoint PtiThisUpCoastStart(nXUpCoastThisStart, nYUpCoastThisStart);

            // Append this new landward cell to the up-coast part-profile
            PtiVUpCoastPartProfileCell.push_back(PtiThisUpCoastStart);

            // And calculate the co-ords of a new landwards cell for the parallel profile
            int
               nXParNew = nXUpCoastThisStart + nXOffset,
               nYParNew = nYUpCoastThisStart + nYOffset;

            // Safety check
            if (! bIsWithinGrid(nXParNew, nYParNew))
            {
//                LogStream << WARN << "02 @@@@ while estimating actual beach erosion on coast " << nCoast << " polygon " << nPoly << " in DOWN-COAST direction (nInlandOffset = " << nInlandOffset << "), hit edge of grid at [" << nParProfEndX << "][" << nParProfEndY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its landward end" << endl;

               KeepWithinGrid(nCoastX, nCoastY, nXParNew, nYParNew);
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
         double const dPower = 2.0 / 3.0;
         double dParProfA = (dParProfCoastElev - dParProfEndElev) /  pow(dParProfileLen, dPower);

         // For the parallel profile, calculate the Dean equilibrium profile of the unconsolidated sediment h(y) = A * y^(2/3) where h(y) is the distance below the highest point in the profile at a distance y from the landward start of the profile
         nParProfLen = PtiVParProfile.size();
         dVParProfileDeanElev.resize(nParProfLen, 0);

         double
            dDistFromParProfStart = 0,
            dInc = dParProfileLen / (nParProfLen-1);
         for (int n = 0; n < nParProfLen; n++)
         {
            double dDistBelowHighest = dParProfA * pow(dDistFromParProfStart, dPower);
            dVParProfileDeanElev[n] = dParProfCoastElev - dDistBelowHighest;
            dDistFromParProfStart += dInc;
         }

         double dParProfTotDiff = 0;

         for (int n = 0; n < nParProfLen; n++)
         {
            CGeom2DIPoint PtiTmp = PtiVParProfile[nParProfLen - n - 1];
            int
               nX = PtiTmp.nGetX(),
               nY = PtiTmp.nGetY();

            // Safety check
            if (! bIsWithinGrid(nX, nY))
            {
//                LogStream << WARN << "03 @@@@ while constructing parallel profile to estimate DOWN-COAST actual beach erosion on coast " << nCoast << " polygon " << nPoly << ", hit edge of grid at [" << nParProfEndX << "][" << nParProfEndY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

               KeepWithinGrid(nCoastX, nCoastY, nX, nY);
               PtiTmp.SetX(nX);
               PtiTmp.SetY(nY);
            }

            double
               dTmpElev = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev(),
               dDiff = dTmpElev - dVParProfileDeanElev[n];

            dParProfTotDiff += dDiff;
         }

//          // DEBUG STUFF -----------------------------------------------------
//          LogStream << "\tFor polygon " << nPoly << " doing DOWN-COAST estimation of actual beach erosion, parallel profile from [" << PtiVParProfile.back().nGetX() << "][" << PtiVParProfile.back().nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiVParProfile.back().nGetX()) << ", " <<  dGridCentroidYToExtCRSY(PtiVParProfile.back().nGetY()) << "} to [" << PtiVParProfile[0].nGetX() << "][" << PtiVParProfile[0].nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiVParProfile[0].nGetX()) << ", " <<  dGridCentroidYToExtCRSY(PtiVParProfile[0].nGetY()) << "}, nParProfLen = " << nParProfLen << " dParProfileLen = " << dParProfileLen << " dParProfCoastElev = " << dParProfCoastElev << " dParProfEndElev = " << dParProfEndElev << " dParProfA = " << dParProfA << endl;
//
//          LogStream << "\tPresent profile = ";
//          for (int n = 0; n < nParProfLen; n++)
//          {
//             CGeom2DIPoint PtiTmp = PtiVParProfile[nParProfLen - n - 1];
//             int
//                nX = PtiTmp.nGetX(),
//                nY = PtiTmp.nGetY();
//
//             // Safety check
//             if (! bIsWithinGrid(nX, nY))
//             {
//                KeepWithinGrid(nX, nY);
//             }
//
//             // Don't do anything to intervention cells
//             if (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION)
//                continue;
// 
//             // Don't do cells twice
//             if (! m_pRasterGrid->m_Cell[nX][nY].bBeachErosionOrDepositionThisTimestep())
//             {
//                LogStream << m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev() << " ";
//             }
//          }
//          LogStream << endl;
//
//          LogStream << "\tParallel Dean profile for erosion = ";
//          for (int n = 0; n < nParProfLen; n++)
//          {
//             LogStream << dVParProfileDeanElev[n] << " ";
//          }
//          LogStream << endl;
//
//          LogStream << "\tnCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << endl << "\tDifference = ";
//          for (int n = 0; n < nParProfLen; n++)
//          {
//             CGeom2DIPoint PtiTmp = PtiVParProfile[nParProfLen - n - 1];
//             int
//                nX = PtiTmp.nGetX(),
//                nY = PtiTmp.nGetY();
//
//             // Safety check
//             if (! bIsWithinGrid(nX, nY))
//             {
//                KeepWithinGrid(nX, nY);
//             }
//
//             // Don't do anything to intervention cells
//             if (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION)
//                continue;
// 
//             // Don't do cells twice
//             if (! m_pRasterGrid->m_Cell[nX][nY].bBeachErosionOrDepositionThisTimestep())
//             {
//                double
//                   dTmpElev = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev(),
//                   dDiff = dTmpElev - dVParProfileDeanElev[n];
//                LogStream << dDiff << " ";
//             }
//          }
//          LogStream << endl << endl;
//          // DEBUG STUFF -----------------------------------------------------

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

      // OK, this value of nInlandOffset gives us enough erosion. So calculate how much fine-, sand-, and coarse sized sediment can be obtained here by working along the parallel profile from the landward end (which is inland from the existing coast, if nInlandOffset > 0)
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
//             LogStream << "\tIn nEstimateActualBeachErosionOnPolygon() going DOWN-COAST, nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << " leaving loop at start of timestep (" << nDistSeawardFromNewCoast << " / " << nParProfLen << ") because enough erosion for polygon, dSedToErodeOnThisPolygon = " << dSedToErodeOnThisPolygon << endl;
            break;
         }

         // Leave the loop if we have done enough erosion for this profile
         if (dSedToErodeOnThisProfile <= 0)
         {
//             LogStream << "\tIn nEstimateActualBeachErosionOnPolygon() going DOWN-COAST, nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << " leaving loop at start of timestep (" << nDistSeawardFromNewCoast << " / " << nParProfLen << ") because enough erosion for profile, dSedToErodeOnThisProfile = " << dSedToErodeOnThisProfile << " dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << endl;

            break;
         }

         CGeom2DIPoint PtiTmp = PtiVParProfile[nParProfLen - nDistSeawardFromNewCoast - 1];
         int
            nX = PtiTmp.nGetX(),
            nY = PtiTmp.nGetY();

         // Safety check
         if (! bIsWithinGrid(nX, nY))
         {
//             LogStream << WARN << "04 @@@@ while estimating sediment size fractions on eroding polygon " << nPoly << " in DOWN-COAST direction, hit edge of grid at [" << nParProfEndX << "][" << nParProfEndY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

            KeepWithinGrid(nCoastX, nCoastY, nX, nY);
            PtiTmp.SetX(nX);
            PtiTmp.SetY(nY);
         }

         // Don't do anything to intervention cells
         if (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION)
            continue;

         // Don't do cells twice
         if (! m_pRasterGrid->m_Cell[nX][nY].bGetActualBeachErosionEstimated())
         {
            // Set flag so we don't estimate erosion/deposition on this cell again
            m_pRasterGrid->m_Cell[nX][nY].SetActualBeachErosionEstimated();

            // Get this cell's current elevation
            double dThisElevNow = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev();

   //          LogStream << "\tnPoly = " << nPoly << " going DOWN-COAST, [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "}  nCoastPoint = " << nCoastPoint << " nDistSeawardFromNewCoast = " << nDistSeawardFromNewCoast << " dThisElevNow = " << dThisElevNow << " Dean Elev = " << dVParProfileDeanElev[nDistSeawardFromNewCoast] << endl;

            // Subtract the two elevations
            double dElevDiff = dThisElevNow - dVParProfileDeanElev[nDistSeawardFromNewCoast];
            if (dElevDiff > SEDIMENT_ELEV_TOLERANCE)
            {
               // The current elevation is higher than the Dean elevation, so we have possible beach erosion (i.e. if not constrained by availability of unconsolidated sediment) here
//                LogStream << "\tnPoly = " << nPoly << " doing DOWN-COAST, possible beach erosion = " << dElevDiff << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "} nCoastPoint = " << nCoastPoint << " nDistSeawardFromNewCoast = " << nDistSeawardFromNewCoast << " dThisElevNow = " << dThisElevNow << " Dean Elev = " << dVParProfileDeanElev[nDistSeawardFromNewCoast] << endl;

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
                  EstimateActualBeachErosionOnCell(nX, nY, nThisLayer, dToErode, dFine, dSand, dCoarse);

                  // Totals for this parallel profile
                  dFineEroded += dFine;
                  dSandEroded += dSand;
                  dCoarseEroded += dCoarse;

                  dTotFineEroded += dFine;
                  dTotSandEroded += dSand;
                  dTotCoarseEroded += dCoarse;

                  double dTmpTot = dFine + dSand + dCoarse;
                  dSedToErodeOnThisProfile -= dTmpTot;
                  dSedToErodeOnThisPolygon -= dTmpTot;

//                   LogStream << "\tIn nEstimateActualBeachErosionOnPolygon(), nPoly = " << nPoly << " going DOWN-COAST, actual beach erosion = " << dTmpTot << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "}  nCoastPoint = " << nCoastPoint << " nDistSeawardFromNewCoast = " << nDistSeawardFromNewCoast << endl;
               }
            }
            else if ((dElevDiff < -SEDIMENT_ELEV_TOLERANCE) && ((m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea()) || (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_DRIFT)))
            {
               // The current elevation is below the Dean elevation, so we have can have beach deposition here provided that we have some previously-eroded unconsolidated sediment (sand and coarse only) to deposit
               double dTotSandAndCoarseEroded = dTotSandEroded + dTotCoarseEroded;
               if (dTotSandAndCoarseEroded > SEDIMENT_ELEV_TOLERANCE)
               {
                  // Assume that the sediment size fractions which are deposited are in the same ratios as the previously-eroded sediment
                  double
                     dToDeposit = tMin(-dElevDiff, dTotSandAndCoarseEroded),
                     dSandProp = dTotSandEroded / dTotSandAndCoarseEroded,
                     dCoarseProp = 1 - dSandProp,
                     dSandToDeposit = tMin(dToDeposit * dSandProp, dTotSandEroded),
                     dCoarseToDeposit = tMin(dToDeposit * dCoarseProp, dTotCoarseEroded);

                  if (dSandToDeposit > SEDIMENT_ELEV_TOLERANCE)
                  {
                     // Totals for this parallel profile
                     dSandEroded -= dSandToDeposit;
                     dSedToErodeOnThisProfile += dSandToDeposit;

                     // Totals for the polygon
                     dTotSandEroded -= dSandToDeposit;
                     dSedToErodeOnThisPolygon += dSandToDeposit;
                  }

                  if (dCoarseToDeposit > SEDIMENT_ELEV_TOLERANCE)
                  {
                     // Totals for this parallel profile
                     dCoarseEroded -= dCoarseToDeposit;
                     dSedToErodeOnThisProfile += dCoarseToDeposit;

                     // Totals for the polygon
                     dTotCoarseEroded -= dCoarseToDeposit;
                     dSedToErodeOnThisPolygon += dCoarseToDeposit;
                  }

//                   LogStream << "\tIn nEstimateActualBeachErosionOnPolygon(), nPoly = " << nPoly << " going DOWN-COAST, beach deposition = " << dSandToDeposit + dCoarseToDeposit << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "}  nCoastPoint = " << nCoastPoint << " nDistSeawardFromNewCoast = " << nDistSeawardFromNewCoast << endl;
               }
            }
         }
      }

//       LogStream << "\tIn nEstimateActualBeachErosionOnPolygon() going DOWN-COAST, nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << " dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << " dSedToErodeOnThisProfile = " <<  dAllSedimentTargetPerProfile << " dSedToErodeOnThisPolygon = " << dSedToErodeOnThisPolygon << " dPotentialErosion = " << dPotentialErosion << endl;

   }

   // Have we eroded the full potential-erosion amount?
   if (dSedToErodeOnThisPolygon > SEDIMENT_ELEV_TOLERANCE)
   {
      // No, so traverse the polygon's existing coastline in an UP-COAST (i.e. decreasing coastpoint indices) sequence, at each coast point fitting a Dean profile which is parallel to the down-coast polygon boundary. Start by finding the seaward end point of the up-coast part-profile, this does not change as the landwards offset changes
      int nIndex = pDownCoastProfile->nGetCellGivenDepth(m_pRasterGrid, m_dDepthOfClosure);
      if (nIndex == INT_NODATA)
      {
         LogStream << m_ulTimestep << ": " << ERR << "in nEstimateActualBeachErosionOnPolygon() for polygon " << nPoly << ", could not find the seaward end point of the down-coast profile (" << nDownCoastProfile << ") for depth of closure = " << m_dDepthOfClosure << ". Lengthen the coastline normals." << endl;

         return RTN_ERR_BAD_BEACH_EROSION_PROFILE;
      }

      // The part-profile length is one greater than nIndex, since pPtiGetCellGivenDepth() returns the index of the cell at depth of closure
      int nDownCoastPartProfileLen = nIndex + 1;

//       assert(bIsWithinGrid(&PtiDownCoastPartProfileSeawardEnd));

      // Store the cell co-ordinates of the boundary part-profile in reverse (sea to coast) order so we can append to the coastward end as we move inland (i.e. as nInlandOffset increases)
      vector<CGeom2DIPoint> PtiVDownCoastPartProfileCell;
      for (int n = 0; n < nDownCoastPartProfileLen; n++)
         PtiVDownCoastPartProfileCell.push_back(*pDownCoastProfile->pPtiGetCellInProfile(nDownCoastPartProfileLen - n - 1));

      int
         nXDownCoastProfileExistingCoastPoint = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nDownCoastProfileCoastPoint)->nGetX(),
         nYDownCoastProfileExistingCoastPoint = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nDownCoastProfileCoastPoint)->nGetY();

      // Store the coast point numbers for this polygon
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

      // Now traverse the polygon's existing coastline, fitting a Dean profile at each coast point
      for (int n = 0; n < nCoastSegLen; n++)
      {
         // Get the coast point
         int nCoastPoint = nVCoastPoint[n];
         CGeom2DIPoint PtiCoastPoint = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint);
         int
            nCoastX = PtiCoastPoint.nGetX(),
            nCoastY = PtiCoastPoint.nGetY();

         // Is this cell an intervention structure?
         if (m_pRasterGrid->m_Cell[nCoastX][nCoastY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION)
            // No erosion possible here, so move on
            continue;

         // Not an intervention structure, so calculate the x-y offset between this coast point, and the coast point of the down-coast normal
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
         if (! bIsWithinGrid(nParProfEndX, nParProfEndY))
         {
//             LogStream << WARN << "05 @@@@ while estimating actual beach erosion on coast " << nCoast << " polygon " << nPoly << " in UP-COAST direction, hit edge of grid at [" << nParProfEndX << "][" << nParProfEndY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

            KeepWithinGrid(nCoastX, nCoastY, nParProfEndX, nParProfEndY);
            PtiVParProfile[0].SetX(nParProfEndX);
            PtiVParProfile[0].SetY(nParProfEndY);
         }

         double
            dParProfCoastElev = m_pRasterGrid->m_Cell[nCoastX][nCoastY].dGetSedimentTopElev(),
            dParProfEndElev = m_pRasterGrid->m_Cell[nParProfEndX][nParProfEndY].dGetSedimentTopElev();

         int
            nParProfLen,
            nInlandOffset = -1;
         vector<double> dVParProfileDeanElev;

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
//                   LogStream << m_ulTimestep << ": reached end of down-coast profile " << nDownCoastProfile << " during up-coast ESTIMATION of actual beach erosion for coast " << nCoast << " polygon " << nPoly << " (nCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << ")" << endl;

                  bEndProfile = true;
                  break;
               }

               // Extend the parallel profile by one cell in the coastward direction. First get the coords of the cell that is nInlandOffset cells seaward from the existing down-coast part-profile start point
               CGeom2DIPoint PtiDownCoastTmp = *pDownCoastProfile->pPtiGetCellInProfile(nInlandOffset);

               // Then get the offset between this PtiDownCoastTmp cell and the existing down-coast part-profile start point, and use the reverse of this offset to get the co-ords of the cell that extends the existing up-coast part-profile landwards
               int
                  nXDownCoastStartOffset = PtiDownCoastTmp.nGetX() - nXDownCoastProfileExistingCoastPoint,
                  nYDownCoastStartOffset = PtiDownCoastTmp.nGetY() - nYDownCoastProfileExistingCoastPoint,
                  nXDownCoastThisStart = nXDownCoastProfileExistingCoastPoint - nXDownCoastStartOffset,
                  nYDownCoastThisStart = nYDownCoastProfileExistingCoastPoint - nYDownCoastStartOffset;

               // Is the new landwards point within the raster grid?
               if (! bIsWithinGrid(nXDownCoastThisStart, nYDownCoastThisStart))
               {
                  // It isn't
//                   LogStream << m_ulTimestep << ": " << WARN << "reached edge of grid at [" << nXDownCoastThisStart << "][" << nYDownCoastThisStart << "] during UP-COAST estimation of actual beach erosion of unconsolidated sediment for coast " << nCoast << " polygon " << nPoly << " (nCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << ")" << endl;

                  // TODO Need to improve this: at present just abandons erosion on this coast point and moves to another coast point
                  bHitEdge = true;
                  break;
               }

               CGeom2DIPoint PtiThisDownCoastStart(nXDownCoastThisStart, nYDownCoastThisStart);

               // Append this new landward cell to the down-coast part-profile
               PtiVDownCoastPartProfileCell.push_back(PtiThisDownCoastStart);

               // And calculate the co-ords of a new landwards cell for the parallel profile
               int
                  nXParNew = nXDownCoastThisStart + nXOffset,
                  nYParNew = nYDownCoastThisStart + nYOffset;

               // Safety check
               if (! bIsWithinGrid(nXParNew, nYParNew))
               {
//                   LogStream << WARN << "06 @@@@ while estimating actual beach erosion on coast " << nCoast << " polygon " << nPoly << " in UP-COAST direction, hit edge of grid at [" << nParProfEndX << "][" << nParProfEndY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its landward end" << endl;

                  KeepWithinGrid(nCoastX, nCoastY, nXParNew, nYParNew);
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
            double const dPower = 2.0 / 3.0;
            double dParProfA = (dParProfCoastElev - dParProfEndElev) /  pow(dParProfileLen, dPower);

            // For the parallel profile, calculate the Dean equilibrium profile of the unconsolidated sediment h(y) = A * y^(2/3) where h(y) is the distance below the highest point in the profile at a distance y from the landward start of the profile
            nParProfLen = PtiVParProfile.size();
            dVParProfileDeanElev.resize(nParProfLen, 0);

            double
               dDistFromParProfStart = 0,
               dInc = dParProfileLen / (nParProfLen-1);
            for (int n = 0; n < nParProfLen; n++)
            {
               double dDistBelowHighest = dParProfA * pow(dDistFromParProfStart, dPower);
               dVParProfileDeanElev[n] = dParProfCoastElev - dDistBelowHighest;
               dDistFromParProfStart += dInc;
            }

            double dParProfTotDiff = 0;
            for (int n = 0; n < nParProfLen; n++)
            {
               CGeom2DIPoint PtiTmp = PtiVParProfile[nParProfLen - n - 1];
               int
                  nX = PtiTmp.nGetX(),
                  nY = PtiTmp.nGetY();

               // Safety check
               if (! bIsWithinGrid(nX, nY))
               {
//                   LogStream << WARN << "07 @@@@ while constructing parallel profile to estimate UP-COAST actual beach erosion on polygon " << nPoly << ", hit edge of grid at [" << nParProfEndX << "][" << nParProfEndY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

                  KeepWithinGrid(nCoastX, nCoastY, nX, nY);
                  PtiTmp.SetX(nX);
                  PtiTmp.SetY(nY);
               }

               double
                  dTmpElev = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev(),
                  dDiff = dTmpElev - dVParProfileDeanElev[n];

               dParProfTotDiff += dDiff;
            }

//             // DEBUG STUFF -----------------------------------------------------
//             LogStream << "\tFor polygon " << nPoly << " doing UP-COAST estimation of actual beach erosion, parallel profile from [" << PtiVParProfile.back().nGetX() << "][" << PtiVParProfile.back().nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiVParProfile.back().nGetX()) << ", " <<  dGridCentroidYToExtCRSY(PtiVParProfile.back().nGetY()) << "} to [" << PtiVParProfile[0].nGetX() << "][" << PtiVParProfile[0].nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiVParProfile[0].nGetX()) << ", " <<  dGridCentroidYToExtCRSY(PtiVParProfile[0].nGetY()) << "}, nParProfLen = " << nParProfLen << " dParProfileLen = " << dParProfileLen << " dParProfCoastElev = " << dParProfCoastElev << " dParProfEndElev = " << dParProfEndElev << " dParProfA = " << dParProfA << endl;
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
//                if (! bIsWithinGrid(nX, nY))
//                {
//                   KeepWithinGrid(nX, nY);
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
//                LogStream << dVParProfileDeanElev[n] << " ";
//             }
//             LogStream << endl;
//
//             LogStream << "\tnCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << endl << "\tDifference = ";
//             for (int n = 0; n < nParProfLen; n++)
//             {
//                CGeom2DIPoint PtiTmp = PtiVParProfile[nParProfLen - n - 1];
//                int
//                   nX = PtiTmp.nGetX(),
//                   nY = PtiTmp.nGetY();
//
//                // Safety check
//                if (! bIsWithinGrid(nX, nY))
//                {
//                   KeepWithinGrid(nX, nY);
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
//                      dDiff = dTmpElev - dVParProfileDeanElev[n];
//                   LogStream << dDiff << " ";
//                }
//             }
//             LogStream << endl << endl;
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
//                LogStream << "\tIn nRouteActualBeachErosionToAdjacentPolygons() estimating actual beach erosion going UP-COAST, nCoast = " << nCoast << " nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << " leaving loop at start of timestep (" << nDistSeawardFromNewCoast << " / " << nParProfLen << ") because enough erosion for polygon, dSedToErodeOnThisPolygon = " << dSedToErodeOnThisPolygon << endl;

               break;
            }

            // Leave the loop if we have done enough erosion for this profile
            if (dSedToErodeOnThisProfile <= 0)
            {
//                LogStream << "\tIn nRouteActualBeachErosionToAdjacentPolygons() estimating actual beach erosion going UP-COAST, nCoast = " << nCoast << " nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << " leaving loop at start of timestep (" << nDistSeawardFromNewCoast << " / " << nParProfLen << ") because enough erosion for profile, dSedToErodeOnThisProfile = " << dSedToErodeOnThisProfile << " dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << endl;

               break;
            }

            CGeom2DIPoint PtiTmp = PtiVParProfile[nParProfLen - nDistSeawardFromNewCoast - 1];
            int
               nX = PtiTmp.nGetX(),
               nY = PtiTmp.nGetY();

            // Safety check
            if (! bIsWithinGrid(nX, nY))
            {
//                LogStream << WARN << "08 @@@@ while constructing parallel profile to estimate UP-COAST actual beach erosion for coast " << nCoast << " polygon " << nPoly << ", hit edge of grid at [" << nParProfEndX << "][" << nParProfEndY << "] for parallel profile from coast point " << nCoastPoint << " at [" << nCoastX << "][" << nCoastY << "]. Constraining this parallel profile at its seaward end" << endl;

               KeepWithinGrid(nCoastX, nCoastY, nX, nY);
               PtiTmp.SetX(nX);
               PtiTmp.SetY(nY);
            }

            // Don't do anything to intervention cells
            if (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION)
               continue;

            // Don't do cells twice
            if (! m_pRasterGrid->m_Cell[nX][nY].bGetActualBeachErosionEstimated())
            {
               // Set flag so we don't do estimate erosion/deposition on this cell a second time
               m_pRasterGrid->m_Cell[nX][nY].SetActualBeachErosionEstimated();

               // Get this cell's current elevation
               double dThisElevNow = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev();

//                LogStream << "\tnPoly = " << nPoly << " going UP-COAST, [" << nX << "][" << nY << "]  = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "} nCoastPoint = " << nCoastPoint << " nDistSeawardFromNewCoast = " << nDistSeawardFromNewCoast << " dThisElevNow = " << dThisElevNow << " Dean Elev = " << dVParProfileDeanElev[nDistSeawardFromNewCoast] << endl;

               // Subtract the two elevations
               double dElevDiff = dThisElevNow - dVParProfileDeanElev[nDistSeawardFromNewCoast];
               if (dElevDiff > SEDIMENT_ELEV_TOLERANCE)
               {
                  // The current elevation is higher than the Dean elevation, so we have possible beach erosion (i.e. if not constrained by availability of unconsolidated sediment) here
//                   LogStream << "\tnPoly = " << nPoly << " doing UP-COAST, possible beach erosion = " << dElevDiff << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "} nCoastPoint = " << nCoastPoint << " nDistSeawardFromNewCoast = " << nDistSeawardFromNewCoast << " dThisElevNow = " << dThisElevNow << " Dean Elev = " << dVParProfileDeanElev[nDistSeawardFromNewCoast] << endl;

                  // Now get the number of the highest layer with non-zero thickness
                  int nThisLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();

                  // Safety check
                  if (nThisLayer == INT_NODATA)
                     return RTN_ERR_NO_TOP_LAYER;

                  if (nThisLayer != NO_NONZERO_THICKNESS_LAYERS)
                  {
                     // We still have at least one layer left with non-zero thickness (i.e. we are not down to basement), and the cell's current elevation is higher than the Dean equilibrium profile elevation. So we can have some beach erosion here
                     double
                        dToErode = tMin(dElevDiff, dSedToErodeOnThisProfile, dSedToErodeOnThisPolygon),
                        dFine = 0,
                        dSand = 0,
                        dCoarse = 0;
      //                assert(dToErode > 0);
                     EstimateActualBeachErosionOnCell(nX, nY, nThisLayer, dToErode, dFine, dSand, dCoarse);

                     // Totals for this parallel profile
                     dFineEroded += dFine;
                     dSandEroded += dSand;
                     dCoarseEroded += dCoarse;

                     // Totals for the polygon
                     dTotFineEroded += dFine;
                     dTotSandEroded += dSand;
                     dTotCoarseEroded += dCoarse;

                     double dTmpTot = dFine + dSand + dCoarse;
                     dSedToErodeOnThisProfile -= dTmpTot;
                     dSedToErodeOnThisPolygon -= dTmpTot;

//                      LogStream << "\tIn nEstimateActualBeachErosionOnPolygon(), nPoly = " << nPoly << " going UP-COAST, actual beach erosion = " << dTmpTot << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "}  nCoastPoint = " << nCoastPoint << " nDistSeawardFromNewCoast = " << nDistSeawardFromNewCoast << endl;
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
                        dToDeposit = tMin(-dElevDiff, dTotSandAndCoarseEroded),
                        dSandProp = dTotSandEroded / dTotSandAndCoarseEroded,
                        dCoarseProp = 1 - dSandProp,
                        dSandToDeposit = tMin(dToDeposit * dSandProp, dTotSandEroded),
                        dCoarseToDeposit = tMin(dToDeposit * dCoarseProp, dTotCoarseEroded);

                     if (dSandToDeposit > SEDIMENT_ELEV_TOLERANCE)
                     {
                        // Totals for this parallel profile
                        dSandEroded -= dSandToDeposit;
                        dSedToErodeOnThisProfile += dSandToDeposit;

                        // Totals for the polygon
                        dTotSandEroded -= dSandToDeposit;
                        dSedToErodeOnThisPolygon += dSandToDeposit;
                     }

                     if (dCoarseToDeposit > SEDIMENT_ELEV_TOLERANCE)
                     {
                        // Totals for this parallel profile
                        dCoarseEroded -= dCoarseToDeposit;
                        dSedToErodeOnThisProfile += dCoarseToDeposit;

                        // Totals for the polygon
                        dTotCoarseEroded -= dCoarseToDeposit;
                        dSedToErodeOnThisPolygon += dCoarseToDeposit;
                     }

//                      LogStream << "\tIn nEstimateActualBeachErosionOnPolygon(), nPoly = " << nPoly << " going UP-COAST, beach deposition = " << dSandToDeposit + dCoarseToDeposit << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " <<  dGridCentroidYToExtCRSY(nY) << "}  nCoastPoint = " << nCoastPoint << " nDistSeawardFromNewCoast = " << nDistSeawardFromNewCoast << endl;
                  }
               }
            }
         }

//          LogStream << "\tIn nRouteActualBeachErosionToAdjacentPolygons() going UP-COAST, nPoly = " << nPoly << " nCoastPoint = " << nCoastPoint << " nInlandOffset = " << nInlandOffset << " dAllSedimentTargetPerProfile = " << dAllSedimentTargetPerProfile << " dSedToErodeOnThisProfile = " <<  dAllSedimentTargetPerProfile << " dSedToErodeOnThisPolygon = " << dSedToErodeOnThisPolygon << " dPotentialErosion = " << dPotentialErosion << endl;
      }
   }

   // How much have we been able to erode?
//    LogStream << "\tIn nEstimateActualBeachErosionOnPolygon() going UP-COAST, nPoly = " << nPoly << " dTotFineEroded = " << dTotFineEroded << " dTotSandEroded = " << dTotSandEroded << " dTotCoarseEroded = " << dTotCoarseEroded << endl;

//    double dEstimatedErosion = dTotFineEroded + dTotSandEroded + dTotCoarseEroded;

   // Is the estimated depth eroded within TOLERANCE of the potential erosion depth-equivalent?
//    if (bFPIsEqual(dEstimatedErosion, dPotentialErosion, TOLERANCE))
//       LogStream << "Polygon " << nPoly << " has estimated actual beach erosion approx equal to potential beach erosion: estimated = " << -dEstimatedErosion << " potential = " << -dPotentialErosion << endl;
//    else
//       LogStream << "Polygon " << nPoly << " has estimated actual beach erosion less than potential beach erosion: estimated = " << -dEstimatedErosion << " potential = " << -dPotentialErosion << endl;

   // And save these values
   pPolygon->SetDeltaEstimatedUnconsFine(-dTotFineEroded);
   pPolygon->SetDeltaEstimatedUnconsSand(-dTotSandEroded);
   pPolygon->SetDeltaEstimatedUnconsCoarse(-dTotCoarseEroded);

   // Save the estimated values
   m_dThisTimestepEstimatedActualFineBeachErosion   += dTotFineEroded;
   m_dThisTimestepEstimatedActualSandBeachErosion   += dTotSandEroded;
   m_dThisTimestepEstimatedActualCoarseBeachErosion += dTotCoarseEroded;

   return RTN_OK;
}


/*===============================================================================================================================

 Estimates actual (supply-limited) erosion of unconsolidated beach sediment on a single cell, returns the depth-equivalents of fine, sand and coarse sediment which could be removed

===============================================================================================================================*/
void CSimulation::EstimateActualBeachErosionOnCell(int const nX, int const nY, int const nThisLayer, double const dPotentialErosion, double& dFine, double& dSand, double& dCoarse)
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
      dTotErodibility = (nFineWeight * m_dFineErodibilityNormalized) + (nSandWeight * m_dSandErodibilityNormalized) + (nCoarseWeight * m_dCoarseErodibilityNormalized);

   if (nFineWeight)
   {
      // Pretend to erode some fine-sized consolidated sediment
      double dFineLowering = (m_dFineErodibilityNormalized * dPotentialErosion) / dTotErodibility;

      // Make sure we don't get -ve amounts left on the cell
      dFine = tMin(dExistingAvailableFine, dFineLowering);
   }

   if (nSandWeight)
   {
      // Pretend to erode some sand-sized consolidated sediment
      double dSandLowering = (m_dSandErodibilityNormalized * dPotentialErosion) / dTotErodibility;

      // Make sure we don't get -ve amounts left on the source cell
      dSand = tMin(dExistingAvailableSand, dSandLowering);
   }

   if (nCoarseWeight)
   {
      // Erode some coarse-sized consolidated sediment
      double dCoarseLowering = (m_dCoarseErodibilityNormalized * dPotentialErosion) / dTotErodibility;

      // Make sure we don't get -ve amounts left on the source cell
      dCoarse = tMin(dExistingAvailableCoarse, dCoarseLowering);
   }

//    LogStream << "\tEstimating actual erosion on [" << nX << "][" << nY << "] dPotentialErosion = " << dPotentialErosion << " dTotActualErosion = " << dTotActualErosion << " dFine = " << dFine << " dSand = " << dSand << " dCoarse = " << dCoarse << endl;
}


