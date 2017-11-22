/*!
 *
 * \file create_polygons.cpp
 * \brief Creates coast polygons for sediment transport calcs
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
#include <assert.h>
#include <iostream>
using std::endl;

#include <stack>
using std::stack;

#include "cme.h"
#include "simulation.h"
#include "coast.h"


/*===============================================================================================================================

 Create polygons, and marks the polygon boundaries on the raster grid

===============================================================================================================================*/
int CSimulation::nCreateAllPolygons(void)
{
   // Global polygon count
   m_nGlobalPolygonID = -1;

   // Do this for each coast
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      // Do this for every point on the coastline
      int
         nPolygon = -1,          // This-coast-only polygon ID
         nNode = -1,
         nPrevProfileCoastPoint = -1,
         nPrevProfile = -1;
      for (int nCoastPoint = 0; nCoastPoint < m_VCoast[nCoast].nGetCoastlineSize(); nCoastPoint++)
      {
         int nThisProfile = m_VCoast[nCoast].nGetProfileNumber(nCoastPoint);
         CGeomProfile* pThisProfile = m_VCoast[nCoast].pGetProfile(nThisProfile);

         if ((nThisProfile != INT_NODATA) && (pThisProfile->bOKIncStartAndEndOfCoast()))
         {
            // There is a valid profile at this coast point
            if (nPrevProfileCoastPoint >= 0)
            {
               // Calculate half the along-coast distance (in coast points) between this profile and the previous (i.e. up-coast) profile
               int nDist = (nCoastPoint - nPrevProfileCoastPoint) / 2;

               // OK, set the node point in the coast object. We do this now, instead of earlier on, since some profiles (i.e. polygon boundaries) may have been marked as invalid
               int nNodePoint = nCoastPoint - nDist;
               m_VCoast[nCoast].SetPolygonNode(nNodePoint, ++nNode);

               // Get the grid CRS co-ords of the coast node
               CGeom2DIPoint PtiNode = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nNodePoint);

               // Get the previous profile, also set some defaults (assuming for now that this polygon is not approximately triangular i.e. both normals do not meet)
               CGeomProfile* pPrevProfile = m_VCoast[nCoast].pGetProfile(nPrevProfile);
               int
                  nPrevProfileEnd = pPrevProfile->nGetProfileSize()-1,
                  nThisProfileEnd = pThisProfile->nGetProfileSize()-1;
               bool bMeetsAtAPoint = false;
               CGeom2DPoint PtCoastwardTip;

//                // DEBUG CODE =============================
//                CGeom2DPoint PtPrevEndTmp = *pPrevProfile->pPtGetPointInProfile(nPrevProfileEnd);
//                double
//                   dXTmp = dExtCRSXToGridX(PtPrevEndTmp.dGetX()),
//                   dYTmp = dExtCRSYToGridY(PtPrevEndTmp.dGetY());
//                assert(bIsWithinValidGrid(dXTmp, dYTmp));
//
//                CGeom2DPoint PtThisEndTmp = *pThisProfile->pPtGetPointInProfile(nThisProfileEnd);
//                dXTmp = dExtCRSXToGridX(PtThisEndTmp.dGetX());
//                dYTmp = dExtCRSYToGridY(PtThisEndTmp.dGetY());
//                assert(bIsWithinValidGrid(dXTmp, dYTmp));
//                // DEBUG CODE =============================

               // Now check to see if the two normals do meet i.e. if they are coincident
               if (pThisProfile->bFindProfileInCoincidentProfiles(nPrevProfile))
               {
                  // Yes they do meet
                  bMeetsAtAPoint = true;

                  // Find the most coastward point at which this normal and the previous normal touch. If they do not touch, the polygon requires a 'joining line'
                  pThisProfile->GetMostCoastwardSharedLineSegment(nPrevProfile, nThisProfileEnd, nPrevProfileEnd);
                  if (nThisProfileEnd == -1)
                  {
                     LogStream << m_ulIteration << ": " << ERR << "profile " << nPrevProfile << " should be coincident with profile " << nThisProfile << " but was not found" << endl;
                     return RTN_ERR_BAD_MULTILINE;
                  }

                  PtCoastwardTip = *pThisProfile->pPtGetPointInProfile(nThisProfileEnd);
               }

               // Next, calculate the polygon's boundary (external CRS), do this in an anti-clockwise sequence
               vector<CGeom2DPoint> PtVBoundary;

               // Start appending points: begin at the node point, then move down-coast as far as the down-coast (this) normal
               for (int i = nNodePoint; i <= nCoastPoint; i++)
                  PtVBoundary.push_back(*m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(i));

               // Use the penultimate coastline point as the start point for the point-in-polygon search later, during flood fill
               int nPointInPolygonStartPoint = PtVBoundary.size() - 2;

               // Append the points in the down-coast normal. Omit the last point of this normal if the the most seaward point of the this normal, and the most seaward point of the up-coast (previous) normal are the same
               int nFinishPoint = nThisProfileEnd;
               if (bMeetsAtAPoint)
                  nFinishPoint--;
               for (int i = 0; i <= nFinishPoint; i++)
               {
                  CGeom2DPoint PtThis = *pThisProfile->pPtGetPointInProfile(i);
                  if (! (PtThis == &PtVBoundary.back()))
                     PtVBoundary.push_back(PtThis);
               }

               // Append the points in the up-coast (previous) normal, in reverse order
               for (int i = nPrevProfileEnd; i >= 0; i--)
               {
                  CGeom2DPoint PtThis = *pPrevProfile->pPtGetPointInProfile(i);
                  if (! (PtThis == &PtVBoundary.back()))
                     PtVBoundary.push_back(PtThis);
               }

               // Append the points from the remaining bit of coast, moving down-coast to finish at the node point. Note that we must include the node point here, in order to obtain a closed polygon
               for (int i = nPrevProfileCoastPoint; i <= nNodePoint; i++)
               {
                  CGeom2DPoint PtThis = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(i);
                  if (! (PtThis == &PtVBoundary.back()))
                     PtVBoundary.push_back(PtThis);
               }

               // Now identify the 'anti-node', this is the seaward point 'opposite' the polygon's coastal node
               CGeom2DIPoint PtiAntiNode;
               if (bMeetsAtAPoint)
                  PtiAntiNode = PtiExtCRSToGrid(&PtCoastwardTip);
               else
               {
                  CGeom2DPoint PtAvg = PtAverage(pThisProfile->pPtGetPointInProfile(nThisProfileEnd), pPrevProfile->pPtGetPointInProfile(nPrevProfileEnd));
                  PtiAntiNode = PtiExtCRSToGrid(&PtAvg);
               }

               // Create the coast's polygon object
               m_VCoast[nCoast].CreatePolygon(++m_nGlobalPolygonID, ++nPolygon, nNodePoint, &PtiNode, &PtiAntiNode, nPrevProfile, nThisProfile, &PtVBoundary, nPrevProfileEnd+1, nThisProfileEnd+1, nPointInPolygonStartPoint);

               // Get a pointer to this polygon object
               CGeomCoastPolygon* pPolygon = m_VCoast[nCoast].pGetPolygon(nPolygon);

               // And store this pointer for simulation-wide access
               m_pVCoastPolygon.push_back(pPolygon);

               // Now rasterize the polygon boundaries: first, the coastline. This is necessary so that sand/coarse sediment derived from platform erosion of the coast cells is correctly added to the containing polygon's unconsolidated sediment
               for (int i = nPrevProfileCoastPoint; i <= nCoastPoint; i++)
               {
                  CGeom2DIPoint PtiToMark = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(i);
                  m_pRasterGrid->m_Cell[PtiToMark.nGetX()][PtiToMark.nGetY()].SetPolygonID(m_nGlobalPolygonID);
               }

               // Do the upcoast normal profile (does the whole length, including any shared line segments. So some cells are marked twice, however this is not a problem)
               int nCellsInProfile = pPrevProfile->nGetNumCellsInProfile();
               for (int i = 0; i < nCellsInProfile; i++)
               {
                  CGeom2DIPoint PtiToMark = *pPrevProfile->pPtiGetCellInProfile(i);
                  m_pRasterGrid->m_Cell[PtiToMark.nGetX()][PtiToMark.nGetY()].SetPolygonID(m_nGlobalPolygonID);
               }

               // Do the downcoast normal profile (again does the whole length, including any shared line segments)
               nCellsInProfile = pThisProfile->nGetNumCellsInProfile();
               for (int i = 0; i < nCellsInProfile; i++)
               {
                  CGeom2DIPoint PtiToMark = *pThisProfile->pPtiGetCellInProfile(i);
                  m_pRasterGrid->m_Cell[PtiToMark.nGetX()][PtiToMark.nGetY()].SetPolygonID(m_nGlobalPolygonID);
               }

               // If the polygon doesn't meet at a point at its seaward end, also need to rasterize the 'joining line'
               if (! bMeetsAtAPoint)
               {
                  CGeom2DPoint
                     PtUpCoastNormalEnd = *pPrevProfile->pPtGetPointInProfile(nPrevProfileEnd),
                     PtDownCoastNormalEnd = *pThisProfile->pPtGetPointInProfile(nThisProfileEnd);

                  RasterizePolygonJoiningLine(&PtUpCoastNormalEnd, &PtDownCoastNormalEnd);

//                   pPolygon->SetNotPointed();
               }
            }

            nPrevProfileCoastPoint = nCoastPoint;
            nPrevProfile = nThisProfile;
         }
      }
   }

   return RTN_OK;
}


/*===============================================================================================================================

 Puts a polygon 'joining line' (the line which is the seaward boundary of the polygon, if the polygon doesn't meet at a point) onto the raster grid

===============================================================================================================================*/
void CSimulation::RasterizePolygonJoiningLine(CGeom2DPoint const* pPt1, CGeom2DPoint const* pPt2)
{
   // The start point of the line, must convert from the external CRS to grid CRS
   double
      dXStart = dExtCRSXToGridX(pPt1->dGetX()),
      dYStart = dExtCRSYToGridY(pPt1->dGetY());

   // The end point of the line, again convert from the external CRS to grid CRS
   double
      dXEnd = dExtCRSXToGridX(pPt2->dGetX()),
      dYEnd = dExtCRSYToGridY(pPt2->dGetY());

   // Safety check, in case the two points are identical (can happen due to rounding errors)
   if ((bFPIsEqual(dXStart, dXEnd, TOLERANCE)) && (bFPIsEqual(dYStart, dYEnd, TOLERANCE)))
      return;

   // Interpolate between cells by a simple DDA line algorithm, see http://en.wikipedia.org/wiki/Digital_differential_analyzer_(graphics_algorithm) Note that Bresenham's algorithm gave occasional gaps
   double
      dXInc = dXEnd - dXStart,
      dYInc = dYEnd - dYStart,
      dLength = tMax(tAbs(dXInc), tAbs(dYInc));

   dXInc /= dLength;
   dYInc /= dLength;

   double
      dX = dXStart,
      dY = dYStart;

   // Process each interpolated point
   for (int m = 0; m <= static_cast<int>(dRound(dLength)); m++)
   {
      int
         nX = static_cast<int>(dX),
         nY = static_cast<int>(dY);

      // Safety check
      if (! bIsWithinValidGrid(nX, nY))
         KeepWithinValidGrid(dXStart, dYStart, nX, nY);

      // Mark this point on the raster grid
      m_pRasterGrid->m_Cell[nX][nY].SetPolygonID(m_nGlobalPolygonID);

      // And increment for next time
      dX += dXInc;
      dY += dYInc;
   }
}


/*===============================================================================================================================

 Marks cells of the raster grid that are within each coastal polygon. The flood fill code used here is adapted from an example by Lode Vandevenne (http://lodev.org/cgtutor/floodfill.html#Scanline_Floodfill_Algorithm_With_Stack)

===============================================================================================================================*/
void CSimulation::MarkPolygonCells(void)
{
   // Do this for each coast
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      // Do this for every coastal polygon
      for (int nPoly = 0; nPoly < m_VCoast[nCoast].nGetNumPolygons(); nPoly++)
      {
         int
            nCellsInPolygon = 0;
//       double dTotDepth = 0;

         CGeomCoastPolygon* pPolygon = m_VCoast[nCoast].pGetPolygon(nPoly);
         int nPolyID = pPolygon->nGetGlobalID();

         // Create an empty stack
         stack<CGeom2DIPoint> PtiStack;

         // Since the polygon's vector boundary does not coincide exactly with the polygon's raster boundary, and the point-in-polygon check gives an indeterminate result if the point is exactly on the polygon's boundary, for safety we must construct a vector 'inner buffer' which is smaller than, and inside, the vector boundary
         int
            nHand = m_VCoast[nCoast].nGetSeaHandedness(),
            nSize = pPolygon->nGetBoundarySize();
         vector<CGeom2DPoint> PtVInnerBuffer;
         for (int i = 0; i < nSize-1; i++)
         {
            int j = i+1;
            if (i == nSize-2)       // We must ignore the duplicated node point
               j = 0;
            CGeom2DPoint
               PtThis = *pPolygon->pPtGetBoundaryPoint(i),
               PtNext = *pPolygon->pPtGetBoundaryPoint(j),
               PtBuffer = PtGetPerpendicular(&PtThis, &PtNext, m_dCellSide, nHand);

            PtVInnerBuffer.push_back(PtBuffer);
         }

//          // DEBUG STUFF
//          LogStream << endl << "Timestep " << m_ulIteration << ": coast " << nCoast << ", polygon " << nPoly << endl;
//          LogStream << "Boundary\t\t\tBuffer" << endl;
//          for (int i = 0; i < pPolygon->nGetBoundarySize()-1; i++)
//             LogStream << "{" << pPolygon->pPtGetBoundaryPoint(i)->dGetX() << ", " << pPolygon->pPtGetBoundaryPoint(i)->dGetY() << "}\t{" << PtVInnerBuffer[i].dGetX() << ", " << PtVInnerBuffer[i].dGetY() << "}" << endl;
//          LogStream << endl;

         // Now look for a point that is within the polygon as a start point for the flood fill. For a first attempt, calculate the polygon's centroid
         CGeom2DPoint PtStart = pPolygon->PtGetCentroid();

         // Is the centroid within the inner buffer?
         if (! bIsWithinPolygon(&PtStart, &PtVInnerBuffer))
         {
            // No, it is not: the polygon must be a concave polygon. So keep looking for a point which is definitely inside the polygon, using an alternative method
            PtStart = PtFindPointInPolygon(&PtVInnerBuffer, pPolygon->nGetPointInPolygonSearchStartPoint());
         }

         // Safety check
         if (PtStart.dGetX() == DBL_NODATA)
         {
            LogStream << m_ulIteration << ": " << ERR << "could not find a flood fill start point for coast " << nCoast << ", polygon " << nPoly << endl;
            break;
         }

         // We have a flood fill start point which is definitely within the polygon so push this point onto the stack
         CGeom2DIPoint PtiStart = PtiExtCRSToGrid(&PtStart);                // Grid CRS
         PtiStack.push(PtiStart);

//          LogStream << m_ulIteration << ": filling polygon " << nPoly << " from [" << PtiStart.nGetX() << "][" << PtiStart.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiStart.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiStart.nGetY()) << "}" << endl;

         // Then do the flood fill: loop until there are no more cell co-ords on the stack
         while (! PtiStack.empty())
         {
            CGeom2DIPoint Pti = PtiStack.top();
            PtiStack.pop();

            int
               nX = Pti.nGetX(),
               nY = Pti.nGetY();

            while ((nX >= 0) && (m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID() == INT_NODATA))
               nX--;

            nX++;

            bool
               bSpanAbove = false,
               bSpanBelow = false;

            while ((nX < m_nXGridMax) && (m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID() == INT_NODATA))
            {
               // Mark the cell as being in this polygon
               m_pRasterGrid->m_Cell[nX][nY].SetPolygonID(nPolyID);
//                LogStream << "[" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;

               // Increment the running totals for this polygon
               nCellsInPolygon++;
//                dTotDepth += m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth();

               if ((! bSpanAbove) && (nY > 0) && (m_pRasterGrid->m_Cell[nX][nY-1].nGetPolygonID() == INT_NODATA))
               {
                  PtiStack.push(CGeom2DIPoint(nX, nY-1));
                  bSpanAbove = true;
               }

               else if (bSpanAbove && (nY > 0) && (m_pRasterGrid->m_Cell[nX][nY-1].nGetPolygonID() != INT_NODATA))
               {
                  bSpanAbove = false;
               }

               if ((! bSpanBelow) && (nY < m_nYGridMax-1) && (m_pRasterGrid->m_Cell[nX][nY+1].nGetPolygonID() == INT_NODATA))
               {
                  PtiStack.push(CGeom2DIPoint(nX, nY+1));
                  bSpanBelow = true;
               }

               else if (bSpanBelow && (nY < m_nYGridMax-1) && (m_pRasterGrid->m_Cell[nX][nY+1].nGetPolygonID() != INT_NODATA))
               {
                  bSpanBelow = false;
               }

               nX++;
            }
         }

//          // Store the number of cells in the interior of the polygon (note that this is an underestimate, it does not include cells in the polygon boundary)
//          pPolygon->SetNumCells(nCellsInPolygon);
//          LogStream << m_ulIteration << ": N cells = " << nCellsInPolygon << " in polygon " << nPoly << endl;

         // Calculate the total volume of seawater on the polygon (m3) and store it
//          double dSeaVolume = dTotDepth * m_dCellSide;
//          pPolygon->SetSeawaterVolume(dSeaVolume);
      }
   }
}


/*===============================================================================================================================

 For between-polygon potential sediment routing: find which are the adjacent polygons, and calc the length of the shared normal between this polygon and the adjacent polygons

 // TODO Will need to change this when length of coastline-normal profiles (and so polygon seaward length) is determined by depth of closure

===============================================================================================================================*/
void CSimulation::DoPolygonSharedBoundaries(void)
{
   // Do this for each coast
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      // Do this for every coastal polygon
      int nNumPolygons = m_VCoast[nCoast].nGetNumPolygons();
      for (int nPoly = 0; nPoly < nNumPolygons; nPoly++)
      {
         CGeomCoastPolygon* pPolygon = m_VCoast[nCoast].pGetPolygon(nPoly);

         double
            dUpCoastTotBoundaryLen = 0,
            dDownCoastTotBoundaryLen = 0;

         vector<int>
            nVUpCoastAdjacentPolygon,
            nVDownCoastAdjacentPolygon;

         vector<double>
            dVUpCoastBoundaryShare,
            dVDownCoastBoundaryShare;

         // Do first for the down-coast profile, then for the up-coast profile
         for (int nDirection = DIRECTION_DOWNCOAST; nDirection <= DIRECTION_UPCOAST; nDirection++)
         {
            // Are we at the start of the coastline, leaving it by going up-coast?
            if ((nPoly == 0) && (nDirection == DIRECTION_UPCOAST))
            {
               // No other polygon is adjacent to the up-coast profile of the start-of-coast polygon
               nVUpCoastAdjacentPolygon.push_back(INT_NODATA);
               dVUpCoastBoundaryShare.push_back(1);

               // Store in the polygon
               pPolygon->SetUpCoastAdjacentPolygons(&nVUpCoastAdjacentPolygon);
               pPolygon->SetUpCoastAdjacentPolygonBoundaryShares(&dVUpCoastBoundaryShare);

               continue;
            }

            // Are we at the end of the coastline, leaving it by going down-coast?
            if ((nPoly == nNumPolygons-1) && (nDirection == DIRECTION_DOWNCOAST))
            {
               // No other polygon is adjacent to the down-coast profile of the end-of-coast polygon
               nVDownCoastAdjacentPolygon.push_back(INT_NODATA);
               dVDownCoastBoundaryShare.push_back(1);

               // Store in the polygon
               pPolygon->SetDownCoastAdjacentPolygons(&nVDownCoastAdjacentPolygon);
               pPolygon->SetDownCoastAdjacentPolygonBoundaryShares(&dVDownCoastBoundaryShare);

               continue;
            }

            // We are not leaving the coastline from one end or other, so there must be at least one other polygon on either side of this one
            int
               nProfile,
               nPointsInProfile;

            if (nDirection == DIRECTION_UPCOAST)
            {
               // Going up-coast
               nProfile = pPolygon->nGetUpCoastProfile();
               nPointsInProfile = pPolygon->nGetUpCoastProfileNumPointsUsed();
            }
            else
            {
               // Going down-coast
               nProfile = pPolygon->nGetDownCoastProfile();
               nPointsInProfile = pPolygon->nGetDownCoastProfileNumPointsUsed();
            }

            CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfile(nProfile);

            for (int nPoint = 0; nPoint < nPointsInProfile-1; nPoint++)
            {
               CGeom2DPoint
                  PtStart = *pProfile->pPtGetPointInProfile(nPoint),
                  PtEnd = *pProfile->pPtGetPointInProfile(nPoint+1);

               // Calculate the length of this segment of the normal profile. Note that it should not be zero, since we checked for duplicate points when creating profiles
               double dDistBetween = dGetDistanceBetween(&PtStart, &PtEnd);

               // Find out what polygons are adjacent
               int nCoincidentProfiles = pProfile->nGetNumCoincidentProfilesInLineSegment(nPoint);
               if (nDirection == DIRECTION_UPCOAST)
               {
                  int nAdj = nPoly - nCoincidentProfiles;

                  // Safety check
                  if (nAdj >= 0)
                     nVUpCoastAdjacentPolygon.push_back(nAdj);
               }
               else
               {
                  int nAdj = nPoly + nCoincidentProfiles;

                  // Safety check
                  if (nAdj < nNumPolygons)
                     nVDownCoastAdjacentPolygon.push_back(nPoly + nCoincidentProfiles);
               }

               // Store the line segment data
               if (nDirection == DIRECTION_UPCOAST)
               {
                  dUpCoastTotBoundaryLen += dDistBetween;
                  dVUpCoastBoundaryShare.push_back(dDistBetween);
               }
               else
               {
                  dDownCoastTotBoundaryLen += dDistBetween;
                  dVDownCoastBoundaryShare.push_back(dDistBetween);
               }
            }

//          assert(dVUpCoastBoundaryShare.size() == nVUpCoastAdjacentPolygon.size());
//          assert(dVDownCoastBoundaryShare.size() == nVDownCoastAdjacentPolygon.size());

//          LogStream << m_ulIteration << ": polygon = " << nPoly << (pPolygon->bIsPointed() ? " IS TRIANGULAR" : "") << endl;
//             LogStream << m_ulIteration << ": polygon = " << nPoly << endl;
//
//          LogStream << "\tUP-COAST Boundary lengths = ";
//          for (unsigned int n = 0; n < dVUpCoastBoundaryShare.size(); n++)
//             LogStream << dVUpCoastBoundaryShare[n] << " ";
//          LogStream << endl;
//          LogStream << "\tTotal UP-COAST boundary length = " << dUpCoastTotBoundaryLen << endl;
//             LogStream << "\tUP-COAST Adjacent polygons = ";
//             for (unsigned int n = 0; n < nVUpCoastAdjacentPolygon.size(); n++)
//                LogStream << nVUpCoastAdjacentPolygon[n] << " ";
//             LogStream << endl;
//
//          LogStream << "\tDOWN-COAST Boundary lengths = ";
//          for (unsigned int n = 0; n < dVDownCoastBoundaryShare.size(); n++)
//             LogStream << dVDownCoastBoundaryShare[n] << " ";
//          LogStream << endl;
//          LogStream << "\tTotal DOWN-COAST boundary length = " << dDownCoastTotBoundaryLen << endl;
//             LogStream << "\tDOWN-COAST Adjacent polygons = ";
//             for (unsigned int n = 0; n < nVDownCoastAdjacentPolygon.size(); n++)
//                LogStream << nVDownCoastAdjacentPolygon[n] << " ";
//             LogStream << endl;

            // Calculate the fraction of the total boundary shared with each adjacent polygon
            if (nDirection == DIRECTION_UPCOAST)
            {
               for (unsigned int n = 0; n < dVUpCoastBoundaryShare.size(); n++)
                  dVUpCoastBoundaryShare[n] /= dUpCoastTotBoundaryLen;

               // Store in the polygon
               pPolygon->SetUpCoastAdjacentPolygons(&nVUpCoastAdjacentPolygon);
               pPolygon->SetUpCoastAdjacentPolygonBoundaryShares(&dVUpCoastBoundaryShare);
            }
            else
            {
               for (unsigned int n = 0; n < dVDownCoastBoundaryShare.size(); n++)
                  dVDownCoastBoundaryShare[n] /= dDownCoastTotBoundaryLen;

               // Store in the polygon
               pPolygon->SetDownCoastAdjacentPolygons(&nVDownCoastAdjacentPolygon);
               pPolygon->SetDownCoastAdjacentPolygonBoundaryShares(&dVDownCoastBoundaryShare);

               // Finally, calculate the distance between the coast node and the antinode of the polygon
               double dPolygonSeawardLen = dGetDistanceBetween(pPolygon->pPtiGetNode(), pPolygon->pPtiGetAntinode());

               // And store it
               m_VCoast[nCoast].AppendPolygonLength(dPolygonSeawardLen);
            }
         }
      }
   }
}


/*===============================================================================================================================

 Determines whether a point is within a polygon: however if the point is exactly on the edge of the polygon, then the result is indeterminate

 Modified from code at http://alienryderflex.com/polygon/, our thanks to Darel Rex Finley (DarelRex@gmail.com)

===============================================================================================================================*/
bool CSimulation::bIsWithinPolygon(CGeom2DPoint const* pPtStart, vector<CGeom2DPoint> const* pPtPoints)
{
   bool bOddNodes = false;

   int
      nPolyCorners = pPtPoints->size(),
      j = nPolyCorners-1;

   double
      dX = pPtStart->dGetX(),
      dY = pPtStart->dGetY();

   for (int i = 0; i < nPolyCorners; i++)
   {
      double
         dCorneriX = pPtPoints->at(i).dGetX(),
         dCorneriY = pPtPoints->at(i).dGetY(),
         dCornerjX = pPtPoints->at(j).dGetX(),
         dCornerjY = pPtPoints->at(j).dGetY();

      if ((dCorneriY < dY && dCornerjY >= dY) || (dCornerjY < dY && dCorneriY >= dY))
      {
         if (dCorneriX + (dY - dCorneriY) / (dCornerjY - dCorneriY) * (dCornerjX - dCorneriX) < dX)
         {
            bOddNodes = ! bOddNodes;
         }
      }

      j = i;
   }

  return bOddNodes;
}


/*===============================================================================================================================

 Finds a point in a polygon: is guaranteed to succeed, as every strictly closed polygon has at least one triangle that is completely contained within the polygon

 Derived from an algorithm at http://stackoverflow.com/questions/9797448/get-a-point-inside-the-polygon

===============================================================================================================================*/
CGeom2DPoint CSimulation::PtFindPointInPolygon(vector<CGeom2DPoint> const* pPtPoints, int const nStartPoint)
{
   int
      nPolySize = pPtPoints->size(),
      nOffSet = 0;
   CGeom2DPoint PtStart;

   do
   {
      // Choose three consecutive points from the polygon
      vector <CGeom2DPoint> nVTestPoints;
      for (int n = 0; n < 3; n++)
      {
         int nIndex = n + nStartPoint + nOffSet;
         if (nIndex > nPolySize-1)
            nIndex -= nPolySize;

         // Safety check
         if (nIndex < 0)
            return CGeom2DPoint(DBL_NODATA, DBL_NODATA);

         nVTestPoints.push_back(pPtPoints->at(nIndex));
      }

      // Increment ready for next time
      nOffSet++;

      // Safety check
      if (nOffSet >= (nPolySize + 3))
         return CGeom2DPoint(DBL_NODATA, DBL_NODATA);

      // Check if the halfway point between the first and the third point is inside the polygon
      PtStart = PtAverage(&nVTestPoints[0], &nVTestPoints[2]);
   }
   while (! bIsWithinPolygon(&PtStart, pPtPoints));

   return PtStart;
}
