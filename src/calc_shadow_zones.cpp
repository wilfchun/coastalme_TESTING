/*!
 *
 * \file calc_shadow_zones.cpp
 * \brief Locates shadow zones, is part of wave propagation calculations
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
#include <cmath>

#include <iostream>
using std::cout;
using std::endl;

#include <stack>
using std::stack;

#include "cme.h"
#include "coast.h"
#include "simulation.h"



/*===============================================================================================================================
 
 Determines whether the wave orientation is onshore or offshore, and up-coast or down-coast
 
===============================================================================================================================*/
bool CSimulation::bOnOrOffShoreAndUpOrDownCoast(double const dCoastAngle, double const dWaveOrientation, int const nSeaHand, bool& bDownCoast)
{
   bool bOnShore;
   double dWaveToCoastAngle = fmod((dWaveOrientation - dCoastAngle + 360), 360);
   
   bDownCoast = ((dWaveToCoastAngle > 270) || (dWaveToCoastAngle < 90)) ? true : false;
   
   if (nSeaHand == RIGHT_HANDED)
   {
      // The sea is on the RHS travelling down-coast
      bOnShore = dWaveToCoastAngle > 180 ? true : false;
   }         
   else
   {
      // The sea is on the LHS travelling down-coast
      bOnShore = dWaveToCoastAngle > 180 ? false : true;
   }
   
   return bOnShore;
}


/*===============================================================================================================================
 
 Given a cell and a wave orientation, finds the 'upwave' cell
 
===============================================================================================================================*/
CGeom2DIPoint CSimulation::PtiFollowWaveOrientation(CGeom2DIPoint const* pPtiLast, double const dWaveOrientationIn, double& dCorrection)
{
   int
      nXLast = pPtiLast->nGetX(),
      nYLast = pPtiLast->nGetY(),
      nXNext = nXLast, 
      nYNext = nYLast;
      
   double dWaveOrientation = dWaveOrientationIn - dCorrection;   
      
   if (dWaveOrientation < 22.5)
   {
      nYNext--;
      dCorrection = 22.5 - dWaveOrientation;
   }
   else if (dWaveOrientation < 67.5)
   {
      nYNext--;
      nXNext++;      
      dCorrection = 67.5 - dWaveOrientation;
   }
   else if (dWaveOrientation < 112.5)
   {
      nXNext++;
      dCorrection = 112.5 - dWaveOrientation;
   }
   else if (dWaveOrientation < 157.5)
   {
      nXNext++;
      nYNext++;
      dCorrection = 157.5 - dWaveOrientation;
   }
   else if (dWaveOrientation < 202.5)
   {
      nYNext++;
      dCorrection = 202.5 - dWaveOrientation;
   }
   else if (dWaveOrientation < 247.5)
   {
      nXNext--;
      nYNext++;
      dCorrection = 247.5 - dWaveOrientation;
   }
   else if (dWaveOrientation < 292.5)
   {
      nXNext--;
      dCorrection = 292.5 - dWaveOrientation;
   }
   else if (dWaveOrientation < 337.5)
   {
      nXNext--;
      nYNext--;
      dCorrection = 337.5 - dWaveOrientation;
   }
   else
   {
      nYNext--;
      dCorrection = 22.5 - dWaveOrientation;
   }
   
   dCorrection = dKeepWithin360(dCorrection);
   
   return CGeom2DIPoint(nXNext, nYNext);   
}


/*===============================================================================================================================

 Finds wave shadow zones and modifies waves in and near them. Note that where up-coast and down-coast shadow zones overlap, the effects on wave values in the overlap area is an additive decrease in wave energy. Changes to wave energy in any down-drift increased-energy zones are also additive.  

===============================================================================================================================*/
int CSimulation::nDoAllShadowZones(void)
{
   // Do this once for each coastline
   for (unsigned int nCoast = 0; nCoast < m_VCoast.size(); nCoast++)
   {
      // =========================================================================================================================
      // The first stage: find coastline start points for possible shadow zone boundaries by sweeping the coastline: first down-coast then up-coast
      int 
         nSeaHand = m_VCoast[nCoast].nGetSeaHandedness(),
         nCoastSize = m_VCoast[nCoast].nGetCoastlineSize();
         
      vector<int> VnShadowZoneCoastPoint;         
      
      for (bool bDownCoast : {true, false})
      {
         if (bDownCoast)
         {
            bool bLastDownCoastAndOnshore = false;
            
            // Work along coast in down-coast direction
            for (int nCoastPoint = 0; nCoastPoint < nCoastSize; nCoastPoint++)
            {
               double dCurvature = m_VCoast[nCoast].dGetDetailedCurvature(nCoastPoint);
               if (dCurvature < 0)
               {
                  // OK, the coast is convex here, now get the flux orientation (a tangent to the coastline)
                  double dFluxOrientation = m_VCoast[nCoast].dGetFluxOrientation(nCoastPoint);
                     
                  // If this coast point is in the active zone, use the breaking wave orientation, otherwise use the deep water wave orientation   
                  double dWaveOrientation;   
                  if (m_VCoast[nCoast].dGetDepthOfBreaking(nCoastPoint) == DBL_NODATA)
                     // Not in active zone
                     dWaveOrientation = m_VCoast[nCoast].dGetDeepWaterWaveOrientation(nCoastPoint);
                  else
                     // In active zone
                     dWaveOrientation = m_VCoast[nCoast].dGetBreakingWaveOrientation(nCoastPoint);
                  
                  // Are waves on- or off-shore, and up- or down-coast?
                  bool bDownCoast = false;   
                  bool bOnShore = bOnOrOffShoreAndUpOrDownCoast(dFluxOrientation, dWaveOrientation, nSeaHand, bDownCoast);
                  
//                   CGeom2DPoint PtTmp1 = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint);
//                   LogStream << m_ulIteration << ": going down-coast, coast point (" << nCoastPoint << ") at {" << PtTmp1.dGetX() << ", " << PtTmp1.dGetY() << "} has " << (bDownCoast ? "down-coast " : "up-coast ") << (bOnShore ? "on-shore" : "off-shore") << " waves, dWaveOrientation = " << dWaveOrientation << " dFluxOrientation = " << dFluxOrientation << endl;
                  
                  if (bDownCoast && (! bOnShore))
                  {
                     // Waves are down-coast and off-shore
//                      CGeom2DPoint PtTmp1 = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint);
//                      LogStream << m_ulIteration << ": going down-coast, waves have off-shore and down-coast component at coast point (" << nCoastPoint << ") at {" << PtTmp1.dGetX() << ", " << PtTmp1.dGetY() << "}" << endl;
                     
                     // If the previous coast point had waves which were down-coast and on-shore, then this could be the boundary of a shadow zone
                     if (bLastDownCoastAndOnshore)
                     {
                        VnShadowZoneCoastPoint.push_back(nCoastPoint);
                        bLastDownCoastAndOnshore = false;
                        
                        CGeom2DPoint PtTmp = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint);
                        LogStream << m_ulIteration << ": possible shadow boundary start found at {" << PtTmp.dGetX() << ", " << PtTmp.dGetY() << "} while going down-coast, this is coast point (" << nCoastPoint << ")" << endl;
                     }
                  }
                  else if (bDownCoast && bOnShore)
                  {
                     bLastDownCoastAndOnshore = true;
                  }
                  else
                  {
                     bLastDownCoastAndOnshore = false;
                  }
               }
            }            
         }
      
         else
         {
            // Moving up-coast
            bool bLastUpCoastAndOnshore = false;

            // Work along coast in up-coast direction
            for (int nCoastPoint = nCoastSize-1; nCoastPoint >= 0; nCoastPoint--)
            {
               double dCurvature = m_VCoast[nCoast].dGetDetailedCurvature(nCoastPoint);
               if (dCurvature < 0)
               {
                  // OK, the coast is convex here, now get the flux orientation (a tangent to the coastline)
                  double dFluxOrientation = m_VCoast[nCoast].dGetFluxOrientation(nCoastPoint);
                  
                  // If this coast point is in the active zone, use the breaking wave orientation, otherwise use the deep water wave orientation   
                  double dWaveOrientation;   
                  if (m_VCoast[nCoast].dGetDepthOfBreaking(nCoastPoint) == DBL_NODATA)
                     // Not in active zone
                     dWaveOrientation = m_VCoast[nCoast].dGetDeepWaterWaveOrientation(nCoastPoint);
                  else
                     // In active zone
                     dWaveOrientation = m_VCoast[nCoast].dGetBreakingWaveOrientation(nCoastPoint);
                  
                  // Are waves on- or off-shore, and up- or down-coast?
                  bool bDownCoast = false;   
                  bool bOnShore = bOnOrOffShoreAndUpOrDownCoast(dFluxOrientation, dWaveOrientation, nSeaHand, bDownCoast);
                  
//                   CGeom2DPoint PtTmp1 = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint);
//                   LogStream << m_ulIteration << ": going up-coast, coast point (" << nCoastPoint << ") at {" << PtTmp1.dGetX() << ", " << PtTmp1.dGetY() << "} has " << (bDownCoast ? "down-coast " : "up-coast ") << (bOnShore ? "on-shore" : "off-shore") << " waves, dWaveOrientation = " << dWaveOrientation << " dFluxOrientation = " << dFluxOrientation << endl;
                  
                  if ((! bDownCoast) && (! bOnShore))
                  {
                     // Waves are up-coast and off-shore
//                      CGeom2DPoint PtTmp1 = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint);
//                      LogStream << m_ulIteration << ": going up-coast, waves have off-shore and down-coast component at coast point (" << nCoastPoint << ") at {" << PtTmp1.dGetX() << ", " << PtTmp1.dGetY() << "}" << endl;
                     
                     // If the previous coast point had waves which were up-coast and on-shore, then this could be the boundary of a shadow zone
                     if (bLastUpCoastAndOnshore)
                     {
                        VnShadowZoneCoastPoint.push_back(nCoastPoint);
                        bLastUpCoastAndOnshore = false;
                        
                        CGeom2DPoint PtTmp = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint);
                        LogStream << m_ulIteration << ": possible shadow boundary start found at {" << PtTmp.dGetX() << ", " << PtTmp.dGetY() << "} while going up-coast, this is coast point (" << nCoastPoint << ")" << endl;
                     }
                  }
                  else if ((! bDownCoast) && bOnShore)
                  {
                     bLastUpCoastAndOnshore = true;
                  }
                  else
                  {
                     bLastUpCoastAndOnshore = false;
                  }
               }
            }            
         }
      }
      
      // =========================================================================================================================
      // The second stage: we have a list of possible shadow zone start points, trace each of these 'up-wave' to identify valid shadow zones
      vector<CGeomILine> VILShadowBoundary;
      vector<int> VnShadowBoundaryEndCoastPoint;
      
      for (unsigned int nStartPoint = 0; nStartPoint < VnShadowZoneCoastPoint.size(); nStartPoint++)
      {
         LogStream << m_ulIteration << ": processing shadow boundary start point " << nStartPoint << " (" << VnShadowZoneCoastPoint.size() << " total)" << endl;

         bool bHitSea = false;
         int nShadowBoundaryCoastPoint = -1;
         
         // From each start point, follow the wave direction
         CGeomILine ILShadowBoundary;
         
         bool
            bHitEdge = false,
            bHitCoast = false;
            
         // If this coast point is in the active zone, start with the breaking wave orientation, otherwise use the deep water wave orientation   
         double dPrevWaveOrientation;   
         if (m_VCoast[nCoast].dGetDepthOfBreaking(nStartPoint) == DBL_NODATA)
         {
            // Not in active zone
            dPrevWaveOrientation = m_VCoast[nCoast].dGetDeepWaterWaveOrientation(nStartPoint);
         }
         else
         {
            // In active zone
            dPrevWaveOrientation = m_VCoast[nCoast].dGetBreakingWaveOrientation(nStartPoint);
         }
            
         CGeom2DIPoint PtiPrev = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(VnShadowZoneCoastPoint[nStartPoint]);
         ILShadowBoundary.Append(&PtiPrev);
         
         int nDist = 0;         
         double 
            dTotOrientation = 0,
            dCorrection = 0;         
         while ((! bHitEdge) && (! bHitCoast))
         {
            if (nDist > 0)
            {
               int
                  nXPrev = PtiPrev.nGetX(),
                  nYPrev = PtiPrev.nGetY();               
               
               if (! m_pRasterGrid->m_Cell[nXPrev][nYPrev].bIsInActiveZone())
               {
                  // The previous cell was outside the active zone, so use its wave orientation value
                  dPrevWaveOrientation = m_pRasterGrid->m_Cell[nXPrev][nYPrev].dGetWaveOrientation();
               }
               else
               {
                  // The previous cell was in the active zone
                  if (bHitSea)
                  {
                     // If this shadow boundary has already hit sea, then we must be getting near a coast: use the average-so-far wave orientation
                     double dAvgOrientationSoFar = dTotOrientation / (nDist + 1);
                     dPrevWaveOrientation = dAvgOrientationSoFar;                     
                  }
                  else
                  {
                     // This shadow boundary has not already hit sea, just use the wave orientation from the previous cell
                     dPrevWaveOrientation = m_pRasterGrid->m_Cell[nXPrev][nYPrev].dGetWaveOrientation();
                     LogStream << m_ulIteration << ": not already hit sea, using previous cell's wave orientation for cell [" << nXPrev << "][" << nYPrev << "] = {" << dGridCentroidXToExtCRSX(nXPrev) << ", " << dGridCentroidYToExtCRSY(nYPrev) << "}" << endl;  
                  }                     
               }
               
               if (dPrevWaveOrientation == DBL_NODATA)
               {
                  LogStream << m_ulIteration << ": dPrevWaveOrientation == DBL_NODATA for cell [" << nXPrev << "][" << nYPrev << "] = {" << dGridCentroidXToExtCRSX(nXPrev) << ", " << dGridCentroidYToExtCRSY(nYPrev) << "}" << endl;                    
                   
                  if (! m_pRasterGrid->m_Cell[nXPrev][nYPrev].bIsInContiguousSea())
                  {
                     // The previous cell was an inland cell, so use the deep water wave orientation
                     dPrevWaveOrientation = m_pRasterGrid->m_Cell[nXPrev][nYPrev].dGetDeepWaterWaveOrientation();
                  }
                  else
                  {
                     double dAvgOrientationSoFar = dTotOrientation / (nDist + 1);
                     dPrevWaveOrientation = dAvgOrientationSoFar;
                  }                  
               }
               
               dTotOrientation += dPrevWaveOrientation;                              
            }
            
            // Go upwave along the previous cell's wave orientation to find the new boundary cell
            CGeom2DIPoint PtiNew = PtiFollowWaveOrientation(&PtiPrev, dPrevWaveOrientation, dCorrection);
            
            // Get the co-ordinates of 'this' cell
            int
               nX = PtiNew.nGetX(),
               nY = PtiNew.nGetY();
            
            // Have we hit the edge of the valid part of the grid?
            if ((m_pRasterGrid->m_Cell[nX][nY].bIsEdgeCell()) || (! bIsWithinValidGrid(&PtiNew)))
            {
               bHitEdge = true;
               
//                LogStream << m_ulIteration << ": shadow boundary " << nStartPoint << " hit edge cell at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;  
               
               continue;               
            }
            
            // Have we hit a sea cell yet?
            if ((nDist > MAX_LAND_LENGTH_OF_SHADOW_ZONE_LINE) && (! bHitSea))
            {
               if (m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea())
                  bHitSea = true;
               else
               {
                  if (nDist >= MAX_LAND_LENGTH_OF_SHADOW_ZONE_LINE)
                     // If we have travelled MAX_LAND_LENGTH_OF_SHADOW_ZONE_LINE cells without hitting sea, then abandon this shadow boundary
                     break;
               }
            }

            // Store the co-ords of every cell which we cross
            ILShadowBoundary.Append(&PtiNew);
            
//             LogStream << m_ulIteration << ": at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;   
            
            // Having hit sea, have we now hit we hit a coast point? Note that two diagonal(ish) raster lines can cross each other without any intersection, so must also test an adjacent cell for intersection (does not matter which adjacent cell)
            if (bHitSea)
            {
               if (m_pRasterGrid->m_Cell[nX][nY].bIsCoastline() || (bIsWithinValidGrid(nX, nY+1) && m_pRasterGrid->m_Cell[nX][nY+1].bIsCoastline()))
               {
                  bHitCoast = true;
                  
//                   LogStream << m_ulIteration << ": shadow boundary " << nStartPoint << " hit coast at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;  
               }
            }
            
            // For next time
            PtiPrev = PtiNew;
            nDist++;
         }

         if (bHitCoast)
         {
            // The shadow zone boundary has hit a coast, but is the shadow zone line trivially short?
            double dShadowLen = dGetDistanceBetween(&ILShadowBoundary[0], &ILShadowBoundary.Back()) * m_dCellSide;
            if (dShadowLen < MIN_LENGTH_OF_SHADOW_ZONE_LINE)
            {
               // Too short, so forget about it
               LogStream << m_ulIteration << ": shadow boundary " << nStartPoint << " is too short, has length " << dShadowLen << " m but the minimum length is " << MIN_LENGTH_OF_SHADOW_ZONE_LINE << " m. It starts at [" << ILShadowBoundary[0].nGetX() << "][" << ILShadowBoundary[0].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILShadowBoundary[0].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILShadowBoundary[0].nGetY()) << "} and hits coast at [" << ILShadowBoundary.Back().nGetX() << "][" << ILShadowBoundary.Back().nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILShadowBoundary.Back().nGetX()) << ", " << dGridCentroidYToExtCRSY(ILShadowBoundary.Back().nGetY()) << "}" << endl;
               
               continue;
            }
                
            // We've found a valid shadow zone. Check the last point in the shadow boundary. Note that occasionally this last cell is not 'above' a cell but one of its neighbouring cells is: in which case, replace the last point in the shadow boundary with the co-ords of this neighbouring cell
            nShadowBoundaryCoastPoint = m_VCoast[nCoast].nGetCoastPointGivenCell(&ILShadowBoundary.Back()); 
            if (nShadowBoundaryCoastPoint == INT_NODATA)
            {
               // Could not find a neighbouring cell which is 'under' the coastline
               LogStream << m_ulIteration << ": no coast point under {" << dGridCentroidXToExtCRSX(ILShadowBoundary.Back().nGetX()) << ", " << dGridCentroidYToExtCRSY(ILShadowBoundary.Back().nGetY()) << "}" << endl;
               
               return RTN_ERR_NO_CELL_UNDER_COASTLINE;               
            }
            
            // Now store the shadow zone boundary information
            VILShadowBoundary.push_back(ILShadowBoundary);            
            VnShadowBoundaryEndCoastPoint.push_back(nShadowBoundaryCoastPoint);
            
            LogStream << m_ulIteration << ": shadow boundary " << nStartPoint << " has created a valid shadow zone. It starts at [" << ILShadowBoundary[0].nGetX() << "][" << ILShadowBoundary[0].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILShadowBoundary[0].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILShadowBoundary[0].nGetY()) << "} and hits coast at [" << ILShadowBoundary.Back().nGetX() << "][" << ILShadowBoundary.Back().nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILShadowBoundary.Back().nGetX()) << ", " << dGridCentroidYToExtCRSY(ILShadowBoundary.Back().nGetY()) << "} which is coast point (" << nShadowBoundaryCoastPoint << "). This will be shadow zone " << VnShadowBoundaryEndCoastPoint.size()-1 << endl;
         }
         
         if (bHitEdge)
         {
            if (CREATE_SHADOW_ZONE_IF_HITS_GRID_EDGE)
            {
               // We are creating shadow zones if we hit the grid edge. But is the shadow zone line trivially short?
               double dShadowLen = dGetDistanceBetween(&ILShadowBoundary[0], &ILShadowBoundary.Back()) * m_dCellSide;
               if (dShadowLen < MIN_LENGTH_OF_SHADOW_ZONE_LINE)
               {
                  // Too short, so forget about it
                  LogStream << m_ulIteration << ": shadow boundary " << nStartPoint << " is too short, has length " << dShadowLen << " m but the minimum length is " << MIN_LENGTH_OF_SHADOW_ZONE_LINE << " m. It starts at [" << ILShadowBoundary[0].nGetX() << "][" << ILShadowBoundary[0].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILShadowBoundary[0].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILShadowBoundary[0].nGetY()) << "} and hits the grid edge at [" << ILShadowBoundary.Back().nGetX() << "][" << ILShadowBoundary.Back().nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILShadowBoundary.Back().nGetX()) << ", " << dGridCentroidYToExtCRSY(ILShadowBoundary.Back().nGetY()) << "}" << endl;
                  
                  break;
               }
                  
               // We've found a valid grid-edge shadow zone, so store it
               VILShadowBoundary.push_back(ILShadowBoundary);
               
               // TODO improve this
               // We need a distance (in cells) between the shadow boundary start and the 'virtual' shadow boundary end: this is the off-grid point where the shadow boundary would have intersected the coastline, if the grid were big enough. This is of course unknowable. So as a best guess, we choose the shorter of the two distances between the point where the shadow boundary hits the valid edge of the grid, and the start or end of the coast
               int nCoastSize = m_VCoast[nCoast].nGetCoastlineSize();
               CGeom2DIPoint PtiCoastStart = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(0);
               CGeom2DIPoint PtiCoastEnd = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastSize-1);
               
               int nDistance = tMin(dGetDistanceBetween(&ILShadowBoundary.Back(), &PtiCoastStart), dGetDistanceBetween(&ILShadowBoundary.Back(), &PtiCoastEnd));
               VnShadowBoundaryEndCoastPoint.push_back(nDistance);               
                              
               LogStream << m_ulIteration << ": shadow boundary " << nStartPoint << " has created a valid shadow zone. It starts at [" << ILShadowBoundary[0].nGetX() << "][" << ILShadowBoundary[0].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILShadowBoundary[0].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILShadowBoundary[0].nGetY()) << "} and hit the grid edge at [" << ILShadowBoundary.Back().nGetX() << "][" << ILShadowBoundary.Back().nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILShadowBoundary.Back().nGetX()) << ", " << dGridCentroidYToExtCRSY(ILShadowBoundary.Back().nGetY()) << "}. The best-guess length of the shadow boundary is " << nDistance << " cells. This will be shadow zone " << VnShadowBoundaryEndCoastPoint.size()-1 << endl;
            }
            else
            {
               // We are not creating shadow zones if we hit the grid edge
               LogStream << m_ulIteration << ": shadow boundary " << nStartPoint << " hits grid edge: ignored. It starts at [" << ILShadowBoundary[0].nGetX() << "][" << ILShadowBoundary[0].nGetY() << "]" << endl;
            }
         }
      }
         
      // =========================================================================================================================
      // The third stage: store the shadow zone boundary, flood fill the shadow zone, then change wave properties by sweeping the shadow zone and the area downdrift from the shadow zone
      for (unsigned int nZone = 0; nZone < VILShadowBoundary.size(); nZone++)
      {
         LogStream << m_ulIteration << ": processing shadow zone " << nZone << " (" << VILShadowBoundary.size() << " total)" << endl;
         
         int nShadowLineLen = VILShadowBoundary[nZone].nGetSize();
         
         // The pre-allocated vector shadow boundary (external CRS), will be in reverse sequence (i.e. start point is last)
         CGeomLine LBoundary(nShadowLineLen-1);
         
         for (int nn = 0; nn < nShadowLineLen; nn++)
         {
            int
               nTmpX = VILShadowBoundary[nZone][nn].nGetX(),
               nTmpY = VILShadowBoundary[nZone][nn].nGetY();
               
            // Mark the cells as shadow zone boundary   
            m_pRasterGrid->m_Cell[nTmpX][nTmpY].SetShadowZoneBoundary();
            
            // Replace value in the vector shadow boundary
            LBoundary[nShadowLineLen - nn - 1] = CGeom2DPoint(dGridCentroidXToExtCRSX(nTmpX), dGridCentroidYToExtCRSY(nTmpY));
            
            // If this is a sea cell, mark the shadow zone boundary cell as being in the shadow zone, but not yet processed (a -ve number)
            if (m_pRasterGrid->m_Cell[nTmpX][nTmpY].bIsInContiguousSea())
               m_pRasterGrid->m_Cell[nTmpX][nTmpY].SetShadowZoneNumber(-(nZone+1));
            
//                LogStream << m_ulIteration << ": shadow zone " << nZone << ", which starts at [" << nTmpX << "][" << nTmpY << "] = {" << dGridCentroidXToExtCRSX(nTmpX) << ", " << dGridCentroidYToExtCRSY(nTmpY) << "} has cell [" << nTmpX << "][" << nTmpY << "] marked as shadow zone boundary" << endl;
         }

         // Store the shadow zone boundary (external CRS), with the end point first            
         m_VCoast[nCoast].AppendShadowZoneBoundary(LBoundary);
         
         int
            nStartX = VILShadowBoundary[nZone][0].nGetX(),
            nStartY = VILShadowBoundary[nZone][0].nGetY(),
            nEndX = VILShadowBoundary[nZone][nShadowLineLen-1].nGetX(),
            nEndY = VILShadowBoundary[nZone][nShadowLineLen-1].nGetY();
            
         // Grid CRS   
         CGeom2DIPoint
            PtiStart(nStartX, nStartY),
            PtiEnd(nEndX, nEndY);
            
         // Flood fill the shadow zone
         int nRet = nFloodFillShadowZone(nCoast, nZone, VnShadowZoneCoastPoint[nZone], &PtiStart, VnShadowBoundaryEndCoastPoint[nZone], &PtiEnd);
         if (nRet != RTN_OK)
         {
            // Could not do a flood fill of the shadow zone
            if (nRet == RTN_ERR_SHADOW_ZONE_FLOOD_START_POINT)
            {
               // Could not find start point for flood fill. How serious this is depends on the length of the shadow zone line
               if (nShadowLineLen < MAX_LEN_SHADOW_LINE_TO_IGNORE)
               {
                  LogStream << m_ulIteration << ": " << WARN << "could not find start point for flood fill of shadow zone " << nZone << " but continuing simulation because this is a small shadow zone (shadow line length = " << nShadowLineLen << " cells)" << endl;
                  
                  continue;
               }
               else
               {
                  LogStream << m_ulIteration << ": " << ERR << "could not find start point for flood fill of shadow zone " << nZone << " (shadow line length = " << nShadowLineLen << " cells)" << endl;
                  return nRet;
               }
            }
            else
               return nRet;
         }
         
         // Sweep the shadow zone, changing wave orientation and height
         int nLengthOfSweep = 0;
         nRet = nSweepShadowZone(nCoast, nZone, VnShadowZoneCoastPoint[nZone], &PtiStart, VnShadowBoundaryEndCoastPoint[nZone], &PtiEnd, nLengthOfSweep);
         if (nRet != RTN_OK)
            return nRet;
         
         // Sweep the coast downdrift of the shadow zone, changing wave height in order to conserve energy
         nRet = nSweepDownDriftFromShadowZone(nCoast, nZone, VnShadowZoneCoastPoint[nZone], &PtiStart, VnShadowBoundaryEndCoastPoint[nZone], nLengthOfSweep);
         if (nRet != RTN_OK)
            return nRet;
      }
   }

   return RTN_OK;
}


/*===============================================================================================================================

 Flood fills a shadow zone

===============================================================================================================================*/
int CSimulation::nFloodFillShadowZone(int const nCoast, int const nZone, int const nShadowBoundaryStartPoint, CGeom2DIPoint const* pPtiStart, int const nCoastPoint, CGeom2DIPoint const* pPtiEnd)
{
   int
      nCoastSeaHand = m_VCoast[nCoast].nGetSeaHandedness(),
      nShadowZoneCoastToCapeSeaHand;

   if (nShadowBoundaryStartPoint > nCoastPoint)
   {
      // The cape point is down-coast from the coast point
      if (nCoastSeaHand == LEFT_HANDED)
         nShadowZoneCoastToCapeSeaHand = LEFT_HANDED;
      else
         nShadowZoneCoastToCapeSeaHand = RIGHT_HANDED;
   }
   else
   {
      // The cape point is up-coast from the coast point
      if (nCoastSeaHand == LEFT_HANDED)
         nShadowZoneCoastToCapeSeaHand = RIGHT_HANDED;
      else
         nShadowZoneCoastToCapeSeaHand = LEFT_HANDED;
   }

   bool bStartPointOK = false;
   CGeom2DIPoint PtiFloodFillStart;
   for (int nOffset = FLOOD_FILL_START_OFFSET; nOffset > 0; nOffset--)
   {
      if (bStartPointOK)
         break;

      double dWeight = 0.05;
      while ((! bStartPointOK) && (dWeight < 1))
      {
         // Find a start point for the flood fill. Because shadow zones are generally triangular, start by choosing a low weighting so that the start point is close to the coast (i.e. not in the narrow apex of the triangle near the cape point). Then go nOffset cells inward, toward the centre of the shadow zone triangle, seems to be about right
         CGeom2DIPoint PtiWeightAvg = PtiWeightedAverage(pPtiEnd, pPtiStart, dWeight);
         
         // Safety check
         if (PtiWeightAvg == *pPtiStart)
         {
            dWeight += 0.05;
            continue;
         }
         
         PtiFloodFillStart = PtiGetPerpendicular(&PtiWeightAvg, pPtiStart, nOffset, nShadowZoneCoastToCapeSeaHand);

         // Safety check
         if (! bIsWithinValidGrid(&PtiFloodFillStart))
         {
            LogStream << m_ulIteration << ": " << ERR << "start point [" << PtiFloodFillStart.nGetX() << "][" << PtiFloodFillStart.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiFloodFillStart.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiFloodFillStart.nGetY()) << "} for flood fill of shadow zone is outside grid" << endl;

            return RTN_ERR_SHADOW_ZONE_FLOOD_FILL_NOGRID;
         }

         if (m_pRasterGrid->m_Cell[PtiFloodFillStart.nGetX()][PtiFloodFillStart.nGetY()].bIsInContiguousSea())
         {
            // Start point is a sea cell, all OK
            LogStream << m_ulIteration << ": shadow zone flood fill start point [" << PtiFloodFillStart.nGetX() << "][" << PtiFloodFillStart.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiFloodFillStart.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiFloodFillStart.nGetY()) << "} is OK for shadow boundary from [" << pPtiStart->nGetX() << "][" << pPtiStart->nGetY() << "] = {" << dGridCentroidXToExtCRSX(pPtiStart->nGetX()) << ", " << dGridCentroidYToExtCRSY(pPtiStart->nGetY()) << "} to [" << pPtiEnd->nGetX() << "][" << pPtiEnd->nGetY() << "] = {" << dGridCentroidXToExtCRSX(pPtiEnd->nGetX()) << ", " << dGridCentroidYToExtCRSY(pPtiEnd->nGetY()) << "}, dWeight = " << dWeight << ", and nOffset = " << nOffset << endl;

            bStartPointOK = true;
         }
         else
         {
            // Start point is not a sea cell
            LogStream << m_ulIteration << ": shadow zone flood fill start point [" << PtiFloodFillStart.nGetX() << "][" << PtiFloodFillStart.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiFloodFillStart.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiFloodFillStart.nGetY()) << "} is NOT a sea cell for shadow boundary from cape point [" << pPtiStart->nGetX() << "][" << pPtiStart->nGetY() << "] = {" << dGridCentroidXToExtCRSX(pPtiStart->nGetX()) << ", " << dGridCentroidYToExtCRSY(pPtiStart->nGetY()) << "} to [" << pPtiEnd->nGetX() << "][" << pPtiEnd->nGetY() << "] = {" << dGridCentroidXToExtCRSX(pPtiEnd->nGetX()) << ", " << dGridCentroidYToExtCRSY(pPtiEnd->nGetY()) << "}, dWeight = " << dWeight << endl;

            dWeight += 0.05;
         }
      }
   }

   if (! bStartPointOK)
   {
      LogStream << m_ulIteration << ": " << ERR << "could not find shadow zone flood fill start point" << endl;

      return RTN_ERR_SHADOW_ZONE_FLOOD_START_POINT;
   }

   // All OK, so create an empty stack
   stack<CGeom2DIPoint> PtiStack;

   // We have a flood fill start point so push this point onto the stack
   PtiStack.push(PtiFloodFillStart);

   // Then do the flood fill: loop until there are no more cell co-ords on the stack
   while (! PtiStack.empty())
   {
      CGeom2DIPoint Pti = PtiStack.top();
      PtiStack.pop();

      int
         nX = Pti.nGetX(),
         nY = Pti.nGetY();

      while ((nX >= 0) && m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea() && (! m_pRasterGrid->m_Cell[nX][nY].bIsinThisShadowZone(-nZone-1)) && (! m_pRasterGrid->m_Cell[nX][nY].bIsShadowZoneBoundary()) && (! m_pRasterGrid->m_Cell[nX][nY].bIsCoastline()))
         nX--;

      nX++;

      bool
         bSpanAbove = false,
         bSpanBelow = false;

      while ((nX < m_nXGridMax) && m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea() && (! m_pRasterGrid->m_Cell[nX][nY].bIsinThisShadowZone(-nZone-1)) && (! m_pRasterGrid->m_Cell[nX][nY].bIsShadowZoneBoundary()) && (! m_pRasterGrid->m_Cell[nX][nY].bIsCoastline()))
      {
         // Mark the cell as being in the shadow zone
         m_pRasterGrid->m_Cell[nX][nY].SetShadowZoneNumber(-nZone-1);

//          LogStream << m_ulIteration << ": [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} marked as shadow zone" << endl;

         if ((! bSpanAbove) && (nY > 0) && m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea() && (! m_pRasterGrid->m_Cell[nX][nY-1].bIsinThisShadowZone(-nZone-1)) && (! m_pRasterGrid->m_Cell[nX][nY-1].bIsShadowZoneBoundary()) && (! m_pRasterGrid->m_Cell[nX][nY-1].bIsCoastline()))
         {
            PtiStack.push(CGeom2DIPoint(nX, nY-1));
            bSpanAbove = true;
         }
         else if (bSpanAbove && (nY > 0) && ((! m_pRasterGrid->m_Cell[nX][nY-1].bIsInContiguousSea()) || m_pRasterGrid->m_Cell[nX][nY-1].bIsinThisShadowZone(-nZone-1) || m_pRasterGrid->m_Cell[nX][nY-1].bIsShadowZoneBoundary() || m_pRasterGrid->m_Cell[nX][nY-1].bIsCoastline()))
         {
            bSpanAbove = false;
         }

         if ((! bSpanBelow) && m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea() && (nY < m_nYGridMax+1) && (! m_pRasterGrid->m_Cell[nX][nY+1].bIsinThisShadowZone(-nZone-1)) && (! m_pRasterGrid->m_Cell[nX][nY+1].bIsShadowZoneBoundary()) && (! m_pRasterGrid->m_Cell[nX][nY+1].bIsCoastline()))
         {
            PtiStack.push(CGeom2DIPoint(nX, nY+1));
            bSpanBelow = true;
         }
         else if (bSpanBelow && (nY < m_nYGridMax-1) && ((! m_pRasterGrid->m_Cell[nX][nY+1].bIsInContiguousSea()) || m_pRasterGrid->m_Cell[nX][nY+1].bIsinThisShadowZone(-nZone-1) || m_pRasterGrid->m_Cell[nX][nY+1].bIsShadowZoneBoundary() || m_pRasterGrid->m_Cell[nX][nY+1].bIsCoastline()))
         {
            bSpanBelow = false;
         }

         nX++;
      }
   }

   return RTN_OK;
}


/*===============================================================================================================================

 Sweep the shadow zone, changing wave orientation and height

===============================================================================================================================*/
int CSimulation::nSweepShadowZone(int const nCoast, int const nZone, int const nShadowBoundaryStartPoint, CGeom2DIPoint const* pPtiStart, int const nCoastPoint, CGeom2DIPoint const* pPtiEnd, int& nCoastSwept)
{
   // TODO merge this and the down drift sweep. We go out (one call at a time) radially from the shadow boundary start point. For each radius, move along the radius and change both shadow zone cells and down drift cells which lie on that radius
   int
      nCoastSize = m_VCoast[nCoast].nGetCoastlineSize(),
      nCoastSeaHand = m_VCoast[nCoast].nGetSeaHandedness(),
      nShadowZoneCoastToCapeSeaHand;

   if (nCoastSeaHand == LEFT_HANDED)
      nShadowZoneCoastToCapeSeaHand = RIGHT_HANDED;
   else
      nShadowZoneCoastToCapeSeaHand = LEFT_HANDED;

   // We are going to sweep the coastline starting from the end point of the shadow zone line, going toward the start point. Which direction is this?
   bool bSweepDownCoast = true;
   if (nCoastPoint > nShadowBoundaryStartPoint)
      bSweepDownCoast = false;

   int nThisEndPoint = nCoastPoint;

   CGeom2DIPoint PtCoastStart = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(0);
   CGeom2DIPoint PtCoastEnd = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastSize-1);

   bool bSomeChanged = false;
   nCoastSwept = 0;
   int nLengthSwept = 0;
   while (true)
   {
      nLengthSwept++;

      if (bSweepDownCoast)
      {
         nThisEndPoint++;
         if (nThisEndPoint > nCoastSize-1)
            break;
      }
      else
      {
         nThisEndPoint--;
         if (nThisEndPoint < 0)
            break;
      }

      if (nThisEndPoint == nShadowBoundaryStartPoint)
         break;

      // Construct a line from the cape to the new end point
      CGeom2DIPoint PtiThisEnd;
      if (nThisEndPoint < 0)
      {
         // The shadow line hit the grid edge from which the coastline begins. We need to move along this edge till we get to the start of the coastline
         if (pPtiEnd->nGetX() == 0)
         {
            // Shadow line hit the W edge of the grid
            PtiThisEnd.SetX(0);

            if (pPtiEnd->nGetY() > PtCoastStart.nGetY())
               PtiThisEnd.SetY(pPtiEnd->nGetY() - nLengthSwept);
            else
               PtiThisEnd.SetY(pPtiEnd->nGetY() + nLengthSwept);
         }
         else if (pPtiEnd->nGetX() == m_nXGridMax-1)
         {
            // Shadow line hit the E edge of the grid
            PtiThisEnd.SetX(m_nXGridMax-1);

            if (pPtiEnd->nGetY() > PtCoastStart.nGetY())
               PtiThisEnd.SetY(pPtiEnd->nGetY() - nLengthSwept);
            else
               PtiThisEnd.SetY(pPtiEnd->nGetY() + nLengthSwept);
         }
         else if (pPtiEnd->nGetY() == 0)
         {
            // Shadow line hit the N edge of the grid
            PtiThisEnd.SetY(0);

            if (pPtiEnd->nGetX() > PtCoastStart.nGetX())
               PtiThisEnd.SetX(pPtiEnd->nGetX() - nLengthSwept);
            else
               PtiThisEnd.SetX(pPtiEnd->nGetX() + nLengthSwept);
         }
         else if (pPtiEnd->nGetY() == m_nYGridMax-1)
         {
            // Shadow line hit the S edge of the grid
            PtiThisEnd.SetY(m_nYGridMax-1);

            if (pPtiEnd->nGetX() > PtCoastStart.nGetX())
               PtiThisEnd.SetX(pPtiEnd->nGetX() - nLengthSwept);
            else
               PtiThisEnd.SetX(pPtiEnd->nGetX() + nLengthSwept);
         }
      }
      else if (nThisEndPoint >= nCoastSize)
      {
         // The shadow line hit the grid edge at which the coastline ends. We need to move along this edge till we get to the end of the coastline
         if (pPtiEnd->nGetX() == 0)
         {
            // Shadow line hit the W edge of the grid
            PtiThisEnd.SetX(0);

            if (pPtiEnd->nGetY() > PtCoastEnd.nGetY())
               PtiThisEnd.SetY(pPtiEnd->nGetY() - nLengthSwept);
            else
               PtiThisEnd.SetY(pPtiEnd->nGetY() + nLengthSwept);
         }
         else if (pPtiEnd->nGetX() == m_nXGridMax-1)
         {
            // Shadow line hit the E edge of the grid
            PtiThisEnd.SetX(m_nXGridMax-1);

            if (pPtiEnd->nGetY() > PtCoastEnd.nGetY())
               PtiThisEnd.SetY(pPtiEnd->nGetY() - nLengthSwept);
            else
               PtiThisEnd.SetY(pPtiEnd->nGetY() + nLengthSwept);
         }
         else if (pPtiEnd->nGetY() == 0)
         {
            // Shadow line hit the N edge of the grid
            PtiThisEnd.SetY(0);

            if (pPtiEnd->nGetX() > PtCoastEnd.nGetX())
               PtiThisEnd.SetX(pPtiEnd->nGetX() - nLengthSwept);
            else
               PtiThisEnd.SetX(pPtiEnd->nGetX() + nLengthSwept);
         }
         else if (pPtiEnd->nGetY() == m_nYGridMax-1)
         {
            // Shadow line hit the S edge of the grid
            PtiThisEnd.SetY(m_nYGridMax-1);

            if (pPtiEnd->nGetX() > PtCoastEnd.nGetX())
               PtiThisEnd.SetX(pPtiEnd->nGetX() - nLengthSwept);
            else
               PtiThisEnd.SetX(pPtiEnd->nGetX() + nLengthSwept);
         }
      }
      else
      {
         // The shadow line hit the coastline, not a grid edge
         PtiThisEnd = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nThisEndPoint);    // grid CRS

         if (bSomeChanged)
            nCoastSwept++;

         bSomeChanged = false;
      }

      int
         nSweepEndX = PtiThisEnd.nGetX(),
         nSweepEndY = PtiThisEnd.nGetY(),
         nXDist = nSweepEndX - pPtiStart->nGetX(),
         nYDist = nSweepEndY - pPtiStart->nGetY();
      double
         dX = pPtiStart->nGetX(),
         dY = pPtiStart->nGetY(),
         dDist = dGetDistanceBetween(pPtiStart, &PtiThisEnd),
         dXIncr = nXDist / dDist,
         dYIncr = nYDist / dDist;

      for (int n = 0; n < static_cast<int>(dRound(dDist)); n++)
      {
         dX += dXIncr;
         dY += dYIncr;

         int
            nX = static_cast<int>(dRound(dX)),
            nY = static_cast<int>(dRound(dY));

         // Safety check
         if (! bIsWithinValidGrid(nX, nY))
            continue;

         int nZoneCode = m_pRasterGrid->m_Cell[nX][nY].nGetShadowZoneNumber();
         if (nZoneCode == (-nZone-1))
         {
            // We are in the shadow zone and have not already processed this cell, so mark it (a +ve number)
            m_pRasterGrid->m_Cell[nX][nY].SetShadowZoneNumber(nZone+1);
            bSomeChanged = true;

            // Next calculate wave angle here: first calculate dOmega, the signed angle subtended between this end point and the start point, and this end point and the end of the shadow boundary
            double dOmega = 180 * dAngleSubtended(pPtiStart, &PtiThisEnd, pPtiEnd) / PI;

            // If dOmega is 90 degrees or more in either direction, set both wave angle and wave height to zero
            if (tAbs(dOmega) >= 90)
            {
               m_pRasterGrid->m_Cell[nX][nY].SetWaveOrientation(0);
               m_pRasterGrid->m_Cell[nX][nY].SetWaveHeight(0);
            }
            else
            {
               // Adapted from equation 12 in Hurst et al.
               double dDeltaShadowWaveAngle = 1.5 * dOmega;

               // Get the pre-existing (i.e. shore-parallel) wave orientation
               double dWaveOrientation = m_pRasterGrid->m_Cell[nX][nY].dGetWaveOrientation();

               double dShadowWaveOrientation;
               if (nShadowZoneCoastToCapeSeaHand == LEFT_HANDED)
                  dShadowWaveOrientation = dWaveOrientation + dDeltaShadowWaveAngle;
               else
                  dShadowWaveOrientation = dWaveOrientation - dDeltaShadowWaveAngle;

               // Set the shadow zone wave orientation
               m_pRasterGrid->m_Cell[nX][nY].SetWaveOrientation(dKeepWithin360(dShadowWaveOrientation));

               // Now calculate wave height within the shadow zone, use equation 13 from Hurst et al.
               double dKp = 0.5 * cos(dOmega * PI / 180);

               // Get the pre-existing (i.e. shore-parallel) wave height
               double dWaveHeight = m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight();

               // Set the shadow zone wave height
               m_pRasterGrid->m_Cell[nX][nY].SetWaveHeight(dKp * dWaveHeight);

               //             LogStream << m_ulIteration << ": nThisEndPoint = " << nThisEndPoint << ", n = " << n << ", nLengthSwept = " << nLengthSwept << ", cape point {" << dGridCentroidXToExtCRSX(pPtiStart->nGetX()) << ", " << dGridCentroidYToExtCRSY(pPtiStart->nGetY()) << "}, end point {" << dGridCentroidXToExtCRSX(pPtiEnd->nGetX()) << ", " << dGridCentroidYToExtCRSY(pPtiEnd->nGetY()) << "}, this point [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}, angle subtended = " << dOmega << " degrees, m_pRasterGrid->m_Cell[" << nX < "][" << nY << "].dGetDeepWaterWaveHeight() = " << m_pRasterGrid->m_Cell[nX][nY].dGetDeepWaterWaveHeight() << " degrees, dDeltaShadowWaveAngle = " << dDeltaShadowWaveAngle << " degrees, dWaveOrientation = " << dWaveOrientation << " degrees, dShadowWaveOrientation = " << dShadowWaveOrientation << " degrees, dWaveHeight = " << dWaveHeight << " m, dKp = " << dKp << ", shadow zone wave height = " << dKp * dWaveHeight << " m" << endl;
            }
         }
      }
   }

   return RTN_OK;
}


/*===============================================================================================================================

 Sweeps the coast downdrift from the shadow zone, changing wave height in order to conserve energy

===============================================================================================================================*/
int CSimulation::nSweepDownDriftFromShadowZone(int const nCoast, int const nZone, int const nShadowBoundaryStartPoint, CGeom2DIPoint const* pPtiStart, int const nShadowBoundaryEndPoint, int const nSweepLength)
{
   int nCoastSize = m_VCoast[nCoast].nGetCoastlineSize();

   // Which direction do we need to sweep?
   bool bSweepDownCoast = false;
   if (nShadowBoundaryStartPoint < nShadowBoundaryEndPoint)
      bSweepDownCoast = true;

   // Start sweeping from the shadow boundary end point
   int 
      nSweepCoastPoint = nShadowBoundaryEndPoint,
      nEdge = NO_DIRECTION;
      
   for (int nSweep = 0; nSweep < nSweepLength; nSweep++)
   {
      int
         nSweepEndX,
         nSweepEndY;
         
      if (nEdge == NO_DIRECTION)
      {
         // Find the cell which is the end point of the sweep line
         CGeom2DIPoint PtiSweepPoint = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nSweepCoastPoint);
         nSweepEndX = PtiSweepPoint.nGetX(),
         nSweepEndY = PtiSweepPoint.nGetY();
            
         if (m_pRasterGrid->m_Cell[nSweepEndX][nSweepEndY].bIsEdgeCell())
         {
            // We have hit a grid-edge cell
            nEdge = m_pRasterGrid->m_Cell[nSweepEndX][nSweepEndY].nGetEdgeCell();
         }
      }
            
      if ((nSweepCoastPoint < 0) || (nSweepCoastPoint >= nCoastSize))
      {
         // The end point of the sweep line is beyond the start or end of the coast
         if (nEdge == NORTH)
            nSweepEndY--;
         else if (nEdge == EAST)
            nSweepEndX++;
         else if (nEdge == SOUTH)
            nSweepEndY++;
         else if (nEdge == EAST)
            nSweepEndX--;
      }
      
      // End of the sweep line, in grid CRS (note that this may be outside the grid)
      CGeom2DIPoint PtiSweepEnd(nSweepEndX, nSweepEndY);

      // End of the sweep line, in external CRS
      CGeom2DPoint PtSweepEnd(dGridCentroidXToExtCRSX(nSweepEndX), dGridCentroidYToExtCRSY(nSweepEndY));

      int
         nSweepStartX = pPtiStart->nGetX(),
         nSweepStartY = pPtiStart->nGetY(),
         nXDist = nSweepEndX - nSweepStartX,
         nYDist = nSweepEndY - nSweepStartY,
         nXLast = -1,
         nYLast = -1;
      double
         dX = nSweepStartX,
         dY = nSweepStartY,
         dDist = dGetDistanceBetween(pPtiStart, &PtiSweepEnd),
         dXIncr = nXDist / dDist,
         dYIncr = nYDist / dDist;

//       LogStream << m_ulIteration << ": sweeping downdrift, nSweepCoastPoint = " << nSweepCoastPoint << ", nSweep = " << nSweep << ", nSweepLength = " << nSweepLength << ", shadow boundary start point {" << dGridCentroidXToExtCRSX(nSweepStartX) << ", " << dGridCentroidYToExtCRSY(nSweepStartY) << "}, sweep line end point {" << dGridCentroidXToExtCRSX(nSweepEndX) << ", " << dGridCentroidYToExtCRSY(nSweepEndY) << "}" << endl;

      for (int n = 0; n < static_cast<int>(dRound(dDist)); n++)
      {
         dX += dXIncr;
         dY += dYIncr;

         int
            nX = static_cast<int>(dRound(dX)),
            nY = static_cast<int>(dRound(dY));

         // Safety check
         if (! bIsWithinValidGrid(nX, nY))
            continue;
         
         if ((nX == nXLast) && (nY == nYLast))
            continue;

         // Get the pre-existing (i.e. shore-parallel) wave height
         double dWaveHeight = m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight();
         if (dWaveHeight == DBL_NODATA)
         {
            // Is not a sea cell            
            LogStream << m_ulIteration << ": ignored, not a sea cell [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;
            
            nXLast = nX;
            nYLast = nY;

            continue;
         }

         int nZoneCode = m_pRasterGrid->m_Cell[nX][nY].nGetShadowZoneNumber();
         if (nZoneCode == (nZone+1))
         {
            // This cell is in 'our' shadow zone, so don't change it
            LogStream << m_ulIteration << ": ignored, in shadow zone " << nZone+1 << " [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;
            
            nXLast = nX;
            nYLast = nY;

            continue;
         }
         
         if (nZoneCode < 0)
         {
            // This cell is in a shadow zone but is not yet processed, so don't change it
            LogStream << m_ulIteration << ": ignored, unprocessed cell in shadow zone " << nZoneCode << " [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;
            
            nXLast = nX;
            nYLast = nY;
            
            continue;            
         }
         
         nZoneCode = m_pRasterGrid->m_Cell[nX][nY].nGetDownDriftZoneNumber();
         if (nZoneCode == (nZone+1))
         {
            // We have already processed this cell for this downdrift zone
            LogStream << m_ulIteration << ": ignored, already done for down drift zone " << nZone+1 << " [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;
            
            nXLast = nX;
            nYLast = nY;

            continue;
         }

         // OK, we are downdrift of the shadow zone area and have not yet processed this cell for this zone, so mark it
         m_pRasterGrid->m_Cell[nX][nY].SetDownDriftZoneNumber(nZone+1);

         // Equation 14 from Hurst et al. NOTE could not get this to work, so used the equation below instead TODO check this with Andres
//         double dKp = 0.5 * (1.0 - sin((PI * 90.0 * nSweep) / (180.0 * nSweepLength)));
//         double dKp = 0.5 + (0.5 * sin((PI * nSweep) / (2.0 * nSweepLength)));
         double dKp = 0.5 + (1.0 * sin((PI * nSweep) / (2.0 * nSweepLength)));
         
         // Set the modified wave height
         m_pRasterGrid->m_Cell[nX][nY].SetWaveHeight(dKp * dWaveHeight);

         LogStream << m_ulIteration << ": nSweepCoastPoint = " << nSweepCoastPoint << ", shadow boundary start point {" << dGridCentroidXToExtCRSX(nSweepStartX) << ", " << dGridCentroidYToExtCRSY(nSweepStartY) << "}, nSweep = " << nSweep << ", sweep line end point {" << dGridCentroidXToExtCRSX(nSweepEndX) << ", " << dGridCentroidYToExtCRSY(nSweepEndY) << "}, n = " << n << ", this point [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}, DOWNDRIFT SWEEP original dWaveHeight = " << dWaveHeight << " m, dKp = " << dKp << ", modified wave height = " << dKp * dWaveHeight << " m" << endl;
         
         nXLast = nX;
         nYLast = nY;
      }
      
      // Move to next sweep end point
      if (bSweepDownCoast)
         nSweepCoastPoint++;
      else
         nSweepCoastPoint--;
   }

   return RTN_OK;
}


