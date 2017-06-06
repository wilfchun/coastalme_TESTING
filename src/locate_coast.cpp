/*!
 *
 * \file locate_coast.cpp
 * \brief Finds the coastline on the raster grid
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
#include <cfloat>

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::ios;

#include <iomanip>
using std::setiosflags;
using std::resetiosflags;
using std::setprecision;
using std::setw;

#include <stack>
using std::stack;

#include "cme.h"
#include "i_line.h"
#include "line.h"
#include "simulation.h"
#include "raster_grid.h"
#include "coast.h"


/*===============================================================================================================================

 First find all connected sea areas, then locate the vector coastline(s), then put these onto the raster grid

===============================================================================================================================*/
int CSimulation::nLocateSeaAndCoasts(void)
{
   // Find all connected sea cells
   FindAllSeaCells();

   // Find every coastline on the raster grid, mark raster cells, then create the vector coastline
   int nRet = nTraceAllCoasts();
   if (nRet != RTN_OK)
      return nRet;
   
   // Have we created any coasts?
   if (m_VCoast.empty())
   {
      cerr << m_ulTimestep << ": " << ERR << "no coastline located" << endl;
      return RTN_ERR_NOCOAST;
   }

   return RTN_OK;
}


/*===============================================================================================================================

 Finds and flags all sea areas which have at least one cell at a grid edge (i.e. does not flag 'inland' seas)

===============================================================================================================================*/
void CSimulation::FindAllSeaCells(void)
{
   // Go along the list of edge cells
   for (unsigned int n = 0; n < m_VEdgeCell.size(); n++)
   {
      if (m_bOmitSearchNorthEdge && m_VEdgeCellEdge[n] == NORTH)
         continue;
      
      if (m_bOmitSearchSouthEdge && m_VEdgeCellEdge[n] == SOUTH)
         continue;
      
      if (m_bOmitSearchWestEdge && m_VEdgeCellEdge[n] == WEST)
         continue;
      
      if (m_bOmitSearchEastEdge && m_VEdgeCellEdge[n] == EAST)
         continue;
      
      int
         nX = m_VEdgeCell[n].nGetX(),
         nY = m_VEdgeCell[n].nGetY();
         
      if ((m_pRasterGrid->m_Cell[nX][nY].bIsInundated()) && (m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth() == 0))
         // This edge cell is below SWL and sea depth remains set to zero
         FloodFillSea(nX, nY);         
   }
}


/*===============================================================================================================================

 Flood-fills all sea cells starting from a given cell. The flood fill code used here is adapted from an example by Lode Vandevenne (http://lodev.org/cgtutor/floodfill.html#Scanline_Floodfill_Algorithm_With_Stack)

===============================================================================================================================*/
void CSimulation::FloodFillSea(int const nXStart, int const nYStart)
{
   // Create an empty stack
   stack<CGeom2DIPoint> PtiStack;

   // Start at the given edge cell, push this onto the stack
   PtiStack.push(CGeom2DIPoint(nXStart, nYStart));

   // Then do the flood fill: loop until there are no more cell co-ords on the stack
   while (! PtiStack.empty())
   {
      CGeom2DIPoint Pti = PtiStack.top();
      PtiStack.pop();

      int
         nX = Pti.nGetX(),
         nY = Pti.nGetY();

      while ((nX >= 0) && (! m_pRasterGrid->m_Cell[nX][nY].bBasementElevIsMissingValue()) && (m_pRasterGrid->m_Cell[nX][nY].bIsInundated()) && (m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth() == 0))
         nX--;

      nX++;
      bool
         bSpanAbove = false,
         bSpanBelow = false;
         
      while ((nX < m_nXGridMax) && (! m_pRasterGrid->m_Cell[nX][nY].bBasementElevIsMissingValue()) && (m_pRasterGrid->m_Cell[nX][nY].bIsInundated()) && (m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth() == 0))
      {
         // Set the sea depth for this cell
         m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();

         // Mark as sea
         m_pRasterGrid->m_Cell[nX][nY].SetInContiguousSea();         
         m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetLFCategory(LF_CAT_SEA);
         
         // Set all sea cells to have deep water (off-shore) wave orientation and height, will change this later for cells closer to the shoreline if we have on-shore waves
         m_pRasterGrid->m_Cell[nX][nY].SetWaveOrientation(m_dDeepWaterWaveOrientation);
         m_pRasterGrid->m_Cell[nX][nY].SetWaveHeight(m_dDeepWaterWaveHeight);
         
         // Now sort out the x-y extremities of the contiguous sea for the bounding box (used later in wave propagation)
         if (nX < m_nXMinBoundingBox)
            m_nXMinBoundingBox = nX;
                                    
         if (nX > m_nXMaxBoundingBox)
            m_nXMaxBoundingBox = nX;
         
         if (nY < m_nYMinBoundingBox)
            m_nYMinBoundingBox = nY;
                                    
         if (nY > m_nYMaxBoundingBox)
            m_nYMaxBoundingBox = nY;
         
         // Update count
         m_ulThisTimestepNumSeaCells++;

         if ((! bSpanAbove) && (nY > 0) && (! m_pRasterGrid->m_Cell[nX][nY-1].bBasementElevIsMissingValue()) && (m_pRasterGrid->m_Cell[nX][nY-1].bIsInundated()) && (m_pRasterGrid->m_Cell[nX][nY-1].dGetSeaDepth() == 0))
         {
            PtiStack.push(CGeom2DIPoint(nX, nY-1));
            bSpanAbove = true;
         }
         else if (bSpanAbove && (nY > 0) && (m_pRasterGrid->m_Cell[nX][nY-1].bBasementElevIsMissingValue()) && (! m_pRasterGrid->m_Cell[nX][nY-1].bIsInundated() || (m_pRasterGrid->m_Cell[nX][nY-1].dGetSeaDepth() != 0)))
         {
            bSpanAbove = false;
         }

         if ((! bSpanBelow) && (nY < m_nYGridMax-1) && (! m_pRasterGrid->m_Cell[nX][nY+1].bBasementElevIsMissingValue()) && (m_pRasterGrid->m_Cell[nX][nY+1].bIsInundated()) && (m_pRasterGrid->m_Cell[nX][nY+1].dGetSeaDepth() == 0))
         {
            PtiStack.push(CGeom2DIPoint(nX, nY+1));
            bSpanBelow = true;
         }
         else if (bSpanBelow && (nY < m_nYGridMax-1) && (m_pRasterGrid->m_Cell[nX][nY+1].bBasementElevIsMissingValue()) && (! m_pRasterGrid->m_Cell[nX][nY+1].bIsInundated() || (m_pRasterGrid->m_Cell[nX][nY+1].dGetSeaDepth() != 0)))
         {
            bSpanBelow = false;
         }

         nX++;
      }
   }
   
   LogStream << m_ulTimestep << ": doing flood fill of sea from [" << nXStart << "][" << nYStart << "] {" << dGridCentroidXToExtCRSX(nXStart) << ", " << dGridCentroidYToExtCRSY(nYStart) << "} with SWL = " << m_dThisTimestepSWL << ", " << m_ulThisTimestepNumSeaCells << " of " << m_ulNumCells << " cells now marked as sea (" <<  setiosflags(ios::fixed) << setprecision(2) << 100.0 * m_ulThisTimestepNumSeaCells / m_ulNumCells << " %)" << endl;
   
//    LogStream << " m_nXMinBoundingBox = " << m_nXMinBoundingBox << " m_nXMaxBoundingBox = " << m_nXMaxBoundingBox << " m_nYMinBoundingBox = " << m_nYMinBoundingBox << " m_nYMaxBoundingBox = " << m_nYMaxBoundingBox << endl;
}


/*===============================================================================================================================

 Locates coastline start points on the edges of the raster grid, and trace the vector coastline(s) from these start points. The vector coastlines are then smoothed

===============================================================================================================================*/
int CSimulation::nTraceAllCoasts(void)
{
   // Go along the list of edge cells
   for (unsigned int n = 0; n < m_VEdgeCell.size()-1; n++)
   {
      if (m_bOmitSearchNorthEdge && (m_VEdgeCellEdge[n] == NORTH || m_VEdgeCellEdge[n+1] == NORTH))
         continue;
      
      if (m_bOmitSearchSouthEdge && (m_VEdgeCellEdge[n] == SOUTH || m_VEdgeCellEdge[n+1] == SOUTH))
         continue;
      
      if (m_bOmitSearchWestEdge && (m_VEdgeCellEdge[n] == WEST || m_VEdgeCellEdge[n+1] == WEST))
         continue;
      
      if (m_bOmitSearchEastEdge && (m_VEdgeCellEdge[n] == EAST || m_VEdgeCellEdge[n+1] == EAST))
         continue;
      
      int
         nXThis = m_VEdgeCell[n].nGetX(),
         nYThis = m_VEdgeCell[n].nGetY(),
         nXNext = m_VEdgeCell[n+1].nGetX(),
         nYNext = m_VEdgeCell[n+1].nGetY();
         
      // Get "Is it sea?" information for 'this' and 'next' cells
      bool
         bThisCellIsSea = m_pRasterGrid->m_Cell[nXThis][nYThis].bIsInContiguousSea(),
         bNextCellIsSea = m_pRasterGrid->m_Cell[nXNext][nYNext].bIsInContiguousSea();

      // Are we at a coast?
      if ((! bThisCellIsSea) && bNextCellIsSea)
      {
         // 'This' cell is just inland, has it already been flagged as a coast cell?
         if (! m_pRasterGrid->m_Cell[nXThis][nYThis].bIsCoastline())
         {
            // It has not, so trace a coastline from 'this' cell
            int nSearchDirection = nGetOppositeDirection(m_VEdgeCellEdge[n]);
            int nRet = nTraceCoastLine(nSearchDirection, LEFT_HANDED, nXThis, nYThis);
            if (nRet != RTN_OK)
               return nRet;
         }
      }
      else if (bThisCellIsSea && (! bNextCellIsSea))
      {
         // The 'next' cell is just inland, has it already been flagged as a coast cell?
         if (! m_pRasterGrid->m_Cell[nXNext][nYNext].bIsCoastline())
         {
            // It has not, so trace a coastline from the 'next' cell
            int nSearchDirection = nGetOppositeDirection(m_VEdgeCellEdge[n+1]);
            int nRet = nTraceCoastLine(nSearchDirection, RIGHT_HANDED, nXNext, nYNext);
            if (nRet != RTN_OK)
               return nRet;
         }
      }
   }
   
   return RTN_OK;
}


/*==============================================================================================================================

 Traces the coastline (which is defined to be just above still water level) on the grid using the 'wall follower' rule for maze traversal (http://en.wikipedia.org/wiki/Maze_solving_algorithm#Wall_follower)

===============================================================================================================================*/
int CSimulation::nTraceCoastLine(int const nStartSearchDirection, int const nHandedness, int const nStartX, int const nStartY)
{
   bool
      bAtCoast = false,
      bHasLeftStartEdge = false;

   int
      nX = nStartX,
      nY = nStartY,
      nSearchDirection = nStartSearchDirection,
      nLastX = INT_NODATA,
      nLastY = INT_NODATA,
      nBeforeLastX = INT_NODATA,
      nBeforeLastY = INT_NODATA,
      nRoundTheLoop = 0,
      nThisLen = 0;
//       nLastLen = 0,
//       nPreLastLen = 0;

   // Temporary coastline as a line of integer points (in grid coordinates)
   CGeomILine LTempGridCRS;

   // Start at this grid-edge point and trace the rest of the coastline using the 'wall follower' rule for maze traversal, keeping next to cells flagged as sea
   do
   {
//       LogStream << "Now at [" << nX << "][" << nY << "] {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;    
//       LogStream << "LTempGridCRS is now:" << endl;
//       for (int n = 0; n < LTempGridCRS.nGetSize(); n++)
//          LogStream << "[" << LTempGridCRS[n].nGetX() << "][" << LTempGridCRS[n].nGetY() << "] {" << dGridCentroidXToExtCRSX(LTempGridCRS[n].nGetX()) << ", " << dGridCentroidYToExtCRSY(LTempGridCRS[n].nGetY()) << "}" << endl;
//       LogStream <<  "=================" << endl;
            
//       // First safety check: oscillation between two points
//       if ((nX == nBeforeLastX) && (nY == nBeforeLastY))
//       {
//          LogStream << m_ulTimestep << ": abandoning tracing coast from start point [" << nStartX << "][" << nStartY << "] {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} since stuck in two-step loop between [" << nX << "][" << nY << "] {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} and [" << nLastX << "][" << nLastY << "] {" << dGridCentroidXToExtCRSX(nLastX) << ", " << dGridCentroidYToExtCRSY(nLastY) << "}" << endl;
//          
//          return RTN_OK;
//       }
      
      nBeforeLastX = nLastX;
      nBeforeLastY = nLastY;
      nLastX = nX;
      nLastY = nY;
      
      // Second safety check
      if (++nRoundTheLoop > m_nCoastMax)
      {
         LogStream << m_ulTimestep << ": abandoning tracing coastline from [" << nStartX << "][" << nStartY << "] {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "}, exceeded maximum length (" << m_nCoastMax << ")" << endl;
         
         for (int n = 0; n < LTempGridCRS.nGetSize(); n++)
            LogStream << "[" << LTempGridCRS[n].nGetX() << "][" << LTempGridCRS[n].nGetY() << "] {" << dGridCentroidXToExtCRSX(LTempGridCRS[n].nGetX()) << ", " << dGridCentroidYToExtCRSY(LTempGridCRS[n].nGetY()) << "}" << endl;
         
         LogStream << endl;
         
         return RTN_OK;        
      }

//       // Third safety check      
//       if (LTempGridCRS.bIsPresent(nX, nY))
//       {
//          // Uh-oh, we are in a loop
//          LogStream << m_ulTimestep << ": abandoning tracing coastline from [" << nStartX << "][" << nStartY << "] {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "}, about to append [" << nX << "][" << nY << "] which is already present in list (size "<< nThisLen << ")" << endl;
//          if (nThisLen > 0)
//          {
//             for (int n = 0; n < LTempGridCRS.nGetSize(); n++)
//                LogStream << "[" << LTempGridCRS[n].nGetX() << "][" << LTempGridCRS[n].nGetY() << "] {" << dGridCentroidXToExtCRSX(LTempGridCRS[n].nGetX()) << ", " << dGridCentroidYToExtCRSY(LTempGridCRS[n].nGetY()) << "}" << endl;
// 
//             LogStream << endl;
//          }
//          
//          return RTN_OK;
//       }

      //       // First safety device
//       nPreLastLen = nLastLen;
//       nLastLen = nThisLen;
//       
//       nThisLen = LTempGridCRS.nGetSize();
//       if ((nRoundTheLoop > 2) && (nThisLen == nPreLastLen))
//       {
//          // List is not growing
//          bSearchedTooLong = true;
//          
//          LogStream << m_ulTimestep << ": abandoning tracing coastline from [" << nStartX << "][" << nStartY << "] {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} since list (size " << nThisLen << ") is not growing" << endl;
//          if (nThisLen > 0)
//          {
//             for (int n = 0; n < LTempGridCRS.nGetSize(); n++)
//                LogStream << "[" << LTempGridCRS[n].nGetX() << "][" << LTempGridCRS[n].nGetY() << "] {" << dGridCentroidXToExtCRSX(LTempGridCRS[n].nGetX()) << ", " << dGridCentroidYToExtCRSY(LTempGridCRS[n].nGetY()) << "}" << endl;
//          }
//          
//          return RTN_OK;
//       }     
      
      

      // OK so far: so have we left the start edge?
      if (! bHasLeftStartEdge)
      {
         if (((nStartSearchDirection == SOUTH) && (nY > nStartY)) || ((nStartSearchDirection == NORTH) && (nY < nStartY)) ||
             ((nStartSearchDirection == EAST) && (nX > nStartX))  || ((nStartSearchDirection == WEST) && (nX < nStartX)))
            bHasLeftStartEdge = true;
      }

      // Leave the loop if the vector coastline has left the start edge, then we find a coast cell which is at an edge (note that this edge could be the same edge from which this coastline started)
      if (bHasLeftStartEdge && bAtCoast)
      {
         if (m_pRasterGrid->m_Cell[nX][nY].bIsEdgeCell())
            break;
      }

      // OK now sort out the next iteration of the search
      bAtCoast = false;

      int
         nXSeaward = 0,
         nYSeaward = 0,
         nSeawardNewDirection = 0,
         nXStraightOn = 0,
         nYStraightOn = 0,
         nXAntiSeaward = 0,
         nYAntiSeaward = 0,
         nAntiSeawardNewDirection = 0,
         nXGoBack = 0,
         nYGoBack = 0,
         nGoBackNewDirection = 0;

      CGeom2DIPoint Pti(nX, nY);
      
      // Set up the variables
      switch (nHandedness)
      {
         case RIGHT_HANDED:
            // The sea is to the right-hand side of the coast as we traverse it. We are just inland, so we need to keep heading right to find the sea
            switch (nSearchDirection)
            {
               case NORTH:
                  // The sea is towards the RHS (E) of the coast, so first try to go right (to the E)
                  nXSeaward = nX+1;
                  nYSeaward = nY;
                  nSeawardNewDirection = EAST;

                  // If can't do this, try to go straight on (to the N)
                  nXStraightOn = nX;
                  nYStraightOn = nY-1;

                  // If can't do either of these, try to go anti-seaward i.e. towards the LHS (W)
                  nXAntiSeaward = nX-1;
                  nYAntiSeaward = nY;
                  nAntiSeawardNewDirection = WEST;

                  // As a last resort, go back (to the S)
                  nXGoBack = nX;
                  nYGoBack = nY+1;
                  nGoBackNewDirection = SOUTH;

                  break;

               case EAST:
                  // The sea is towards the RHS (S) of the coast, so first try to go right (to the S)
                  nXSeaward = nX;
                  nYSeaward = nY+1;
                  nSeawardNewDirection = SOUTH;

                  // If can't do this, try to go straight on (to the E)
                  nXStraightOn = nX+1;
                  nYStraightOn = nY;

                  // If can't do either of these, try to go anti-seaward i.e. towards the LHS (N)
                  nXAntiSeaward = nX;
                  nYAntiSeaward = nY-1;
                  nAntiSeawardNewDirection = NORTH;

                  // As a last resort, go back (to the W)
                  nXGoBack = nX-1;
                  nYGoBack = nY;
                  nGoBackNewDirection = WEST;

                  break;

               case SOUTH:
                  // The sea is towards the RHS (W) of the coast, so first try to go right (to the W)
                  nXSeaward = nX-1;
                  nYSeaward = nY;
                  nSeawardNewDirection = WEST;

                  // If can't do this, try to go straight on (to the S)
                  nXStraightOn = nX;
                  nYStraightOn = nY+1;

                  // If can't do either of these, try to go anti-seaward i.e. towards the LHS (E)
                  nXAntiSeaward = nX+1;
                  nYAntiSeaward = nY;
                  nAntiSeawardNewDirection = EAST;

                  // As a last resort, go back (to the N)
                  nXGoBack = nX;
                  nYGoBack = nY-1;
                  nGoBackNewDirection = NORTH;

                  break;

               case WEST :
                  // The sea is towards the RHS (N) of the coast, so first try to go right (to the N)
                  nXSeaward = nX;
                  nYSeaward = nY-1;
                  nSeawardNewDirection = NORTH;

                  // If can't do this, try to go straight on (to the W)
                  nXStraightOn = nX-1;
                  nYStraightOn = nY;

                  // If can't do either of these, try to go anti-seaward i.e. towards the LHS (S)
                  nXAntiSeaward = nX;
                  nYAntiSeaward = nY+1;
                  nAntiSeawardNewDirection = SOUTH;

                  // As a last resort, go back (to the E)
                  nXGoBack = nX+1;
                  nYGoBack = nY;
                  nGoBackNewDirection = EAST;

                  break;
            }
            break;

         case LEFT_HANDED:
            // The sea is to the left-hand side of the coast as we traverse it. We are just inland, so we need to keep heading left to find the sea
            switch (nSearchDirection)
            {
               case NORTH:
                  // The sea is towards the LHS (W) of the coast, so first try to go left (to the W)
                  nXSeaward = nX-1;
                  nYSeaward = nY;
                  nSeawardNewDirection = WEST;

                  // If can't do this, try to go straight on (to the N)
                  nXStraightOn = nX;
                  nYStraightOn = nY-1;

                  // If can't do either of these, try to go anti-seaward i.e. towards the RHS (E)
                  nXAntiSeaward = nX+1;
                  nYAntiSeaward = nY;
                  nAntiSeawardNewDirection = EAST;

                  // As a last resort, go back (to the S)
                  nXGoBack = nX;
                  nYGoBack = nY+1;
                  nGoBackNewDirection = SOUTH;

                  break;

               case EAST :
                  // The sea is towards the LHS (N) of the coast, so first try to go left (to the N)
                  nXSeaward = nX;
                  nYSeaward = nY-1;
                  nSeawardNewDirection = NORTH;

                  // If can't do this, try to go straight on (to the E)
                  nXStraightOn = nX+1;
                  nYStraightOn = nY;

                  // If can't do either of these, try to go anti-seaward i.e. towards the RHS (S)
                  nXAntiSeaward = nX;
                  nYAntiSeaward = nY+1;
                  nAntiSeawardNewDirection = SOUTH;

                  // As a last resort, go back (to the W)
                  nXGoBack = nX-1;
                  nYGoBack = nY;
                  nGoBackNewDirection = WEST;

                  break;

               case SOUTH:
                  // The sea is towards the LHS (E) of the coast, so first try to go left (to the E)
                  nXSeaward = nX+1;
                  nYSeaward = nY;
                  nSeawardNewDirection = EAST;

                  // If can't do this, try to go straight on (to the S)
                  nXStraightOn = nX;
                  nYStraightOn = nY+1;

                  // If can't do either of these, try to go anti-seaward i.e. towards the RHS (W)
                  nXAntiSeaward = nX-1;
                  nYAntiSeaward = nY;
                  nAntiSeawardNewDirection = WEST;

                  // As a last resort, go back (to the N)
                  nXGoBack = nX;
                  nYGoBack = nY-1;
                  nGoBackNewDirection = NORTH;

                  break;

               case WEST :
                  // The sea is towards the LHS (S) of the coast, so first try to go left (to the S)
                  nXSeaward = nX;
                  nYSeaward = nY+1;
                  nSeawardNewDirection = SOUTH;

                  // If can't do this, try to go straight on (to the W)
                  nXStraightOn = nX-1;
                  nYStraightOn = nY;

                  // If can't do either of these, try to go anti-seaward i.e. towards the RHS (N)
                  nXAntiSeaward = nX;
                  nYAntiSeaward = nY-1;
                  nAntiSeawardNewDirection = NORTH;

                  // As a last resort, go back (to the E)
                  nXGoBack = nX+1;
                  nYGoBack = nY;
                  nGoBackNewDirection = EAST;

                  break;
            }
            break;
      }

      // Now do the actual search for this timestep: first try going in the direction of the sea. Is this seaward cell still within the grid?
      if (bIsWithinGrid(nXSeaward, nYSeaward))
      {
         // It is, so check if the cell in the seaward direction is a sea cell
         if (m_pRasterGrid->m_Cell[nXSeaward][nYSeaward].bIsInContiguousSea())
         {
            // There is sea in this seaward direction, so we are on the coast
            bAtCoast = true;
            
            // Has the current cell already marked been marked as a coast cell?
            if (! m_pRasterGrid->m_Cell[nX][nY].bIsCoastline())
            {
               // Not already marked, is this an intervention cell with the top above SWL?
               if ((m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION) && (m_pRasterGrid->m_Cell[nX][nY].dGetInterventionTopElev() > m_dThisTimestepSWL))
               {
                  // It is, so mark as coast and add it to the vector object
                  m_pRasterGrid->m_Cell[nX][nY].SetAsCoastline(true);
                  LTempGridCRS.Append(&Pti);                     
               }
               else if (m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev() > m_dThisTimestepSWL)
               {   
                  // The sediment top is above SWL so mark as coast and add it to the vector object
                  m_pRasterGrid->m_Cell[nX][nY].SetAsCoastline(true);
                  LTempGridCRS.Append(&Pti);
               }
            }
         }
         else
         {
            // The seaward cell is not a sea cell, so we will move to it next time
            nX = nXSeaward;
            nY = nYSeaward;

            // And set a new search direction, to keep turning seaward
            nSearchDirection = nSeawardNewDirection;
            continue;
         }
      }

      // OK, we couldn't move seaward (but we may have marked the current cell as coast) so next try to move straight on. Is this straight-ahead cell still within the grid?
      if (bIsWithinGrid(nXStraightOn, nYStraightOn))
      {
         // It is, so check if there is sea immediately in front
         if (m_pRasterGrid->m_Cell[nXStraightOn][nYStraightOn].bIsInContiguousSea())
         {
            // Sea is in front, so we are on the coast
            bAtCoast = true;            
            
            // Has the current cell already marked been marked as a coast cell?
            if (! m_pRasterGrid->m_Cell[nX][nY].bIsCoastline())
            {
               // Not already marked, is this an intervention cell with the top above SWL?
               if ((m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION) && (m_pRasterGrid->m_Cell[nX][nY].dGetInterventionTopElev() > m_dThisTimestepSWL))
               {
                  // It is, so mark as coast and add it to the vector object
                  m_pRasterGrid->m_Cell[nX][nY].SetAsCoastline(true);
                  LTempGridCRS.Append(&Pti);                     
               }
               else if (m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev() > m_dThisTimestepSWL)
               {   
                  // The sediment top is above SWL so mark as coast and add it to the vector object
                  m_pRasterGrid->m_Cell[nX][nY].SetAsCoastline(true);
                  LTempGridCRS.Append(&Pti);
               }
            }
         }
         else
         {
            // The straight-ahead cell is not a sea cell, so we will move to it next time
            nX = nXStraightOn;
            nY = nYStraightOn;

            // The search direction remains unchanged
            continue;
         }
      }

      // Couldn't move either seaward or straight on (but we may have marked the current cell as coast) so next try to move in the anti-seaward direction. Is this anti-seaward cell still within the grid?
      if (bIsWithinGrid(nXAntiSeaward, nYAntiSeaward))
      {
         // It is, so check if there is sea in this anti-seaward cell
         if (m_pRasterGrid->m_Cell[nXAntiSeaward][nYAntiSeaward].bIsInContiguousSea())
         {
            // There is sea on the anti-seaward side, so we are on the coast
            bAtCoast = true;
            
            // Has the current cell already marked been marked as a coast cell?
            if (! m_pRasterGrid->m_Cell[nX][nY].bIsCoastline())
            {
               // Not already marked, is this an intervention cell with the top above SWL?
               if ((m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION) && (m_pRasterGrid->m_Cell[nX][nY].dGetInterventionTopElev() > m_dThisTimestepSWL))
               {
                  // It is, so mark as coast and add it to the vector object
                  m_pRasterGrid->m_Cell[nX][nY].SetAsCoastline(true);
                  LTempGridCRS.Append(&Pti);                     
               }
               else if (m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev() > m_dThisTimestepSWL)
               {   
                  // The sediment top is above SWL so mark as coast and add it to the vector object
                  m_pRasterGrid->m_Cell[nX][nY].SetAsCoastline(true);
                  LTempGridCRS.Append(&Pti);
               }
            }
         }
         else
         {
            // The anti-seaward cell is not a sea cell, so we will move to it next time
            nX = nXAntiSeaward;
            nY = nYAntiSeaward;

            // And set a new search direction, to keep turning seaward
            nSearchDirection = nAntiSeawardNewDirection;
            continue;
         }
      }

      // Could not move to the seaward side, move straight ahead, or move to the anti-seaward side, so we must be in a single-cell dead end! As a last resort, turn round and move back to where we just came from
      nX = nXGoBack;
      nY = nYGoBack;

      // And change the search direction
      nSearchDirection = nGoBackNewDirection;
   }
   while (true);

   int nCoastSize = LTempGridCRS.nGetSize();
   
   if (nCoastSize == 0)
   {
      // Zero-length coastline
      LogStream << m_ulTimestep << ": " << WARN << "abandoning zero-length coastline from [" << nStartX << "][" << nStartY << "] {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "}" << endl;

      return RTN_OK;      
   }
   
   // OK, we have finished tracing this coastline on the grid, so check to see if the coastline is too short
   if (nCoastSize < m_nCoastMin)
   {
      // The vector coastline is unreasonably small, so abandon it
      LogStream << m_ulTimestep << ": ignoring temporary coastline from [" << nStartX << "][" << nStartY << "] {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} to [" << LTempGridCRS[nCoastSize-1].nGetX() << "][" << LTempGridCRS[nCoastSize-1].nGetY() << "] {" << dGridCentroidXToExtCRSX(LTempGridCRS[nCoastSize-1].nGetX()) << ", " << dGridCentroidYToExtCRSY(LTempGridCRS[nCoastSize-1].nGetY()) << "} since size (" << nCoastSize << ") is less than minimum (" << m_nCoastMin << ")" << endl;

      // Unmark these cells
      for (int n = 0; n < nCoastSize; n++)
         m_pRasterGrid->m_Cell[LTempGridCRS[n].nGetX()][LTempGridCRS[n].nGetY()].SetAsCoastline(false);

      return RTN_OK;
   }
   
   int 
      nEndX = nX,
      nEndY = nY,
      nCoastEndX = LTempGridCRS[nCoastSize-1].nGetX(),
      nCoastEndY = LTempGridCRS[nCoastSize-1].nGetY();
      
   if ((nCoastEndX != nEndX) || (nCoastEndY != nEndY))
   {
      // The grid-edge cell at nEndX, nEndY is not already at end of LTempGridCRS. But is the final cell in LTempGridCRS already at the edge of the grid?
      if (! m_pRasterGrid->m_Cell[nCoastEndX][nCoastEndY].bIsEdgeCell())
      {
         // The final cell in LTempGridCRS is not a grid-edge cell, so add the grid-edge cell and mark the cell as coastline
         LTempGridCRS.Append(nEndX, nEndY);
         nCoastSize++;
         
         m_pRasterGrid->m_Cell[nEndX][nEndY].SetAsCoastline(true);         
      }
   }

   // Need to specify start edge and end edge for smoothing routines
   int
      nStartEdge = m_pRasterGrid->m_Cell[nStartX][nStartY].nGetEdgeCell(),
      nEndEdge = m_pRasterGrid->m_Cell[nEndX][nEndY].nGetEdgeCell();

   // Next, convert the grid coordinates in LTempGridCRS (integer values stored as doubles) to external CRS coordinates (which will probably be non-integer, again stored as doubles). This is done now, so that smoothing is more effective
   CGeomLine LTempExtCRS;
   for (int j = 0; j < nCoastSize; j++)
      LTempExtCRS.Append(dGridCentroidXToExtCRSX(LTempGridCRS[j].nGetX()), dGridCentroidYToExtCRSY(LTempGridCRS[j].nGetY()));

   // Now do some smoothing of the vector output, if desired
   if (m_nCoastSmooth == SMOOTH_RUNNING_MEAN)
      LTempExtCRS = LSmoothCoastRunningMean(&LTempExtCRS, nStartEdge, nEndEdge);
   else if (m_nCoastSmooth == SMOOTH_SAVITZKY_GOLAY)
      LTempExtCRS = LSmoothCoastSavitzkyGolay(&LTempExtCRS, nStartEdge, nEndEdge);

   // Create a new coastline object and append to it the vector of coastline objects
   CRWCoast CoastTmp;
   m_VCoast.push_back(CoastTmp);
   int nCoast = m_VCoast.size()-1;

   CGeom2DPoint PtLast(DBL_MIN, DBL_MIN);
   for (int j = 0; j < nCoastSize; j++)
   {
      // Store the smoothed points (in external CRS) in the coast's m_LCoastline object, also append dummy values to the other attribute vectors
      if (PtLast != &LTempExtCRS[j])        // Avoid duplicate points
      {
         m_VCoast[nCoast].AppendToCoastline(LTempExtCRS[j].dGetX(), LTempExtCRS[j].dGetY());

         // Also store the locations of the corresponding unsmoothed points (in raster-grid CRS) in the coast's m_VCellsMarkedAsCoastline vector
         m_VCoast[nCoast].AppendCellMarkedAsCoastline(&LTempGridCRS[j]);
      }

      PtLast = LTempExtCRS[j];
   }

   // Next, set values for the coast's other attributes. First set the coast's handedness, and start and end edges
   m_VCoast[nCoast].SetSeaHandedness(nHandedness);
   m_VCoast[nCoast].SetStartEdge(nStartEdge);
   m_VCoast[nCoast].SetEndEdge(nEndEdge);

   LogStream << m_ulTimestep << ": coastline " << nCoast << " created, from [" << nStartX << "][" << nStartY << "] to [" << nEndX << "][" << nEndY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} to {" << dGridCentroidXToExtCRSX(nEndX) << ", " << dGridCentroidYToExtCRSY(nEndY) << "} with " << nCoastSize << " points, handedness = " << (nHandedness == LEFT_HANDED ? "left" : "right") << endl;
   LogStream << m_ulTimestep << ": smoothed coastline " << nCoast << " runs from {" << LTempExtCRS[0].dGetX() << ", " << LTempExtCRS[0].dGetY() << "} to {" << LTempExtCRS[nCoastSize-1].dGetX() << ", " << LTempExtCRS[nCoastSize-1].dGetY() << "} i.e. from the ";
   if (nStartEdge == NORTH)
      LogStream << "north";
   else if (nStartEdge == SOUTH)
      LogStream << "south";
   else if (nStartEdge == WEST)
      LogStream << "west";
   else if (nStartEdge == EAST)
      LogStream << "east";
   LogStream << " edge to the ";
   if (nEndEdge == NORTH)
      LogStream << "north";
   else if (nEndEdge == SOUTH)
      LogStream << "south";
   else if (nEndEdge == WEST)
      LogStream << "west";
   else if (nEndEdge == EAST)
      LogStream << "east";
   LogStream << " edge" << endl;   
//    LogStream << "-----------------" << endl;
//    for (int kk = 0; kk < m_VCoast.back().nGetCoastlineSize(); kk++)
//       LogStream << kk << " [" << m_VCoast.back().pPtiGetCellMarkedAsCoastline(kk)->nGetX() << "][" << m_VCoast.back().pPtiGetCellMarkedAsCoastline(kk)->nGetY() << "] {" << dGridCentroidXToExtCRSX(m_VCoast.back().pPtiGetCellMarkedAsCoastline(kk)->nGetX()) << ", " << dGridCentroidYToExtCRSY(m_VCoast.back().pPtiGetCellMarkedAsCoastline(kk)->nGetY()) << "}" << endl;
//    LogStream << "-----------------" << endl;  

   // Next calculate the curvature of the vector coastline
   DoCoastCurvature(nCoast, nHandedness);

   // Calculate values for the coast's flux orientation vector
   CalcCoastTangents(nCoast);
   
   return RTN_OK;
}

