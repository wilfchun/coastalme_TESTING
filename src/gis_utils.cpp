/*!
 *
 * \file gis_utils.cpp
 * \brief Various GIS-related functions. This version will build with GDAL version 2
 * \details TODO A more detailed description of these routines.
 * \author David Favis-Mortlock
 * \author Andres Payo

 * \date 2017
 * \copyright GNU General Public License
 *
 */

/*===============================================================================================================================

 This file is part of CoastalME, the Coastal Modelling Environment.

 CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

===============================================================================================================================*/
#include <assert.h>
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::ios;

#include <cmath>
#include <cfloat>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include "cme.h"
#include "simulation.h"
#include "simulation.h"
#include "raster_grid.h"


/*
 Notes re. co-ordinate systems used

 1. In the raster CRS, cell[0][0] is at the top left (NW) corner of the grid. Raster grid co-oordinate [0][0] is actually the top left (NW) corner of this cell.

 2. We assume that the grid CRS and external CRS have parallel axes. If they have not, see http://www.gdal.org/classGDALDataset.html which says that:

   To convert between pixel/line (P,L) raster space, and projection coordinates (Xp,Yp) space
      Xp = padfTransform[0] + P*padfTransform[1] + L*padfTransform[2];
      Yp = padfTransform[3] + P*padfTransform[4] + L*padfTransform[5];

   In a north-up image, padfTransform[1] is the pixel width, and padfTransform[5] is the pixel height. The upper left corner of the upper left pixel is at position
      (padfTransform[0], padfTransform[3]).

 3. Usually, raster grid CRS values are integer, i.e. they refer to a point which is at the centroid of a cell. They may also be -ve or greater than m_nXGridMax-1 i.e. may refer to a point which lies outside any cell of the raster grid.

*/

/*==============================================================================================================================

 Given the integer X-axis ordinate of a cell in the raster-grid CRS, returns the external CRS X-axis ordinate of the cell's centroid

===============================================================================================================================*/
double CSimulation::dGridCentroidXToExtCRSX(int const nGridX) const
{
   return (m_dGeoTransform[0] + (nGridX * m_dGeoTransform[1]) + (m_dGeoTransform[1] / 2));
}


/*==============================================================================================================================

 Given the integer Y-axis ordinate of a cell in the raster-grid CRS, returns the external CRS Y-axis ordinate of the cell's centroid

===============================================================================================================================*/
double CSimulation::dGridCentroidYToExtCRSY(int const nGridY) const
{
   return (m_dGeoTransform[3] + (nGridY * m_dGeoTransform[5]) + (m_dGeoTransform[5] / 2));
}

/*==============================================================================================================================

 Transforms a pointer to a CGeom2DIPoint in the raster-grid CRS (assumed to be the centroid of a cell) to the equivalent CGeom2DPoint in the external CRS

===============================================================================================================================*/
CGeom2DPoint CSimulation::PtGridCentroidToExt(CGeom2DIPoint const* pPtiIn) const
{
    int
      nGridX = pPtiIn->nGetX(),
      nGridY = pPtiIn->nGetY();

   double
      dX = m_dGeoTransform[0] + (nGridX * m_dGeoTransform[1]) + (m_dGeoTransform[1] / 2),
      dY = m_dGeoTransform[3] + (nGridY * m_dGeoTransform[5]) + (m_dGeoTransform[5] / 2);

   return CGeom2DPoint(dX, dY);
}


/*==============================================================================================================================

 Given a real-valued X-axis ordinate in the raster-grid CRS (i.e. not the centroid of a cell), returns the external CRS X-axis ordinate

===============================================================================================================================*/
double CSimulation::dGridXToExtCRSX(double const dGridX) const
{
   return (m_dGeoTransform[0] + (dGridX * m_dGeoTransform[1]));
}


/*==============================================================================================================================

 Given a real-valued Y-axis ordinate in the raster-grid CRS (i.e. not the centroid of a cell), returns the external CRS Y-axis ordinate

===============================================================================================================================*/
double CSimulation::dGridYToExtCRSY(double const dGridY) const
{
   return (m_dGeoTransform[3] + (dGridY * m_dGeoTransform[5]));
}


/*==============================================================================================================================

 Transforms an X-axis ordinate in the external CRS to the equivalent X-axis ordinate in the raster-grid CRS (the result may not be integer, and may be outside the grid)

===============================================================================================================================*/
double CSimulation::dExtCRSXToGridX(double const dExtCRSX) const
{
   return ((dExtCRSX - m_dGeoTransform[0]) / m_dGeoTransform[1]);
}


/*==============================================================================================================================

 Transforms a Y-axis ordinate in the external CRS to the equivalent Y-axis ordinate in the raster-grid CRS (the result may not be integer, and may be outside the grid)

===============================================================================================================================*/
double CSimulation::dExtCRSYToGridY(double const dExtCRSY) const
{
   return ((dExtCRSY - m_dGeoTransform[3]) / m_dGeoTransform[5]);
}


/*==============================================================================================================================

 Transforms a pointer to a CGeom2DPoint in the external CRS to the equivalent CGeom2DIPoint in the raster-grid CRS (both values rounded)

===============================================================================================================================*/
CGeom2DIPoint CSimulation::PtiExtCRSToGrid(CGeom2DPoint const* pPtIn) const
{
   double
      dX = pPtIn->dGetX(),
      dY = pPtIn->dGetY();

   int
      nX = dRound((dX - m_dGeoTransform[0]) / m_dGeoTransform[1]),
      nY = dRound((dY - m_dGeoTransform[3]) / m_dGeoTransform[5]);

   return CGeom2DIPoint(nX, nY);
}


/*==============================================================================================================================

 Returns the distance (in external CRS) between two points

===============================================================================================================================*/
double CSimulation::dGetDistanceBetween(CGeom2DPoint const* Pt1, CGeom2DPoint const* Pt2)
{
   double
      dXDist = Pt1->dGetX() - Pt2->dGetX(),
      dYDist = Pt1->dGetY() - Pt2->dGetY();

   return hypot(dXDist, dYDist);
}


/*==============================================================================================================================

 Returns the distance (in grid units) between two grid cell points

===============================================================================================================================*/
double CSimulation::dGetDistanceBetween(CGeom2DIPoint const* Pti1, CGeom2DIPoint const* Pti2)
{
   double
      dXDist = Pti1->nGetX() - Pti2->nGetX(),
      dYDist = Pti1->nGetY() - Pti2->nGetY();

   return hypot(dXDist, dYDist);
}


/*===============================================================================================================================

 Returns twice the signed area of a triangle

 ===============================================================================================================================*/
double CSimulation::dTriangleAreax2(CGeom2DPoint const* pPtA, CGeom2DPoint const* pPtB, CGeom2DPoint const* pPtC)
{
   return (pPtB->dGetX() - pPtA->dGetX()) * (pPtC->dGetY() - pPtA->dGetY()) - (pPtB->dGetY() - pPtA->dGetY()) * (pPtC->dGetX() - pPtA->dGetX());
}


/*==============================================================================================================================

 Checks whether the supplied point (an x-y pair, in the grid CRS) is within the raster grid, and is a valid cell (i.e. the basement DEM is not NODATA)

===============================================================================================================================*/
bool CSimulation::bIsWithinValidGrid(int const nX, int const nY) const
{
   if ((nX < 0) || (nX >= m_nXGridMax) || (nY < 0) || (nY >= m_nYGridMax) || m_pRasterGrid->m_Cell[nX][nY].bBasementElevIsMissingValue())
      return false;

   return true;
}


/*==============================================================================================================================

 Checks whether the supplied point (a reference to a CGeom2DIPoint, in the grid CRS) is within the raster grid, and is a valid cell (i.e. the basement DEM is not NODATA)

===============================================================================================================================*/
bool CSimulation::bIsWithinValidGrid(CGeom2DIPoint const* Pti) const
{
   int nX = Pti->nGetX();

   if ((nX < 0) || (nX >= m_nXGridMax))
      return false;

   int nY = Pti->nGetY();

   if ((nY < 0) || (nY >= m_nYGridMax))
      return false;

   if (m_pRasterGrid->m_Cell[nX][nY].bBasementElevIsMissingValue())
      return false;

   return true;
}


/*==============================================================================================================================

 Constrains the second supplied point (both are CGeom2DIPoints, in the grid CRS) to be an a valid cell within the raster grid

===============================================================================================================================*/
void CSimulation::KeepWithinValidGrid(CGeom2DIPoint const* Pti0, CGeom2DIPoint* Pti1) const
{
   KeepWithinValidGrid(Pti0->nGetX(), Pti0->nGetY(), *Pti1->pnGetX(), *Pti1->pnGetY());
}


/*==============================================================================================================================

 Given two points in the grid CRS (the first point asumed to be within the valid part of the raster grid, and the points assumed not to be coincident), this routine modifies the value of the second point so that it is on a line joining the original points and is a valid cell within the raster grid

 However in some cases (e.g. if the first point is at the edge of the valid part of the raster grid) then the second cell will be coincident with the first cell, and the line joining them is thus of zero length. The calling routine has to be able to handle this

===============================================================================================================================*/
void CSimulation::KeepWithinValidGrid(int const nX0, int const nY0, int& nX1, int& nY1) const
{
   int
      nDiffX = nX0 - nX1,
      nDiffY = nY0 - nY1;

   if (nDiffX == 0)
   {
      // The two points have the same x co-ords, so we just need to constrain the y co-ord
      if (nY1 < nY0)
      {
         nY1 = -1;

         do
         {
            nY1++;
         }
         while (m_pRasterGrid->m_Cell[nX1][nY1].bBasementElevIsMissingValue());

         return;
      }
      else
      {
         nY1 = m_nYGridMax;

         do
         {
            nY1--;
         }
         while (m_pRasterGrid->m_Cell[nX1][nY1].bBasementElevIsMissingValue());

         return;
      }
   }
   else if (nDiffY == 0)
   {
      // The two points have the same y co-ords, so we just need to constrain the x co-ord
      if (nX1 < nX0)
      {
         nX1 = -1;

         do
         {
            nX1++;
         }
         while (m_pRasterGrid->m_Cell[nX1][nY1].bBasementElevIsMissingValue());

         return;
      }
      else
      {
         nX1 = m_nXGridMax;

         do
         {
            nX1--;
         }
         while (m_pRasterGrid->m_Cell[nX1][nY1].bBasementElevIsMissingValue());

         return;
      }
   }
   else
   {
      // The two points have different x co-ords and different y co-ords, so we have to work harder. First find which of the co-ords is the greatest distance outside the grid, and constrain that co-ord for efficiency (since this will reduce the number of times round the loop). Note that both may be inside the grid, if the incorrect co-ord is in the invalid margin, in which case arbitrarily contrain the x co-ord
      int
         nXDistanceOutside = 0,
         nYDistanceOutside = 0;

      if (nX1 < 0)
         nXDistanceOutside = -nX1;
      else if (nX1 >= m_nXGridMax)
         nXDistanceOutside = nX1 - m_nXGridMax + 1;

      if (nY1 < 0)
         nYDistanceOutside = -nY1;
      else if (nY1 >= m_nYGridMax)
         nXDistanceOutside = nY1 - m_nYGridMax + 1;

      if (nXDistanceOutside >= nYDistanceOutside)
      {
         // Constrain the x co-ord
         if (nX1 < nX0)
         {
            // The incorrect x co-ord is less than the correct x co-ord: constrain it and find the y co-ord
            nX1 = -1;

            do
            {
               nX1++;

               nY1 = nY0 + dRound(((nX1 - nX0) * nDiffY) / static_cast<double>(nDiffX));
            }
            while ((nY1 < 0) || (nY1 >= m_nYGridMax) || (m_pRasterGrid->m_Cell[nX1][nY1].bBasementElevIsMissingValue()));

            return;
         }

         else
         {
            // The incorrect x co-ord is greater than the correct x-co-ord: constrain it and find the y co-ord
            nX1 = m_nXGridMax;

            do
            {
               nX1--;

               nY1 = nY0 + dRound(((nX1 - nX0) * nDiffY) / static_cast<double>(nDiffX));
            }
            while ((nY1 < 0) || (nY1 >= m_nYGridMax) || (m_pRasterGrid->m_Cell[nX1][nY1].bBasementElevIsMissingValue()));

            return;
         }
      }
      else
      {
         // Constrain the y co-ord
         if (nY1 < nY0)
         {
            // The incorrect y co-ord is less than the correct y-co-ord: constrain it and find the x co-ord
            nY1 = -1;

            do
            {
               nY1++;

               nX1 = nX0 + dRound(((nY1 - nY0) * nDiffX) / static_cast<double>(nDiffY));
            }
            while ((nX1 < 0) || (nX1 >= m_nXGridMax) || (m_pRasterGrid->m_Cell[nX1][nY1].bBasementElevIsMissingValue()));

            return;
         }

         else
         {
            // The incorrect y co-ord is greater than the correct y co-ord: constrain it and find the x co-ord
            nY1 = m_nYGridMax;

            do
            {
               nY1--;

               nX1 = nX0 + dRound(((nY1 - nY0) * nDiffX) / static_cast<double>(nDiffY));
            }
            while ((nX1 < 0) || (nX1 >= m_nXGridMax) || (m_pRasterGrid->m_Cell[nX1][nY1].bBasementElevIsMissingValue()));

            return;
         }
      }
   }
}


/*==============================================================================================================================

 Constrains the supplied angle to be within 0 and 360 degrees

===============================================================================================================================*/
double CSimulation::dKeepWithin360(double const dAngle)
{
   double dNewAngle = dAngle;

   // Sort out -ve angles
   while (dNewAngle < 0)
      dNewAngle += 360;

   // Sort out angles > 360
   while (dNewAngle > 360)
      dNewAngle -= 360;

   return dNewAngle;
}


/*==============================================================================================================================

 Returns a point (external CRS) which is the average of (i.e. is midway between) two other external CRS points

===============================================================================================================================*/
CGeom2DPoint CSimulation::PtAverage(CGeom2DPoint const* pPt1, CGeom2DPoint const* pPt2)
{
   double
      dPt1X = pPt1->dGetX(),
      dPt1Y = pPt1->dGetY(),
      dPt2X = pPt2->dGetX(),
      dPt2Y = pPt2->dGetY(),
      dPtAvgX = (dPt1X + dPt2X) / 2,
      dPtAvgY = (dPt1Y + dPt2Y) / 2;

   return CGeom2DPoint(dPtAvgX, dPtAvgY);
}


/*==============================================================================================================================

 Returns an integer point (grid CRS) which is the approximate average of (i.e. is midway between) two other grid CRS integer points

===============================================================================================================================*/
CGeom2DIPoint CSimulation::PtiAverage(CGeom2DIPoint const* pPti1, CGeom2DIPoint const* pPti2)
{
   int
      nPti1X = pPti1->nGetX(),
      nPti1Y = pPti1->nGetY(),
      nPti2X = pPti2->nGetX(),
      nPti2Y = pPti2->nGetY(),
      nPtiAvgX = dRound((nPti1X + nPti2X) / 2.0),
      nPtiAvgY = dRound((nPti1Y + nPti2Y) / 2.0);

   return CGeom2DIPoint(nPtiAvgX, nPtiAvgY);
}


/*==============================================================================================================================

 Returns an integer point (grid CRS) which is the weighted average of two other grid CRS integer points. The weight must be <= 1, if the weight is < 0.5 then the output point is closer to the first point, if the weight is > 0.5 then the output point is closer to the second point

===============================================================================================================================*/
CGeom2DIPoint CSimulation::PtiWeightedAverage(CGeom2DIPoint const* pPti1, CGeom2DIPoint const* pPti2, double const dWeight)
{
   int
      nPti1X = pPti1->nGetX(),
      nPti1Y = pPti1->nGetY(),
      nPti2X = pPti2->nGetX(),
      nPti2Y = pPti2->nGetY();
   double dOtherWeight = 1.0 - dWeight;

   int
      nPtiWeightAvgX = dRound(dWeight * nPti2X) + (dOtherWeight * nPti1X),
      nPtiWeightAvgY = dRound((dWeight * nPti2Y) + (dOtherWeight * nPti1Y));

   return CGeom2DIPoint(nPtiWeightAvgX, nPtiWeightAvgY);
}


/*==============================================================================================================================

 Returns a point (external CRS) which is the average of a vector of external CRS points

===============================================================================================================================*/
CGeom2DPoint CSimulation::PtAverage(vector<CGeom2DPoint>* pVIn)
{
   int nSize = pVIn->size();
   if (nSize == 0)
      return CGeom2DPoint(DBL_NODATA, DBL_NODATA);

   double
      dAvgX = 0,
      dAvgY = 0;

   for (int n = 0; n < nSize; n++)
   {
      dAvgX += pVIn->at(n).dGetX();
      dAvgY += pVIn->at(n).dGetY();
   }

   dAvgX /= nSize;
   dAvgY /= nSize;

   return CGeom2DPoint(dAvgX, dAvgY);
}


/*==============================================================================================================================

 Returns a vector which is perpendicular to an existing vector

===============================================================================================================================*/
// vector<CGeom2DPoint> CSimulation::VGetPerpendicular(CGeom2DPoint const* PtStart, CGeom2DPoint const* PtNext, double const dDesiredLength, int const nHandedness)
// {
//    // Returns a two-point vector which passes through PtStart with a scaled length
//    double dXLen = PtNext->dGetX() - PtStart->dGetX();
//    double dYLen = PtNext->dGetY() - PtStart->dGetY();
//
//    double dLength = hypot(dXLen, dYLen);
//    double dScaleFactor = dDesiredLength / dLength;
//
//    // The difference vector is (dXLen, dYLen), so the perpendicular difference vector is (-dYLen, dXLen) or (dYLen, -dXLen)
//    CGeom2DPoint EndPt;
//    if (nHandedness == RIGHT_HANDED)
//    {
//       EndPt.SetX(PtStart->dGetX() + (dScaleFactor * dYLen));
//       EndPt.SetY(PtStart->dGetY() - (dScaleFactor * dXLen));
//    }
//    else
//    {
//       EndPt.SetX(PtStart->dGetX() - (dScaleFactor * dYLen));
//       EndPt.SetY(PtStart->dGetY() + (dScaleFactor * dXLen));
//    }
//
//    vector<CGeom2DPoint> VNew;
//    VNew.push_back(*PtStart);
//    VNew.push_back(EndPt);
//    return VNew;
// }



/*==============================================================================================================================

 Returns a CGeom2DPoint which is the 'other' point of a two-point vector passing through PtStart, and which is perpendicular to the two-point vector from PtStart to PtNext

===============================================================================================================================*/
CGeom2DPoint CSimulation::PtGetPerpendicular(CGeom2DPoint const* PtStart, CGeom2DPoint const* PtNext, double const dDesiredLength, int const nHandedness)
{
   double
      dXLen = PtNext->dGetX() - PtStart->dGetX(),
      dYLen = PtNext->dGetY() - PtStart->dGetY(),
      dLength;

   if (dXLen == 0)
      dLength = dYLen;
   else if (dYLen == 0)
      dLength = dXLen;
   else
      dLength = hypot(dXLen, dYLen);

   double dScaleFactor = dDesiredLength / dLength;

   // The difference vector is (dXLen, dYLen), so the perpendicular difference vector is (-dYLen, dXLen) or (dYLen, -dXLen)
   CGeom2DPoint EndPt;
   if (nHandedness == RIGHT_HANDED)
   {
      EndPt.SetX(PtStart->dGetX() + (dScaleFactor * dYLen));
      EndPt.SetY(PtStart->dGetY() - (dScaleFactor * dXLen));
   }
   else
   {
      EndPt.SetX(PtStart->dGetX() - (dScaleFactor * dYLen));
      EndPt.SetY(PtStart->dGetY() + (dScaleFactor * dXLen));
   }

   return EndPt;
}


/*==============================================================================================================================

 Returns a CGeom2DIPoint (grid CRS) which is the 'other' point of a two-point vector passing through PtiStart, and which is perpendicular to the two-point vector from PtiStart to PtiNext

*==============================================================================================================================*/
CGeom2DIPoint CSimulation::PtiGetPerpendicular(CGeom2DIPoint const* PtiStart, CGeom2DIPoint const* PtiNext, double const dDesiredLength, int const nHandedness)
{
   double
      dXLen = PtiNext->nGetX() - PtiStart->nGetX(),
      dYLen = PtiNext->nGetY() - PtiStart->nGetY(),
      dLength;

   if (dXLen == 0)
      dLength = dYLen;
   else if (dYLen == 0)
      dLength = dXLen;
   else
      dLength = hypot(dXLen, dYLen);

   double dScaleFactor = dDesiredLength / dLength;

   // The difference vector is (dXLen, dYLen), so the perpendicular difference vector is (-dYLen, dXLen) or (dYLen, -dXLen)
   CGeom2DIPoint EndPti;
   if (nHandedness == RIGHT_HANDED)
   {
      EndPti.SetX(PtiStart->nGetX() + (dScaleFactor * dYLen));
      EndPti.SetY(PtiStart->nGetY() - (dScaleFactor * dXLen));
   }
   else
   {
      EndPti.SetX(PtiStart->nGetX() - (dScaleFactor * dYLen));
      EndPti.SetY(PtiStart->nGetY() + (dScaleFactor * dXLen));
   }

   return EndPti;
}


/*==============================================================================================================================

 Returns the signed angle BAC (in radians) subtended between three CGeom2DIPoints B A C. From http://stackoverflow.com/questions/3057448/angle-between-3-vertices

*==============================================================================================================================*/
double CSimulation::dAngleSubtended(CGeom2DIPoint const* pPtiA, CGeom2DIPoint const* pPtiB, CGeom2DIPoint const* pPtiC)
{
   int
      nXDistBtoA = pPtiB->nGetX() - pPtiA->nGetX(),
      nYDistBtoA = pPtiB->nGetY() - pPtiA->nGetY(),
      nXDistCtoA = pPtiC->nGetX() - pPtiA->nGetX(),
      nYDistCtoA = pPtiC->nGetY() - pPtiA->nGetY();

   double
      dDotProduct = nXDistBtoA * nXDistCtoA + nYDistBtoA * nYDistCtoA,
      dPseudoCrossProduct = nXDistBtoA * nYDistCtoA - nYDistBtoA * nXDistCtoA,
      dAngle = atan2(dPseudoCrossProduct, dDotProduct);

   return dAngle;
}


/*==============================================================================================================================

 Checks whether the selected raster GDAL driver supports file creation, 32-bit doubles, etc.

===============================================================================================================================*/
bool CSimulation::bCheckRasterGISOutputFormat(void)
{
   // Register all available GDAL raster and vector drivers (GDAL 2)
   GDALAllRegister();

   // If the user hasn't specified a GIS output format, assume that we will use the same GIS format as the input basement DEM
   if (m_strRasterGISOutFormat.empty())
      m_strRasterGISOutFormat = m_strGDALBasementDEMDriverCode;

   // Load the raster GDAL driver
   GDALDriver* pDriver = GetGDALDriverManager()->GetDriverByName(m_strRasterGISOutFormat.c_str());
   if (NULL == pDriver)
   {
      // Can't load raster GDAL driver. Incorrectly specified?
      cerr << ERR << "Unknown raster GIS output format '" << m_strRasterGISOutFormat << "'." << endl;
      return false;
   }

   // Get the metadata for this raster driver
   char** papszMetadata = pDriver->GetMetadata();

//    for (int i = 0; papszMetadata[i] != NULL; i++)
//       cout << papszMetadata[i] << endl;
//    cout << endl;

   // For GDAL2, need to test if this is a raster driver
   if (! CSLFetchBoolean(papszMetadata, GDAL_DCAP_RASTER, FALSE))
   {
      // This is not a raster driver
      cerr << ERR << "GDAL driver '" << m_strRasterGISOutFormat << "' is not a raster driver. Choose another format." << endl;
      return false;
   }

   // This driver is OK, so store its longname and the default file extension
   m_strGDALRasterOutputDriverLongname  = CSLFetchNameValue(papszMetadata, "DMD_LONGNAME");
   m_strGDALRasterOutputDriverExtension = CSLFetchNameValue(papszMetadata, "DMD_EXTENSION");

   // Set up any defaults for raster files that are created using this driver
   SetRasterFileCreationDefaults();

   // Now do various tests of the driver's capabilities
   if (! CSLFetchBoolean(papszMetadata, GDAL_DCAP_CREATE, FALSE))
   {
      // This raster driver does not support the Create() method, does it support CreateCopy()?
      if (! CSLFetchBoolean(papszMetadata, GDAL_DCAP_CREATECOPY, FALSE))
      {
         cerr << ERR << "Cannot write using raster GDAL driver '" << m_strRasterGISOutFormat << " since neither Create() or CreateCopy() are supported'. Choose another GDAL raster format." << endl;
         return false;
      }

      // Can't use Create() but can use CreateCopy()
      m_bGDALCanCreate = false;
   }

   // Next, test to see what data types the driver can write and from this, work out the largest int and float we can write
   if (strstr(CSLFetchNameValue(papszMetadata, "DMD_CREATIONDATATYPES"), "Float"))
   {
      m_bGDALCanWriteFloat = true;
      m_GDALWriteFloatDataType = GDT_Float32;
   }

   if (strstr(CSLFetchNameValue(papszMetadata, "DMD_CREATIONDATATYPES"), "UInt32"))
   {
      m_bGDALCanWriteInt32 = true;

      m_GDALWriteIntDataType = GDT_UInt32;
      m_lGDALMaxCanWrite = UINT32_MAX;
      m_lGDALMinCanWrite = 0;

      if (! m_bGDALCanWriteFloat)
         m_GDALWriteFloatDataType = GDT_UInt32;

      return true;
   }

   if (strstr(CSLFetchNameValue(papszMetadata, "DMD_CREATIONDATATYPES"), "Int32"))
   {
      m_bGDALCanWriteInt32 = true;

      m_GDALWriteIntDataType = GDT_Int32;
      m_lGDALMaxCanWrite = INT32_MAX;
      m_lGDALMinCanWrite = INT32_MIN;

      if (! m_bGDALCanWriteFloat)
         m_GDALWriteFloatDataType = GDT_Int32;

      return true;
   }

   if (strstr(CSLFetchNameValue(papszMetadata, "DMD_CREATIONDATATYPES"), "UInt16"))
   {
      m_bGDALCanWriteInt32 = false;

      m_GDALWriteIntDataType = GDT_UInt16;
      m_lGDALMaxCanWrite = UINT16_MAX;
      m_lGDALMinCanWrite = 0;

      if (! m_bGDALCanWriteFloat)
         m_GDALWriteFloatDataType = GDT_UInt16;

      return true;
   }

   if (strstr(CSLFetchNameValue(papszMetadata, "DMD_CREATIONDATATYPES"), "Int16"))
   {
      m_bGDALCanWriteInt32 = false;

      m_GDALWriteIntDataType = GDT_Int16;
      m_lGDALMaxCanWrite = INT16_MAX;
      m_lGDALMinCanWrite = INT16_MIN;

      if (! m_bGDALCanWriteFloat)
         m_GDALWriteFloatDataType = GDT_Int16;

      return true;
   }

   if (strstr(CSLFetchNameValue(papszMetadata, "DMD_CREATIONDATATYPES"), "Byte"))
   {
      m_bGDALCanWriteInt32 = false;

      m_GDALWriteIntDataType = GDT_Byte;
      m_lGDALMaxCanWrite = UINT8_MAX;
      m_lGDALMinCanWrite = 0;

      if (! m_bGDALCanWriteFloat)
         m_GDALWriteFloatDataType = GDT_Byte;

      return true;
   }

   // This driver does not even support byte output
   cerr << ERR << "Cannot write using raster GDAL driver '" << m_strRasterGISOutFormat << ", not even byte output is supported'. Choose another GIS raster format." << endl;
   return false;
}


/*==============================================================================================================================

 Checks whether the selected vector OGR driver supports file creation etc.

===============================================================================================================================*/
bool CSimulation::bCheckVectorGISOutputFormat(void)
{
   // Load the vector GDAL driver (NOTE this assumes that GDALAllRegister() has already been called)
   GDALDriver* pDriver = GetGDALDriverManager()->GetDriverByName(m_strVectorGISOutFormat.c_str());
   if (NULL == pDriver)
   {
      // Can't load vector GDAL driver. Incorrectly specified?
      cerr << ERR << "Unknown vector GIS output format '" << m_strVectorGISOutFormat << "'." << endl;
      return false;
   }

   // Get the metadata for this vector driver
   char** papszMetadata = pDriver->GetMetadata();

   // For GDAL2, need to test if this is a vector driver
   if (! CSLFetchBoolean(papszMetadata, GDAL_DCAP_VECTOR, FALSE))
   {
      // This is not a vector driver
      cerr << ERR << "GDAL driver '" << m_strVectorGISOutFormat << "' is not a vector driver. Choose another format." << endl;
      return false;
   }

   if (! CSLFetchBoolean(papszMetadata, GDAL_DCAP_CREATE, FALSE))
   {
      // Driver does not support create() method
      cerr << ERR << "Cannot write vector GIS files using GDAL driver '" << m_strRasterGISOutFormat << "'. Choose another format." << endl;
      return false;
   }

   // Driver is OK, now set some options for individual drivers
   if (m_strVectorGISOutFormat == "ESRI Shapefile")
   {
      // Set this, so that just a single dataset-with-one-layer shapefile is created, rather than a directory
      // (see http://www.gdal.org/ogr/drv_shapefile.html)
      m_strOGRVectorOutputExtension = ".shp";

   }
   // TODO Others

   return true;
}


/*==============================================================================================================================

 The bSaveAllRasterGISFiles member function saves the raster GIS files using values from the RasterGrid array

==============================================================================================================================*/
bool CSimulation::bSaveAllRasterGISFiles(void)
{
   // Increment file number
   m_nGISSave++;

   // Set for next save
   if (m_bSaveRegular)
      m_dRSaveTime += m_dRSaveInterval;
   else
      m_nThisSave = tMin(++m_nThisSave, m_nUSave);

   // These are always written
   if (! bWriteRasterGISFloat(PLOT_SEDIMENT_TOP_ELEV, &PLOT_SEDIMENT_TOP_ELEV_TITLE))
      return false;

   // These are always written
   if (! bWriteRasterGISFloat(PLOT_TOP_ELEV, &PLOT_TOP_ELEV_TITLE))
      return false;

   if (! bWriteRasterGISFloat(PLOT_LOCAL_CONS_SLOPE, &PLOT_LOCAL_CONS_SLOPE_TITLE))
      return false;

   if (! bWriteRasterGISFloat(PLOT_SEA_DEPTH, &PLOT_SEA_DEPTH_TITLE))
      return false;

   if (! bWriteRasterGISFloat(PLOT_WAVE_HEIGHT, &PLOT_WAVE_HEIGHT_TITLE))
      return false;

   if (! bWriteRasterGISFloat(PLOT_BEACH_PROTECTION, &PLOT_BEACH_PROTECTION_TITLE))
      return false;

   if (! bWriteRasterGISFloat(PLOT_POTENTIAL_PLATFORM_EROSION, &PLOT_POTENTIAL_PLATFORM_EROSION_TITLE))
      return false;

   if (! bWriteRasterGISFloat(PLOT_ACTUAL_PLATFORM_EROSION, &PLOT_ACTUAL_PLATFORM_EROSION_TITLE))
      return false;

   if (! bWriteRasterGISFloat(PLOT_TOTAL_POTENTIAL_PLATFORM_EROSION, &PLOT_TOTAL_POTENTIAL_PLATFORM_EROSION_TITLE))
      return false;

   if (! bWriteRasterGISFloat(PLOT_TOTAL_ACTUAL_PLATFORM_EROSION, &PLOT_TOTAL_ACTUAL_PLATFORM_EROSION_TITLE))
      return false;

   if (! bWriteRasterGISFloat(PLOT_POTENTIAL_BEACH_EROSION, &PLOT_POTENTIAL_BEACH_EROSION_TITLE))
      return false;

   if (! bWriteRasterGISFloat(PLOT_ACTUAL_BEACH_EROSION, &PLOT_ACTUAL_BEACH_EROSION_TITLE))
      return false;

   if (! bWriteRasterGISFloat(PLOT_TOTAL_POTENTIAL_BEACH_EROSION, &PLOT_TOTAL_POTENTIAL_BEACH_EROSION_TITLE))
      return false;

   if (! bWriteRasterGISFloat(PLOT_TOTAL_ACTUAL_BEACH_EROSION, &PLOT_TOTAL_ACTUAL_BEACH_EROSION_TITLE))
      return false;

   if (! bWriteRasterGISFloat(PLOT_BEACH_DEPOSITION, &PLOT_BEACH_DEPOSITION_TITLE))
      return false;

   if (! bWriteRasterGISFloat(PLOT_TOTAL_BEACH_DEPOSITION, &PLOT_TOTAL_BEACH_DEPOSITION_TITLE))
      return false;

   if (! bWriteRasterGISInt(PLOT_LANDFORM, &PLOT_LANDFORM_TITLE))
      return false;

   // These are optional
   if (m_bAvgWaveHeightSave)
      if (! bWriteRasterGISFloat(PLOT_AVG_WAVE_HEIGHT, &PLOT_AVG_WAVE_HEIGHT_TITLE))
         return false;

   if (m_bAvgSeaDepthSave)
      if (! bWriteRasterGISFloat(PLOT_AVG_SEA_DEPTH, &PLOT_AVG_SEA_DEPTH_TITLE))
      return false;

   if (m_bSuspSedSave)
      if (! bWriteRasterGISFloat(PLOT_SUSPENDED_SEDIMENT, &PLOT_SUSPENDED_SEDIMENT_TITLE))
         return false;

   if (m_bAvgSuspSedSave)
       if (! bWriteRasterGISFloat(PLOT_AVG_SUSPENDED_SEDIMENT, &PLOT_AVG_SUSPENDED_SEDIMENT_TITLE))
         return false;

   if (m_bBasementElevSave)
      if (! bWriteRasterGISFloat(PLOT_BASEMENT_ELEV, &PLOT_BASEMENT_ELEV_TITLE))
         return false;

   for (int nLayer = 0; nLayer < m_nLayers; nLayer++)
   {
      if (m_bFineUnconsSedSave)
      {
         if (! bWriteRasterGISFloat(PLOT_FINEUNCONSSED, &PLOT_FINEUNCONSSED_TITLE, nLayer))
            return false;
      }

      if (m_bSandUnconsSedSave)
      {
         if (! bWriteRasterGISFloat(PLOT_SANDUNCONSSED, &PLOT_SANDUNCONSSED_TITLE, nLayer))
            return false;
      }

      if (m_bCoarseUnconsSedSave)
      {
         if (! bWriteRasterGISFloat(PLOT_COARSEUNCONSSED, &PLOT_COARSEUNCONSSED_TITLE, nLayer))
            return false;
      }

      if (m_bFineConsSedSave)
      {
         if (! bWriteRasterGISFloat(PLOT_FINECONSSED, &PLOT_FINECONSSED_TITLE, nLayer))
            return false;
      }

      if (m_bSandConsSedSave)
      {
         if (! bWriteRasterGISFloat(PLOT_SANDCONSSED, &PLOT_SANDCONSSED_TITLE, nLayer))
            return false;
      }

      if (m_bCoarseConsSedSave)
      {
         if (! bWriteRasterGISFloat(PLOT_COARSECONSSED, &PLOT_COARSECONSSED_TITLE, nLayer))
            return false;
      }
   }

   if (m_bSliceSave)
   {
      for (int i = 0; i < static_cast<int>(m_VdSliceElev.size()); i++)
      {
         if (! bWriteRasterGISInt(PLOT_SLICE, &PLOT_SLICE_TITLE, m_VdSliceElev[i]))
            return false;
      }
   }

   if (m_bRasterCoastlineSave)
   {
      if (! bWriteRasterGISInt(PLOT_RASTER_COAST, &PLOT_RASTER_COAST_TITLE))
         return false;
   }

   if (m_bRasterNormalSave)
   {
      if (! bWriteRasterGISInt(PLOT_RASTER_NORMAL, &PLOT_RASTER_NORMAL_TITLE))
         return false;
   }

   if (m_bActiveZoneSave)
   {
      if (! bWriteRasterGISInt(PLOT_ACTIVE_ZONE, &PLOT_ACTIVE_ZONE_TITLE))
         return false;
   }

   if (m_bCliffCollapseSave)
   {
      if (! bWriteRasterGISFloat(PLOT_CLIFF_COLLAPSE, &PLOT_CLIFF_COLLAPSE_TITLE))
         return false;
   }

   if (m_bTotCliffCollapseSave)
   {
      if (! bWriteRasterGISFloat(PLOT_TOTAL_CLIFF_COLLAPSE, &PLOT_TOTAL_CLIFF_COLLAPSE_TITLE))
         return false;
   }

   if (m_bCliffCollapseDepositionSave)
   {
      if (! bWriteRasterGISFloat(PLOT_CLIFF_COLLAPSE_DEPOSIT, &PLOT_CLIFF_COLLAPSE_DEPOSIT_TITLE))
         return false;
   }

   if (m_bTotCliffCollapseDepositionSave)
   {
      if (! bWriteRasterGISFloat(PLOT_TOTAL_CLIFF_COLLAPSE_DEPOSIT, &PLOT_TOTAL_CLIFF_COLLAPSE_DEPOSIT_TITLE))
         return false;
   }

   if (m_bRasterPolygonSave)
   {
      if (! bWriteRasterGISInt(PLOT_RASTER_POLYGON, &PLOT_RASTER_POLYGON_TITLE))
         return false;
   }

   if (m_bPotentialPlatformErosionMaskSave)
   {
      if (! bWriteRasterGISInt(PLOT_POTENTIAL_PLATFORM_EROSION_MASK, &PLOT_POTENTIAL_PLATFORM_EROSION_MASK_TITLE))
      return false;
   }

   if (m_bSeaMaskSave)
   {
      if (! bWriteRasterGISInt(PLOT_INUNDATION_MASK, &PLOT_INUNDATION_MASK_TITLE))
      return false;
   }

   if (m_bBeachMaskSave)
   {
      if (! bWriteRasterGISInt(PLOT_BEACH_MASK, &PLOT_BEACH_MASK_TITLE))
      return false;
   }

   if (m_bInterventionClassSave)
   {
      if (! bWriteRasterGISInt(PLOT_INTERVENTION_CLASS, &PLOT_INTERVENTION_CLASS_TITLE))
         return false;
   }

   if (m_bInterventionHeightSave)
   {
      if (! bWriteRasterGISFloat(PLOT_INTERVENTION_HEIGHT, &PLOT_INTERVENTION_HEIGHT_TITLE))
         return false;
   }

   if (m_bShadowZoneCodesSave)
   {
      if (! bWriteRasterGISInt(PLOT_SHADOW_ZONE_CODES, &PLOT_SHADOW_ZONE_CODES_TITLE))
         return false;
   }
   
   if (m_bDeepWaterWaveOrientationSave)
   {
      if (! bWriteRasterGISFloat(PLOT_DEEP_WATER_WAVE_ORIENTATION, &PLOT_DEEP_WATER_WAVE_ORIENTATION_TITLE))
         return false;
   }
   
   if (m_bDeepWaterWaveHeightSave)
   {
      if (! bWriteRasterGISFloat(PLOT_DEEP_WATER_WAVE_HEIGHT, &PLOT_DEEP_WATER_WAVE_HEIGHT_TITLE))
         return false;
   }
   
   return true;
}


/*==============================================================================================================================

 The bSaveAllvectorGISFiles member function saves the vector GIS files

==============================================================================================================================*/
bool CSimulation::bSaveAllVectorGISFiles(void)
{
   // Always written
   if (! bWriteVectorGIS(PLOT_COAST, &PLOT_COAST_TITLE))
      return false;

   if (! bWriteVectorGIS(PLOT_NORMALS, &PLOT_NORMALS_TITLE))
      return false;

   if (m_bInvalidNormalsSave)
   {
      if (! bWriteVectorGIS(PLOT_INVALID_NORMALS, &PLOT_INVALID_NORMALS_TITLE))
         return false;
   }

   if (m_bCoastCurvatureSave)
   {
      if (! bWriteVectorGIS(PLOT_COAST_CURVATURE, &PLOT_COAST_CURVATURE_TITLE))
         return false;
   }

   if (m_bWaveAngleSave)
   {
      if (! bWriteVectorGIS(PLOT_WAVE_AND_HEIGHT, &PLOT_WAVE_AND_HEIGHT_TITLE))
         return false;
   }

   if (m_bAvgWaveAngleSave)
   {
      if (! bWriteVectorGIS(PLOT_AVG_WAVE_AND_HEIGHT, &PLOT_AVG_WAVE_AND_HEIGHT_TITLE))
         return false;
   }

   if (m_bWaveEnergySinceCollapseSave)
   {
      if (! bWriteVectorGIS(PLOT_WAVE_ENERGY_SINCE_COLLAPSE, &PLOT_WAVE_ENERGY_SINCE_COLLAPSE_TITLE))
         return false;
   }

   if (m_bMeanWaveEnergySave)
   {
      if (! bWriteVectorGIS(PLOT_MEAN_WAVE_ENERGY, &PLOT_MEAN_WAVE_ENERGY_TITLE))
         return false;
   }

   if (m_bBreakingWaveHeightSave)
   {
      if (! bWriteVectorGIS(PLOT_BREAKING_WAVE_HEIGHT, &PLOT_BREAKING_WAVE_HEIGHT_TITLE))
         return false;
   }

   if (m_bBeachProtectionSave)
   {
      if (! bWriteVectorGIS(PLOT_BEACH_PROTECTION, &PLOT_BEACH_PROTECTION_TITLE))
         return false;
   }
   if (m_bPolygonNodeSave)
   {
      if (! bWriteVectorGIS(PLOT_POLYGON_NODES, &PLOT_POLYGON_NODES_TITLE))
         return false;
   }

   if (m_bPolygonBoundarySave)
   {
      if (! bWriteVectorGIS(PLOT_POLYGON_BOUNDARY, &PLOT_POLYGON_BOUNDARY_TITLE))
         return false;
   }

   if (m_bCliffNotchSave)
   {
      if (! bWriteVectorGIS(PLOT_CLIFF_NOTCH_SIZE, &PLOT_CLIFF_NOTCH_SIZE_TITLE))
         return false;
   }

   if (m_bShadowZoneLineSave)
   {
      if (! bWriteVectorGIS(PLOT_SHADOW_ZONE_BOUNDARY, &PLOT_SHADOW_ZONE_BOUNDARY_TITLE))
         return false;
   }

   return true;
}


/*==============================================================================================================================

 Finds the max and min values in order to scale raster output if we cannot write floats or 32-bit integers

==============================================================================================================================*/
void CSimulation::GetRasterOutputMinMax(int const nDataItem, double& dMin, double& dMax, int const nLayer, double const dElev)
{
   // If this is a binary mask layer, we already know the max and min values
   if ((nDataItem == PLOT_POTENTIAL_PLATFORM_EROSION_MASK) ||
       (nDataItem == PLOT_INUNDATION_MASK) ||
       (nDataItem == PLOT_BEACH_MASK) ||
       (nDataItem == PLOT_RASTER_COAST) ||
       (nDataItem == PLOT_RASTER_NORMAL) ||
       (nDataItem == PLOT_ACTIVE_ZONE))
   {
      dMin = 0;
      dMax = 1;

      return;
   }

   // Not a binary mask layer, so we must find the max and min values
   dMin = DBL_MAX;
   dMax = DBL_MIN;

   double dTmp;
   for (int nY = 0; nY < m_nYGridMax; nY++)
   {
      for (int nX = 0; nX < m_nXGridMax; nX++)
      {
         switch (nDataItem)
         {
            case (PLOT_SLICE):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].nGetLayerAtElev(dElev);
               break;
            }

            case (PLOT_LANDFORM):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory();
               break;
            }

            case (PLOT_INTERVENTION_CLASS):
            {
               dTmp = INT_NODATA;
               if (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_INTERVENTION)
                  dTmp = m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFSubCategory();
               break;
            }

            case (PLOT_INTERVENTION_HEIGHT):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetInterventionHeight();
               break;
            }

            case (PLOT_RASTER_POLYGON):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();
               break;
            }

            case (PLOT_BASEMENT_ELEV):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetBasementElev();
               break;
            }

            case (PLOT_SEDIMENT_TOP_ELEV):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev();
               break;
            }

            case (PLOT_TOP_ELEV):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].nGetShadowZoneCode();
               break;
            }

            case (PLOT_LOCAL_CONS_SLOPE):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetLocalConsSlope();
               break;
            }

            case (PLOT_SEA_DEPTH):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth();
               break;
            }

            case (PLOT_AVG_SEA_DEPTH):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTotSeaDepth() / m_ulTimestep;
               break;
            }

            case (PLOT_WAVE_HEIGHT):
            {
               if (! m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea())
                  dTmp = DBL_NODATA;
               else
                  dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight();
               break;
            }

            case (PLOT_AVG_WAVE_HEIGHT):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTotWaveHeight() / m_ulTimestep;
               break;
            }

            case (PLOT_BEACH_PROTECTION):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetBeachProtectionFactor();
               if (dTmp != DBL_NODATA)
                  dTmp = 1 - dTmp;                 // Output the inverse, seems more intuitive
               break;
            }

            case (PLOT_POTENTIAL_PLATFORM_EROSION):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetPotentialPlatformErosion();
               break;
            }

            case (PLOT_ACTUAL_PLATFORM_EROSION):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetActualPlatformErosion();
               break;
            }

            case (PLOT_TOTAL_POTENTIAL_PLATFORM_EROSION):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTotPotentialPlatformErosion();
               break;
            }

            case (PLOT_TOTAL_ACTUAL_PLATFORM_EROSION):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTotActualPlatformErosion();
               break;
            }

            case (PLOT_POTENTIAL_BEACH_EROSION):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetPotentialBeachErosion();
               break;
            }

            case (PLOT_ACTUAL_BEACH_EROSION):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetActualBeachErosion();
               break;
            }

            case (PLOT_TOTAL_POTENTIAL_BEACH_EROSION):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTotPotentialBeachErosion();
               break;
            }

            case (PLOT_TOTAL_ACTUAL_BEACH_EROSION):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTotActualBeachErosion();
               break;
            }

            case (PLOT_BEACH_DEPOSITION):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetBeachDeposition();
               break;
            }

            case (PLOT_TOTAL_BEACH_DEPOSITION):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTotBeachDeposition();
               break;
            }

            case (PLOT_SUSPENDED_SEDIMENT):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetSuspendedSediment();
               break;
            }

            case (PLOT_AVG_SUSPENDED_SEDIMENT):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTotSuspendedSediment() / m_ulTimestep;
               break;
            }

            case (PLOT_FINEUNCONSSED):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nLayer)->pGetUnconsolidatedSediment()->dGetFine();
               break;
            }

            case (PLOT_SANDUNCONSSED):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nLayer)->pGetUnconsolidatedSediment()->dGetSand();
               break;
            }

            case (PLOT_COARSEUNCONSSED):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nLayer)->pGetUnconsolidatedSediment()->dGetCoarse();
               break;
            }

            case (PLOT_FINECONSSED):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nLayer)->pGetConsolidatedSediment()->dGetFine();
               break;
            }

            case (PLOT_SANDCONSSED):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nLayer)->pGetConsolidatedSediment()->dGetSand();
               break;
            }

            case (PLOT_COARSECONSSED):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nLayer)->pGetConsolidatedSediment()->dGetCoarse();
               break;
            }

            case (PLOT_CLIFF_COLLAPSE):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetCliffCollapse();
               break;
            }

            case (PLOT_TOTAL_CLIFF_COLLAPSE):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTotCliffCollapse();
               break;
            }

            case (PLOT_CLIFF_COLLAPSE_DEPOSIT):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetCliffCollapseDeposition();
               break;
            }

            case (PLOT_TOTAL_CLIFF_COLLAPSE_DEPOSIT):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTotCliffCollapseDeposition();
               break;
            }

            case (PLOT_SHADOW_ZONE_CODES):
            {
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTotCliffCollapseDeposition();
               break;
            }
         }

         if (dTmp != DBL_NODATA)
         {
            if (dTmp > dMax)
               dMax = dTmp;

            if (dTmp < dMin)
               dMin = dTmp;
         }
      }
   }
}


/*==============================================================================================================================

 Sets per-driver defaults for raster files created using GDAL

===============================================================================================================================*/
void CSimulation::SetRasterFileCreationDefaults(void)
{
   string
      strDriver = strToLower(&m_strRasterGISOutFormat),
      strComment =  "Created by " + PROGNAME + " for " + PLATFORM + " " + strGetBuild() + " running on " + strGetComputerName();

   // TODO Do these for all commonly-used file types
   if (strDriver == "aaigrid")
   {

   }

   else if (strDriver == "bmp")
   {

   }

   else if (strDriver == "gtiff")
   {
      if (m_bWorldFile)
         m_papszGDALRasterOptions = CSLSetNameValue(m_papszGDALRasterOptions, "TFW", "YES");

//       if (m_bCompressGTIFF)
//       {
//          m_papszGDALRasterOptions = CSLSetNameValue(m_papszGDALRasterOptions, "NUM_THREADS", "ALL_CPUS");
//          m_papszGDALRasterOptions = CSLSetNameValue(m_papszGDALRasterOptions, "COMPRESS", "LZW");
//       }
   }

   else if (strDriver == "hfa")
   {
      m_papszGDALRasterOptions = CSLSetNameValue(m_papszGDALRasterOptions, "NBITS", "4");
   }

   else if (strDriver == "jpeg")
   {
      m_papszGDALRasterOptions = CSLSetNameValue(m_papszGDALRasterOptions, "COMMENT", strComment.c_str());
      m_papszGDALRasterOptions = CSLSetNameValue(m_papszGDALRasterOptions, "QUALITY", "95");
   }

   else if (strDriver == "png")
   {
      if (m_bWorldFile)
         m_papszGDALRasterOptions = CSLSetNameValue(m_papszGDALRasterOptions, "WORLDFILE", "YES");

//       m_papszGDALRasterOptions = CSLSetNameValue(m_papszGDALRasterOptions, "TITLE", "This is the title");
//       m_papszGDALRasterOptions = CSLSetNameValue(m_papszGDALRasterOptions, "DESCRIPTION", "This is a description");
//       m_papszGDALRasterOptions = CSLSetNameValue(m_papszGDALRasterOptions, "COPYRIGHT", "This is some copyright statement");
      m_papszGDALRasterOptions = CSLSetNameValue(m_papszGDALRasterOptions, "COMMENT", strComment.c_str());
      m_papszGDALRasterOptions = CSLSetNameValue(m_papszGDALRasterOptions, "NBITS", "4");
   }

   else if (strDriver == "rst")
   {

   }

//    for (int i = 0; m_papszGDALRasterOptions[i] != NULL; i++)
//       cout << m_papszGDALRasterOptions[i] << endl;
//    cout << endl;

}


/*==============================================================================================================================

 Returns the opposite diection

===============================================================================================================================*/
int CSimulation::nGetOppositeDirection(int const nDirection)
{
   switch (nDirection)
   {
      case NORTH:
         return SOUTH;

      case NORTH_EAST:
         return SOUTH-WEST;

      case EAST:
         return WEST;

      case SOUTH_EAST:
         return NORTH-WEST;

      case SOUTH:
         return NORTH;

      case SOUTH_WEST:
         return NORTH-EAST;

      case WEST:
         return EAST;

      case NORTH_WEST:
         return SOUTH-EAST;
   }

   // Should never get here
   return NO_DIRECTION;
}
