/*!
 *
 * \file coast.cpp
 * \brief CRWCoast routines
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
//#include <assert.h>

#include "cme.h"
#include "coast.h"
#include "line.h"
#include "i_line.h"


CRWCoast::CRWCoast(void)
:  m_nSeaHandedness(NULL_HANDED),
   m_nStartEdge(INT_NODATA),
   m_nEndEdge(INT_NODATA),
   m_dCurvatureDetailedMean(0),
   m_dCurvatureDetailedSTD(0),
   m_dCurvatureSmoothMean(0),
   m_dCurvatureSmoothSTD(0)
{
}

CRWCoast::~CRWCoast(void)
{
   for (unsigned int i = 0; i < m_pVLandforms.size(); i++)
      delete m_pVLandforms[i];

   for (unsigned int i = 0; i < m_pVPolygon.size(); i++)
      delete m_pVPolygon[i];
}


void CRWCoast::SetSeaHandedness(int const nNewHandedness)
{
   m_nSeaHandedness = nNewHandedness;
}

int CRWCoast::nGetSeaHandedness(void) const
{
   return m_nSeaHandedness;
}


void CRWCoast::SetStartEdge(int const nEdge)
{
   m_nStartEdge = nEdge;
}

int CRWCoast::nGetStartEdge(void) const
{
   return m_nStartEdge;
}


void CRWCoast::SetEndEdge(int const nEdge)
{
   m_nEndEdge = nEdge;
}

int CRWCoast::nGetEndEdge(void) const
{
   return m_nEndEdge;
}


void CRWCoast::SetCoastlineExtCRS(CGeomLine const* pLCoast)
{
   m_LCoastlineExtCRS = *pLCoast;
   
   int nLen = m_LCoastlineExtCRS.nGetSize();
   
   m_VnProfileNumber = vector<int>(nLen, INT_NODATA);   
   m_VnPolygonNode = vector<int>(nLen, INT_NODATA);
   m_VnBreakingDistance = vector<int>(nLen, INT_NODATA);

   m_VdCurvatureDetailed = vector<double>(nLen, DBL_NODATA);
   m_VdCurvatureSmooth = vector<double>(nLen, DBL_NODATA);
   m_VdBreakingWaveHeight = vector<double>(nLen, DBL_NODATA);
   m_VdBreakingWaveAngle = vector<double>(nLen, DBL_NODATA);
   m_VdDepthOfBreaking = vector<double>(nLen, DBL_NODATA);
   m_VdFluxOrientation = vector<double>(nLen, DBL_NODATA);
   m_VdWaveEnergy = vector<double>(nLen, DBL_NODATA);  
}



void CRWCoast::AppendPointToCoastlineExtCRS(double const dX, double const dY)
{
   // Appends a coastline point (in external CRS), also appends dummy values for curvature, breaking wave height, wave angle, and flux orientation
   m_LCoastlineExtCRS.Append(dX, dY);

   m_VnProfileNumber.push_back(INT_NODATA);
   m_VnPolygonNode.push_back(INT_NODATA);
   m_VnBreakingDistance.push_back(INT_NODATA);

   m_VdCurvatureDetailed.push_back(DBL_NODATA);
   m_VdCurvatureSmooth.push_back(DBL_NODATA);
   m_VdBreakingWaveHeight.push_back(DBL_NODATA);
   m_VdBreakingWaveAngle.push_back(DBL_NODATA);
   m_VdDepthOfBreaking.push_back(DBL_NODATA);
   m_VdFluxOrientation.push_back(DBL_NODATA);
   m_VdWaveEnergy.push_back(DBL_NODATA);
}

CGeomLine* CRWCoast::pLGetCoastlineExtCRS(void)
{
   return &m_LCoastlineExtCRS;
}

CGeom2DPoint* CRWCoast::pPtGetCoastlinePointExtCRS(int const n)
{
   // Point is in external CRS NOTE no check to see that n is < m_LCoastlineExtCRS.Size()
   return &m_LCoastlineExtCRS[n];
}

int CRWCoast::nGetCoastlineSize(void) const
{
   return m_LCoastlineExtCRS.nGetSize();
}

// void CRWCoast::DisplayCoastline(void)
// {
//    m_LCoastlineExtCRS.Display();
// }


void CRWCoast::SetCoastlineGridCRS(CGeomILine const* pILCoastCells)
{
   m_ILCellsMarkedAsCoastline = *pILCoastCells;
}

// void CRWCoast::AppendCellMarkedAsCoastline(CGeom2DIPoint const* pPti)
// {
//    m_ILCellsMarkedAsCoastline.Append(*pPti);
// }
// 
// void CRWCoast::AppendCellMarkedAsCoastline(int const nX, int const nY)
// {
//    m_ILCellsMarkedAsCoastline.Append(CGeom2DIPoint(nX, nY));
// }

CGeom2DIPoint* CRWCoast::pPtiGetCellMarkedAsCoastline(int const n)
{
   // NOTE No check to see if n < size()
   return &m_ILCellsMarkedAsCoastline[n];
}

// int CRWCoast::nGetNCellsMarkedAsCoastline(void) const
// {
//    return m_ILCellsMarkedAsCoastline.size();
// }

// double CRWCoast::dGetCoastlineSegmentLength(int const m, int const n)
// {
//    // NOTE no check to see that m is < m_LCoastlineExtCRS.Size(), same for n
//    if (m == n)
//       return 0;
//
//    return hypot(m_LCoastlineExtCRS[n].dGetX() - m_LCoastlineExtCRS[m].dGetX(), m_LCoastlineExtCRS[n].dGetY() - m_LCoastlineExtCRS[m].dGetY());
// }

// double CRWCoast::dGetCoastlineLengthSoFar(int const n)
// {
//    // NOTE no check to see that n is < m_LCoastlineExtCRS.Size()
//    double dLen = 0;
//    for (int m = 0; m < n; m++)
//       dLen += dGetCoastlineSegmentLength(m, m+1);
//    return dLen;
// }

//! Returns the coastline number given a cell, or INT_NODATA if the cell is not 'under' a coastline 
int CRWCoast::nGetCoastPointGivenCell(CGeom2DIPoint const* pPtiCell)
{
   int nIndex = INT_NODATA;
   
   for (int n = 0; n < m_ILCellsMarkedAsCoastline.nGetSize(); n++)
   {
      if (m_ILCellsMarkedAsCoastline[n] == *pPtiCell)
         nIndex = n;
   }
   
   return nIndex;
}


double CRWCoast::dGetDetailedCurvature(int const nCoastPoint) const
{
   // NOTE no sanity check for nCoastPoint < m_VdCurvatureDetailed.Size()
   return m_VdCurvatureDetailed[nCoastPoint];
}

void CRWCoast::SetDetailedCurvature(int const nCoastPoint, double const dCurvature)
{
   // NOTE no check to see if nCoastPoint < m_VdCurvatureDetailed.size()
   m_VdCurvatureDetailed[nCoastPoint] = dCurvature;
}

vector<double>* CRWCoast::pVGetDetailedCurvature(void)
{
   return &m_VdCurvatureDetailed;
}

double CRWCoast::dGetSmoothCurvature(int const nCoastPoint) const
{
   // NOTE no sanity check for nCoastPoint < m_VdCurvatureSmooth.Size()
   return m_VdCurvatureSmooth[nCoastPoint];
}

void CRWCoast::SetSmoothCurvature(int const nCoastPoint, double const dCurvature)
{
   // NOTE no check to see if nCoastPoint < m_VdCurvatureSmooth.size()
   m_VdCurvatureSmooth[nCoastPoint] = dCurvature;
}

vector<double>* CRWCoast::pVGetSmoothCurvature(void)
{
   return &m_VdCurvatureSmooth;
}

void CRWCoast::SetDetailedCurvatureMean(double const dMean)
{
   m_dCurvatureDetailedMean = dMean;
}

double CRWCoast::dGetDetailedCurvatureMean(void) const
{
   return m_dCurvatureDetailedMean;
}

void CRWCoast::SetDetailedCurvatureSTD(double const dSTD)
{
   m_dCurvatureDetailedSTD = dSTD;
}

double CRWCoast::dGetDetailedCurvatureSTD(void) const
{
   return m_dCurvatureDetailedSTD;
}

void CRWCoast::SetSmoothCurvatureMean(double const dMean)
{
   m_dCurvatureSmoothMean = dMean;
}

double CRWCoast::dGetSmoothCurvatureMean(void) const
{
   return m_dCurvatureSmoothMean;
}


void CRWCoast::SetSmoothCurvatureSTD(double const dSTD)
{
   m_dCurvatureSmoothSTD = dSTD;
}

double CRWCoast::dGetSmoothCurvatureSTD(void) const
{
   return m_dCurvatureSmoothSTD;
}


CGeomProfile* CRWCoast::pGetProfile(int const nProfile)
{
   // NOTE No safety check that nProfile < m_VProfile.size()
   return &m_VProfile[nProfile];
}

void CRWCoast::AppendProfile(int const nCoastPoint, int const nProfile)
{
   CGeomProfile Profile(nCoastPoint);
   m_VProfile.push_back(Profile);

   m_VnProfileNumber[nCoastPoint] = nProfile;
}

// void CRWCoast::ReplaceProfile(int const nProfile, vector<CGeom2DPoint> const* pPtVProfileNew)
// {
//    // NOTE No safety check that nProfile < m_VProfile.size()
//    m_VProfile[nProfile].SetAllPointsInProfile(pPtVProfileNew);
// }

int CRWCoast::nGetNumProfiles(void) const
{
   return m_VProfile.size();
}

bool CRWCoast::bIsNormalProfileStartPoint(int const nCoastPoint) const
{
   // NOTE no sanity check for nCoastPoint < m_VnProfileNumber.Size()
   if (m_VnProfileNumber[nCoastPoint] != INT_NODATA)
      return true;

   return false;
}

int CRWCoast::nGetProfileNumber(int const nCoastPoint) const
{
   // Returns INT_NODATA if no profile at this point
   return m_VnProfileNumber[nCoastPoint];
}


void CRWCoast::CreateAlongCoastlineProfileIndex(void)
{
   // Creates an index containing the numbers of the coastline-normal profiles in along-coast sequence
   for (int nCoastPoint = 0; nCoastPoint < m_LCoastlineExtCRS.nGetSize(); nCoastPoint++)
      if (m_VnProfileNumber[nCoastPoint] != INT_NODATA)
         m_VnProfileCoastIndex.push_back(m_VnProfileNumber[nCoastPoint]);
}

int CRWCoast::nGetProfileAtAlongCoastlinePosition(int const n) const
{
   // Returns the number of the coastline-normal profile which is at position n in the along-coast sequence
   return m_VnProfileCoastIndex[n];
}

// int CRWCoast::nGetAlongCoastlineIndexOfProfile(int const nProfile)
// {
//    // Returns the along-coastline position of a coastline-normal profile
//    for (unsigned int n = 0; n < m_VnProfileCoastIndex.size(); n++)
//       if (m_VnProfileCoastIndex[n] == nProfile)
//          return n;
//    return -1;
// }


void CRWCoast::SetBreakingWaveHeight(int const nCoastPoint, double const dHeight)
{
   // NOTE no check to see if nCoastPoint < m_VdBreakingWaveHeight.size()
   m_VdBreakingWaveHeight[nCoastPoint] = dHeight;
}

double CRWCoast::dGetBreakingWaveHeight(int const nCoastPoint) const
{
   // NOTE no check to see if nCoastPoint < m_VdBreakingWaveHeight.size()
   return m_VdBreakingWaveHeight[nCoastPoint];
}

void CRWCoast::SetBreakingWaveOrientation(int const nCoastPoint, double const dOrientation)
{
   // NOTE no check to see if nCoastPoint < m_VdBreakingWaveAngle.size()
   m_VdBreakingWaveAngle[nCoastPoint] = dOrientation;
}

double CRWCoast::dGetBreakingWaveOrientation(int const nCoastPoint) const
{
   // NOTE no check to see if nCoastPoint < m_VdBreakingWaveAngle.size()
   return m_VdBreakingWaveAngle[nCoastPoint];
}

void CRWCoast::SetDepthOfBreaking(int const nCoastPoint, double const dDepth)
{
   // NOTE no check to see if nCoastPoint < m_VdDepthOfBreaking.size()
   m_VdDepthOfBreaking[nCoastPoint] = dDepth;
}

double CRWCoast::dGetDepthOfBreaking(int const nCoastPoint) const
{
   // NOTE no check to see if nCoastPoint < m_VdDepthOfBreaking.size()
   return m_VdDepthOfBreaking[nCoastPoint];
}

void CRWCoast::SetBreakingDistance(int const nCoastPoint, int const nDist)
{
   // NOTE no check to see if nCoastPoint < m_VnBreakingDistance.size()
   m_VnBreakingDistance[nCoastPoint] = nDist;
}

int CRWCoast::nGetBreakingDistance(int const nCoastPoint) const
{
   // NOTE no check to see if nCoastPoint < m_VnBreakingDistance.size()
   return m_VnBreakingDistance[nCoastPoint];
}

void CRWCoast::SetFluxOrientation(int const nCoastPoint, double const dOrientation)
{
   // NOTE no check to see if nCoastPoint < m_VdFluxOrientation.size()
   m_VdFluxOrientation[nCoastPoint] = dOrientation;
}

double CRWCoast::dGetFluxOrientation(int const nCoastPoint) const
{
   // NOTE no check to see if nCoastPoint < m_VdFluxOrientation.size()
   return m_VdFluxOrientation[nCoastPoint];
}

void CRWCoast::SetWaveEnergy(int const nCoastPoint, double const dEnergy)
{
   // NOTE no check to see if nCoastPoint < m_VdWaveEnergy.size()
//    assert(bIsFinite(dEnergy));
   m_VdWaveEnergy[nCoastPoint] = dEnergy;
}

double CRWCoast::dGetWaveEnergy(int const nCoastPoint) const
{
   // NOTE no check to see if nCoastPoint < m_VdWaveEnergy.size()
//    assert(bIsFinite(m_VdWaveEnergy[nCoastPoint]));

   return m_VdWaveEnergy[nCoastPoint];
}


void CRWCoast::AppendCoastLandform(CACoastLandform* pCoastLandform)
{
   m_pVLandforms.push_back(pCoastLandform);
}

CACoastLandform* CRWCoast::pGetCoastLandform(int const nCoastPoint)
{
   // NOTE no check to see if nCoastPoint < m_ILCellsMarkedAsCoastline.size()
   return m_pVLandforms[nCoastPoint];
}


void CRWCoast::SetPolygonNode(int const nPoint, int const nNode)
{
   // NOTE no check to see if nPoint < m_VnPolygonNode.size()
   m_VnPolygonNode[nPoint] = nNode;
}

int CRWCoast::nGetPolygonNode(int const nPoint) const
{
   // NOTE no check to see if nPoint < m_VnPolygonNode.size()
   return m_VnPolygonNode[nPoint];
}

void CRWCoast::CreatePolygon(int const nGlobalID, int const nCoastID, int const nCoastPoint, CGeom2DIPoint const* PtiNode, CGeom2DIPoint const* PtiAntiNode, int const nProfileUpCoast, int const nProfileDownCoast, vector<CGeom2DPoint> const* pVIn, int const nPointsUpCoastProfile, int const nPointsDownCoastProfile, int const nPointInPolygonStartPoint)
{
   CGeomCoastPolygon* pPolygon = new CGeomCoastPolygon(nGlobalID, nCoastID, nCoastPoint, nProfileUpCoast, nProfileDownCoast, pVIn, nPointsUpCoastProfile, nPointsDownCoastProfile, PtiNode, PtiAntiNode, nPointInPolygonStartPoint);

   m_pVPolygon.push_back(pPolygon);
}

int CRWCoast::nGetNumPolygons(void) const
{
   return m_pVPolygon.size();
}

CGeomCoastPolygon* CRWCoast::pGetPolygon(int const nPoly) const
{
   // NOTE no check to see if nPoint < m_VnPolygonNode.size()
   return m_pVPolygon[nPoly];
}


void CRWCoast::AppendPolygonLength(const double dLength)
{
   m_VdPolygonLength.push_back(dLength);
}

double CRWCoast::dGetPolygonLength(int const nIndex) const
{
   // NOTE no check to see if nIndex < m_VdPolygonLength.size()
   return m_VdPolygonLength[nIndex];
}


int CRWCoast::nGetNumShadowZoneBoundaries(void)
{
   return m_LShadowZoneBoundary.size();
}

void CRWCoast::AppendShadowZoneBoundary(const CGeomLine LBoundary)
{
   m_LShadowZoneBoundary.push_back(LBoundary);
}

CGeomLine* CRWCoast::pGetShadowZoneBoundary(int const n)
{
   // NOTE: no check to see if n < m_LShadowZoneBoundary.size()
   return &m_LShadowZoneBoundary[n];
}
