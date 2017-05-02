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


void CRWCoast::AppendToCoastline(double const dX, double const dY)
{
   // Appends a coastline point (in external CRS), also appends dummy values for curvature, breaking wave height, wave angle, and flux orientation
   m_LCoastline.Append(dX, dY);

   m_nVProfileNumber.push_back(INT_NODATA);
   m_nVPolygonNode.push_back(INT_NODATA);

   m_dVCurvatureDetailed.push_back(DBL_NODATA);
   m_dVCurvatureSmooth.push_back(DBL_NODATA);
   m_dVBreakingWaveHeight.push_back(DBL_NODATA);
   m_dVBreakingWaveAngle.push_back(DBL_NODATA);
   m_dVDepthOfBreaking.push_back(DBL_NODATA);
   m_dVFluxOrientation.push_back(DBL_NODATA);
   m_dVWaveEnergy.push_back(DBL_NODATA);
   m_nVBreakingDistance.push_back(INT_NODATA);
}

CGeomLine* CRWCoast::pLGetCoastline(void)
{
   return &m_LCoastline;
}

CGeom2DPoint* CRWCoast::pPtGetVectorCoastlinePoint(int const n)
{
   // Point is in external CRS NOTE no check to see that n is < m_LCoastline.Size()
   return &m_LCoastline[n];
}

int CRWCoast::nGetCoastlineSize(void) const
{
   return m_LCoastline.nGetSize();
}

// void CRWCoast::DisplayCoastline(void)
// {
//    m_LCoastline.Display();
// }

void CRWCoast::AppendCellMarkedAsCoastline(CGeom2DIPoint* Pti)
{
   m_VCellsMarkedAsCoastline.push_back(*Pti);
}

void CRWCoast::AppendCellMarkedAsCoastline(int const nX, int const nY)
{
   m_VCellsMarkedAsCoastline.push_back(CGeom2DIPoint(nX, nY));
}

// void CRWCoast::SetCellsMarkedAsCoastline(vector<CGeom2DIPoint>* VNewPoints)
// {
//    m_VCellsMarkedAsCoastline = *VNewPoints;
// }

CGeom2DIPoint* CRWCoast::pPtiGetCellMarkedAsCoastline(int const n)
{
   // NOTE No check to see if n < size()
   return &m_VCellsMarkedAsCoastline[n];
}

// int CRWCoast::nGetNCellsMarkedAsCoastline(void) const
// {
//    return m_VCellsMarkedAsCoastline.size();
// }

// double CRWCoast::dGetCoastlineSegmentLength(int const m, int const n)
// {
//    // NOTE no check to see that m is < m_LCoastline.Size(), same for n
//    if (m == n)
//       return 0;
//
//    return hypot(m_LCoastline[n].dGetX() - m_LCoastline[m].dGetX(), m_LCoastline[n].dGetY() - m_LCoastline[m].dGetY());
// }

// double CRWCoast::dGetCoastlineLengthSoFar(int const n)
// {
//    // NOTE no check to see that n is < m_LCoastline.Size()
//    double dLen = 0;
//    for (int m = 0; m < n; m++)
//       dLen += dGetCoastlineSegmentLength(m, m+1);
//    return dLen;
// }

//! Returns the coastline number given a cell, or INT_NODATA if the cell is not 'under' a coastline 
int CRWCoast::nGetCoastPointGivenCell(CGeom2DIPoint const* pPtiCell)
{
   int nIndex = INT_NODATA;
   
   for (unsigned int n = 0; n < m_VCellsMarkedAsCoastline.size(); n++)
   {
      if (m_VCellsMarkedAsCoastline[n] == *pPtiCell)
         nIndex = n;
   }
   
   return nIndex;
}


double CRWCoast::dGetDetailedCurvature(int const nCoastPoint) const
{
   // NOTE no sanity check for nCoastPoint < m_dVCurvatureDetailed.Size()
   return m_dVCurvatureDetailed[nCoastPoint];
}

void CRWCoast::SetDetailedCurvature(int const nCoastPoint, double const dCurvature)
{
   // NOTE no check to see if nCoastPoint < m_dVCurvatureDetailed.size()
   m_dVCurvatureDetailed[nCoastPoint] = dCurvature;
}

vector<double>* CRWCoast::pVGetDetailedCurvature(void)
{
   return &m_dVCurvatureDetailed;
}

double CRWCoast::dGetSmoothCurvature(int const nCoastPoint) const
{
   // NOTE no sanity check for nCoastPoint < m_dVCurvatureSmooth.Size()
   return m_dVCurvatureSmooth[nCoastPoint];
}

void CRWCoast::SetSmoothCurvature(int const nCoastPoint, double const dCurvature)
{
   // NOTE no check to see if nCoastPoint < m_dVCurvatureSmooth.size()
   m_dVCurvatureSmooth[nCoastPoint] = dCurvature;
}

vector<double>* CRWCoast::pVGetSmoothCurvature(void)
{
   return &m_dVCurvatureSmooth;
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

   m_nVProfileNumber[nCoastPoint] = nProfile;
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
   // NOTE no sanity check for nCoastPoint < m_nVProfileNumber.Size()
   if (m_nVProfileNumber[nCoastPoint] != INT_NODATA)
      return true;

   return false;
}

int CRWCoast::nGetProfileNumber(int const nCoastPoint) const
{
   // Returns INT_NODATA if no profile at this point
   return m_nVProfileNumber[nCoastPoint];
}


void CRWCoast::CreateAlongCoastlineProfileIndex(void)
{
   // Creates an index containing the numbers of the coastline-normal profiles in along-coast sequence
   for (int nCoastPoint = 0; nCoastPoint < m_LCoastline.nGetSize(); nCoastPoint++)
      if (m_nVProfileNumber[nCoastPoint] != INT_NODATA)
         m_nVProfileCoastIndex.push_back(m_nVProfileNumber[nCoastPoint]);
}

int CRWCoast::nGetProfileAtAlongCoastlinePosition(int const n) const
{
   // Returns the number of the coastline-normal profile which is at position n in the along-coast sequence
   return m_nVProfileCoastIndex[n];
}

// int CRWCoast::nGetAlongCoastlineIndexOfProfile(int const nProfile)
// {
//    // Returns the along-coastline position of a coastline-normal profile
//    for (unsigned int n = 0; n < m_nVProfileCoastIndex.size(); n++)
//       if (m_nVProfileCoastIndex[n] == nProfile)
//          return n;
//    return -1;
// }


void CRWCoast::SetBreakingWaveHeight(int const nCoastPoint, double const dHeight)
{
   // NOTE no check to see if nCoastPoint < m_dVBreakingWaveHeight.size()
   m_dVBreakingWaveHeight[nCoastPoint] = dHeight;
}

double CRWCoast::dGetBreakingWaveHeight(int const nCoastPoint) const
{
   // NOTE no check to see if nCoastPoint < m_dVBreakingWaveHeight.size()
   return m_dVBreakingWaveHeight[nCoastPoint];
}

void CRWCoast::SetBreakingWaveOrientation(int const nCoastPoint, double const dOrientation)
{
   // NOTE no check to see if nCoastPoint < m_dVBreakingWaveAngle.size()
   m_dVBreakingWaveAngle[nCoastPoint] = dOrientation;
}

double CRWCoast::dGetBreakingWaveOrientation(int const nCoastPoint) const
{
   // NOTE no check to see if nCoastPoint < m_dVBreakingWaveAngle.size()
   return m_dVBreakingWaveAngle[nCoastPoint];
}

void CRWCoast::SetDepthOfBreaking(int const nCoastPoint, double const dDepth)
{
   // NOTE no check to see if nCoastPoint < m_dVDepthOfBreaking.size()
   m_dVDepthOfBreaking[nCoastPoint] = dDepth;
}

double CRWCoast::dGetDepthOfBreaking(int const nCoastPoint) const
{
   // NOTE no check to see if nCoastPoint < m_dVDepthOfBreaking.size()
   return m_dVDepthOfBreaking[nCoastPoint];
}

void CRWCoast::SetBreakingDistance(int const nCoastPoint, int const nDist)
{
   // NOTE no check to see if nCoastPoint < m_nVBreakingDistance.size()
   m_nVBreakingDistance[nCoastPoint] = nDist;
}

int CRWCoast::nGetBreakingDistance(int const nCoastPoint) const
{
   // NOTE no check to see if nCoastPoint < m_nVBreakingDistance.size()
   return m_nVBreakingDistance[nCoastPoint];
}

void CRWCoast::SetFluxOrientation(int const nCoastPoint, double const dOrientation)
{
   // NOTE no check to see if nCoastPoint < m_dVFluxOrientation.size()
   m_dVFluxOrientation[nCoastPoint] = dOrientation;
}

double CRWCoast::dGetFluxOrientation(int const nCoastPoint) const
{
   // NOTE no check to see if nCoastPoint < m_dVFluxOrientation.size()
   return m_dVFluxOrientation[nCoastPoint];
}

void CRWCoast::SetWaveEnergy(int const nCoastPoint, double const dEnergy)
{
   // NOTE no check to see if nCoastPoint < m_dVWaveEnergy.size()
//    assert(bIsFinite(dEnergy));
   m_dVWaveEnergy[nCoastPoint] = dEnergy;
}

double CRWCoast::dGetWaveEnergy(int const nCoastPoint) const
{
   // NOTE no check to see if nCoastPoint < m_dVWaveEnergy.size()
//    assert(bIsFinite(m_dVWaveEnergy[nCoastPoint]));

   return m_dVWaveEnergy[nCoastPoint];
}


void CRWCoast::AppendCoastLandform(CACoastLandform* pCoastLandform)
{
   m_pVLandforms.push_back(pCoastLandform);
}

CACoastLandform* CRWCoast::pGetCoastLandform(int const nCoastPoint)
{
   // NOTE no check to see if nCoastPoint < m_VCellsMarkedAsCoastline.size()
   return m_pVLandforms[nCoastPoint];
}


void CRWCoast::SetPolygonNode(int const nPoint, int const nNode)
{
   // NOTE no check to see if nPoint < m_nVPolygonNode.size()
   m_nVPolygonNode[nPoint] = nNode;
}

int CRWCoast::nGetPolygonNode(int const nPoint) const
{
   // NOTE no check to see if nPoint < m_nVPolygonNode.size()
   return m_nVPolygonNode[nPoint];
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
   // NOTE no check to see if nPoint < m_nVPolygonNode.size()
   return m_pVPolygon[nPoly];
}


void CRWCoast::AppendPolygonLength(const double dLength)
{
   m_dVPolygonLength.push_back(dLength);
}

double CRWCoast::dGetPolygonLength(int const nIndex) const
{
   // NOTE no check to see if nIndex < m_dVPolygonLength.size()
   return m_dVPolygonLength[nIndex];
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
