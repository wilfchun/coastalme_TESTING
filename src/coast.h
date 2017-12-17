/*!
 * \class CRWCoast
 * \brief Real-world class used to represent coastline objects
 * \details TODO This is a more detailed description of the CRWCoast class.
 * \author David Favis-Mortlock
 * \author Andres Payo

 * \date 2017
 * \copyright GNU General Public License
 *
 * \file coast.h
 * \brief Contains CRWCoast definitions
 *
 */

#ifndef COAST_H
#define COAST_H
/*===============================================================================================================================

 This file is part of CoastalME, the Coastal Modelling Environment.

 CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

===============================================================================================================================*/
#include "cme.h"
#include "profile.h"
#include "cell.h"
#include "coast_landform.h"
#include "coast_polygon.h"


class CGeomProfile;
class CACoastLandform;
class CGeomCoastPolygon;

class CRWCoast
{
private:
   int
      m_nSeaHandedness,          // Direction of the sea from the coastline, travelling down-coast (i.e. in direction of increasing coast point indices)
      m_nStartEdge,
      m_nEndEdge;

   double
      m_dCurvatureDetailedMean,
      m_dCurvatureDetailedSTD,
      m_dCurvatureSmoothMean,
      m_dCurvatureSmoothSTD;
      
   CGeomLine
      m_LCoastlineExtCRS;           // Smoothed line of points (external CRS) giving the plan view of the vector coast

   // The following have the same length as m_LCoastlineExtCRS (which may be different each timestep)
   CGeomILine
      m_ILCellsMarkedAsCoastline;   // Unsmoothed integer x-y co-ords (grid CRS) of the cell marked as coastline at each point on the vector coastline.
                                    // Note that this is the same as point zero in profile coords
   vector<int>
      m_VnProfileNumber,            // At each point on m_LCoastlineExtCRS: INT_NODATA if no profile there, otherwise the profile number
      m_VnBreakingDistance,         // Distance of breaking (in cells), at each point on m_LCoastlineExtCRS
      m_VnPolygonNode;              // At every point on m_LCoastlineExtCRS: INT_NODATA if no nodepoint there, otherwise the node (point of greatest concave curvature) number for a coast polygon
   vector<double>
      m_VdCurvatureDetailed,        // Detailed curvature at each point on m_LCoastlineExtCRS
      m_VdCurvatureSmooth,          // Smoothed curvature at each point on m_LCoastlineExtCRS
      m_VdDeepWaterWaveHeight,      // The deep water wave height at the end of a normal drawn from each point on m_LCoastlineExtCRS
      m_VdDeepWaterWaveOrientation, // The deep water wave orientation at the end of a normal drawn from each point on m_LCoastlineExtCRS
      m_VdBreakingWaveHeight,       // The breaking wave height on a normal drawn from each point on m_LCoastlineExtCRS
      m_VdBreakingWaveOrientation,  // The breaking wave orientation on a normal drawn from each point on m_LCoastlineExtCRS
      m_VdDepthOfBreaking,          // The depth of breaking on a normal drawn from each point on m_LCoastlineExtCRS
      m_VdFluxOrientation,          // As in the COVE model, is the orientation alongshore energy/sediment movement; a +ve flux is in direction of increasing indices along coast. At each point on m_LCoastlineExtCRS
      m_VdWaveEnergy;               // Wave energy at each point on m_LCoastlineExtCRS
      
   vector<CACoastLandform*>
      m_pVLandforms;                // Pointer to a coastal landform object, at each point on m_LCoastlineExtCRS

   // These do not have the same length as m_LCoastlineExtCRS
   vector<CGeomProfile>
      m_VProfile;                   // Coast profile objects, in the sequence in which they were created (concave coastline curvature)
   vector<int>
      m_VnProfileCoastIndex;        // Indices of coast profiles sorted into along-coastline sequence, size = number of profiles
   vector<CGeomCoastPolygon*>
      m_pVPolygon;                  // Coast polygons, size = number of polygons
   vector<double>
      m_VdPolygonLength;            // Lengths of coast polygons, size = number of polygons
   vector<CGeomLine>
      m_LShadowBoundary;            // Lines which delineate the edge of a shadow zone, ext CRS
   vector<CGeomLine>
      m_LDowndriftBoundary;         // Lines which delineate the edge of a downdrift zone, ext CRS
      
public:
   CRWCoast(void);
   ~CRWCoast(void);

   void SetSeaHandedness(int const);
   int nGetSeaHandedness(void) const;

   void SetStartEdge(int const);
   int nGetStartEdge(void) const;

   void SetEndEdge(int const);
   int nGetEndEdge(void) const;

   void SetCoastlineExtCRS(CGeomLine const*);
   void AppendPointToCoastlineExtCRS(double const, double const);
   CGeomLine* pLGetCoastlineExtCRS(void);
   CGeom2DPoint* pPtGetCoastlinePointExtCRS(int const);
   
   int nGetCoastlineSize(void) const;
//    double dGetCoastlineSegmentLength(int const, int const);
//    double dGetCoastlineLengthSoFar(int const);
//    void DisplayCoastline(void);
   
   void SetCoastlineGridCRS(CGeomILine const*);
//    void AppendCellMarkedAsCoastline(CGeom2DIPoint const*);
//    void AppendCellMarkedAsCoastline(int const, int const);
   CGeom2DIPoint* pPtiGetCellMarkedAsCoastline(int const);
//    int nGetNCellsMarkedAsCoastline(void) const;
   int nGetCoastPointGivenCell(CGeom2DIPoint*);

   double dGetDetailedCurvature(int const) const;
   void SetDetailedCurvature(int const, double const);
   vector<double>* pVGetDetailedCurvature(void);
   double dGetSmoothCurvature(int const) const;
   void SetSmoothCurvature(int const, double const);
   vector<double>* pVGetSmoothCurvature(void);
   void SetDetailedCurvatureMean(double const);
   double dGetDetailedCurvatureMean(void) const;
   void SetDetailedCurvatureSTD(double const);
   double dGetDetailedCurvatureSTD(void) const;
   void SetSmoothCurvatureMean(double const);
   double dGetSmoothCurvatureMean(void) const;
   void SetSmoothCurvatureSTD(double const);
   double dGetSmoothCurvatureSTD(void) const;
   
   CGeomProfile* pGetProfile(int const);
   void AppendProfile(int const, int const);
//    void ReplaceProfile(int const, vector<CGeom2DPoint> const*);
   int nGetNumProfiles(void) const;
   bool bIsProfileStartPoint(int const) const;
   int nGetProfileNumber(int const) const;

   void CreateAlongCoastlineProfileIndex(void);
   int nGetProfileAtAlongCoastlinePosition(int const) const;
   int nGetDownCoastProfileNumber(int const nProfile) const;   
//    int nGetAlongCoastlineIndexOfProfile(int const);

   void SetDeepWaterWaveHeight(int const, double const);
   double dGetDeepWaterWaveHeight(int const) const;
   
   void SetDeepWaterWaveOrientation(int const, double const);
   double dGetDeepWaterWaveOrientation(int const) const;

   void SetBreakingWaveHeight(int const, double const);
   double dGetBreakingWaveHeight(int const) const;
   
   void SetBreakingWaveOrientation(int const, double const);
   double dGetBreakingWaveOrientation(int const) const;
   
   void SetDepthOfBreaking(int const, double const);
   double dGetDepthOfBreaking(int const) const;
   
   void SetBreakingDistance(int const, int const);
   int nGetBreakingDistance(int const) const;
   
   void SetFluxOrientation(int const, double const);
   double dGetFluxOrientation(int const) const;
   
   void SetWaveEnergy(int const, double const);
   double dGetWaveEnergy(int const) const;

   void AppendCoastLandform(CACoastLandform*);
   CACoastLandform* pGetCoastLandform(int const);

   void SetPolygonNode(int const, int const);
   int nGetPolygonNode(int const) const;
   void CreatePolygon(int const, int const, int const, CGeom2DIPoint const*, CGeom2DIPoint const*, int const, int const, vector<CGeom2DPoint> const*, int const, int const, int const);
   int nGetNumPolygons(void) const;
   CGeomCoastPolygon* pGetPolygon(int const) const;

   void AppendPolygonLength(const double);
   double dGetPolygonLength(int const) const;
   
   int nGetNumShadowBoundaries(void);
   void AppendShadowBoundary(const CGeomLine);
   CGeomLine* pGetShadowBoundary(int const);
   
   int nGetNumDowndriftBoundaries(void);
   void AppendDowndriftBoundary(const CGeomLine);
   CGeomLine* pGetDowndriftBoundary(int const);
};
#endif //COAST_H

