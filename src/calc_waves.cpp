/*!
 *
 * \file calc_waves.cpp
 * \brief Simulates wave propagation
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
// #include <assert.h>
#include <cmath>

#include <string>
using std::stoi;

#include <iostream>
using std::ifstream;
using std::cout;
using std::endl;
using std::ios;

#include <iomanip>
using std::setprecision;
using std::setiosflags;
using std::resetiosflags;
using std::setw;

#include <algorithm>
using std::sort;
using std::remove;
using std::reverse;

#include <stack>
using std::stack;

#include "cme.h"
#include "coast.h"
#include "simulation.h"

#include "hermite_cubic.h"
#include "linearinterp.h"



/*===============================================================================================================================

 Give every coast point a value for deep water wave height and direction TODO this may not be realistic, may need to use end-of-profile value instead (how?)

===============================================================================================================================*/
int CSimulation::nSetAllCoastpointDeepWaterWaveValues(void)
{   
   // For each coastline, put a value for deep water wave height and direction at each coastline point
   int
      nNextProfile = -1,
      nDistFromPrevProfile = 0,
      nDistToNextProfile = 0;
      
   double
      dPrevProfileDeepWaterWaveHeight = 0,
      dPrevProfileDeepWaterWaveOrientation = 0,
      dNextProfileDeepWaterWaveHeight = 0,
      dNextProfileDeepWaterWaveOrientation = 0,
      dDist = 0;
      
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      for (int nPoint = 0; nPoint < m_VCoast[nCoast].nGetCoastlineSize(); nPoint++)
      {
         // We are going down-coast
         if (m_VCoast[nCoast].bIsProfileStartPoint(nPoint))
         {
            // OK, a coastline-normal profile begins at this coastline point, so set the deep water wave values at this coastline point to be the values at the seaward end of the coastline normal
            int nProfile = m_VCoast[nCoast].nGetProfileNumber(nPoint);
            CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfile(nProfile);
            
            double
               dThisDeepWaterWaveHeight = pProfile->dGetDeepWaterWaveHeight(),
               dThisDeepWaterWaveOrientation = pProfile->dGetDeepWaterWaveOrientation();
               
            m_VCoast[nCoast].SetDeepWaterWaveHeight(nPoint, dThisDeepWaterWaveHeight);
            m_VCoast[nCoast].SetDeepWaterWaveOrientation(nPoint, dThisDeepWaterWaveOrientation);
            
            // Reset for next time
            nDistFromPrevProfile = 0;
            dPrevProfileDeepWaterWaveHeight = dThisDeepWaterWaveHeight;
            dPrevProfileDeepWaterWaveOrientation = dThisDeepWaterWaveOrientation;            
            
            // Find the next profile
            nNextProfile = m_VCoast[nCoast].nGetDownCoastProfileNumber(nProfile);
            if (nNextProfile == INT_NODATA)
            {
               // We are at the end of the coast
               break;               
            }
            
            CGeomProfile* pNextProfile = m_VCoast[nCoast].pGetProfile(nNextProfile);

            // And the distance (in along-coast points) to the next profile
            nDistToNextProfile = pNextProfile->nGetNumCoastPoint() - nPoint;
            dDist = nDistToNextProfile;
            
            // And the next profile's deep water wave values
            dNextProfileDeepWaterWaveHeight = pNextProfile->dGetDeepWaterWaveHeight();
            dNextProfileDeepWaterWaveOrientation = pNextProfile->dGetDeepWaterWaveOrientation();            
                        
//             LogStream << m_ulIteration << ": coast point = " << nPoint << " IS PROFILE START, dThisDeepWaterWaveHeight = " << dThisDeepWaterWaveHeight << ", dThisDeepWaterWaveOrientation = " << dThisDeepWaterWaveOrientation << endl;
         }
         
         else
         {
            // This coast point is not the start of a coastline normal, so set the deep water wave values to a weighted average of those from the up-coast and down-coast profiles
            nDistFromPrevProfile++;
            nDistToNextProfile--;
            
            double
               dPrevWeight = (dDist - nDistFromPrevProfile) / dDist,
               dNextWeight = (dDist - nDistToNextProfile) / dDist,
               dThisDeepWaterWaveHeight = (dPrevWeight * dPrevProfileDeepWaterWaveHeight) + (dNextWeight * dNextProfileDeepWaterWaveHeight),
               dThisDeepWaterWaveOrientation = dKeepWithin360((dPrevWeight * dPrevProfileDeepWaterWaveOrientation) + (dNextWeight * dNextProfileDeepWaterWaveOrientation));
               
            m_VCoast[nCoast].SetDeepWaterWaveHeight(nPoint, dThisDeepWaterWaveHeight);
            m_VCoast[nCoast].SetDeepWaterWaveOrientation(nPoint, dThisDeepWaterWaveOrientation);
            
//             LogStream << m_ulIteration << ": coast point = " << nPoint << " dThisDeepWaterWaveHeight = " << dThisDeepWaterWaveHeight << " dThisDeepWaterWaveOrientation = " << dThisDeepWaterWaveOrientation << endl;
         }
      }
   }
   
   return RTN_OK;
}


/*===============================================================================================================================

 Simulates wave propagation along all coastline-normal profiles

===============================================================================================================================*/
int CSimulation::nDoAllPropagateWaves(void)
{
   // Set up all-profile vectors to hold the wave attribute data at every profile point on all profiles
   vector<bool>VbBreakingAll;

   vector<int>
      VnXAll,
      VnYAll;
      
   vector<double>
      VdHeightXAll,
      VdHeightYAll;
      
   // Calculate wave properties for every coast
   bool bSomeNonStartOrEndOfCoastProfiles = false;
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      int
         nCoastSize = m_VCoast[nCoast].nGetCoastlineSize(),
         nNumProfiles = m_VCoast[nCoast].nGetNumProfiles();

      // Calculate wave properties at every point along each valid profile, and for the cells under the profiles. Do this in the original (curvature-related) profile sequence
      for (int nProfile = 0; nProfile < nNumProfiles; nProfile++)
      {
         vector<bool> VbBreaking;
         vector<int>
            VnX,
            VnY;
         vector<double>
            VdHeightX,
            VdHeightY;
            
         int nRet = nCalcWavePropertiesOnProfile(nCoast, nCoastSize, nProfile, &VnX, &VnY, &VdHeightX, &VdHeightY, &VbBreaking);
         if (nRet != RTN_OK)
            return nRet;
         
         // Are the waves off-shore? If so, do nothing more with this profile. The wave values for cells have already been given the off-shore value
         if (VbBreaking.empty())
            continue;
         
         // Is this a start of coast or end of coast profile?
         if ((! m_VCoast[nCoast].pGetProfile(nProfile)->bStartOfCoast()) && (! m_VCoast[nCoast].pGetProfile(nProfile)->bEndOfCoast()))
         {
            // It is neither a start of coast or an end of coast profile, so set switch
            bSomeNonStartOrEndOfCoastProfiles = true;            
         }
         
         // TEST
//          for (int nn = 0; nn < VnX.size(); nn++)
//          {
//             LogStream << "nProfile = " << nProfile << " nn = " << nn << " VnX[nn] = " << VnX[nn] << " VnY[nn] = " << VnY[nn] << " VdHeightX[nn] = " << VdHeightX[nn] << " VdHeightY[nn] = " << VdHeightY[nn] << " VbBreaking[nn] = " << VbBreaking[nn] << endl;
//          }
//          LogStream << endl;
         // TEST
         
         // Append to the all-profile vectors
         VnXAll.insert(VnXAll.end(), VnX.begin(), VnX.end());
         VnYAll.insert(VnYAll.end(), VnY.begin(), VnY.end());
         VdHeightXAll.insert(VdHeightXAll.end(), VdHeightX.begin(), VdHeightX.end());
         VdHeightYAll.insert(VdHeightYAll.end(), VdHeightY.begin(), VdHeightY.end());
         VbBreakingAll.insert(VbBreakingAll.end(), VbBreaking.begin(), VbBreaking.end());
      }
   }
   
   // OK, do we have some profiles other than start of coast or end of coast profiles in the all-profile vectors? We need to check this, because GDALGridCreate() in nInterpolateWavePropertiesToWithinPolygonCells() does not work if we give it only a start of coast or an end of cell profile to work with
   if (! bSomeNonStartOrEndOfCoastProfiles)
   {
      LogStream << m_ulIteration << ": waves are on-shore only for start and/or end of coast profiles" << endl;
      
      return RTN_OK;      
   }

//    // DEBUG CODE ============================================
//    LogStream << "Out of loop" << endl;
//    for (int nn = 0; nn < VnXAll.size(); nn++)
//    {
//       LogStream << "nn = " << nn << " VnXAll[nn] = " << VnXAll[nn] << " VnYAll[nn] = " << VnYAll[nn] << " VdHeightXAll[nn] = " << VdHeightXAll[nn] << " VdHeightYAll[nn] = " << VdHeightYAll[nn] << " VbBreakingAll[nn] = " << VbBreakingAll[nn] << endl;
//    }
//    LogStream << endl;
//    // DEBUG CODE ============================================

   // Are the waves off-shore for every profile? If so, do nothing more
   if (VbBreakingAll.empty())
   {
      LogStream << m_ulIteration << ": waves off-shore for all profiles" << endl;
      
      return RTN_OK;      
   }

   // Some waves are on-shore, so interpolate the wave attributes from all profile points to all within-polygon sea cells
   int nRet = nInterpolateWavePropertiesToWithinPolygonCells(&VnXAll, &VnYAll, &VdHeightXAll, &VdHeightYAll);
   if (nRet != RTN_OK)
      return nRet;
 
   // Interpolate the wave attributes from all profile points to all sea cells that are in the active zone
   //nRet = nInterpolateWavePropertiesToActiveZoneCells(&VnXAll, &VnYAll, &VbBreakingAll);
   nRet = nInterpolateWavePropertiesToActiveZoneCells();
   if (nRet != RTN_OK)
      return nRet;

   // Find the shadow zones and then modify waves in and adjacent to them
   nRet = nDoAllShadowZones();
   if (nRet != RTN_OK)
      return nRet;

   // Fill in artefactual 'holes' in active zone and wave property patterns
   CalcD50AndFillWaveCalcHoles();

   // Modify the wave breaking properties (wave height, wave dir, breaking depth, breaking distance) for coastline points within the shadow zone
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      int nNumProfiles = m_VCoast[nCoast].nGetNumProfiles();

      for (int nProfile = 0; nProfile < nNumProfiles; nProfile++)
         ModifyBreakingWavePropertiesWithinShadowZoneToCoastline(nCoast, nProfile);
   }

   // Interpolate these wave properties for all remaining coastline points. Do this in along-coastline sequence
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      int
         nCoastSize = m_VCoast[nCoast].nGetCoastlineSize(),
         nNumProfiles = m_VCoast[nCoast].nGetNumProfiles();

      // Interpolate these wave properties for all remaining coastline points. Do this in along-coastline sequence, but do not do this for the end-of-coastline profile (which is the final one)
      for (int n = 0; n < nNumProfiles-1; n++)
         InterpolateWavePropertiesToCoastline(nCoast, n, nNumProfiles);

      // Calculate wave energy at every point on the coastline
      for (int nCoastPoint = 0; nCoastPoint < nCoastSize; nCoastPoint++)
      {
         // Equation 4 from Walkden & Hall, 2005
         double dBreakingWaveHeight = m_VCoast[nCoast].dGetBreakingWaveHeight(nCoastPoint);
         if (dBreakingWaveHeight != DBL_NODATA)
         {
            double dErosiveWaveForce = pow(dBreakingWaveHeight, WALKDEN_HALL_PARAM_1) * pow(m_dWavePeriod, WALKDEN_HALL_PARAM_2);

            // Calculate total wave energy at each coast point during this timestep
            double dWaveEnergy = dErosiveWaveForce * m_dTimeStep * 3600;
            
            m_VCoast[nCoast].SetWaveEnergyAtBreaking(nCoastPoint, dWaveEnergy);
         }
      }
   }

   return RTN_OK;
}


/*===============================================================================================================================

 Calculates the angle between the wave direction and a normal to the coastline tangent. If wave direction has a component which is down-coast (i.e. in the direction with increasing coast point numbers), then the angle returned is +ve. If wave direction has a component which is up-coast (i.e. in the direction with decreasing coast point numbers), then the angle returned is -ve. If waves are in an off-shore direction, DBL_NODATA is returned

===============================================================================================================================*/
double CSimulation::dCalcWaveAngleToCoastNormal(double const dCoastAngle, double const dWaveOrientation, int const nSeaHand)
{
   double dWaveToNormalAngle = 0;
   
   if (nSeaHand == LEFT_HANDED)
      // Left-handed coast
      dWaveToNormalAngle = fmod((dWaveOrientation - dCoastAngle + 360), 360) - 90;
   else
      // Right-handed coast
      dWaveToNormalAngle = fmod((dWaveOrientation - dCoastAngle + 360), 360) - 270;

   if ((dWaveToNormalAngle >= 90) || (dWaveToNormalAngle <= -90))
      dWaveToNormalAngle = DBL_NODATA;

   return dWaveToNormalAngle;
}


/*===============================================================================================================================
 
 Calculates wave properties along a coastline-normal profile using either the COVE linear wave theory approach or the external CShore model
 
===============================================================================================================================*/
int CSimulation::nCalcWavePropertiesOnProfile(int const nCoast, int const nCoastSize, int const nProfile, vector<int>* pVnX, vector<int>* pVnY, vector<double>* pVdHeightX, vector<double>* pVdHeightY, vector<bool>* pVbBreaking)
{
   CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfile(nProfile);
   
   // Calculate some wave properties based on the wave period following Airy wave theory
   m_dC_0 = (m_dG * m_dWavePeriod) / (2 * PI);           // Deep water (offshore) wave celerity (m/s)
   m_dL_0 = m_dC_0 * m_dWavePeriod;                      // Deep water (offshore) wave length (m)
   
   // Only do this for profiles without problems. Still do start- and end-of-coast profiles however
   if (! pProfile->bOKIncStartAndEndOfCoast())
   {
      LogStream << m_ulIteration << ": invalid profile " << nProfile << endl;
      
      return RTN_OK;
   }
   
   int
      nSeaHand = m_VCoast[nCoast].nGetSeaHandedness(),
      nCoastPoint = pProfile->nGetNumCoastPoint();
   
   // Get the flux orientation (the orientation of a line which is tangential to the coast) at adjacent coastline points. Note special treatment for the coastline end points
   double
      dFluxOrientationThis = m_VCoast[nCoast].dGetFluxOrientation(nCoastPoint),
      dFluxOrientationPrev = 0,
      dFluxOrientationNext = 0;
   
   if (nCoastPoint == 0)
   {
      dFluxOrientationPrev = dFluxOrientationThis;
      dFluxOrientationNext = m_VCoast[nCoast].dGetFluxOrientation(1);
   }
   else if (nCoastPoint == nCoastSize-1)
   {
      dFluxOrientationPrev = m_VCoast[nCoast].dGetFluxOrientation(nCoastPoint-2);
      dFluxOrientationNext = dFluxOrientationThis;
   }
   else
   {
      dFluxOrientationPrev = m_VCoast[nCoast].dGetFluxOrientation(nCoastPoint-1);
      dFluxOrientationNext = m_VCoast[nCoast].dGetFluxOrientation(nCoastPoint+1);
   }
   
   // Get the deep water wave orientation for this profile
   double dDeepWaterWaveOrientation = pProfile->dGetDeepWaterWaveOrientation();
   
   // Calculate the angle between the deep water wave direction and a normal to the coast tangent
   double dWaveToNormalAngle = dCalcWaveAngleToCoastNormal(dFluxOrientationThis, dDeepWaterWaveOrientation, nSeaHand);
   
   // Are the waves off-shore?
   if (dWaveToNormalAngle == DBL_NODATA)
   {
      // They are so, do nothing (each cell under the profile has already been initialised with deep water wave height and wave direction)
//       LogStream << m_ulIteration << ": profile " << nProfile << " has sea to " << (m_VCoast[nCoast].nGetSeaHandedness() == RIGHT_HANDED ? "right" : "left") << " dWaveToNormalAngle = " << dWaveToNormalAngle << " which is off-shore" << endl;
      
      return RTN_OK;
   }
   
//    LogStream << m_ulIteration << ": profile = " << nProfile << " has sea to " << (m_VCoast[nCoast].nGetSeaHandedness() == RIGHT_HANDED ? "right" : "left") << " dWaveToNormalAngle = " << dWaveToNormalAngle << " which is " << (dWaveToNormalAngle < 0 ? "DOWN" : "UP") << "-coast" << endl;
   
   // Calculate the angle between the deep water wave direction and a normal to the coast tangent for the previous coast point
   double dWaveToNormalAnglePrev;
   if (nCoastPoint > 0)
   {
      // Get the deep water wave orientation for the up-coast point
      double dPrevDeepWaterWaveOrientation = m_VCoast[nCoast].dGetDeepWaterWaveOrientation(nCoastPoint-1);
      
      dWaveToNormalAnglePrev = dCalcWaveAngleToCoastNormal(dFluxOrientationPrev, dPrevDeepWaterWaveOrientation, nSeaHand);
   }
   else
      dWaveToNormalAnglePrev = dWaveToNormalAngle;
   
//    if (dWaveToNormalAnglePrev == DBL_NODATA)
//       LogStream << "\tPrevious profile, dWaveToNormalAnglePrev = " << dWaveToNormalAnglePrev << " which is off-shore" << endl;
//    else
//       LogStream << "\tPrevious profile, dWaveToNormalAnglePrev = " << dWaveToNormalAnglePrev << " which is " << (dWaveToNormalAnglePrev < 0 ? "DOWN" : "UP") << "-coast" << endl;
   
   // Calculate the angle between the deep water wave direction and a normal to the coast tangent for the next coast point
   double dWaveToNormalAngleNext;
   if (nCoastPoint < nCoastSize-1)
   {
      // Get the deep water wave orientation for the down-coast point
      double dNextDeepWaterWaveOrientation = m_VCoast[nCoast].dGetDeepWaterWaveOrientation(nCoastPoint+1);
      
      dWaveToNormalAngleNext = dCalcWaveAngleToCoastNormal(dFluxOrientationNext, dNextDeepWaterWaveOrientation, nSeaHand);
   }
   else
      dWaveToNormalAngleNext = dWaveToNormalAngle;
   
//    if (dWaveToNormalAngleNext == DBL_NODATA)
//       LogStream << "\tNext profile, dWaveToNormalAngleNext = " << dWaveToNormalAngleNext << " which is off-shore" << endl;
//    else
//       LogStream << "\tNext profile, dWaveToNormalAngleNext = " << dWaveToNormalAngleNext << " which is " << (dWaveToNormalAngleNext < 0 ? "DOWN" : "UP") << "-coast" << endl;
   
   // Following Ashton and Murray (2006), if we have high-angle waves then use the flux orientation of the previous (up-coast) profile, if transitioning from diffusive to antidiffusive use flux maximizing angle (45 degrees)
   if ((dWaveToNormalAngle > 0) && (dWaveToNormalAnglePrev != DBL_NODATA) && (dWaveToNormalAnglePrev > 0))
   {
      if (dWaveToNormalAngle > 45)
      {
         if (dWaveToNormalAnglePrev < 45)
         {
            dWaveToNormalAngle = 45;
            //            LogStream << "\tA1" << endl;
         }
         else
         {
            dWaveToNormalAngle = dWaveToNormalAnglePrev;
            //            LogStream << "\tA2" << endl;
         }
      }
   }
   else if ((dWaveToNormalAngle < 0) && (dWaveToNormalAngleNext != DBL_NODATA) && (dWaveToNormalAngleNext < 0))
   {
      if (dWaveToNormalAngle < -45)
      {
         if (dWaveToNormalAngleNext > -45)
         {
            dWaveToNormalAngle = -45;
            //            LogStream << "\tB1" << endl;
         }
         else
         {
            dWaveToNormalAngle = dWaveToNormalAngleNext;
            //            LogStream << "\tB2" << endl;
         }
      }
   }
   else if ((dWaveToNormalAngle > 45) && (dWaveToNormalAnglePrev != DBL_NODATA) && (dWaveToNormalAnglePrev > 0))
   {
      // The wave direction here has an up-coast (decreasing indices) component: so for high-angle waves use the orientation from the up-coast (previous) profile
      //      LogStream << "\tCCC" << endl;
      
      dWaveToNormalAngle = dFluxOrientationPrev;
   }
   else if ((dWaveToNormalAngle < -45)  && (dWaveToNormalAngleNext != DBL_NODATA) && (dWaveToNormalAngleNext < 0))
   {
      // The wave direction here has a down-coast (increasing indices) component: so for high-angle waves use the orientation from the down-coast (next) profile
      //      LogStream << "\tDDD" << endl;
      
      dWaveToNormalAngle = dFluxOrientationNext;
   }
   
   // Initialize the wave properties at breaking for this profile
   bool bBreaking = false;
   int
      nProfileSize = pProfile->nGetNumCellsInProfile(),
      nProfileBreakingDist = 0;
   double
      dProfileBreakingWaveHeight = 0,
      dProfileBreakingWaveOrientation = 0,
      dProfileBreakingDepth = 0,
      dProfileWaveHeight = DBL_NODATA,
      dProfileDeepWaterWaveHeight = pProfile->dGetDeepWaterWaveHeight(),
      dProfileWaveOrientation = DBL_NODATA,
      dProfileDeepWaterWaveOrientation = pProfile->dGetDeepWaterWaveOrientation();
   vector<bool>
      VbWaveIsBreaking(nProfileSize,0);
   vector<double>
      VdWaveHeight(nProfileSize,0),
      VdWaveDirection(nProfileSize,0);
   
   if (m_nWavePropagationModel == MODEL_CSHORE)
   {
      // We are using CShore to propagate the waves, so create an input file for this profile
      double
         dCShoreTimeStep = 3600,     // In seconds, not important because we are not using CShore to erode the profile, just to get the hydrodynamics
         dSurgeLevel = 0.0,          // Not used, but in the future we might include surge in the calculations
         dWaveFriction = CSHORE_FRICTION_FACTOR;
      
      // Set up vectors for the coastline-normal profile elevations. The length of this vector line is given by the number of cells 'under' the profile. Thus each point on the vector relates to a single cell in the grid. This assumes that all points on the profile vector are equally spaced (not quite true, depends on the orientation of the line segments which comprise the profile)
      vector<double>
         VdProfileZ,                  // Initial (pre-erosion) elevation of both consolidated and unconsolidated sediment for cells 'under' the profile, in CShore units
         VdProfileDistXY;             // Along-profile distance measured from the seaward limit, in CShore units
      
      // The elevation of each of these profile points is the elevation of the centroid of the cell that is 'under' the point. However we cannot always be confident that this is the 'true' elevation of the point on the vector since (unless the profile runs planview N-S or W-E) the vector does not always run exactly through the centroid of the cell
      int nRet = nGetThisProfileElevationVectorsForCShore(nCoast, nProfile, nProfileSize, &VdProfileDistXY, &VdProfileZ);
      if (nRet != RTN_OK)
         return nRet;
      
      if (VdProfileDistXY.empty())
      {
         // VdProfileDistXY has not been populated
         LogStream << m_ulIteration << ": VdProfileDistXY is empty for profile " << nProfile << endl;
         
         return RTN_ERR_CSHORE_EMPTY_PROFILE;
      }
      
      // Move to the CShore folder
      nRet = chdir(CSHOREDIR.c_str());
      
      char szBuf[BUF_SIZE] = "";
      string strCWD = getcwd(szBuf, BUF_SIZE);
      
      // TODO Andres check this re. CShore input requirements. Constrain the wave to normal angle to be between -80 and 80 degrees, this is a requirement of CShore
      dWaveToNormalAngle = tMax(dWaveToNormalAngle, -80.0);
      dWaveToNormalAngle = tMin(dWaveToNormalAngle, 80.0);
      
      // Create the file which will be read by CShore
      nRet = nCreateCShoreInfile(dCShoreTimeStep, m_dWavePeriod, dProfileDeepWaterWaveHeight, dWaveToNormalAngle, dSurgeLevel, dWaveFriction, &VdProfileDistXY, &VdProfileZ);
      if (nRet != RTN_OK)
         return nRet;
      
      // Set the error flag: this will be changed to 0 within CShore if CShore returns correctly
      nRet = -1;
      
      // Run CShore for this profile
      cshore(&nRet);
      
      // Check for error
      if (nRet != 0)
         return RTN_ERR_CSHORE_ERROR;
      
      // Fetch the CShore results by reading files written by CShore
      vector<double>
         VdFreeSurfaceStd(VdProfileDistXY.size(), 0),          // This is converted to Hrms by Hrms = sqr(8)*FreeSurfaceStd
         VdSinWaveAngleRadians(VdProfileDistXY.size(), 0),     // This is converted to deg by asin(VdSinWaveAngleRadians)*(180/pi)
         VdFractionBreakingWaves(VdProfileDistXY.size(), 0);   // Is 0 if no wave breaking, and 1 if all waves breaking
      
      string
         strOSETUP = "OSETUP",
         strOYVELO = "OYVELO",
         strOPARAM = "OPARAM";
      
      nRet = nReadCShoreOutput(&strOSETUP, 4, 4, &VdProfileDistXY, &VdFreeSurfaceStd);
      if (nRet != RTN_OK)
         return nRet;
      
      nRet = nReadCShoreOutput(&strOYVELO, 4, 2, &VdProfileDistXY, &VdSinWaveAngleRadians);
      if (nRet != RTN_OK)
         return nRet;
      
      nRet = nReadCShoreOutput(&strOPARAM, 4, 3, &VdProfileDistXY, &VdFractionBreakingWaves);
      if (nRet != RTN_OK)
         return nRet;
      
      // Clean up the CShore outputs
#ifdef _WIN32
      nRet = system("./clean.bat")
#else
      nRet = system("./clean.sh");
#endif
      
      // And return to the CoastalME folder
      nRet = chdir(m_strCMEDir.c_str());
      
      // Convert CShore outputs to wave height and wave direction and update wave profile attributes
      for (int nProfilePoint = (nProfileSize-1); nProfilePoint >= 0; nProfilePoint--)
      {
         int
            nX = pProfile->pPtiGetCellInProfile(nProfilePoint)->nGetX(),
            nY = pProfile->pPtiGetCellInProfile(nProfilePoint)->nGetY();
         
         VdWaveHeight[nProfilePoint] = sqrt(8) * VdFreeSurfaceStd[nProfilePoint];
         
         double dAlpha = asin(VdSinWaveAngleRadians[nProfilePoint]) * (180/PI);
         if (nSeaHand == LEFT_HANDED)
            VdWaveDirection[nProfilePoint] = dKeepWithin360(dAlpha + 90 + dFluxOrientationThis);
         else
            VdWaveDirection[nProfilePoint] = dKeepWithin360(dAlpha + 270 + dFluxOrientationThis);
         
         if ((VdFractionBreakingWaves[nProfilePoint] >= 0.99) & (! bBreaking))
         {
            bBreaking = true;
            dProfileBreakingWaveHeight = VdWaveHeight[nProfilePoint];
            dProfileBreakingWaveOrientation = VdWaveDirection[nProfilePoint];
            dProfileBreakingDepth = m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth();  // Water depth for the cell 'under' this point in the profile
            nProfileBreakingDist = nProfilePoint;
            
            //             LogStream << m_ulIteration << ": CShore breaking at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} nProfile = " << nProfile << ", nProfilePoint = " << nProfilePoint << ", dBreakingWaveHeight = " << dBreakingWaveHeight << ", dBreakingWaveOrientation = " << dBreakingWaveOrientation << ", dProfileBreakingDepth = " << dProfileBreakingDepth << ", nProfileBreakingDist = " << nProfileBreakingDist << endl;
         }
         
         VbWaveIsBreaking[nProfilePoint] = bBreaking;
      }
   }
   
   else if (m_nWavePropagationModel == MODEL_COVE)
   {
      // We are using COVE's linear wave theory to propagate the waves
      double dDepthLookupMax = m_dWaveDepthRatioForWaveCalcs * dProfileDeepWaterWaveHeight;
      
      // Go landwards along the profile, calculating wave height and wave angle for every inundated point on the profile (don't do point zero, this is on the coastline) until the waves start to break  after breaking wave height is assumed to decrease linearly to zero at the shoreline and wave angle is equalt to wave angle at breaking
      for (int nProfilePoint = (nProfileSize-1); nProfilePoint > 0; nProfilePoint--)
      {
         int
            nX = pProfile->pPtiGetCellInProfile(nProfilePoint)->nGetX(),
            nY = pProfile->pPtiGetCellInProfile(nProfilePoint)->nGetY();
         
         // Safety check
         if (! m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea())
            continue;
         
         double dSeaDepth = m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth();  // Water depth for the cell 'under' this point in the profile
         
         if (dSeaDepth > dDepthLookupMax)
         {
            // Sea depth is too large relative to wave height to feel the bottom, so do nothing since each cell under the profile has already been initialised with deep water wave height and wave direction
            dProfileWaveHeight = dProfileDeepWaterWaveHeight;
            dProfileWaveOrientation = dProfileDeepWaterWaveOrientation;
         }
         else
         {
            if (! bBreaking)
            {
               // Start calculating wave properties using linear wave theory
               double dL = m_dL_0 * sqrt(tanh((2 * PI * dSeaDepth) / m_dL_0));               // Wavelength (m) in intermediate-shallow waters
               double dC = m_dC_0 * tanh((2 * PI * dSeaDepth) / dL);                         // Wave speed (m/s) set by dSeaDepth, dL and m_dC_0
               double dk = 2 * PI / dL;                                                      // Wave number (1/m)
               double dn = ((2 * dSeaDepth * dk) / (sinh(2 * dSeaDepth * dk)) + 1) / 2;      // Shoaling factor
               double dKs = sqrt(m_dC_0 / (dn * dC * 2));                                    // Shoaling coefficient
               double dAlpha = (180 / PI) * asin((dC / m_dC_0) * sin((PI / 180) * dWaveToNormalAngle));  // Calculate angle between wave direction and the normal to the coast tangent
               double dKr = sqrt(cos((PI / 180) * dWaveToNormalAngle) / cos((PI / 180) * dAlpha));       // Refraction coefficient
               dProfileWaveHeight = dProfileDeepWaterWaveHeight * dKs * dKr;                               // Calculate wave height, based on the previous (more seaward) wave height
               if (nSeaHand == LEFT_HANDED)
                  dProfileWaveOrientation = dKeepWithin360(dAlpha + 90 + dFluxOrientationThis);
               else
                  dProfileWaveOrientation = dKeepWithin360(dAlpha + 270 + dFluxOrientationThis);
               
               // Test to see if the wave breaks at this depth
               if (dProfileWaveHeight > (dSeaDepth * m_dBreakingWaveHeightDeptRatio))
               {
                  // It does
                  bBreaking = true;
                  dProfileBreakingWaveHeight = dProfileWaveHeight;
                  dProfileBreakingWaveOrientation = dProfileWaveOrientation;
                  dProfileBreakingDepth = dSeaDepth;
                  nProfileBreakingDist = nProfilePoint;
               }
            }
            else
            {
               // Wave has already broken
               dProfileWaveOrientation = dProfileBreakingWaveOrientation;                          // Wave orientation remains equal to wave orientation at breaking
               
               //dProfileWaveHeight = dProfileBreakingWaveHeight * (nProfilePoint / nProfileBreakingDist);    // Wave height decreases linearly to zero at shoreline
               dProfileWaveHeight = dSeaDepth * m_dBreakingWaveHeightDeptRatio;                    // Wave height is limited by depth
            }
         }
         
         // Save current wave attributes
         VdWaveDirection[nProfilePoint] = dProfileWaveOrientation;
         VdWaveHeight[nProfilePoint] = dProfileWaveHeight;
         VbWaveIsBreaking[nProfilePoint] = bBreaking;
      }
   }
   
   // Go landwards along the profile, fetching the calculated wave height and wave angle for every inundated point on this profile
   for (int nProfilePoint = (nProfileSize-1); nProfilePoint > 0; nProfilePoint--)
   {
      int
         nX = pProfile->pPtiGetCellInProfile(nProfilePoint)->nGetX(),
         nY = pProfile->pPtiGetCellInProfile(nProfilePoint)->nGetY();
      
      // Safety check
      if (! m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea())
         continue;
      
      // Fetch the wave attributes calculated for this profile: wave height, wave angle, and whether is in the active zone
      double 
         dWaveHeight = VdWaveHeight[nProfilePoint],
         dWaveOrientation = VdWaveDirection[nProfilePoint];
      bBreaking = VbWaveIsBreaking[nProfilePoint];
      
      // Update RasterGrid wave properties
      m_pRasterGrid->m_Cell[nX][nY].SetInActiveZone(bBreaking);
      m_pRasterGrid->m_Cell[nX][nY].SetWaveHeight(dWaveHeight);
      m_pRasterGrid->m_Cell[nX][nY].SetWaveOrientation(dWaveOrientation);
      
      // And store the wave properties for this point in the all-profiles vectors
      pVnX->push_back(nX);
      pVnY->push_back(nY);
      pVdHeightX->push_back(dWaveHeight * sin(dWaveOrientation * PI/180));
      pVdHeightY->push_back(dWaveHeight * cos(dWaveOrientation * PI/180));
      pVbBreaking->push_back(bBreaking);
   }
   
   // Update wave attributes along the coastline object
   // Wave height at the coast is calculated irrespectively if waves are breaking or not
   
   //cout << "Wave Height at the coast is " << VdWaveHeight[0] << endl;
      m_VCoast[nCoast].SetCoastWaveHeight(nCoastPoint, VdWaveHeight[0]);
   
   if (nProfileBreakingDist > 0)
   {
      // This coast point is in the active zone, so set breaking wave height, breaking wave angle, and depth of breaking for the coast point
      m_VCoast[nCoast].SetBreakingWaveHeight(nCoastPoint, dProfileBreakingWaveHeight);
      m_VCoast[nCoast].SetBreakingWaveOrientation(nCoastPoint, dProfileBreakingWaveOrientation);
      m_VCoast[nCoast].SetDepthOfBreaking(nCoastPoint, dProfileBreakingDepth);
      m_VCoast[nCoast].SetBreakingDistance(nCoastPoint, nProfileBreakingDist);
      
//       LogStream << m_ulIteration << ": nProfile = " << nProfile << ", nCoastPoint = " << nCoastPoint << " in active zone, dBreakingWaveHeight = " << dBreakingWaveHeight << endl;
   }
   else
   {
      // This coast point is not in the active zone, so breaking wave height, breaking wave angle, and depth of breaking are all meaningless
      m_VCoast[nCoast].SetBreakingWaveHeight(nCoastPoint, DBL_NODATA);
      m_VCoast[nCoast].SetBreakingWaveOrientation(nCoastPoint, DBL_NODATA);
      m_VCoast[nCoast].SetDepthOfBreaking(nCoastPoint, DBL_NODATA);
      m_VCoast[nCoast].SetBreakingDistance(nCoastPoint, INT_NODATA);
      
//       LogStream << m_ulIteration << ": nProfile = " << nProfile << ", nCoastPoint = " << nCoastPoint << " NOT in active zone" << endl;
   }
   
   return RTN_OK;
}


/*===============================================================================================================================
 
 Create the CShore input file
 
===============================================================================================================================*/
int CSimulation::nCreateCShoreInfile(double dTimestep, double dWavePeriod, double dHrms, double dWaveAngle , double dSurgeLevel, double dWaveFriction, vector<double> const* pVdXdist, vector<double> const* pVdBottomElevation)
{
   // Initialize inifile from infileTemplate
   int nRet = system("cp infileTemplate infile");       // The infileTemplate must be in the working directory
   if (nRet == -1)
      return RTN_ERR_CSHORE_OUTPUT_FILE;
   
   string strFName = "infile";
   
   // We have all the inputs in the CShore format, so we can create the input file
   ofstream OutStream;
   OutStream.open(strFName.c_str(), ios::out | ios::app);
   if (OutStream.fail())
   {
      // Error, cannot open file for writing
      LogStream << m_ulIteration << ": " << ERR << "cannot open " << strFName << " for output" << endl;
      return RTN_ERR_CSHORE_OUTPUT_FILE;
   }
   
   // OK, write to the file
   OutStream << setiosflags(ios::fixed) << setprecision(4); OutStream << setw(11) << m_dBreakingWaveHeightDeptRatio << "                        -> GAMMA" << endl;
   OutStream << setiosflags(ios::fixed) << setprecision(0); OutStream << setw(11) << 0 << "                        -> ILAB" << endl;
   OutStream << setiosflags(ios::fixed) << setprecision(0); OutStream << setw(11) << 1 << "                        -> NWAVE" << endl;
   OutStream << setiosflags(ios::fixed) << setprecision(0); OutStream << setw(11) << 1 << "                        -> NSURGE" << endl;
   
   OutStream << setiosflags(ios::fixed) << setprecision(2); OutStream << setw(11) << 0.0;
   OutStream << setiosflags(ios::fixed) << setprecision(4); OutStream << setw(11) << dWavePeriod << setw(11) << dHrms << setw(11) << dWaveAngle << endl;
   OutStream << setiosflags(ios::fixed) << setprecision(2); OutStream << setw(11) << dTimestep;
   OutStream << setiosflags(ios::fixed) << setprecision(4); OutStream << setw(11) << dWavePeriod << setw(11) << dHrms << setw(11) << dWaveAngle << endl;
   
   OutStream << setiosflags(ios::fixed) << setprecision(2); OutStream << setw(11) << 0.0;
   OutStream << setiosflags(ios::fixed) << setprecision(4); OutStream << setw(11) << dSurgeLevel << endl;
   OutStream << setiosflags(ios::fixed) << setprecision(2); OutStream << setw(11) << dTimestep;
   OutStream << setiosflags(ios::fixed) << setprecision(4); OutStream << setw(11) << dSurgeLevel << endl;
   
   OutStream << setw(8) << pVdXdist->size() << "                        -> NBINP" << endl;
   OutStream << setiosflags(ios::fixed) << setprecision(4);
   for (unsigned int i = 0; i < pVdXdist->size(); i++)
      OutStream << setw(11) << pVdXdist->at(i) << setw(11) << pVdBottomElevation->at(i) << setw(11) << dWaveFriction << endl;
   
   return RTN_OK;
}


/*===============================================================================================================================
 
 Get profile horizontal distance and bottom elevation vectors in CShore units
 
===============================================================================================================================*/
int CSimulation::nGetThisProfileElevationVectorsForCShore(int const nCoast, int const nProfile, int const nProfSize, vector<double>* VdDistXY, vector<double>* VdVZ)
{
   int
      nX1 = 0,
      nY1 = 0;
   
   double
      dXDist,
      dYDist,
      dProfileDistXY = 0;
   
   CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfile(nProfile);
   
   for (int i = nProfSize-1; i >= 0; i--)
   {
      int
         nX = pProfile->pPtiVGetCellsInProfile()->at(i).nGetX(),
         nY = pProfile->pPtiVGetCellsInProfile()->at(i).nGetY();
      
      // Calculate the horizontal distance relative to the most seaward point
      if (i == nProfSize-1)
         dProfileDistXY = 0;
      else
      {
         dXDist = dGridCentroidXToExtCRSX(nX1) - dGridCentroidXToExtCRSX(nX),
         dYDist = dGridCentroidYToExtCRSY(nY1) - dGridCentroidYToExtCRSY(nY),
         dProfileDistXY = dProfileDistXY + hypot(dXDist, dYDist);
      }
      
      // Update the cell indexes, the initial cell is now the previous one
      nX1 = nX;
      nY1 = nY;
      
      // Get the number of the highest layer with non-zero thickness
      int const nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();
      
      // Safety check
      if (nTopLayer == INT_NODATA)
         return RTN_ERR_NO_TOP_LAYER;
      
      if (nTopLayer == NO_NONZERO_THICKNESS_LAYERS)
         // TODO We are down to basement, decide what to do
         return RTN_OK;
      
      // Get the elevation for both consolidated and unconsolidated sediment on this cell
      double VdProfileZ = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev() - m_dThisTimestepSWL;
      VdVZ->push_back(VdProfileZ);
      
      // And store the X-Y plane distance from the start of the profile
      VdDistXY->push_back(dProfileDistXY);
   }
   
   return RTN_OK;
}


/*==============================================================================================================================
 
 The CShore lookup reads a CShore output file and creates a vector holding interpolated values. The interpolation may be simple linear or a more advanced hermite cubic method
 
==============================================================================================================================*/
int CSimulation::nReadCShoreOutput(string const* strCShoreFilename, int const nExpectedColumns, int const nCShorecolumn, vector<double> const* pVdDistXY, vector<double>* pVdMyInterpolatedValues)
{
   // TODO Make this a user input
   // Select the interpolation method to be used: 0 for simple linear or 1 for hermite cubic
   int InterpMethodOption = CSHORE_INTERPOLATION_LINEAR;
//       int InterpMethodOption = CSHORE_INTERPOLATION_HERMITE_CUBIC;
   
   // Read in the first column (contains XY distance relative to seaward limit) and CShore column from the CShore output file
   ifstream InStream;
   InStream.open(strCShoreFilename->c_str(), ios::in);
   
   // Did it open OK?
   if (! InStream.is_open())
   {
      // Error: cannot open CShore file for input
      LogStream << m_ulIteration << ": " << ERR << "cannot open " << *strCShoreFilename << " for input" << endl;
      
      return RTN_ERR_CSHORE_OUTPUT_FILE;
   }
   
   // Opened OK, so set up the vectors to hold the CShore output data
   vector<double>
      VdXYDistCShore,
      VdValuesCShore;
   
   // And read in the data
   int 
      n = -1,
      nExpectedRows = 0;
   string strLineIn;
   while (getline(InStream, strLineIn))
   {
      n++;
      if (n == 0)
      {
         // The header line
         vector<string> VstrItems = strSplit(&strLineIn, SPACE);
         nExpectedRows = stoi(VstrItems[1].c_str());
      }
      else
      {
         // The data
         vector<string> VstrItems = strSplit(&strLineIn, SPACE);
         
         int nCols = VstrItems.size();         
         if (nCols != nExpectedColumns)
         {
            // Error: did not get nExpectedColumns CShore output columns
            LogStream << m_ulIteration << ": " << ERR << "expected " << nExpectedColumns << " CShore output columns but read " << nCols << " columns from header section of file " << strCShoreFilename << endl;
            
            return RTN_ERR_CSHORE_OUTPUT_FILE;        
         }
         // Number of columns is OK
         VdXYDistCShore.push_back(atof(VstrItems[0].c_str()));         
         VdValuesCShore.push_back(atof(VstrItems[nCShorecolumn-1].c_str()));
      }         
   }
   
   // Check that we have read nExpectedRows from the file
   int nReadRows = VdXYDistCShore.size();
   if (nReadRows != nExpectedRows)
   {
      // Error: did not get nExpectedRows CShore output rows
      LogStream << m_ulIteration << ": " << ERR << "expected " << nExpectedRows << " CShore output rows, but read " << nReadRows << " rows from file " << strCShoreFilename << endl;
      
      return RTN_ERR_CSHORE_OUTPUT_FILE;
   }
   
   if (nReadRows < 2)
   {
      // Error: cannot interpolate values if we have only one value
      LogStream << m_ulIteration << ": " << ERR << "only " << nReadRows << " CShore output rows in file " << strCShoreFilename << ". Try lengthening the coastline normals." << endl;
      
      return RTN_ERR_CSHORE_OUTPUT_FILE;
   }   
   
   // The output is OK, so change the origin of the across-shore distance from the CShore convention to the one used here (i.e. with the origin at the shoreline)
   vector<double>  VdXYDistCME(nReadRows, 0);
   for (int i = 0; i < nReadRows; i++)
      VdXYDistCME[i] = VdXYDistCShore[nReadRows-1] - VdXYDistCShore[i];
   
   // Reverse the CShore XYdistance and value vectors (i.e. first point is at the shoreline and must be in strictly ascending order)
   reverse(VdXYDistCME.begin(), VdXYDistCME.end());
   reverse(VdValuesCShore.begin(), VdValuesCShore.end());
   
   // Now we have everything ready to do the interpolation
   if (InterpMethodOption == CSHORE_INTERPOLATION_HERMITE_CUBIC)
   {
      // Using the hermite cubic approach: calculate the first derivative of CShore values (needed for the hermite interpolant)
      vector<double> VdValuesCShoreDeriv(nReadRows, 0);
      for (int i = 1; i < nReadRows-1; i++)
      {
         // Calculate the horizontal distance increment between two adjacent points (not always the same distance because it depend on profile-cells centroid location)
         double dX = VdXYDistCME[i+1] - VdXYDistCME[i-1];    // This is always positive
         VdValuesCShoreDeriv[i] = (VdValuesCShore[i+1] - VdValuesCShore[i-1]) / (2 * dX);
      }
      
      VdValuesCShoreDeriv[0] = VdValuesCShoreDeriv[1];
      VdValuesCShoreDeriv[nReadRows-1] = VdValuesCShoreDeriv[nReadRows-2];
      
      // Interpolate the CShore values
      int nSize = pVdDistXY->size();
      vector<double>
         VdDistXYCopy(pVdDistXY->begin(), pVdDistXY->end()),
         //dVInter(nSize, 0.),
         VdDeriv(nSize, 0),         // First derivative at the sample points: calculated by the spline function but not subsequently used
         VdDeriv2(nSize, 0),        // Second derivative at the sample points, ditto
         VdDeriv3(nSize, 0);        // Third derivative at the sample points, ditto
      
      // Calculate the value of erosion potential (is a -ve value) for each of the sample values of DepthOverDB, and store it for use in the look-up function
      hermite_cubic_spline_value(nReadRows, &(VdXYDistCME.at(0)), &(VdValuesCShore.at(0)), &(VdValuesCShoreDeriv.at(0)), nSize, &(VdDistXYCopy[0]), &(pVdMyInterpolatedValues->at(0)), &(VdDeriv[0]), &(VdDeriv2[0]), &(VdDeriv3[0]));
   }
   else
   {
      // Using the simple linear approach
      
      // TODO We get an error if pVdDistXY->size() == 1, Andres to check why we get this
//       if (pVdDistXY->size() == 1)
      
      
      vector<double> VdDistXYCopy(pVdDistXY->begin(), pVdDistXY->end());
      *pVdMyInterpolatedValues = VdInterp1(&VdXYDistCME, &VdValuesCShore, &VdDistXYCopy);
   }
   
   return RTN_OK;   
}


/*===============================================================================================================================

 Modifies the wave breaking properties at coastline points of profiles within the shadow zone
 
===============================================================================================================================*/
void CSimulation::ModifyBreakingWavePropertiesWithinShadowZoneToCoastline(int const nCoast, int const nProfIndex)
{
   int nProfile = m_VCoast[nCoast].nGetProfileFromAlongCoastProfileIndex(nProfIndex);
   CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfile(nProfile);
   
   // Only do this for profiles without problems, including the start and end-of-coast profile
   if (! pProfile->bOKIncStartAndEndOfCoast())
      return;
   
   int
      nThisCoastPoint = pProfile->nGetNumCoastPoint(),
      nProfileSize = pProfile->nGetNumCellsInProfile();
   
   int nThisBreakingDist = m_VCoast[nCoast].nGetBreakingDistance(nThisCoastPoint);
   double
      dThisBreakingWaveHeight = m_VCoast[nCoast].dGetBreakingWaveHeight(nThisCoastPoint),       // This could be DBL_NODATA
      dThisBreakingWaveOrientation = m_VCoast[nCoast].dGetBreakingWaveOrientation(nThisCoastPoint),
      dThisBreakingDepth = m_VCoast[nCoast].dGetDepthOfBreaking(nThisCoastPoint);
   bool
      bModfiedWaveHeightisBreaking = false,
      bProfileIsinShadowZone = false;
   
   // Traverse the profile landwards, checking if any profile cell is within the shadow zone
   for (int nProfilePoint = (nProfileSize-1); nProfilePoint >= 0; nProfilePoint--)
   {
      int
         nX = pProfile->pPtiGetCellInProfile(nProfilePoint)->nGetX(),
         nY = pProfile->pPtiGetCellInProfile(nProfilePoint)->nGetY();
      
      // If there is any cell profile  within the shadow zone and waves are breaking then modify wave breaking properties otherwise continue
      if (m_pRasterGrid->m_Cell[nX][nY].bIsinAnyShadowZone())
      {
         bProfileIsinShadowZone = true;
         
         // Check if the new wave height is breaking
         double
            dSeaDepth = m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth(),
            dWaveHeight = m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight(),
            dWaveOrientation = m_pRasterGrid->m_Cell[nX][nY].dGetWaveOrientation();
         
         if (dWaveHeight > (dSeaDepth * m_dBreakingWaveHeightDeptRatio) && (! bModfiedWaveHeightisBreaking))
         {
            // It is breaking
            bModfiedWaveHeightisBreaking = true;
            
            dThisBreakingWaveHeight = dWaveHeight;
            dThisBreakingWaveOrientation = dWaveOrientation;
            dThisBreakingDepth = dSeaDepth;
            nThisBreakingDist = nProfilePoint;
         }
      }
   }
   
   // Update breaking wave properties along coastal line object (Wave height, dir, distance). TODO update the active zone cells
   if (bProfileIsinShadowZone && bModfiedWaveHeightisBreaking) // Modified wave height is still breaking
   {
      // This coast point is in the active zone, so set breaking wave height, breaking wave angle, and depth of breaking for the coast point
      m_VCoast[nCoast].SetBreakingWaveHeight(nThisCoastPoint, dThisBreakingWaveHeight);
      m_VCoast[nCoast].SetBreakingWaveOrientation(nThisCoastPoint, dThisBreakingWaveOrientation);
      m_VCoast[nCoast].SetDepthOfBreaking(nThisCoastPoint, dThisBreakingDepth);
      m_VCoast[nCoast].SetBreakingDistance(nThisCoastPoint, nThisBreakingDist);
      
      //       LogStream << m_ulIteration << ": nProfile = " << nProfile << ", nCoastPoint = " << nCoastPoint << " in active zone, dBreakingWaveHeight = " << dBreakingWaveHeight << endl;
   }
   /*else if (bProfileIsinShadowZone && !bModfiedWaveHeightisBreaking)
    *   {
    *      // This coast point is no longer in the active zone
    *      m_VCoast[nCoast].SetBreakingWaveHeight(nThisCoastPoint, DBL_NODATA);
    *      m_VCoast[nCoast].SetBreakingWaveOrientation(nThisCoastPoint, DBL_NODATA);
    *      m_VCoast[nCoast].SetDepthOfBreaking(nThisCoastPoint, DBL_NODATA);
    *      m_VCoast[nCoast].SetBreakingDistance(nThisCoastPoint, INT_NODATA);
    * 
    * //       LogStream << m_ulIteration << ": nProfile = " << nProfile << ", nCoastPoint = " << nCoastPoint << " NOT in active zone" << endl;
}*/
   
   return;
}


/*===============================================================================================================================

 Interpolates wave properties from profiles to the coastline points between two profiles. Do this by weighting the wave properties so that they change smoothly between the two profiles

===============================================================================================================================*/
void CSimulation::InterpolateWavePropertiesToCoastline(int const nCoast, int const nProfIndex, int const nNumProfiles)
{
   int nProfile = m_VCoast[nCoast].nGetProfileFromAlongCoastProfileIndex(nProfIndex);
   CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfile(nProfile);

   // Only do this for profiles without problems, including the start-of-coast profile (but not the end-of-coast profile)
   if (! pProfile->bOKIncStartOfCoast())
      return;

   int nThisCoastPoint = pProfile->nGetNumCoastPoint();

   // For the breaking wave stuff, to go into the in-between coastline points
   int const nThisBreakingDist = m_VCoast[nCoast].nGetBreakingDistance(nThisCoastPoint);
   double
      dThisBreakingWaveHeight = m_VCoast[nCoast].dGetBreakingWaveHeight(nThisCoastPoint),       // This could be DBL_NODATA
      dThisCoastWaveHeight = m_VCoast[nCoast].dGetCoastWaveHeight(nThisCoastPoint),       
      dThisBreakingWaveOrientation = m_VCoast[nCoast].dGetBreakingWaveOrientation(nThisCoastPoint),
      dThisBreakingDepth = m_VCoast[nCoast].dGetDepthOfBreaking(nThisCoastPoint);

   // Get the number of the next profile along the coast. If this next profile has a problem, go to the one after that, etc
   int nNextProfile = -1;
   for (int nNextProfIndex = nProfIndex+1; nNextProfIndex < nNumProfiles; nNextProfIndex++)
   {
      nNextProfile = m_VCoast[nCoast].nGetProfileFromAlongCoastProfileIndex(nNextProfIndex);

      if (m_VCoast[nCoast].pGetProfile(nNextProfile)->bOKIncStartAndEndOfCoast())
         break;
   }

   int const
      nNextCoastPoint = m_VCoast[nCoast].pGetProfile(nNextProfile)->nGetNumCoastPoint(),
      nDistBetween = nNextCoastPoint - nThisCoastPoint - 1;

   if (nDistBetween <= 0)
      // Nothing to do
      return;

   int nNextBreakingDist = m_VCoast[nCoast].nGetBreakingDistance(nNextCoastPoint);
   double
      dNextBreakingWaveHeight = m_VCoast[nCoast].dGetBreakingWaveHeight(nNextCoastPoint),          // This could be DBL_NODATA
      dNextCoastWaveHeight = m_VCoast[nCoast].dGetCoastWaveHeight(nNextCoastPoint),          
      dNextBreakingWaveOrientation = m_VCoast[nCoast].dGetBreakingWaveOrientation(nNextCoastPoint),
      dNextBreakingDepth = m_VCoast[nCoast].dGetDepthOfBreaking(nNextCoastPoint);

   // If both this profile and the next profile are not in the active zone, then do no more
   if ((dThisBreakingWaveHeight == DBL_NODATA) && (dNextBreakingWaveHeight == DBL_NODATA))
   {
//       LogStream << m_ulIteration << ": both profile " << nProfile << " at coast point " << nThisCoastPoint << ", and profile " << nNextProfile << " at coast point " << nNextCoastPoint << ", are not in the active zone" << endl;
      return;
   }

   // OK, at least one of the two profiles is in the active zone
   if (dThisBreakingWaveHeight == DBL_NODATA)
   {
      // The next profile must be in the active zone, so use values from the next profile
      for (int n = nThisCoastPoint; n < nNextCoastPoint; n++)
      {
         // Set the breaking wave height, breaking wave angle, and depth of breaking for this coast point
         m_VCoast[nCoast].SetBreakingWaveHeight(n, dNextBreakingWaveHeight);
	 m_VCoast[nCoast].SetCoastWaveHeight(n, dNextCoastWaveHeight);
         m_VCoast[nCoast].SetBreakingWaveOrientation(n, dNextBreakingWaveOrientation);
         m_VCoast[nCoast].SetDepthOfBreaking(n, dNextBreakingDepth);
         m_VCoast[nCoast].SetBreakingDistance(n, nNextBreakingDist);
      }

      return;
   }

   if (dNextBreakingWaveHeight == DBL_NODATA)
   {
      // This profile must be in the active zone, so use values from this profile
      for (int n = nThisCoastPoint+1; n <= nNextCoastPoint; n++)
      {
         // Set the breaking wave height, breaking wave angle, and depth of breaking for this coast point
         m_VCoast[nCoast].SetBreakingWaveHeight(n, dThisBreakingWaveHeight);
	 m_VCoast[nCoast].SetCoastWaveHeight(n, dThisCoastWaveHeight);
         m_VCoast[nCoast].SetBreakingWaveOrientation(n, dThisBreakingWaveOrientation);
         m_VCoast[nCoast].SetDepthOfBreaking(n, dThisBreakingDepth);
         m_VCoast[nCoast].SetBreakingDistance(n, nThisBreakingDist);
      }

      return;
   }

   // The values for both this profile point and the next profile point are fine, so do a weighted interpolation between this profile and the next profile
   for (int n = nThisCoastPoint+1; n < nNextCoastPoint; n++)
   {
      int nDist = n - nThisCoastPoint;

      double
         dBreakingWaveHeight = 0,
	 dCoastWaveHeight = 0,
         dBreakingWaveOrientation = 0,
         dBreakingDepth = 0,
         dBreakingDist = 0;

      if ((dNextBreakingDepth > 0) && (dThisBreakingDepth > 0))
      {
         double
            dThisWeight = (nDistBetween - nDist) / static_cast<double>(nDistBetween),
            dNextWeight = 1 - dThisWeight;

//          assert(dThisWeight >= 0);
//          assert(dNextWeight >= 0);

         dBreakingWaveHeight = (dThisWeight * dThisBreakingWaveHeight) + (dNextWeight * dNextBreakingWaveHeight),
         dBreakingWaveOrientation = (dThisWeight * dThisBreakingWaveOrientation) + (dNextWeight * dNextBreakingWaveOrientation),
         dBreakingDepth = (dThisWeight * dThisBreakingDepth) + (dNextWeight * dNextBreakingDepth),
         dBreakingDist = (dThisWeight * nThisBreakingDist) + (dNextWeight * nNextBreakingDist);
      }
      else if (dThisBreakingDepth > 0)
      {
         dBreakingWaveHeight = dNextBreakingWaveHeight,
         dBreakingWaveOrientation = dNextBreakingWaveOrientation,
         dBreakingDepth = dNextBreakingDepth,
         dBreakingDist = nNextBreakingDist;
      }
      else if (dNextBreakingDepth > 0)
      {
         dBreakingWaveHeight = dThisBreakingWaveHeight,
         dBreakingWaveOrientation = dThisBreakingWaveOrientation,
         dBreakingDepth = dThisBreakingDepth,
         dBreakingDist = nThisBreakingDist;
      }

      // Always interpolate the wave height at coast
      double
            dThisWeight = (nDistBetween - nDist) / static_cast<double>(nDistBetween),
            dNextWeight = 1 - dThisWeight;
	dCoastWaveHeight = (dThisWeight * dThisCoastWaveHeight) + (dNextWeight * dNextCoastWaveHeight);   
	
      // Set the breaking wave height, breaking wave angle, and depth of breaking for this coast point
      m_VCoast[nCoast].SetBreakingWaveHeight(n, dBreakingWaveHeight);
      m_VCoast[nCoast].SetCoastWaveHeight(n, dCoastWaveHeight);
      m_VCoast[nCoast].SetBreakingWaveOrientation(n, dBreakingWaveOrientation);
      m_VCoast[nCoast].SetDepthOfBreaking(n, dBreakingDepth);
      m_VCoast[nCoast].SetBreakingDistance(n, static_cast<int>(dRound(dBreakingDist)));
   }
}


/*===============================================================================================================================

 Calculates tangents to a coastline: the tangent is assumed to be the orientation of energy/sediment flux along a coast. The tangent is specified as an angle (in degrees) measured clockwise from north. Based on a routine by Martin Hurst

===============================================================================================================================*/
void CSimulation::CalcCoastTangents(int const nCoast)
{
   int nCoastSize = m_VCoast[nCoast].nGetCoastlineSize();

   for (int nCoastPoint = 0; nCoastPoint < nCoastSize; nCoastPoint++)
   {
      double
         dXDiff,
         dYDiff;

      if (nCoastPoint == 0)
      {
         // For the point at the start of the coastline: use the straight line from 'this' point to the next point
         CGeom2DPoint PtThis = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint);              // In external CRS
         CGeom2DPoint PtAfter = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint+1);           // In external CRS

         dXDiff = PtAfter.dGetX() - PtThis.dGetX();
         dYDiff = PtAfter.dGetY() - PtThis.dGetY();
      }
      else if (nCoastPoint == nCoastSize-1)
      {
         // For the point at the end of the coastline: use the straight line from the point before to 'this' point
         CGeom2DPoint PtBefore = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint-1);          // In external CRS
         CGeom2DPoint PtThis = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint);              // In external CRS

         dXDiff = PtThis.dGetX() - PtBefore.dGetX();
         dYDiff = PtThis.dGetY() - PtBefore.dGetY();
      }
      else
      {
         // For coastline points not at the start or end of the coast: start with a straight line which links the coastline points before and after 'this' coastline point
         CGeom2DPoint PtBefore = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint-1);           // In external CRS
         CGeom2DPoint PtAfter = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint+1);            // In external CRS

         dXDiff = PtAfter.dGetX() - PtBefore.dGetX();
         dYDiff = PtAfter.dGetY() - PtBefore.dGetY();
      }

      // Calculate angle between line and north point, measured clockwise (the azimuth)
      if (bFPIsEqual(dYDiff, 0, TOLERANCE))
      {
         // The linking line runs either W-E or E-W
         if (dXDiff > 0)
            // It runs W to E
            m_VCoast[nCoast].SetFluxOrientation(nCoastPoint, 90);
         else
            // It runs E to W
            m_VCoast[nCoast].SetFluxOrientation(nCoastPoint, 270);
      }
      else if (bFPIsEqual(dXDiff, 0, TOLERANCE))
      {
         // The linking line runs N-S or S-N
         if (dYDiff > 0)
            // It runs S to N
            m_VCoast[nCoast].SetFluxOrientation(nCoastPoint, 0);
         else
            // It runs N to S
            m_VCoast[nCoast].SetFluxOrientation(nCoastPoint, 180);
      }
      else
      {
         // The linking line runs neither W-E nor N-S so we have to work a bit harder to find the angle between it and the azimuth
         double dAzimuthAngle;
         if (dXDiff > 0)
            dAzimuthAngle= (180 / PI) * (PI * 0.5 - atan(dYDiff / dXDiff));
         else
            dAzimuthAngle = (180 / PI) * (PI * 1.5 - atan(dYDiff / dXDiff));

         m_VCoast[nCoast].SetFluxOrientation(nCoastPoint, dAzimuthAngle);
      }

//      LogStream << m_ulIteration << ": nCoastPoint = " << nCoastPoint << " FluxOrientation = " << m_VCoast[nCoast].dGetFluxOrientation(nCoastPoint) << endl;
   }
}


/*===============================================================================================================================

 Calculates an average d50 for each polygon, and also fills in 'holes' in active zone and wave calcs i.e. orphan cells which should have been included in the active zone but which have been omitted because of rounding problems

===============================================================================================================================*/
void CSimulation::CalcD50AndFillWaveCalcHoles(void)
{
   vector<int> VnPolygonD50Count(m_nGlobalPolygonID+1, 0);
   vector<double> VdPolygonD50(m_nGlobalPolygonID+1, 0);

   for (int nX = 0; nX < m_nXGridMax; nX++)
   {
      for (int nY = 0; nY < m_nYGridMax; nY++)
      {
         if (m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea())
         {
            // This is a sea cell, is it in the active zone?
            if (m_pRasterGrid->m_Cell[nX][nY].bIsInActiveZone())
            {
               // It is, so does it have unconsolidated sediment on it?
               double dTmpd50 = m_pRasterGrid->m_Cell[nX][nY].dGetUnconsD50();
               if (dTmpd50 != DBL_NODATA)
               {
                  // It does, so which polygon is it in?
                  int nID = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();
                  if (nID != INT_NODATA)
                  {
                     VnPolygonD50Count[nID]++;
                     VdPolygonD50[nID] += dTmpd50;
                  }
               }
            }

/*            
            // Now fill in wave calc holes, start by looking at the cell's N-S and W-E neighbours
            int
               nXTmp,
               nYTmp,
               nActive = 0,
               nShadow = 0,
               nShadowNum = 0,
               nDownDrift = 0,
               nDownDriftNum = 0,
               nCoast = 0,
               nRead = 0;
            double
               dWaveHeight = 0,
               dWaveAngle = 0;

            // North
            nXTmp = nX;
            nYTmp = nY-1;
            if (bIsWithinValidGrid(nXTmp, nYTmp))
            {
               if (m_pRasterGrid->m_Cell[nXTmp][nYTmp].bIsInContiguousSea())
               {
                  nRead++;
                  dWaveHeight += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetWaveHeight();
                  dWaveAngle += (m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetWaveOrientation());

                  if (m_pRasterGrid->m_Cell[nXTmp][nYTmp].bIsInActiveZone())
                     nActive++;

                  int nTmp = m_pRasterGrid->m_Cell[nXTmp][nYTmp].nGetShadowZoneNumber();
                  if (nTmp != 0)
                  {
                     nShadow++;
                     nShadowNum = nTmp;
                  }
                  
                  nTmp = m_pRasterGrid->m_Cell[nXTmp][nYTmp].nGetDownDriftZoneNumber();
                  if ( nTmp > 0)
                  {
                     nDownDrift++;
                     nDownDriftNum = nTmp;
                  }          
               }
               else
                  nCoast++;
            }
               
            // East
            nXTmp = nX+1;
            nYTmp = nY;
            if (bIsWithinValidGrid(nXTmp, nYTmp))
            {
               if (m_pRasterGrid->m_Cell[nXTmp][nYTmp].bIsInContiguousSea())
               {
                  nRead++;
                  dWaveHeight += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetWaveHeight();
                  dWaveAngle += (m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetWaveOrientation());

                  if (m_pRasterGrid->m_Cell[nXTmp][nYTmp].bIsInActiveZone())
                     nActive++;

                  int nTmp = m_pRasterGrid->m_Cell[nXTmp][nYTmp].nGetShadowZoneNumber();
                  if (nTmp != 0)
                  {
                     nShadow++;
                     nShadowNum = nTmp;
                  }
                  
                  nTmp = m_pRasterGrid->m_Cell[nXTmp][nYTmp].nGetDownDriftZoneNumber();
                  if ( nTmp > 0)
                  {
                     nDownDrift++;
                     nDownDriftNum = nTmp;
                  }
               }
               else
                  nCoast++;
            }
               
            // South
            nXTmp = nX;
            nYTmp = nY+1;
            if (bIsWithinValidGrid(nXTmp, nYTmp))
            {
               if (m_pRasterGrid->m_Cell[nXTmp][nYTmp].bIsInContiguousSea())
               {
                  nRead++;
                  dWaveHeight += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetWaveHeight();
                  dWaveAngle += (m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetWaveOrientation());

                  if (m_pRasterGrid->m_Cell[nXTmp][nYTmp].bIsInActiveZone())
                     nActive++;

                  int nTmp = m_pRasterGrid->m_Cell[nXTmp][nYTmp].nGetShadowZoneNumber();
                  if (nTmp != 0)
                  {
                     nShadow++;
                     nShadowNum = nTmp;
                  }
                  
                  nTmp = m_pRasterGrid->m_Cell[nXTmp][nYTmp].nGetDownDriftZoneNumber();
                  if ( nTmp > 0)
                  {
                     nDownDrift++;
                     nDownDriftNum = nTmp;
                  }
               }
               else
                  nCoast++;
            }
               
            // West
            nXTmp = nX-1;
            nYTmp = nY;
            if (bIsWithinValidGrid(nXTmp, nYTmp))
            {
               if (m_pRasterGrid->m_Cell[nXTmp][nYTmp].bIsInContiguousSea())
               {
                  nRead++;
                  dWaveHeight += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetWaveHeight();
                  dWaveAngle += (m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetWaveOrientation());

                  if (m_pRasterGrid->m_Cell[nXTmp][nYTmp].bIsInActiveZone())
                     nActive++;

                  int nTmp = m_pRasterGrid->m_Cell[nXTmp][nYTmp].nGetShadowZoneNumber();
                  if (nTmp != 0)
                  {
                     nShadow++;
                     nShadowNum = nTmp;
                  }
                  
                  nTmp = m_pRasterGrid->m_Cell[nXTmp][nYTmp].nGetDownDriftZoneNumber();
                  if ( nTmp > 0)
                  {
                     nDownDrift++;
                     nDownDriftNum = nTmp;
                  } 
               }
               else
                  nCoast++;
            }
            
            if (nRead > 0)
            {
               // Calculate the average of neighbours
               dWaveHeight /= nRead;
               dWaveAngle /= nRead;
               dWaveAngle = dKeepWithin360(dWaveAngle);

               // If this sea cell has four active-zone neighbours, then it must also be in the active zone: give it wave height and orientation which is the average of its neighbours
               if (nActive == 4)
               {
                  m_pRasterGrid->m_Cell[nX][nY].SetInActiveZone(true);
                  m_pRasterGrid->m_Cell[nX][nY].SetWaveHeight(dWaveHeight);
                  m_pRasterGrid->m_Cell[nX][nY].SetWaveOrientation(dWaveAngle);                  
               }

               // If this sea cell has a wave height which is the same as its deep-water wave height, but its neighbours have a different average wave height, then give it the average of its neighbours
               double dDeepWaterWaveHeight = m_pRasterGrid->m_Cell[nX][nY].dGetDeepWaterWaveHeight();
               if ((m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight() == dDeepWaterWaveHeight) && (! bFPIsEqual(dDeepWaterWaveHeight, dWaveHeight, TOLERANCE)))
               {
                  m_pRasterGrid->m_Cell[nX][nY].SetWaveHeight(dWaveHeight);
               }

               // If this sea cell has a wave orientation which is the same as its deep-water wave orientation, but its neighbours have a different average wave orientation, then give it the average of its neighbours
               double dDeepWaterWaveOrientation = m_pRasterGrid->m_Cell[nX][nY].dGetDeepWaterWaveOrientation();
               if ((m_pRasterGrid->m_Cell[nX][nY].dGetWaveOrientation() == dDeepWaterWaveOrientation) && (! bFPIsEqual(dDeepWaterWaveOrientation, dWaveAngle, TOLERANCE)))
               {
                  m_pRasterGrid->m_Cell[nX][nY].SetWaveOrientation(dWaveAngle);
               }
               
               // Is this sea cell is not already marked as in a shadow zone (note could be marked as in a shadow zone but not yet processed: a -ve number)?
               int nShadowZoneCode = m_pRasterGrid->m_Cell[nX][nY].nGetShadowZoneNumber();
               if (nShadowZoneCode <= 0) 
               {
                  // If the cell has four neighbours which are all in a shadow zone, or four neighbours some of which are shadow zone and the remainder downdrift zone, or four neighbours some of which are shadow zone and the remainder coast; then it should also be in the shadow zone: give it the average of its neighbours
                  if (nShadow == 4)
                  {
                     m_pRasterGrid->m_Cell[nX][nY].SetShadowZoneNumber(nShadowNum);
                     m_pRasterGrid->m_Cell[nX][nY].SetWaveHeight(dWaveHeight);
                     m_pRasterGrid->m_Cell[nX][nY].SetWaveOrientation(dWaveAngle);
                  }
                  else if (nShadow + nDownDrift == 4)
                  {
                     m_pRasterGrid->m_Cell[nX][nY].SetShadowZoneNumber(nShadowNum);
                     m_pRasterGrid->m_Cell[nX][nY].SetWaveHeight(dWaveHeight);
                     m_pRasterGrid->m_Cell[nX][nY].SetWaveOrientation(dWaveAngle);
                  }
                  else if (nShadow + nCoast == 4)
                  {
                     m_pRasterGrid->m_Cell[nX][nY].SetShadowZoneNumber(nShadowNum);
                     m_pRasterGrid->m_Cell[nX][nY].SetWaveHeight(dWaveHeight);
                     m_pRasterGrid->m_Cell[nX][nY].SetWaveOrientation(dWaveAngle);
                  }
                  
               }

               // If this sea cell is not marked as in a downdrift zone but has four neighbours which are in a downdrift zone, then it should also be in the downdrift zone: give it the average of its neighbours
               int nDownDriftZoneCode = m_pRasterGrid->m_Cell[nX][nY].nGetDownDriftZoneNumber();
               if ((nDownDriftZoneCode == 0) && (nDownDrift == 4))
               {
                  m_pRasterGrid->m_Cell[nX][nY].SetDownDriftZoneNumber(nDownDriftNum);
                  m_pRasterGrid->m_Cell[nX][nY].SetWaveHeight(dWaveHeight);
                  m_pRasterGrid->m_Cell[nX][nY].SetWaveOrientation(dWaveAngle);
               }
            }
*/
          }
      }
   }

   // Calculate the average d50 for every polygon
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      for (int nPoly = 0; nPoly < m_VCoast[nCoast].nGetNumPolygons(); nPoly++)
      {
         CGeomCoastPolygon* pPolygon = m_VCoast[nCoast].pGetPolygon(nPoly);
         int nID = pPolygon->nGetGlobalID();

         if (VnPolygonD50Count[nID] > 0)
            VdPolygonD50[nID] /= VnPolygonD50Count[nID];

         pPolygon->SetAvgUnconsD50(VdPolygonD50[nID]);
      }
   }
}

