/*!
 *
 * \file read_input.cpp
 * \brief Reads non-GIS input files
 * \details TODO A more detailed description of these routines.
 * \author David Favis-Mortlock
 * \author Andres Payo

 * \date 2018
 * \copyright GNU General Public License
 *
 */

/*==============================================================================================================================

 This file is part of CoastalME, the Coastal Modelling Environment.

 CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public  License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

==============================================================================================================================*/
#include <stdlib.h>                 // for atof()
#include <fstream>
using std::ifstream;

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::ios;

#include "cme.h"
#include "simulation.h"


/*==============================================================================================================================

 The bReadIni member function reads the initialization file

==============================================================================================================================*/
bool CSimulation::bReadIni(void)
{
   m_strCMEIni = m_strCMEDir;
   m_strCMEIni.append(CME_INI);

   // The .ini file is assumed to be in the CoastalME executable's directory
   string strFilePathName(m_strCMEIni);

   // Tell the user what is happening
   cout << READFILELOC << strFilePathName << endl;

   // Create an ifstream object
   ifstream InStream;

   // Try to open .ini file for input
   InStream.open(strFilePathName.c_str(), ios::in);

   // Did it open OK?
   if (! InStream.is_open())
   {
      // Error: cannot open .ini file for input
      cerr << ERR << "cannot open " << strFilePathName << " for input" << endl;
      return false;
   }

   int i = 0;
   string strRec, strErr;

   while (getline(InStream, strRec))
   {
      // Trim off leading and trailing whitespace
      strRec = strTrimLeft(&strRec);
      strRec = strTrimRight(&strRec);

      // If it is a blank line or a comment then ignore it
      if ((! strRec.empty()) && (strRec[0] != QUOTE1) && (strRec[0] != QUOTE2))
      {
         // It isn't so increment counter
         i++;

         // Find the colon: note that lines MUST have a colon separating data from leading description portion
         size_t nPos = strRec.find(':');
         if (nPos == string::npos)
         {
            // Error: badly formatted line (no colon)
            cerr << ERR << "badly formatted line (no ':') in " << strFilePathName << endl << "'" << strRec << "'" << endl;
            return false;
         }

         if (nPos == strRec.size()-1)
         {
            // Error: badly formatted line (colon with nothing following)
            cerr << ERR << "badly formatted line (nothing following ':') in " << strFilePathName << endl << "'" << strRec << "'" << endl;
            return false;
         }

         // Strip off leading portion (the bit up to and including the colon)
         string strRH = strRec.substr(nPos+1);

         // Remove leading whitespace
         strRH = strTrimLeft(&strRH);

         // Look for a trailing comment, if found then terminate string at that point and trim off any trailing whitespace
         nPos = strRH.rfind(QUOTE1);
         if (nPos != string::npos)
            strRH = strRH.substr(0, nPos+1);

         nPos = strRH.rfind(QUOTE2);
         if (nPos != string::npos)
            strRH = strRH.substr(0, nPos+1);

         // Remove trailing whitespace
         strRH = strTrimRight(&strRH);

         switch (i)
         {
         case 1:
            // The main input run-data filename
            if (strRH.empty())
               strErr = "path and name of main datafile";
            else
            {
               // First check that we don't already have an input run-data filename, e.g. one entered on the command-line
               if (m_strDataPathName.empty())
               {
                  // We don't: so first check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
                  if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == '~') || (strRH[1] == ':'))
                     // It has an absolute path, so use it 'as is'
                     m_strDataPathName = strRH;
                  else
                  {
                     // It has a relative path, so prepend the CoastalME dir
                     m_strDataPathName = m_strCMEDir;
                     m_strDataPathName.append(strRH);
                  }
               }
            }
            break;

         case 2:
            // Path for CoastalME output
            if (strRH.empty())
               strErr = "path for CoastalME output";
            else
            {
               // Check for trailing slash on CoastalME output directory name (is vital)
               if (strRH[strRH.size()-1] != PATH_SEPARATOR)
                  strRH.append(&PATH_SEPARATOR);

               // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
               if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == '~') || (strRH[1] == ':'))
                  // It is an absolute path, so use it 'as is'
                  m_strOutPath = strRH;
               else
               {
                  // It is a relative path, so prepend the CoastalME dir
                  m_strOutPath = m_strCMEDir;
                  m_strOutPath.append(strRH);
               }
            }
            break;

         case 3:
            // Email address, only useful if running under Linux/Unix
            if (! strRH.empty())
            {
               // Something was entered, do rudimentary check for valid email address
               if (strRH.find('@') == string::npos)
                  strErr = "email address for messages";
               else
                  m_strMailAddress = strRH;
            }
            break;
         }

         // Did an error occur?
         if (! strErr.empty())
         {
            // Error in input to initialisation file
            cerr << ERR << "reading " << strErr << " in " << strFilePathName << endl << "'" << strRec << "'" << endl;
            InStream.close();

            return false;
         }
      }
   }

   InStream.close();
   return true;
}


/*==============================================================================================================================

 Reads the run details input file and does some initialization

==============================================================================================================================*/
bool CSimulation::bReadRunData(void)
{
   // Create an ifstream object
   ifstream InStream;

   // Try to open run details file for input
   InStream.open(m_strDataPathName.c_str(), ios::in);

   // Did it open OK?
   if (! InStream.is_open())
   {
      // Error: cannot open run details file for input
      cerr << ERR << "cannot open " << m_strDataPathName << " for input" << endl;
      return false;
   }

   int
      i = 0,
      nMult = 0;
      size_t nPos = 0;
   string strRec, strErr;
   
   while (getline(InStream, strRec))
   {
      // Trim off leading and trailing whitespace
      strRec = strTrimLeft(&strRec);
      strRec = strTrimRight(&strRec);

      // If it is a blank line or a comment then ignore it
      if ((! strRec.empty()) && (strRec[0] != QUOTE1) && (strRec[0] != QUOTE2))
      {
         // It isn't so increment counter
         i++;

         // Find the colon: note that lines MUST have a colon separating data from leading description portion
         nPos = strRec.find(':');
         if (nPos == string::npos)
         {
            // Error: badly formatted line (no colon)
            cerr << ERR << "badly formatted line (no ':') in " << m_strDataPathName << endl << strRec << endl;
            return false;
         }

         // Strip off leading portion (the bit up to and including the colon)
         string strRH = strRec.substr(nPos+1);

         // Remove leading whitespace after the colon
         strRH = strTrimLeft(&strRH);

         // Look for trailing comments, if found then terminate string at that point and trim off any trailing whitespace
         bool bFound = true;
         while (bFound)
         {
            bFound = false;

            nPos = strRH.rfind(QUOTE1);
            if (nPos != string::npos)
            {
               strRH = strRH.substr(0, nPos);
               bFound = true;
            }

            nPos = strRH.rfind(QUOTE2);
            if (nPos != string::npos)
            {
               strRH = strRH.substr(0, nPos);
               bFound = true;
            }

            // Trim trailing spaces
            strRH = strTrimRight(&strRH);
         }

#ifdef _WIN32
            // For Windows, make sure has backslashes, not Unix-style slashes
            strRH = pstrChangeToBackslash(&strRH);
#endif
         bool bFirst = true;
         int nRet = 0;
         double dMult = 0;
         string strTmp;

         switch (i)
         {
         // ---------------------------------------------- Run Information -----------------------------------------------------
         case 1:
            // Text output file names, don't change case
            if (strRH.empty())
               strErr = "output file names";
            else
            {
               m_strRunName = strRH;

               m_strOutFile = m_strOutPath;
               m_strOutFile.append(strRH);
               m_strOutFile.append(OUTEXT);

               m_strLogFile = m_strOutPath;
               m_strLogFile.append(strRH);
               m_strLogFile.append(LOGEXT);
            }
            break;

         case 2:
            // Duration of simulation (in hours, days, months, or years): sort out multiplier and user units, as used in the per-timestep output
            strRH = strToLower(&strRH);

            nRet = nDoSimulationTimeMultiplier(&strRH);
            if (nRet != RTN_OK)
            {
               strErr = "units for duration of simulation";
               break;
            }

            // And now calculate the duration of the simulation in hours: first find whitespace between the number and the unit
            nPos = strRH.rfind(SPACE);
            if (nPos == string::npos)
            {
               strErr = "format of duration simulation line";
               break;
            }

            // cut off rh bit of string
            strRH = strRH.substr(0, nPos+1);

            // remove trailing spaces
            strRH = strTrimRight(&strRH);

            // Calculate the duration of the simulation in hours
            m_dSimDuration = atof(strRH.c_str()) * m_dDurationUnitsMult;

            if (m_dSimDuration <= 0)
               strErr = "duration of simulation must be greater than zero";
            break;

         case 3:
            // Timestep of simulation (in hours or days)
            strRH = strToLower(&strRH);

            dMult = dGetTimeMultiplier(&strRH);
            if (static_cast<int>(dMult) == TIME_UNKNOWN)
            {
               strErr = "units for simulation timestep";
               break;
            }

            // we have the multiplier, now calculate the timestep in hours: look for the whitespace between the number and unit
            nPos = strRH.rfind(SPACE);
            if (nPos == string::npos)
            {
               strErr = "format of simulation timestep line";
               break;
            }

            // cut off rh bit of string
            strRH = strRH.substr(0, nPos+1);

            // remove trailing spaces
            strRH = strTrimRight(&strRH);

            m_dTimeStep = atof(strRH.c_str()) * dMult;         // in hours

            if (m_dTimeStep <= 0)
               strErr = "timestep of simulation must be greater than zero";

            if (m_dTimeStep >= m_dSimDuration)
               strErr = "timestep of simulation must be less than the duration of the simulation";

            break;

         case 4:
            // Save interval(s)
            strRH = strToLower(&strRH);

            // First get the multiplier
            dMult = dGetTimeMultiplier(&strRH);
            if (static_cast<int>(dMult) == TIME_UNKNOWN)
            {
               strErr = "units for save intervals";
               break;
            }

            // search for a space from the right-hand end of the string
            nPos = strRH.rfind(SPACE);
            if (nPos == string::npos)
            {
               strErr = "format of save times/intervals line";
               break;
            }

            // cut off rh bit of string
            strRH = strRH.substr(0, nPos+1);

            // remove trailing spaces
            strRH = strTrimRight(&strRH);

            // Now find out whether we're dealing with a single regular save interval or a vector of save times: again check for a space
            nPos = strRH.find(SPACE);
            if (nPos != string::npos)
            {
               // There's another space, so we must have more than one number
               m_bSaveRegular = false;

               strRH += SPACE;

               do
               {
                  // Put the number into the array
                  if (m_nUSave > static_cast<int>(SAVEMAX)-1)
                  {
                     strErr = "too many save intervals";
                     break;
                  }
                  m_dUSaveTime[m_nUSave++] = atof(strRH.substr(0, nPos).c_str()) * dMult;        // convert to hours

                  // Trim off the number and remove leading whitespace
                  strRH = strRH.substr(nPos, strRH.size()-nPos);
                  strRH = strTrimLeft(&strRH);

                  // Now look for another space
                  nPos = strRH.find(SPACE);
               }
//               while (nPos != string::npos);
               while (strRH.size() > 1);

               // Check that we have at least 2 numbers
               if (m_nUSave < 1)
               {
                  strErr = "must have at least two save times";
                  break;
               }

               if (m_dUSaveTime[0] < m_dTimeStep)
               {
                  strErr = "first save time cannot be less than timestep";
                  break;
               }

               // Put a dummy save interval as the last entry in the array: this is needed to stop problems at end of run
               m_dUSaveTime[m_nUSave] = m_dSimDuration + 1;
            }
            else
            {
               // There isn't a space, so we have only one number
               m_bSaveRegular = true;
               m_dRSaveInterval = atof(strRH.c_str()) * nMult;                   // convert to hours
               if (m_dRSaveTime <= m_dTimeStep)
                  strErr = "save interval cannot be less than timestep";
               else
                  // Set up for first save
                  m_dRSaveTime = m_dRSaveInterval;
            }
            break;

         case 5:
            // Random number seed(s)
            m_ulRandSeed[0] = atol(strRH.c_str());
            if (0 == m_ulRandSeed[0])
            {
               strErr = "random number seed";
               break;
            }

            // TODO rewrite this, similar to reading raster slice elevations
            // OK, so find out whether we're dealing with a single seed or more than one: check for a space
            nPos = strRH.find(SPACE);
            if (nPos != string::npos)
            {
               // There's a space, so we must have more than one number
               int n = 0;
               do
               {
                  // Trim off the part before the first space then remove leading whitespace
                  strRH = strRH.substr(nPos, strRH.size()-nPos);
                  strRH = strTrimLeft(&strRH);

                  // Put the number into the array
                  m_ulRandSeed[++n] = atol(strRH.c_str());

                  // Now look for another space
                  nPos = strRH.find(SPACE);
               }
               while ((n < NRNG) && (nPos != string::npos));
            }
            else
            {
               // Only one seed specified, so make all seeds the same
               for (int n = 1; n < NRNG; n++)
                  m_ulRandSeed[n] = m_ulRandSeed[n-1];
            }
            break;

         case 6:
            // Raster GIS files to output
            strRH = strToLower(&strRH);

            // These are always saved
            m_bSedimentTopSurfSave                 =
            m_bTopSurfSave                         =
            m_bSeaDepthSave                        =
            m_bWaveHeightSave                      =
            m_bPotentialPlatformErosionSave        =
            m_bActualPlatformErosionSave           =
            m_bTotalPotentialPlatformErosionSave   =
            m_bTotalActualPlatformErosionSave      =
            m_bPotentialBeachErosionSave           =
            m_bActualBeachErosionSave              =
            m_bTotalPotentialBeachErosionSave      =
            m_bTotalActualBeachErosionSave         =
            m_bLandformSave                        = true;

            // Now look for "all"
            if (strRH.find(RASTER_ALL_OUTPUT_CODE) != string::npos)
            {
               m_bAvgSeaDepthSave                  =
               m_bAvgWaveHeightSave                =
               m_bAvgWaveOrientationSave           =
               m_bBeachProtectionSave              =
               m_bBasementElevSave                 =
               m_bSuspSedSave                      =
               m_bAvgSuspSedSave                   =
               m_bFineUnconsSedSave                =
               m_bSandUnconsSedSave                =
               m_bCoarseUnconsSedSave              =
               m_bFineConsSedSave                  =
               m_bSandConsSedSave                  =
               m_bCoarseConsSedSave                =
               m_bRasterCoastlineSave              =
               m_bRasterNormalSave                 =
               m_bDistWeightSave                   =
               m_bActiveZoneSave                   =
               m_bCliffCollapseSave                =
               m_bTotCliffCollapseSave             =
               m_bCliffCollapseDepositionSave      =
               m_bTotCliffCollapseDepositionSave   =
               m_bRasterPolygonSave                =
               m_bPotentialPlatformErosionMaskSave =
               m_bSeaMaskSave                      =
               m_bBeachMaskSave                    =
               m_bInterventionClassSave            =
               m_bInterventionHeightSave           =
               m_bShadowZoneCodesSave              = 
               m_bDeepWaterWaveOrientationSave     = 
               m_bDeepWaterWaveHeightSave          = 
               m_bPolygonUnconsSedUpOrDownDrift    = 
               m_bPolygonUnconssedGainOrLoss       = true;
            }
            else
            {
               // These are only saved if the user specified the code
               if (strRH.find(RASTER_AVG_SEA_DEPTH_NAME) != string::npos)
               {
                  m_bAvgSeaDepthSave = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_AVG_SEA_DEPTH_NAME);
               }

               if (strRH.find(RASTER_AVG_WAVE_HEIGHT_NAME) != string::npos)
               {
                  m_bAvgWaveHeightSave = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_AVG_WAVE_HEIGHT_NAME);
               }

               if (strRH.find(RASTER_AVG_WAVE_ORIENTATION_NAME) != string::npos)
               {
                  m_bAvgWaveOrientationSave = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_AVG_WAVE_ORIENTATION_NAME);
               }
               
               if (strRH.find(RASTER_BEACH_PROTECTION_NAME) != string::npos)
               {
                  m_bBeachProtectionSave = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_BEACH_PROTECTION_NAME);
               }

               if (strRH.find(RASTER_BASEMENT_ELEVATION_NAME) != string::npos)
               {
                  m_bBasementElevSave = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_BASEMENT_ELEVATION_NAME);
               }

               if (strRH.find(RASTER_SUSP_SED_NAME) != string::npos)
               {
                  m_bSuspSedSave = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_SUSP_SED_NAME);
               }

               if (strRH.find(RASTER_AVG_SUSP_SED_NAME) != string::npos)
               {
                  m_bAvgSuspSedSave = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_AVG_SUSP_SED_NAME);
               }

               if (strRH.find(RASTER_FINE_UNCONS_NAME) != string::npos)
               {
                  m_bFineUnconsSedSave = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_FINE_UNCONS_NAME);
               }

               if (strRH.find(RASTER_SAND_UNCONS_NAME) != string::npos)
               {
                  m_bSandUnconsSedSave = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_SAND_UNCONS_NAME);
               }

               if (strRH.find(RASTER_COARSE_UNCONS_NAME) != string::npos)
               {
                  m_bCoarseUnconsSedSave = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_COARSE_UNCONS_NAME);
               }

               if (strRH.find(RASTER_FINE_CONS_NAME) != string::npos)
               {
                  m_bFineConsSedSave = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_FINE_CONS_NAME);
               }

               if (strRH.find(RASTER_SAND_CONS_NAME) != string::npos)
               {
                  m_bSandConsSedSave = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_SAND_CONS_NAME);
               }

               if (strRH.find(RASTER_COARSE_CONS_NAME) != string::npos)
               {
                  m_bCoarseConsSedSave = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_COARSE_CONS_NAME);
               }

               if (strRH.find(RASTER_COAST_NAME) != string::npos)
               {
                  m_bRasterCoastlineSave = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_COAST_NAME);
               }

               if (strRH.find(RASTER_COAST_NORMAL_NAME) != string::npos)
               {
                  m_bRasterNormalSave = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_COAST_NORMAL_NAME);
               }

               if (strRH.find(RASTER_ACTIVE_ZONE_NAME) != string::npos)
               {
                  m_bActiveZoneSave = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_ACTIVE_ZONE_NAME);
               }

               if (strRH.find(RASTER_CLIFF_COLLAPSE_NAME) != string::npos)
               {
                  m_bCliffCollapseSave = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_CLIFF_COLLAPSE_NAME);
               }

               if (strRH.find(RASTER_TOTAL_CLIFF_COLLAPSE_NAME) != string::npos)
               {
                  m_bTotCliffCollapseSave = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_TOTAL_CLIFF_COLLAPSE_NAME);
               }

               if (strRH.find(RASTER_CLIFF_COLLAPSE_DEPOSITION_NAME) != string::npos)
               {
                  m_bCliffCollapseDepositionSave = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_CLIFF_COLLAPSE_DEPOSITION_NAME);
               }

               if (strRH.find(RASTER_TOTAL_CLIFF_COLLAPSE_DEPOSITION_NAME) != string::npos)
               {
                  m_bTotCliffCollapseDepositionSave = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_TOTAL_CLIFF_COLLAPSE_DEPOSITION_NAME);
               }

               if (strRH.find(RASTER_POLYGON_NAME) != string::npos)
               {
                  m_bRasterPolygonSave = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_POLYGON_NAME);
               }

               if (strRH.find(RASTER_POTENTIAL_PLATFORM_EROSION_MASK_NAME) != string::npos)
               {
                  m_bPotentialPlatformErosionMaskSave = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_POTENTIAL_PLATFORM_EROSION_MASK_NAME);
               }

               if (strRH.find(RASTER_INUNDATION_MASK_NAME) != string::npos)
               {
                  m_bSeaMaskSave = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_INUNDATION_MASK_NAME);
               }

               if (strRH.find(RASTER_BEACH_MASK_NAME) != string::npos)
               {
                  m_bBeachMaskSave = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_BEACH_MASK_NAME);
               }

               if (strRH.find(RASTER_INTERVENTION_CLASS_NAME) != string::npos)
               {
                  m_bInterventionClassSave = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_INTERVENTION_CLASS_NAME);
               }

               if (strRH.find(RASTER_INTERVENTION_HEIGHT_NAME) != string::npos)
               {
                  m_bInterventionHeightSave = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_INTERVENTION_HEIGHT_NAME);
               }

               if (strRH.find(RASTER_SHADOW_ZONE_NAME) != string::npos)
               {
                  m_bShadowZoneCodesSave = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_SHADOW_ZONE_NAME);
               }

               if (strRH.find(RASTER_DEEP_WATER_WAVE_ORIENTATION_NAME) != string::npos)
               {
                  m_bDeepWaterWaveOrientationSave = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_DEEP_WATER_WAVE_ORIENTATION_NAME);
               }
               
               if (strRH.find(RASTER_DEEP_WATER_WAVE_HEIGHT_NAME) != string::npos)
               {
                  m_bDeepWaterWaveHeightSave = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_DEEP_WATER_WAVE_HEIGHT_NAME);
               }
               
               if (strRH.find(RASTER_POLYGON_UPDRIFT_OR_DOWNDRIFT_NAME) != string::npos)
               {
                  m_bPolygonUnconsSedUpOrDownDrift = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_POLYGON_UPDRIFT_OR_DOWNDRIFT_NAME);
               }
               
               if (strRH.find(RASTER_POLYGON_GAIN_OR_LOSS_NAME) != string::npos)
               {
                  m_bPolygonUnconssedGainOrLoss = true;
                  strRH = strRemoveSubstr(&strRH, &RASTER_POLYGON_GAIN_OR_LOSS_NAME);
               }
                  
               // Check to see if all codes have been removed
               strRH = strTrimLeft(&strRH);
               if (! strRH.empty())
                  strErr = "raster GIS output file list";
            }
            break;

         case 7:
            // Raster GIS output format (note must retain original case). Blank means use same format as input DEM file (if possible)
            m_strRasterGISOutFormat = strTrimLeft(&strRH);
            break;

         case 8:
            // If needed, also output GIS raster world file?
            strRH = strToLower(&strRH);

            m_bScaleRasterOutput = false;
            if (strRH.find("y") != string::npos)
               m_bScaleRasterOutput = true;
            break;

         case 9:
            // If needed, scale GIS raster output values?
            strRH = strToLower(&strRH);

            m_bWorldFile = false;
            if (strRH.find("y") != string::npos)
               m_bWorldFile = true;
            break;

         case 10:
            // Elevations for raster slice output, if desired
            if (! strRH.empty())
            {
               m_bSliceSave = true;

               // OK, so find out whether we're dealing with a single seed or more than one: check for a space
               nPos = strRH.find(SPACE);
               if (nPos != string::npos)
               {
                  // There's a space, so we must have more than one number
                  do
                  {
                     // Get LH bit
                     strTmp = strRH.substr(0, nPos);
                     m_VdSliceElev.push_back(atof(strTmp.c_str()));

                     // Get the RH bit
                     strRH = strRH.substr(nPos, strRH.size()-nPos);
                     strRH = strTrimLeft(&strRH);

                     // Now look for another space
                     nPos = strRH.find(SPACE);
                  }
                  while (nPos != string::npos);
               }
               // Read either the single number, or the left-over number
               m_VdSliceElev.push_back(atof(strTmp.c_str()));
            }
            break;

         case 11:
            // Vector GIS files to output
            strRH = strToLower(&strRH);

            // These are always saved
            m_bCoastSave =
            m_bWaveAngleAndHeightSave = true;

               // First look for "all"
            if (strRH.find(VECTOR_ALL_OUTPUT_CODE) != string::npos)
            {
               m_bNormalsSave                 =
               m_bInvalidNormalsSave          =
               m_bAvgWaveAngleAndHeightSave            =
               m_bWaveEnergySinceCollapseSave =
               m_bMeanWaveEnergySave          =
               m_bBreakingWaveHeightSave      =
               m_bCoastCurvatureSave          =
               m_bPolygonNodeSave             =
               m_bPolygonBoundarySave         =
               m_bCliffNotchSave              =
               m_bShadowBoundarySave          = 
               m_bShadowDowndriftBoundarySave =
               m_bDeepWaterWaveAngleAndHeightSave      = true;
            }
            else
            {
               // These are only saved if the user specified the code
               if (strRH.find(VECTOR_NORMALS_CODE) != string::npos)
               {
                  m_bNormalsSave = true;
                  strRH = strRemoveSubstr(&strRH, &VECTOR_NORMALS_CODE);
               }

               if (strRH.find(VECTOR_INVALID_NORMALS_CODE) != string::npos)
               {
                  m_bInvalidNormalsSave = true;
                  strRH = strRemoveSubstr(&strRH, &VECTOR_INVALID_NORMALS_CODE);
               }

               if (strRH.find(VECTOR_AVG_WAVE_ANGLE_AND_HEIGHT_CODE) != string::npos)
               {
                  m_bAvgWaveAngleAndHeightSave = true;
                  strRH = strRemoveSubstr(&strRH, &VECTOR_AVG_WAVE_ANGLE_AND_HEIGHT_CODE);
               }

               if (strRH.find(VECTOR_COAST_CURVATURE_CODE) != string::npos)
               {
                  m_bCoastCurvatureSave = true;
                  strRH = strRemoveSubstr(&strRH, &VECTOR_COAST_CURVATURE_CODE);
               }

               if (strRH.find(VECTOR_WAVE_ENERGY_SINCE_COLLAPSE_CODE) != string::npos)
               {
                  m_bWaveEnergySinceCollapseSave = true;
                  strRH = strRemoveSubstr(&strRH, &VECTOR_WAVE_ENERGY_SINCE_COLLAPSE_CODE);
               }

               if (strRH.find(VECTOR_MEAN_WAVE_ENERGY_CODE) != string::npos)
               {
                  m_bMeanWaveEnergySave = true;
                  strRH = strRemoveSubstr(&strRH, &VECTOR_MEAN_WAVE_ENERGY_CODE);
               }

               if (strRH.find(VECTOR_BREAKING_WAVE_HEIGHT_CODE) != string::npos)
               {
                  m_bBreakingWaveHeightSave = true;
                  strRH = strRemoveSubstr(&strRH, &VECTOR_BREAKING_WAVE_HEIGHT_CODE);
               }

               if (strRH.find(VECTOR_POLYGON_NODE_SAVE_CODE) != string::npos)
               {
                  m_bPolygonNodeSave = true;
                  strRH = strRemoveSubstr(&strRH, &VECTOR_POLYGON_NODE_SAVE_CODE);
               }

               if (strRH.find(VECTOR_POLYGON_BOUNDARY_SAVE_CODE) != string::npos)
               {
                  m_bPolygonBoundarySave = true;
                  strRH = strRemoveSubstr(&strRH, &VECTOR_POLYGON_BOUNDARY_SAVE_CODE);
               }

               if (strRH.find(VECTOR_CLIFF_NOTCH_SIZE_CODE) != string::npos)
               {
                  m_bCliffNotchSave = true;
                  strRH = strRemoveSubstr(&strRH, &VECTOR_CLIFF_NOTCH_SIZE_CODE);
               }

               if (strRH.find(VECTOR_SHADOW_BOUNDARY_CODE) != string::npos)
               {
                  m_bShadowBoundarySave = true;
                  strRH = strRemoveSubstr(&strRH, &VECTOR_SHADOW_BOUNDARY_CODE);
               }

               if (strRH.find(VECTOR_DOWNDRIFT_BOUNDARY_CODE) != string::npos)
               {
                  m_bShadowDowndriftBoundarySave = true;
                  strRH = strRemoveSubstr(&strRH, &VECTOR_DOWNDRIFT_BOUNDARY_CODE);
               }
               
               if (strRH.find(VECTOR_DEEP_WATER_WAVE_ANGLE_AND_HEIGHT_CODE) != string::npos)
               {
                  m_bDeepWaterWaveAngleAndHeightSave = true;
                  strRH = strRemoveSubstr(&strRH, &VECTOR_DEEP_WATER_WAVE_ANGLE_AND_HEIGHT_CODE);
               }
               
               // Check to see if all codes have been removed
               strRH = strTrimLeft(&strRH);
               if (! strRH.empty())
                  strErr = "vector GIS output file list";
            }
            break;

         case 12:
            // Vector GIS output format (note must retain original case)
            m_strVectorGISOutFormat = strRH;

            if (strRH.empty())
               strErr = "vector GIS output format";
            break;

         case 13:
            // Time series files to output
            strRH = strToLower(&strRH);

            // First check for "all"
            if (strRH.find(RASTER_ALL_OUTPUT_CODE) != string::npos)
            {
               m_bSeaAreaTS                  =
               m_bStillWaterLevelTS          =
               m_bActualPlatformErosionTS    =
               m_bCliffCollapseErosionTS     =
               m_bCliffCollapseDepositionTS  =
               m_bCliffCollapseNetTS         =
               m_bBeachErosionTS             =
               m_bBeachDepositionTS          =
               m_bBeachSedimentChangeNetTS   =
               m_bSuspSedTS                  = true;
            }
            else
            {
               if (strRH.find(TIME_SERIES_SEA_AREA_CODE) != string::npos)
               {
                  m_bSeaAreaTS = true;
                  strRH = strRemoveSubstr(&strRH, &TIME_SERIES_SEA_AREA_CODE);
               }

               if (strRH.find(TIME_SERIES_STILL_WATER_LEVEL_CODE) != string::npos)
               {
                  m_bStillWaterLevelTS = true;
                  strRH = strRemoveSubstr(&strRH, &TIME_SERIES_STILL_WATER_LEVEL_CODE);
               }

               if (strRH.find(TIME_SERIES_PLATFORM_EROSION_CODE) != string::npos)
               {
                  m_bActualPlatformErosionTS = true;
                  strRH = strRemoveSubstr(&strRH, &TIME_SERIES_PLATFORM_EROSION_CODE);
               }

               if (strRH.find(TIME_SERIES_CLIFF_COLLAPSE_EROSION_CODE) != string::npos)
               {
                  m_bCliffCollapseErosionTS = true;
                  strRH = strRemoveSubstr(&strRH, &TIME_SERIES_CLIFF_COLLAPSE_EROSION_CODE);
               }

               if (strRH.find(TIME_SERIES_CLIFF_COLLAPSE_DEPOSITION_CODE) != string::npos)
               {
                  m_bCliffCollapseDepositionTS = true;
                  strRH = strRemoveSubstr(&strRH, &TIME_SERIES_CLIFF_COLLAPSE_DEPOSITION_CODE);
               }

               if (strRH.find(TIME_SERIES_CLIFF_COLLAPSE_NET_CODE) != string::npos)
               {
                  m_bCliffCollapseNetTS = true;
                  strRH = strRemoveSubstr(&strRH, &TIME_SERIES_CLIFF_COLLAPSE_NET_CODE);
               }
               
               if (strRH.find(TIME_SERIES_SUSPENDED_SEDIMENT_CODE) != string::npos)
               {
                  m_bSuspSedTS = true;
                  strRH = strRemoveSubstr(&strRH, &TIME_SERIES_SUSPENDED_SEDIMENT_CODE);
               }

               // Check to see if all codes have been removed
               strRH = strTrimLeft(&strRH);
               if (! strRH.empty())
                  strErr = "time-series output file list";
            }
            break;

         case 14:
            // Vector coastline smoothing algorithm: 0 = none, 1 = running mean, 2 = Savitsky-Golay
            m_nCoastSmooth = atoi(strRH.c_str());
            if ((m_nCoastSmooth < SMOOTH_NONE) || (m_nCoastSmooth > SMOOTH_SAVITZKY_GOLAY))
               strErr = "coastline vector smoothing algorithm";
            break;

         case 15:
            // Size of coastline smoothing window: must be odd
            m_nCoastSmoothWindow = atoi(strRH.c_str());
            if ((m_nCoastSmoothWindow <= 0) || !(m_nCoastSmoothWindow % 2))
               strErr = "size of coastline vector smoothing window (must be > 0 and odd)";
            break;

         case 16:
            // Order of coastline profile smoothing polynomial for Savitsky-Golay: usually 2 or 4, max is 6
            m_nSavGolCoastPoly = atoi(strRH.c_str());
            if ((m_nSavGolCoastPoly <= 0) || (m_nSavGolCoastPoly > 6))
               strErr = "value of Savitsky-Golay polynomial for coastline smoothing (must be <= 6)";
            break;

         case 17:
            // Grid edge(s) to omit when searching for coastline [NSWE]
            strRH = strToLower(&strRH);

            if (strRH.find("n") != string::npos)
               m_bOmitSearchNorthEdge = true;
            if (strRH.find("s") != string::npos)
               m_bOmitSearchSouthEdge = true;
            if (strRH.find("w") != string::npos)
               m_bOmitSearchWestEdge = true;
            if (strRH.find("e") != string::npos)
               m_bOmitSearchEastEdge = true;
            break;

         case 18:
            // Profile slope running-mean smoothing window size: must be odd
            m_nProfileSmoothWindow = atoi(strRH.c_str());
            if ((m_nProfileSmoothWindow <= 0) || !(m_nProfileSmoothWindow % 2))
               strErr = "size of profile vector smoothing window (must be > 0 and odd)";
            break;

         case 19:
            // Max local slope (m/m)
            m_dProfileMaxSlope = atof(strRH.c_str());
            if (m_dProfileMaxSlope <= 0)
               strErr = "max local slope must be greater than zero";
            break;

         case 20:
            // Weighting for simple smoothing of sediment layers
            m_dSimpleSmoothWeight = atof(strRH.c_str());
            if (m_dSimpleSmoothWeight <= 0)
               strErr = "weighting for simple smoothing of sediment layers must be greater than zero";
            if (m_dSimpleSmoothWeight > 1)
               strErr = "weighting for simple smoothing of sediment layers must be less than or equal to one";
            break;

         case 21:
            // Vertical tolerance to include beach cells in smoothing
            m_dBeachSmoothingVertTolerance = atof(strRH.c_str());
            if (m_dBeachSmoothingVertTolerance < 0)
               strErr = "vertical tolerance for beach cells to be included in smoothing must be greater than or equal to zero";
            break;


         // ------------------------------------------------- Raster GIS layers ------------------------------------------------
         case 22:
            // Number of sediment layers
            m_nLayers = atoi(strRH.c_str());
            if (m_nLayers < 1)
            {
               strErr = "must be at least one sediment layer";
               break;
            }

            // OK we know the number of layers, so add elements to these vectors
            for (int i = 0; i < m_nLayers; i++)
            {
               m_VstrInitialFineUnconsSedimentFile.push_back("");
               m_VstrInitialSandUnconsSedimentFile.push_back("");
               m_VstrInitialCoarseUnconsSedimentFile.push_back("");
               m_VstrInitialFineConsSedimentFile.push_back("");
               m_VstrInitialSandConsSedimentFile.push_back("");
               m_VstrInitialCoarseConsSedimentFile.push_back("");
               m_VstrGDALIUFDriverCode.push_back("");
               m_VstrGDALIUFDriverDesc.push_back("");
               m_VstrGDALIUFProjection.push_back("");
               m_VstrGDALIUFDataType.push_back("");
               m_VstrGDALIUSDriverCode.push_back("");
               m_VstrGDALIUSDriverDesc.push_back("");
               m_VstrGDALIUSProjection.push_back("");
               m_VstrGDALIUSDataType.push_back("");
               m_VstrGDALIUCDriverCode.push_back("");
               m_VstrGDALIUCDriverDesc.push_back("");
               m_VstrGDALIUCProjection.push_back("");
               m_VstrGDALIUCDataType.push_back("");
               m_VstrGDALICFDriverCode.push_back("");
               m_VstrGDALICFDriverDesc.push_back("");
               m_VstrGDALICFProjection.push_back("");
               m_VstrGDALICFDataType.push_back("");
               m_VstrGDALICSDriverCode.push_back("");
               m_VstrGDALICSDriverDesc.push_back("");
               m_VstrGDALICSProjection.push_back("");
               m_VstrGDALICSDataType.push_back("");
               m_VstrGDALICCDriverCode.push_back("");
               m_VstrGDALICCDriverDesc.push_back("");
               m_VstrGDALICCProjection.push_back("");
               m_VstrGDALICCDataType.push_back("");
            }
            break;

         case 23:
            // Basement DEM file (can be blank)
            if (! strRH.empty())
            {
#ifdef _WIN32
               // For Windows, make sure has backslashes, not Unix-style slashes
               strRH = pstrChangeToBackslash(&strRH);
#endif
               // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
               if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == '~') || (strRH[1] == ':'))
                  // It has an absolute path, so use it 'as is'
                  m_strInitialBasementDEMFile = strRH;
               else
               {
                  // It has a relative path, so prepend the CoastalME dir
                  m_strInitialBasementDEMFile = m_strCMEDir;
                  m_strInitialBasementDEMFile.append(strRH);
               }
            }
            break;

         case 24:
            // Read 6 x sediment files for each layer
            for (int nLayer = 0; nLayer < m_nLayers; nLayer++)
            {
               for (int j = 1; j <= 6; j++)
               {
                  if (! bFirst)
                  {
                     do
                     {
                        if (! (getline(InStream, strRec)))
                        {
                           cerr << ERR << "premature end of file in " << m_strDataPathName << endl;
                           return false;
                        }

                        // Trim off leading and trailing whitespace
                        strRec = strTrimLeft(&strRec);
                        strRec = strTrimRight(&strRec);
                     }
                     // If it is a blank line or a comment then ignore it
                     while (strRec.empty() || (strRec[0] == QUOTE1) || (strRec[0] == QUOTE2));

                     // Not blank or a comment, so find the colon: lines MUST have a colon separating data from leading description portion
                     nPos = strRec.find(':');
                     if (nPos == string::npos)
                     {
                        // Error: badly formatted line (no colon)
                        cerr << ERR << "badly formatted line (no ':') in " << m_strDataPathName << endl << strRec << endl;
                        return false;
                     }

                     // Strip off leading portion (the bit up to and including the colon)
                     strRH = strRec.substr(nPos+1);

                     // Remove leading whitespace after the colon
                     strRH = strTrimLeft(&strRH);

                     // Look for a trailing comment, if found then terminate string at that point and trim off any trailing whitespace
                     nPos = strRH.rfind(QUOTE1);
                     if (nPos != string::npos)
                        strRH = strRH.substr(0, nPos);

                     nPos = strRH.rfind(QUOTE2);
                     if (nPos != string::npos)
                        strRH = strRH.substr(0, nPos);

                     // Trim trailing spaces
                     strRH = strTrimRight(&strRH);

#ifdef _WIN32
                     // For Windows, make sure has backslashes, not Unix-style slashes
                     strRH = pstrChangeToBackslash(&strRH);
#endif
                  }

                  bFirst = false;

                  switch (j)
                  {
                  case 1:
                     // Initial fine unconsolidated sediment depth GIS file (can be blank)
                     if (! strRH.empty())
                     {
#ifdef _WIN32
                        // For Windows, make sure has backslashes, not Unix-style slashes
                        strRH = pstrChangeToBackslash(&strRH);
#endif

                        // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
                        if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == '~') || (strRH[1] == ':'))
                        {
                           // It has an absolute path, so use it 'as is'
                           m_VstrInitialFineUnconsSedimentFile[nLayer] = strRH;
                        }
                        else
                        {
                           // It has a relative path, so prepend the CoastalME dir
                           m_VstrInitialFineUnconsSedimentFile[nLayer] = m_strCMEDir;
                           m_VstrInitialFineUnconsSedimentFile[nLayer].append(strRH);
                        }
                     }
                     break;

                  case 2:
                     // Initial sand unconsolidated sediment depth GIS file (can be blank)
                     if (! strRH.empty())
                     {
#ifdef _WIN32
                        // For Windows, make sure has backslashes, not Unix-style slashes
                        strRH = pstrChangeToBackslash(&strRH);
#endif
                        // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
                        if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == '~') || (strRH[1] == ':'))
                        {
                           // It has an absolute path, so use it 'as is'
                           m_VstrInitialSandUnconsSedimentFile[nLayer] = strRH;
                        }
                        else
                        {
                           // It has a relative path, so prepend the CoastalME dir
                           m_VstrInitialSandUnconsSedimentFile[nLayer] = m_strCMEDir;
                           m_VstrInitialSandUnconsSedimentFile[nLayer].append(strRH);
                        }
                     }
                     break;

                  case 3:
                     // Initial coarse unconsolidated sediment depth GIS file (can be blank)
                     if (! strRH.empty())
                     {
#ifdef _WIN32
                        // For Windows, make sure has backslashes, not Unix-style slashes
                        strRH = pstrChangeToBackslash(&strRH);
#endif
                        // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
                        if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == '~') || (strRH[1] == ':'))
                        {
                           // It has an absolute path, so use it 'as is'
                           m_VstrInitialCoarseUnconsSedimentFile[nLayer] = strRH;
                        }
                        else
                        {
                           // It has a relative path, so prepend the CoastalME dir
                           m_VstrInitialCoarseUnconsSedimentFile[nLayer] = m_strCMEDir;
                           m_VstrInitialCoarseUnconsSedimentFile[nLayer].append(strRH);
                        }
                     }
                     break;

                  case 4:
                     // Initial fine consolidated sediment depth GIS file (can be blank)
                     if (! strRH.empty())
                     {
#ifdef _WIN32
                        // For Windows, make sure has backslashes, not Unix-style slashes
                        strRH = pstrChangeToBackslash(&strRH);
#endif
                        // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
                        if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == '~') || (strRH[1] == ':'))
                        {
                           // It has an absolute path, so use it 'as is'
                           m_VstrInitialFineConsSedimentFile[nLayer] = strRH;
                        }
                        else
                        {
                           // It has a relative path, so prepend the CoastalME dir
                           m_VstrInitialFineConsSedimentFile[nLayer] = m_strCMEDir;
                           m_VstrInitialFineConsSedimentFile[nLayer].append(strRH);
                        }
                     }
                     break;

                  case 5:
                     // Initial sand consolidated sediment depth GIS file (can be blank)
                     if (! strRH.empty())
                     {
#ifdef _WIN32
                        // For Windows, make sure has backslashes, not Unix-style slashes
                        strRH = pstrChangeToBackslash(&strRH);
#endif
                        // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
                        if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == '~') || (strRH[1] == ':'))
                        {
                           // It has an absolute path, so use it 'as is'
                           m_VstrInitialSandConsSedimentFile[nLayer] = strRH;
                        }
                        else
                        {
                           // It has a relative path, so prepend the CoastalME dir
                           m_VstrInitialSandConsSedimentFile[nLayer] = m_strCMEDir;
                           m_VstrInitialSandConsSedimentFile[nLayer].append(strRH);
                        }
                     }
                     break;

                  case 6:
                     // Initial coarse consolidated sediment depth GIS file (can be blank)
                     if (! strRH.empty())
                     {
#ifdef _WIN32
                        // For Windows, make sure has backslashes, not Unix-style slashes
                        strRH = pstrChangeToBackslash(&strRH);
#endif
                        // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
                        if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == '~') || (strRH[1] == ':'))
                        {
                           // It has an absolute path, so use it 'as is'
                           m_VstrInitialCoarseConsSedimentFile[nLayer] = strRH;
                        }
                        else
                        {
                           // It has a relative path, so prepend the CoastalME dir
                           m_VstrInitialCoarseConsSedimentFile[nLayer] = m_strCMEDir;
                           m_VstrInitialCoarseConsSedimentFile[nLayer].append(strRH);
                        }
                     }
                     break;
                  }

                  // Did an error occur?
                  if (! strErr.empty())
                  {
                     // Error in input to run details file
                     cerr << ERR << "reading " << strErr << " in " << m_strDataPathName << endl << "'" << strRec << "'" << endl;
                     InStream.close();
                     return false;
                  }
               }
            }
            break;

         case 25:
            // Initial suspended sediment depth GIS file (can be blank)
            if (! strRH.empty())
            {
#ifdef _WIN32
               // For Windows, make sure has backslashes, not Unix-style slashes
               strRH = pstrChangeToBackslash(&strRH);
#endif
               // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
               if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == '~') || (strRH[1] == ':'))
               {
                  // It has an absolute path, so use it 'as is'
                  m_strInitialSuspSedimentFile = strRH;
               }
               else
               {
                  // It has a relative path, so prepend the CoastalME dir
                  m_strInitialSuspSedimentFile = m_strCMEDir;
                  m_strInitialSuspSedimentFile.append(strRH);
               }
            }
            break;

         case 26:
            // Initial Landform class GIS file (can be blank)
            if (! strRH.empty())
            {
#ifdef _WIN32
               // For Windows, make sure has backslashes, not Unix-style slashes
               strRH = pstrChangeToBackslash(&strRH);
#endif
               // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
               if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == '~') || (strRH[1] == ':'))
               {
                  // It has an absolute path, so use it 'as is'
                  m_strInitialLandformFile = strRH;
               }
               else
               {
                  // It has a relative path, so prepend the CoastalME dir
                  m_strInitialLandformFile = m_strCMEDir;
                  m_strInitialLandformFile.append(strRH);
               }
            }
            break;

         case 27:
            // Initial Intervention class GIS file (can be blank: if so then intervention height file must also be blank)
            if (! strRH.empty())
            {
#ifdef _WIN32
               // For Windows, make sure has backslashes, not Unix-style slashes
               strRH = pstrChangeToBackslash(&strRH);
#endif
               // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
               if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == '~') || (strRH[1] == ':'))
               {
                  // It has an absolute path, so use it 'as is'
                  m_strInterventionClassFile = strRH;
               }
               else
               {
                  // It has a relative path, so prepend the CoastalME dir
                  m_strInterventionClassFile = m_strCMEDir;
                  m_strInterventionClassFile.append(strRH);
               }
            }
            break;

         case 28:
            // Initial Intervention height GIS file (can be blank: if so then intervention class file must also be blank)
            if (strRH.empty())
            {
               if (! m_strInterventionClassFile.empty())
                  strErr = "must specify both intervention class and intervention height files";
               break;
            }
            else
            {
               if (m_strInterventionClassFile.empty())
               {
                  strErr = "must specify both intervention class and intervention height files";
                  break;
               }

               #ifdef _WIN32
               // For Windows, make sure has backslashes, not Unix-style slashes
               strRH = pstrChangeToBackslash(&strRH);
               #endif

               // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
               if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == '~') || (strRH[1] == ':'))
               {
                  // It has an absolute path, so use it 'as is'
                  m_strInterventionHeightFile = strRH;
               }
               else
               {
                  // It has a relative path, so prepend the CoastalME dir
                  m_strInterventionHeightFile = m_strCMEDir;
                  m_strInterventionHeightFile.append(strRH);
               }
            }
            break;

            // ---------------------------------------------------- Hydrology data ------------------------------------------------
         case 29:
            // Wave propagation model [0 = COVE, 1 = CShore]
            m_nWavePropagationModel = atoi(strRH.c_str());
            if ((m_nWavePropagationModel != MODEL_COVE) && (m_nWavePropagationModel != MODEL_CSHORE))
               strErr = "switch for wave propagation model must be 0 or 1";
            break;

         case 30:
            // Density of sea water (kg/m3)
            m_dSeaWaterDensity = atof(strRH.c_str());
            if (m_dSeaWaterDensity <= 0)
               strErr = "sea water density must be greater than zero";
            break;

         case 31:
            // Initial still water level (m), or per-timestep SWL file TODO
            m_dOrigSWL = atof(strRH.c_str());
            break;

         case 32:
            // Final still water level (m) [blank = same as initial SWL]
            if (strRH.empty())
               m_dFinalSWL = m_dOrigSWL;
            else
               m_dFinalSWL = atof(strRH.c_str());
            break;

         case 33:
            // Wave period (sec)
            m_dWavePeriod = atof(strRH.c_str());
            if (m_dWavePeriod <= 0)
               strErr = "wave period must be greater than zero";
            break;

         case 34:
            // Deep water wave height (m) or a file of point vectors giving deep water wave height (m) and orientation (for units, see below)
            // TODO need option for multiple files
            if (isdigit(strRH.at(0)))    // if start with a number is a single value, if a file, filename must NOT start with number
            {
               // Just one value of wave height for all deep water cells
               m_bSingleDeepWaterWaveValues = true;

               m_dAllCellsDeepWaterWaveHeight = atof(strRH.c_str());
               if (m_dAllCellsDeepWaterWaveHeight <= 0)
                  strErr = "deep water wave height must be greater than zero";
            }
            else
            {
               // We are reading deep water wave height and deep water wave orientation from a file
               m_bSingleDeepWaterWaveValues = false;
               
               // TODO for the moment, calculate these values in only at the first timestep
               m_VulDeepWaterWaveValuesAtTimestep.push_back(1);

               if (strRH.empty())
               {
                  strErr = "deep water wave height must be either a number or a filename";
                  break;
               }

               #ifdef _WIN32
               // For Windows, make sure has backslashes, not Unix-style slashes
               strRH = pstrChangeToBackslash(&strRH);
               #endif

               // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
               if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == '~') || (strRH[1] == ':'))
               {
                  // It has an absolute path, so use it 'as is'
                  m_strDeepWaterWaveValuesFile = strRH;
               }
               else
               {
                  // It has a relative path, so prepend the CoastalME dir
                  m_strDeepWaterWaveValuesFile = m_strCMEDir;
                  m_strDeepWaterWaveValuesFile.append(strRH);
               }
            }
            break;

         case 35:
            // Deep water wave orientation in input CRS: this is the oceanographic convention i.e. direction TOWARDS which the waves move (in degrees clockwise from north)
            if (m_bSingleDeepWaterWaveValues)
            {
               // Only read this if we also have just a single value of wave height for all deep water cells
               m_dAllCellsDeepWaterWaveOrientation = atof(strRH.c_str());
               if (m_dAllCellsDeepWaterWaveOrientation < 0)
                  strErr = "deep water wave orientation must be zero degrees or more";
               else if (m_dAllCellsDeepWaterWaveOrientation >= 360)
                  strErr = "deep water wave orientation must be less than 360 degrees";
            }
            break;

//          case 35:
//             // Tide data file (can be blank)
//             if (! strRH.empty())
//             {
// #ifdef _WIN32
//                // For Windows, make sure has backslashes, not Unix-style slashes
//                strRH = pstrChangeToBackslash(&strRH);
// #endif
//                // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
//                if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == '~') || (strRH[1] == ':'))
//                   // It has an absolute path, so use it 'as is'
//                   m_strTideDataFile = strRH;
//                else
//                {
//                   // It has a relative path, so prepend the CoastalME dir
//                   m_strTideDataFile = m_strCMEDir;
//                   m_strTideDataFile.append(strRH);
//                }
//             }
//             break;

	 case 36:
            // Breaking wave height-to-depth ratio 
            m_dBreakingWaveHeightDeptRatio = atof(strRH.c_str());
            if (m_dBreakingWaveHeightDeptRatio <= 0)
               strErr = "breaking wave height to depth ratio must be greater than zero";
            break;
         // ----------------------------------------------------- Sediment data ------------------------------------------------
         case 37:
            // Simulate coast platform erosion?
            strRH = strToLower(&strRH);

            m_bDoCoastPlatformErosion = false;
            if (strRH.find("y") != string::npos)
               m_bDoCoastPlatformErosion = true;
            break;

         case 38:
            // R (resistance to erosion) values along profile, see Walkden & Hall, 2011
            m_dR = atof(strRH.c_str());
            if (m_dR <= 0)
               strErr = "R values must be greater than zero";
            break;

         case 39:
            // Beach sediment transport at grid edges [0 = closed, 1 = open, 2 = re-circulate]
            m_nUnconsSedimentHandlingAtGridEdges = atoi(strRH.c_str());
            if ((m_nUnconsSedimentHandlingAtGridEdges < GRID_EDGE_CLOSED) || (m_nUnconsSedimentHandlingAtGridEdges > GRID_EDGE_RECIRCULATE))
               strErr = "switch for handling of beach sediment at grid edges must be 0, 1, or 2";
            break;

         case 40:
            // Beach erosion/deposition equation [0 = CERC, 1 = Kamphuis]
            m_nBeachErosionDepositionEquation = atoi(strRH.c_str());
            if ((m_nBeachErosionDepositionEquation != EQUATION_CERC) && (m_nBeachErosionDepositionEquation != EQUATION_KAMPHUIS))
               strErr = "switch for beach erosion/deposition equation must be 0 or 1";
            break;

         case 41:
            // Median size of fine sediment (mm)      [0 = default, only for Kamphuis eqn]
            m_dD50Fine = atof(strRH.c_str());
            if (m_dD50Fine < 0)
               strErr = "median particle size of fine sediment must be greater than zero";
            else if (m_dD50Fine == 0)
               // Use default value
               m_dD50Fine = 0.0625;
            break;

         case 42:
            // Median size of sand sediment (mm)      [0 = default, only for Kamphuis eqn]
            m_dD50Sand = atof(strRH.c_str());
            if (m_dD50Sand < 0)
               strErr = "median particle size of sand sediment must be greater than zero";
            else if (m_dD50Sand == 0)
               // Use default value
               m_dD50Sand = 0.42;
            break;

         case 43:
            // Median size of coarse sediment (mm)      [0 = default, only for Kamphuis eqn]
            m_dD50Coarse = atof(strRH.c_str());
            if (m_dD50Coarse < 0)
               strErr = "median particle size of coarse sediment must be greater than zero";
            else if (m_dD50Coarse == 0)
               // Use default value
               m_dD50Coarse = 19.0;
            break;

         case 44:
            // Density of beach sediment (kg/m3)
            m_dBeachSedimentDensity = atof(strRH.c_str());
            if (m_dBeachSedimentDensity <= 0)
               strErr = "density of beach sediment must be greater than zero";
            break;

         case 45:
            // Beach sediment porosity
            m_dBeachSedimentPorosity = atof(strRH.c_str());
            if (m_dBeachSedimentPorosity <= 0)
               strErr = "porosity of beach sediment must be greater than zero";
            break;

         case 46:
            // Erodibility of fine-sized sediment
            m_dFineErodibility = atof(strRH.c_str());
            if (m_dFineErodibility < 0)
               strErr = "cannot have negative erodibility values";
            break;

         case 47:
            // Erodibility of sand-sized sediment
            m_dSandErodibility = atof(strRH.c_str());
            if (m_dSandErodibility < 0)
               strErr = "cannot have negative erodibility values";
            break;

         case 48:
            // Erodibility of coarse-sized sediment
            m_dCoarseErodibility = atof(strRH.c_str());
            if (m_dCoarseErodibility < 0)
            {
               strErr = "cannot have negative erodibility values";
               break;
            }

            if ((m_dFineErodibility + m_dSandErodibility + m_dCoarseErodibility) <= 0)
               strErr = "must have at least one non-zero erodibility value";
            break;

/*
     ; LAYER 1, LOWEST
   dErodFine                                                      : 0.8
   dErodSand                                                      : 0.7
   dErodCoarse                                                    : 0.9
   ; LAYER 2
   dErodFine                                                      : 1.0
   dErodSand                                                      : 0.9
   dErodCoarse                                                    : 0.6
   ; LAYER 3, HIGHEST
   dErodFine                                                      : 0.95
   dErodSand                                                      : 0.4
   dErodCoarse                                                    : 0.3
*/

         case 49:
            // Transport parameter KLS in CERC equation
            m_dKLS = atof(strRH.c_str());
            if (m_dKLS <= 0)
               strErr = "transport parameter KLS for CERC equation must be greater than zero";
            break;

         case 50:
            // Transport parameter for Kamphuis equation
            m_dKamphuis = atof(strRH.c_str());
            if (m_dKamphuis <= 0)
               strErr = "transport parameter for Kamphuis equation must be greater than zero";
            break;

         case 51:
            // Dean profile start, height above still water level (m)
            m_dDeanProfileStartAboveSWL = atof(strRH.c_str());
            if (m_dDeanProfileStartAboveSWL < 0)
               strErr = "height above SWL for Dean profile start must be greater than or equal to zero";
            break;

         // ------------------------------------------------ Cliff collapse data -----------------------------------------------
         case 52:
            // Simulate cliff collapse?
            strRH = strToLower(&strRH);

            m_bDoCliffCollapse = false;
            if (strRH.find("y") != string::npos)
               m_bDoCliffCollapse = true;
            break;

         case 53:
            // Cliff resistance to erosion
            m_dCliffErosionResistance = atof(strRH.c_str());
            if (m_dCliffErosionResistance <= 0)
               strErr = "cliff resistance to erosion must be greater than 0";
            break;

         case 54:
            // Notch overhang at collapse (m)
            m_dNotchOverhangAtCollapse = atof(strRH.c_str());
            if (m_dNotchOverhangAtCollapse <= 0)
               strErr = "cliff notch overhang at collapse must be greater than 0";
            break;

         case 55:
            // Notch base below still water level (m)
            m_dNotchBaseBelowSWL = atof(strRH.c_str());
            if (m_dNotchBaseBelowSWL < 0)
               strErr = "cliff notch base below still water level must be greater than 0";
            break;

         case 56:
            // Scale parameter A for cliff deposition (m^(1/3)) [0 = auto]
            m_dCliffDepositionA = atof(strRH.c_str());
            if (m_dCliffDepositionA < 0)
               strErr = "scale parameter A for cliff deposition must be 0 [= auto] or greater";
            break;

         case 57:
            // Planview width of cliff deposition talus (in cells) [must be odd]
            m_nCliffDepositionPlanviewWidth = atoi(strRH.c_str());
            if ((m_nCliffDepositionPlanviewWidth % 2) == 0)
               strErr = "planview width of cliff deposition must be odd";
            if (m_nCliffDepositionPlanviewWidth <= 0)
               strErr = "planview width of cliff deposition must be greater than 0";
            break;

         case 58:
            // Planview length of cliff deposition talus (m)
            m_dCliffDepositionPlanviewLength = atof(strRH.c_str());
            if (m_dCliffDepositionPlanviewLength <= 0)
               strErr = "planview length of cliff deposition must be greater than 0";
            break;

         case 59:
            // Height of landward end of talus, as a fraction of cliff elevation
            m_dCliffDepositionHeightFrac = atof(strRH.c_str());
            if (m_dCliffDepositionHeightFrac < 0)
               strErr = "height of cliff collapse (as a fraction of cliff elevation) must be 0 or greater";
            break;

         // ------------------------------------------------------ Other data --------------------------------------------------
         case 60:
            // Gravitational acceleration (m2/s)
            m_dG = atof(strRH.c_str());
            if (m_dG <= 0)
               strErr = "gravitational acceleration must be greater than zero";
            break;

         case 61:
            // Spacing of coastline normals (m)
            m_dCoastNormalAvgSpacing = atof(strRH.c_str());
            if (m_dCoastNormalAvgSpacing == 0)
               m_nCoastNormalAvgSpacing = MIN_PROFILE_SPACING;    // In cells, we will set m_dCoastNormalAvgSpacing later when we know m_dCellSide
            else if (m_dCoastNormalAvgSpacing < 0)
               strErr = "spacing of coastline normals must be greater than zero";
            break;

         case 62:
            // Random factor for spacing of normals  [0 to 1, 0 = deterministic]
            m_dCoastNormalRandSpaceFact = atof(strRH.c_str());
            if (m_dCoastNormalRandSpaceFact < 0)
               strErr = "spacing of coastline normals must be greater than zero";
            else if (m_dCoastNormalRandSpaceFact > 1)
               strErr = "spacing of coastline normals must be less than one";
            break;

         case 63:
            // Length of coastline normals (m)
            m_dCoastNormalLength = atof(strRH.c_str());
            if (m_dCoastNormalLength <= 0)
               strErr = "length of coastline normals must be greater than zero";
            break;

         case 64:
            // Maximum number of 'cape' normals
            m_nNaturalCapeNormals = atoi(strRH.c_str());
            if (m_nNaturalCapeNormals < 0)
               strErr = "number of 'cape' normals must be zero or greater";
            break;

         case 65:
            // Start depth for wave calcs (ratio to deep water wave height)     : 5
            m_dWaveDepthRatioForWaveCalcs = atof(strRH.c_str());

            if (m_dWaveDepthRatioForWaveCalcs <= 0)
               strErr = "start depth for wave calcs must be greater than zero";
            break;

         // ----------------------------------------------------- Testing only -------------------------------------------------
         case 66:
            // Output profile data?
            strRH = strToLower(&strRH);

            m_bOutputProfileData = false;
            if (strRH.find("y") != string::npos)
               m_bOutputProfileData = true;

            // TODO WHAT ABOUT RANDOMNESS OF PROFILE SPACING now that profile location is determined by curvature??

            break;

         case 67:
            // Numbers of profiles to be saved
            if (m_bOutputProfileData)
            {
               vector<string> strTmp = strSplit(&strRH, SPACE);
               for (unsigned int j = 0; j < strTmp.size(); j++)
               {
                  strTmp[j] = strTrim(&strTmp[j]);
                  int nTmp = atoi(strTmp[j].c_str());
                  if (nTmp < 0)
                  {
                     strErr = "Profile number for saving must be zero or greater";
                     break;
                  }

                  m_VnProfileToSave.push_back(nTmp);
               }
            }
            break;

         case 68:
           // Timesteps to save profile for output
            if (m_bOutputProfileData)
            {
               vector<string> strTmp = strSplit(&strRH, SPACE);
               for (unsigned int j = 0; j < strTmp.size(); j++)
               {
                  strTmp[j] = strTrim(&strTmp[j]);
                  unsigned long ulTmp = atol(strTmp[j].c_str());
                  if (ulTmp < 1)
                  {
                     strErr = "Timestep for profile saves must be one or greater";
                     break;
                  }

                  m_VulProfileTimestep.push_back(ulTmp);
               }
            }
            break;

         case 69:
            // Output parallel profile data?
            strRH = strToLower(&strRH);

            m_bOutputParallelProfileData = false;
            if (strRH.find("y") != string::npos)
               m_bOutputParallelProfileData = true;
            break;

         case 70:
            // Output erosion potential look-up data?
            strRH = strToLower(&strRH);

            m_bOutputLookUpData = false;
            if (strRH.find("y") != string::npos)
               m_bOutputLookUpData = true;
            break;

         case 71:
            // Erode shore platform in alternate direction each timestep?
            strRH = strToLower(&strRH);

            m_bErodeShorePlatformAlternateDirection = false;
            if (strRH.find("y") != string::npos)
            m_bErodeShorePlatformAlternateDirection = true;
               break;
         }

         // Did an error occur?
         if (! strErr.empty())
         {
            // Error in input to run details file
            cerr << endl << ERR << strErr << ".\nPlease edit " << m_strDataPathName << " and change this line:" << endl << "'" << strRec << "'" << endl << endl;
            InStream.close();
            return false;
         }
      }
   }

   // Close file
   InStream.close();

   // Finally, need to check that we have at least one raster file, so that we know the grid size and units (and preferably also the projection)
   bool bNoRasterFiles = true;
   if ((! m_strInitialBasementDEMFile.empty()) || (! m_strInitialSuspSedimentFile.empty()) || (! m_strInitialLandformFile.empty()) || (! m_strInterventionHeightFile.empty()))
      bNoRasterFiles = false;

   for (int i = 0; i < m_nLayers; i++)
   {
      if ((! m_VstrInitialFineUnconsSedimentFile[i].empty()) || (! m_VstrInitialSandUnconsSedimentFile[i].empty()) || (! m_VstrInitialCoarseUnconsSedimentFile[i].empty()) || (! m_VstrInitialFineConsSedimentFile[i].empty()) || (! m_VstrInitialSandConsSedimentFile[i].empty()) || (! m_VstrInitialCoarseConsSedimentFile[i].empty()))
         bNoRasterFiles = false;
   }

   if (bNoRasterFiles)
   {
      // No raster files
      cerr << ERR << "at least one raster GIS file is needed" << endl;
      return false;
   }

   return true;
}

/*==============================================================================================================================

 Reads the tide data

==============================================================================================================================*/
// int CSimulation::nReadTideData()
// {
//    // Create an ifstream object
//    ifstream InStream;
//
//    // Try to open the file for input
//    InStream.open(m_strTideDataFile.c_str(), ios::in);
//
//    // Did it open OK?
//    if (! InStream.is_open())
//    {
//       // Error: cannot open tide data file for input
//       cerr << ERR << "cannot open " << m_strTideDataFile << " for input" << endl;
//       return RTN_ERR_TIDEDATAFILE;
//    }
//
//    // Opened OK
//    string strRec;
//
//    // Now read the data from the file
//    while (getline(InStream, strRec))
//    {
//       // Trim off leading and trailing whitespace
//       strRec = strTrimLeft(&strRec);
//       strRec = strTrimRight(&strRec);
//
//       // If it is a blank line or a comment then ignore it
//       if ((strRec.empty()) || (strRec[0] == QUOTE1) || (strRec[0] == QUOTE2))
//          continue;
//
//       // Convert to double then append the value to the vector TODO check for floating point validity
//       m_VdTideData.push_back(strtod(strRec.c_str(), NULL));
//    }
//
//    // Close file
//    InStream.close();
//
//    return RTN_OK;
// }


/*==============================================================================================================================

 Reads the shape of the erosion potential distribution (see shape function in Walkden & Hall, 2005)

==============================================================================================================================*/
int CSimulation::nReadShapeFunction()
{
   // Sort out the path and filename
   m_strShapeFunctionFile = m_strCMEDir;
   m_strShapeFunctionFile.append(SCAPEDIR);
   m_strShapeFunctionFile.append(SCAPESHAPEFUNCTIONFILE);

   // Create an ifstream object
   ifstream InStream;

   // Try to open the file for input
   InStream.open(m_strShapeFunctionFile.c_str(), ios::in);

   // Did it open OK?
   if (! InStream.is_open())
   {
      // Error: cannot open shape function file for input
      cerr << ERR << "cannot open " << m_strShapeFunctionFile << " for input" << endl;
      return RTN_ERR_SCAPESHAPEFUNCTIONFILE;
   }

   // Opened OK
   int nExpected = 0, nRead = 0;
   string strRec;

   // Read in the number of data lines expected
   InStream >> nExpected;

   // Set up the vectors to hold the input data
   vector<double>
      VdDepthOverDB,
      VdErosionPotential,
      VdErosionPotentialFirstDeriv;

   // Now read the rest of the data from the file to get the Erosion Potential Shape function
   while (getline(InStream, strRec))
   {
      // Trim off leading and trailing whitespace
      strRec = strTrimLeft(&strRec);
      strRec = strTrimRight(&strRec);

      // If it is a blank line or a comment then ignore it
      if ((strRec.empty()) || (strRec[0] == QUOTE1) || (strRec[0] == QUOTE2))
         continue;

      // It isn't so increment counter
      nRead++;

      // Split the string, and remove whitespace
      vector<string> strTmp = strSplit(&strRec, SPACE);
      for (unsigned int i = 0; i < strTmp.size(); i++)
         strTmp[i] = strTrim(&strTmp[i]);

      // Convert to doubles then append the values to the vectors TODO check for floating point validity
      VdDepthOverDB.push_back(strtod(strTmp[0].c_str(), NULL));
      VdErosionPotential.push_back(strtod(strTmp[1].c_str(), NULL));
      VdErosionPotentialFirstDeriv.push_back(strtod(strTmp[2].c_str(), NULL));
   }

   // Close file
   InStream.close();

   // Did we read in what we expected?
   if (nExpected != nRead)
   {
      cout << ERR << "read in " << nRead << " lines from " << m_strShapeFunctionFile << " but " << nExpected << " lines expected" << endl;
      return RTN_ERR_SCAPESHAPEFUNCTIONFILE;
   }

   // OK, now use this data to create a look-up table to be used for the rest of the simulation
   if (! bCreateErosionPotentialLookUp(&VdDepthOverDB, &VdErosionPotential, &VdErosionPotentialFirstDeriv))
   {
      cout << ERR << " in " << m_strShapeFunctionFile << ", erosion potential function is unbounded for high values of depth over DB" << endl;
      return RTN_ERR_SCAPESHAPEFUNCTIONFILE;
   }

   return RTN_OK;
}
