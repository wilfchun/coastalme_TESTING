/*!
 *
 * \file write_output.cpp
 * \brief Writes non-GIS output files
 * \details TODO A more detailed description of this routine.
 * \author David Favis-Mortlock
 * \author Andres Payo

 * \date 2020
 * \copyright GNU General Public License
 *
 */

/*==============================================================================================================================

 This file is part of CoastalME, the Coastal Modelling Environment.

 CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public  License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

==============================================================================================================================*/
#include <assert.h>

#include <ctime>
using std::localtime;

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::showpos;
using std::noshowpos;

#include <iomanip>
using std::setiosflags;
using std::resetiosflags;
using std::setprecision;
using std::setw;
using std::setfill;
using std::put_time;

#include <sstream>
using std::stringstream;

#include "cme.h"
#include "simulation.h"
#include "interpolate.h"


/*==============================================================================================================================

 Writes beginning-of-run information to Out and Log files

==============================================================================================================================*/
void CSimulation::WriteStartRunDetails(void)
{
   // Set the Out file output format to fixed point
   OutStream << std::fixed;

   // Start outputting stuff
   OutStream << PROGRAM_NAME << " for " << PLATFORM << " " << strGetBuild() << " on " << strGetComputerName() << endl << endl;

   LogStream << PROGRAM_NAME << " for " << PLATFORM << " " << strGetBuild() << " on " << strGetComputerName() << endl << endl;

   // ----------------------------------------------- Run Information ----------------------------------------------------------
   OutStream << "RUN DETAILS" << endl;
   OutStream << " Name                                                      \t: " << m_strRunName << endl;
   OutStream << " Run started                                               \t: " << put_time(localtime(&m_tSysStartTime), "%T %A %d %B %Y") << endl;

   // Same info. for Log file
   LogStream << m_strRunName << " run started at " << put_time(localtime(&m_tSysStartTime), "%T on %A %d %B %Y") << endl << endl;

   // Contine with Out file
   OutStream << " Initialization file                                       \t: "
#ifdef _WIN32
      << pstrChangeToForwardSlash(&m_strCMEIni) << endl;
#else
      << m_strCMEIni << endl;
#endif

   OutStream << " Input data read from                                      \t: "
#ifdef _WIN32
      << pstrChangeToForwardSlash(&m_strDataPathName) << endl;
#else
      << m_strDataPathName << endl;
#endif

   OutStream << " Simulation start date/time                                \t: ";
   // hh:mm:ss dd/mm/yyyy
   char cPrev = OutStream.fill('0');
   OutStream << setw(2) << m_nSimStartHour << COLON << setw(2) << m_nSimStartMin << COLON << setw(2) << m_nSimStartSec << SPACE << setw(2) << m_nSimStartDay << SLASH << setw(2) << m_nSimStartMonth << SLASH << setw(2) << m_nSimStartYear << endl;
   OutStream.fill(cPrev);

   OutStream << " Duration of simulation                                    \t: ";
   OutStream << strDispSimTime(m_dSimDuration) << endl;
   if (m_bSaveRegular)
   {
      // Saves at regular intervals
      OutStream << " Time between saves                                        \t: ";
      OutStream << strDispSimTime(m_dRSaveInterval) << endl;
   }
   else
   {
      // Saves at user-defined intervals
      OutStream << " Saves at                                                  \t: ";
      string strTmp;
      for (int i = 0; i < m_nUSave; i++)
      {
         strTmp.append(strDispSimTime(m_dUSaveTime[i]));
         strTmp.append(", ");
      }

      // Also at end of run
      strTmp.append(strDispSimTime(m_dSimDuration));
      OutStream << strTmp << endl;
   }

   OutStream << " Random number seeds                                       \t: ";
   {
      for (int i = 0; i < NRNG; i++)
         OutStream << m_ulRandSeed[i] << '\t';
   }
   OutStream << endl;

   OutStream << "*First random numbers generated                            \t: " << ulGetRand0() << '\t' << ulGetRand1() << endl;
   OutStream << " Raster GIS output format                                  \t: " << m_strGDALRasterOutputDriverLongname << endl;
   OutStream << " Raster output values scaled (if needed)                   \t: " << (m_bScaleRasterOutput ? "Y": "N") << endl;
   OutStream << " Raster world files created (if needed)                    \t: " << (m_bWorldFile ? "Y": "N") << endl;
   OutStream << " Raster GIS files saved                                    \t: " << strListRasterFiles() << endl;
   if (m_bSliceSave)
   {
      OutStream << std::fixed << setprecision(2);
      OutStream << " Elevations for 'slice' raster output files                \t: ";
      for (int i = 0; i < static_cast<int>(m_VdSliceElev.size()); i++)
         OutStream << m_VdSliceElev[i] << " ";
      OutStream << endl;
   }

   OutStream << " Vector GIS output format                                  \t: " << m_strVectorGISOutFormat << endl;
   OutStream << " Vector GIS files saved                                    \t: " << strListVectorFiles() << endl;
   OutStream << " Output file (this file)                                   \t: "
#ifdef _WIN32
      << pstrChangeToForwardSlash(&m_strOutFile) << endl;
#else
      << m_strOutFile << endl;
#endif
   OutStream << " Log file                                                  \t: "
#ifdef _WIN32
      << pstrChangeToForwardSlash(&m_strLogFile) << endl;
#else
      << m_strLogFile << endl;
#endif

   OutStream << " Optional time series files saved                          \t: " << strListTSFiles() << endl;

   OutStream << " Coastline vector smoothing algorithm                      \t: ";
   switch (m_nCoastSmooth)
   {
      case SMOOTH_NONE:
      {
         OutStream << "none";
         break;
      }

      case SMOOTH_RUNNING_MEAN:
      {
         OutStream << "running mean";
         break;
      }

      case SMOOTH_SAVITZKY_GOLAY:
      {
         OutStream << "Savitzky-Golay";
         break;
      }
   }
   OutStream << endl;

   OutStream << " Grid edge(s) to omit when searching for coastline         \t: " << (m_bOmitSearchNorthEdge ? "N" : "") << (m_bOmitSearchSouthEdge ? "S" : "") << (m_bOmitSearchWestEdge ? "W" : "") << (m_bOmitSearchEastEdge ? "E" : "") << endl;

   if (m_nCoastSmooth != SMOOTH_NONE)
   {
      OutStream << " Size of coastline vector smoothing window                 \t: " << m_nCoastSmoothWindow << endl;

      if (m_nCoastSmooth == SMOOTH_SAVITZKY_GOLAY)
         OutStream << " Savitzky-Golay coastline smoothing polynomial order       \t: " << m_nSavGolCoastPoly << endl;
   }
   OutStream << " Size of profile slope smoothing window                    \t: " << m_nProfileSmoothWindow << endl;
   OutStream << " Max local slope on profile (m/m)                          \t: " << m_dProfileMaxSlope << endl;
   OutStream << " Vertical tolerance for beach to be included in smoothing  \t: " << m_dMaxBeachElevAboveSWL << " m" << endl;
   OutStream << endl;


   // --------------------------------------------------- Raster GIS stuff -------------------------------------------------------
   OutStream << "Raster GIS Input Files" << endl;
   OutStream << " Basement DEM file                                         \t: "
#ifdef _WIN32
      << pstrChangeToForwardSlash(&m_strInitialBasementDEMFile) << endl;
#else
      << m_strInitialBasementDEMFile << endl;
#endif
   OutStream << " Basement DEM driver code                                  \t: " << m_strGDALBasementDEMDriverCode << endl;
   OutStream << " GDAL basement DEM driver description                      \t: " << m_strGDALBasementDEMDriverDesc << endl;
   OutStream << " GDAL basement DEM projection                              \t: " << m_strGDALBasementDEMProjection << endl;
   OutStream << " GDAL basement DEM data type                               \t: " << m_strGDALBasementDEMDataType << endl;
   OutStream << " Grid size (X by Y)                                        \t: " << m_nXGridMax << " by " << m_nYGridMax << endl;
   OutStream << resetiosflags(ios::floatfield);
   OutStream << std::fixed << setprecision(1);
   OutStream << "*Coordinates of NW corner of grid (external CRS)           \t: " << m_dNorthWestXExtCRS << ", " << m_dNorthWestYExtCRS << endl;
   OutStream << "*Coordinates of SE corner of grid (external CRS)           \t: " << m_dSouthEastXExtCRS << ", " << m_dSouthEastYExtCRS << endl;
   OutStream << "*Cell size                                                 \t: " << m_dCellSide << " m" << endl;
   OutStream << "*Grid area                                                 \t: " << m_dExtCRSGridArea << " m^2" << endl;
   OutStream << std::fixed << setprecision(2);
   OutStream << "*Grid area                                                 \t: " << m_dExtCRSGridArea * 1e-6 << " km^2" << endl;

   if (! m_strInitialLandformFile.empty())
   {
      OutStream << " Initial Landform Class file                               \t: " << m_strInitialLandformFile << endl;
      OutStream << " GDAL Initial Landform Class file driver code              \t: " << m_strGDALLDriverCode << endl;
      OutStream << " GDAL Initial Landform Class file driver description       \t: " << m_strGDALLDriverDesc << endl;
      OutStream << " GDAL Initial Landform Class file projection               \t: " << m_strGDALLProjection << endl;
      OutStream << " GDAL Initial Landform Class file data type                \t: " << m_strGDALLDataType << endl;
      OutStream << endl;
   }

   if (! m_strInterventionClassFile.empty())
   {
      OutStream << " Intervention Class file                                   \t: " << m_strInterventionClassFile << endl;
      OutStream << " GDAL Intervention Class file driver code                  \t: " << m_strGDALICDriverCode << endl;
      OutStream << " GDAL Intervention Class file driver description           \t: " << m_strGDALICDriverDesc << endl;
      OutStream << " GDAL Intervention Class file projection                   \t: " << m_strGDALICProjection << endl;
      OutStream << " GDAL Intervention Class file data type                    \t: " << m_strGDALICDataType << endl;
      OutStream << endl;
   }

   if (! m_strInterventionHeightFile.empty())
   {
      OutStream << " Intervention Height file                                  \t: " << m_strInterventionHeightFile << endl;
      OutStream << " GDAL Intervention Height file driver code                 \t: " << m_strGDALIHDriverCode << endl;
      OutStream << " GDAL Intervention Height file driver description          \t: " << m_strGDALIHDriverDesc << endl;
      OutStream << " GDAL Intervention Height file projection                  \t: " << m_strGDALIHProjection << endl;
      OutStream << " GDAL Intervention Height file data type                   \t: " << m_strGDALIHDataType << endl;
      OutStream << endl;
   }

   if (! m_strInitialSuspSedimentFile.empty())
   {
      OutStream << " Initial Susp Sediment file                                \t: " << m_strInitialSuspSedimentFile << endl;
      OutStream << " GDAL Initial Susp Sediment file driver code               \t: " << m_strGDALISSDriverCode << endl;
      OutStream << " GDAL Initial Susp Sediment file driver description        \t: " << m_strGDALISSDriverDesc << endl;
      OutStream << " GDAL Initial Susp Sediment file projection                \t: " << m_strGDALISSProjection << endl;
      OutStream << " GDAL Initial Susp Sediment file data type                 \t: " << m_strGDALISSDataType << endl;
      OutStream << endl;
   }

   for (int i = 0; i < m_nLayers; i++)
   {
      if (m_nLayers == 1)
         OutStream << " Only one layer" << endl;
      else
         OutStream << " Layer " << i << (i == 0 ? "(Top)" : "") << (i == m_nLayers-1 ? "(Bottom)" : "") << endl;

      if (! m_VstrInitialFineUnconsSedimentFile[i].empty())
      {
         OutStream << "    Initial Fine Uncons Sediment file                      \t: " << m_VstrInitialFineUnconsSedimentFile[i] << endl;
         OutStream << "    GDAL Initial Fine Uncons Sediment file driver code     \t: " << m_VstrGDALIUFDriverCode[i] << endl;
         OutStream << "    GDAL Initial Fine Uncons Sediment file driver desc     \t: " << m_VstrGDALIUFDriverDesc[i] << endl;
         OutStream << "    GDAL Initial Fine Uncons Sediment file projection      \t: " << m_VstrGDALIUFProjection[i] << endl;
         OutStream << "    GDAL Initial Fine Uncons Sediment file data type       \t: " << m_VstrGDALIUFDataType[i] << endl;
         OutStream << endl;
      }

      if (! m_VstrInitialSandUnconsSedimentFile[i].empty())
      {
         OutStream << "    Initial Sand Uncons Sediment file                      \t: " << m_VstrInitialSandUnconsSedimentFile[i] << endl;
         OutStream << "    GDAL Initial Sand Uncons Sediment file driver code     \t: " << m_VstrGDALIUSDriverCode[i] << endl;
         OutStream << "    GDAL Initial Sand Uncons Sediment file driver desc     \t: " << m_VstrGDALIUSDriverDesc[i] << endl;
         OutStream << "    GDAL Initial Sand Uncons Sediment file projection      \t: " << m_VstrGDALIUSProjection[i] << endl;
         OutStream << "    GDAL Initial Sand Uncons Sediment file data type       \t: " << m_VstrGDALIUSDataType[i] << endl;
         OutStream << endl;
      }

      if (! m_VstrInitialCoarseUnconsSedimentFile[i].empty())
      {
         OutStream << "    Initial Coarse Uncons Sediment file                    \t: " << m_VstrInitialCoarseUnconsSedimentFile[i] << endl;
         OutStream << "    GDAL Initial Coarse Uncons Sediment file driver code   \t: " << m_VstrGDALIUCDriverCode[i] << endl;
         OutStream << "    GDAL Initial Coarse Uncons Sediment file driver desc   \t: " << m_VstrGDALIUCDriverDesc[i] << endl;
         OutStream << "    GDAL Initial Coarse Uncons Sediment file projection    \t: " << m_VstrGDALIUCProjection[i] << endl;
         OutStream << "    GDAL Initial Coarse Uncons Sediment file data type     \t: " << m_VstrGDALIUCDataType[i] << endl;
         OutStream << endl;
      }

      if (! m_VstrInitialFineConsSedimentFile[i].empty())
      {
         OutStream << "    Initial Fine Cons Sediment file                        \t: " << m_VstrInitialFineConsSedimentFile[i] << endl;
         OutStream << "    GDAL Initial Fine Cons Sediment file driver code       \t: " << m_VstrGDALICFDriverCode[i] << endl;
         OutStream << "    GDAL Initial Fine Cons Sediment file driver desc       \t: " << m_VstrGDALICFDriverDesc[i] << endl;
         OutStream << "    GDAL Initial Fine Cons Sediment file projection        \t: " << m_VstrGDALICFProjection[i] << endl;
         OutStream << "    GDAL Initial Fine Cons Sediment file data type         \t: " << m_VstrGDALICFDataType[i] << endl;
         OutStream << endl;
      }

      if (! m_VstrInitialSandConsSedimentFile[i].empty())
      {
         OutStream << "    Initial Sand Cons Sediment file                        \t: " << m_VstrInitialSandConsSedimentFile[i] << endl;
         OutStream << "    GDAL Initial Sand Cons Sediment file driver code       \t: " << m_VstrGDALICSDriverCode[i] << endl;
         OutStream << "    GDAL Initial Sand Cons Sediment file driver desc       \t: " << m_VstrGDALICSDriverDesc[i] << endl;
         OutStream << "    GDAL Initial Sand Cons Sediment file projection        \t: " << m_VstrGDALICSProjection[i] << endl;
         OutStream << "    GDAL Initial Sand Cons Sediment file data type         \t: " << m_VstrGDALICSDataType[i] << endl;
         OutStream << endl;
      }

      if (! m_VstrInitialCoarseConsSedimentFile[i].empty())
      {
         OutStream << "    Initial Coarse Cons Sediment file                      \t: " << m_VstrInitialCoarseConsSedimentFile[i] << endl;
         OutStream << "    GDAL Initial Coarse Cons Sediment file driver code     \t: " << m_VstrGDALICCDriverCode[i] << endl;
         OutStream << "    GDAL Initial Coarse Cons Sediment file driver desc     \t: " << m_VstrGDALICCDriverDesc[i] << endl;
         OutStream << "    GDAL Initial Coarse Cons Sediment file projection      \t: " << m_VstrGDALICCProjection[i] << endl;
         OutStream << "    GDAL Initial Coarse Cons Sediment file data type       \t: " << m_VstrGDALICCDataType[i] << endl;
         OutStream << endl;
      }
   }
//   OutStream << endl;

   // ---------------------------------------------------- Vector GIS stuff ------------------------------------------------------
   OutStream << "Vector GIS Input Files" << endl;

   if (m_bSingleDeepWaterWaveValues)
      OutStream << " None" << endl;
   else
   {

      OutStream << " Deep water wave stations shapefile                        \t: " << m_strDeepWaterWaveStationsShapefile << endl;
      OutStream << " GDAL/OGR deep water wave stations shapefile driver code   \t: " << m_strOGRDWWVDriverCode << endl;
      OutStream << " GDAL/OGR deep water wave stations shapefile data type     \t: " << m_strOGRDWWVDataType << endl;
      OutStream << " GDAL/OGR deep water wave stations shapefile geometry      \t: " << m_strOGRDWWVGeometry << endl;
      OutStream << " Deep water wave values file                               \t: " << m_strDeepWaterWavesTimeSeriesFile << endl;

      if (m_dWaveDataWrapHours > 0)
         OutStream << " Deep water wave values will wrap every " << m_dWaveDataWrapHours << " hours" << endl;
   }
   OutStream << endl;

   // -------------------------------------------------------- Other data --------------------------------------------------------
   OutStream << "Other Input Data" << endl;

   OutStream << " Wave propagation model                                    \t: ";
   if (m_nWavePropagationModel == WAVE_MODEL_COVE)
      OutStream << "COVE";
   else if (m_nWavePropagationModel == WAVE_MODEL_CSHORE)
      OutStream << "CShore";
   OutStream << endl;
   OutStream << " Density of sea water                                     \t: " << resetiosflags(ios::floatfield) << std::fixed << setprecision(0) << m_dSeaWaterDensity << " kg/m^3" << endl;
   OutStream << " Initial still water level                                 \t: " << resetiosflags(ios::floatfield) << std::fixed << setprecision(1) << m_dOrigSWL << " m" << endl;
   OutStream << " Final still water level                                   \t: " << resetiosflags(ios::floatfield) << std::fixed << setprecision(1) << m_dFinalSWL << " m" << endl;
   if (m_bSingleDeepWaterWaveValues)
   {
      OutStream << " Deep water wave height                                    \t: " << m_dAllCellsDeepWaterWaveHeight << " m" << endl;
      OutStream << " Deep water wave orientation                               \t: " << m_dAllCellsDeepWaterWaveAngle << " degrees" << endl;
      OutStream << " Wave period                                               \t: " << m_dAllCellsDeepWaterWavePeriod << " s" << endl;

   }
   else
   {
      OutStream << " Maximum User input Deep water wave height                 \t: " << m_dMaxUserInputWaveHeight << " m" << endl;
      OutStream << " Maximum User input Deep waterWave period                  \t: " << m_dMaxUserInputWavePeriod << " s" << endl;
   }
   OutStream << " Start depth for wave calcs (*deep water wave height)      \t: " << m_dWaveDepthRatioForWaveCalcs << endl;
   OutStream << "*Depth of closure                                          \t: " << resetiosflags(ios::floatfield) << std::fixed << setprecision(3) << m_dDepthOfClosure << " m" << endl;
//    OutStream << " Tide data file                                            \t: " << m_strTideDataFile << endl;
   OutStream << " Do coast platform erosion?                                \t: " << (m_bDoCoastPlatformErosion ? "Y": "N") << endl;
   OutStream << " R value                                                   \t: " << resetiosflags(ios::floatfield) << std::fixed << setprecision(2) << m_dR << endl;
   OutStream << " Handling of beach sediment at grid edges                  \t: ";
   if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_CLOSED)
      OutStream << "closed";
   else if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_OPEN)
      OutStream << "open";
   else if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_RECIRCULATE)
      OutStream << "recirculate";
   OutStream << endl;
   OutStream << " Beach potential erosion/deposition equation               \t: ";
   if (m_nBeachErosionDepositionEquation == UNCONS_SEDIMENT_EQUATION_CERC)
      OutStream << "CERC";
   else if (m_nBeachErosionDepositionEquation == UNCONS_SEDIMENT_EQUATION_KAMPHUIS)
      OutStream << "Kamphuis";
   OutStream << endl;
   OutStream << " Median particle size of fine sediment                     \t: " << resetiosflags(ios::floatfield) << std::fixed << m_dD50Fine << " mm" << endl;
   OutStream << " Median particle size of sand sediment                     \t: " << resetiosflags(ios::floatfield) << std::fixed << m_dD50Sand << " mm" << endl;
   OutStream << " Median particle size of coarse sediment                   \t: " << resetiosflags(ios::floatfield) << std::fixed << m_dD50Coarse << " mm" << endl;
   OutStream << " Beach sediment density                                    \t: " << resetiosflags(ios::floatfield) << std::fixed << m_dBeachSedimentDensity << " kg/m^3" << endl;
   OutStream << " Beach sediment porosity                                   \t: " << resetiosflags(ios::floatfield) << std::fixed << m_dBeachSedimentPorosity << endl;
   OutStream << " Fine-sized sediment relative erodibility                  \t: " << resetiosflags(ios::floatfield) << std::fixed << setprecision(1) << m_dFineErodibility << endl;
   OutStream << " Sand-sized sediment relative erodibility                  \t: " << resetiosflags(ios::floatfield) << m_dSandErodibility << endl;
   OutStream << " Coarse-sized sediment relative erodibility                \t: " << m_dCoarseErodibility << endl;
   if (m_nBeachErosionDepositionEquation == UNCONS_SEDIMENT_EQUATION_CERC)
      OutStream << " Transport parameter KLS for CERC equation                 \t: " << resetiosflags(ios::floatfield) << std::fixed << setprecision(3) << m_dKLS << endl;
   if (m_nBeachErosionDepositionEquation == UNCONS_SEDIMENT_EQUATION_KAMPHUIS)
      OutStream << " Transport parameter for Kamphuis equation                 \t: " << resetiosflags(ios::floatfield) << std::fixed << setprecision(3) << m_dKamphuis << endl;
   OutStream << " Height of Dean profile start above SWL                    \t: " << resetiosflags(ios::floatfield) << std::fixed << setprecision(1) << m_dDeanProfileStartAboveSWL  << " m" << endl;
   OutStream << " Do cliff collapse?                                        \t: " << (m_bDoCliffCollapse ? "Y": "N") << endl;
   OutStream << " Cliff erodibility                                         \t: " << m_dCliffErosionResistance << endl;
   OutStream << " Notch overhang to initiate collapse                       \t: " << m_dNotchOverhangAtCollapse << " m" << endl;
   OutStream << " Notch base below SWL                                      \t: " << m_dNotchBaseBelowSWL << " m" << endl;
   OutStream << " Scale parameter A for cliff deposition                    \t: ";
   if (m_dCliffDepositionA == 0)
      OutStream << "auto";
   else
      OutStream << m_dCliffDepositionA << "  m^(1/3)";
   OutStream << endl;
   OutStream << " Planview width of cliff deposition talus                  \t: " << resetiosflags(ios::floatfield) << std::fixed << m_dCliffDepositionPlanviewWidth << " m" << endl;
   OutStream << " Planview length of cliff deposition talus                 \t: " << m_dCliffDepositionPlanviewLength << " m" << endl;
   OutStream << " Height of talus at land end (fraction of cliff elevation) \t: " << m_dCliffDepositionHeightFrac << endl;
   OutStream << " Gravitational acceleration                                \t: " << resetiosflags(ios::floatfield) << std::fixed << m_dG << " m^2/s" << endl;
   OutStream << " Minimum spacing of coastline normals                      \t: " << resetiosflags(ios::floatfield) << std::fixed << m_dCoastNormalAvgSpacing << " m" << endl;
   OutStream << " Random factor for spacing of normals                      \t: " << resetiosflags(ios::floatfield) << std::fixed << m_dCoastNormalRandSpaceFact << endl;
   OutStream << " Length of coastline normals                               \t: " << m_dCoastNormalLength << " m" << endl;
   OutStream << " Maximum number of 'cape' normals                          \t: " << m_nNaturalCapeNormals << endl;
   OutStream << endl;
/*
   OutStream << std::fixed << setprecision(8);
   OutStream << " Erosion potential shape function:" << endl;
   OutStream << "\tDepth over DB\tErosion potential\tFirst derivative of erosion potential" << endl;
   for (int i = 0; i < m_VdDepthOverDB.size(); i++)
      OutStream << "\t" << m_VdDepthOverDB[i] << "\t\t" << m_VdErosionPotential[i] << "\t\t" << m_VdErosionPotentialFirstDeriv[i] << endl;
   OutStream << endl;
*/
   // ------------------------------------------------------ Testing only --------------------------------------------------------
   OutStream << "Testing only" << endl;

   OutStream << " Output profile data?                                      \t: " << (m_bOutputProfileData ? "Y": "N") << endl;
   OutStream << " Profile numbers to be saved                               \t: ";
   for (unsigned int i = 0; i < m_VnProfileToSave.size(); i++)
      OutStream << m_VnProfileToSave[i] << SPACE;
   OutStream << endl;
   OutStream << " Timesteps when profiles are saved                         \t: ";
   for (unsigned int i = 0; i < m_VlProfileTimestep.size(); i++)
      OutStream << m_VlProfileTimestep[i] << SPACE;
   OutStream << endl;
   OutStream << " Output parallel profile data?                             \t: " << (m_bOutputParallelProfileData ? "Y": "N") << endl;
   OutStream << " Output erosion potential look-up data?                    \t: " << (m_bOutputLookUpData ? "Y": "N");
   if (m_bOutputLookUpData)
      OutStream << " (see " << m_strOutPath << EROSION_POTENTIAL_LOOKUP_FILE << ")";
   OutStream << endl;
   OutStream << " Erode coast in alternate directions?                      \t: " << (m_bErodeShorePlatformAlternateDirection ? "Y": "N") << endl;

   OutStream << endl << endl;

   // -------------------------------------------------- Per-timestep output ----------------------------------------------------
   OutStream << std::fixed << setprecision(2);

   // Write per-timestep headers to .out file
   OutStream << PERITERHEAD << endl;
   OutStream << "Sea depth in metres. All erosion and deposition values in millimetres" << endl;
   OutStream << "GISn = GIS files saved as <filename>n." << endl;
   OutStream << endl;

   OutStream << PERITERHEAD1 << endl;
   OutStream << PERITERHEAD2 << endl;
   OutStream << PERITERHEAD3 << endl;
   OutStream << PERITERHEAD4 << endl;
   OutStream << PERITERHEAD5 << endl;
}


/*===============================================================================================================================

 Write the results for this timestep to the .out file

===============================================================================================================================*/
bool CSimulation::bWritePerTimestepResults(void)
{
   OutStream << resetiosflags(ios::floatfield);
   OutStream << std::fixed << setprecision(0);

   // Output timestep and simulated time info ===================================================================================
   OutStream << setw(4) << m_ulIter;
   OutStream << setw(7) << m_dSimElapsed;                            // In hours
   OutStream << resetiosflags(ios::floatfield);
   OutStream << std::fixed << setprecision(0);
   OutStream << setw(7) << m_dSimElapsed / (24 * 365.25);            // In years

   // Output average sea depth (m) per sea cell =================================================================================
   OutStream << resetiosflags(ios::floatfield);
   OutStream << std::fixed << setprecision(2);
   double dAvgSeaDepth = m_dThisTimestepTotSeaDepth / static_cast<double>(m_ulThisTimestepNumSeaCells);
   OutStream << setw(6) << dAvgSeaDepth;
   OutStream << " ";

   // Output the this-timestep % of sea cells with potential shore platform erosion =============================================
   OutStream << std::fixed << setprecision(1);
   OutStream << setw(6) << 100 * static_cast<double>(m_ulThisTimestepNumPotentialPlatformErosionCells) / static_cast<double>(m_ulThisTimestepNumSeaCells);

   // Output per-timestep potential shore platform erosion in mm (average for all sea cells)
   OutStream << setw(6) << 1000 * m_dThisTimestepPotentialPlatformErosion / static_cast<double>(m_ulThisTimestepNumSeaCells);

   // Output per-timestep potential shore platform erosion in mm (average for all cells with potential shore platform erosion)
   OutStream << std::fixed << setprecision(0);
   if (m_ulThisTimestepNumPotentialPlatformErosionCells > 0)
      OutStream << setw(6) << 1000 * m_dThisTimestepPotentialPlatformErosion / static_cast<double>(m_ulThisTimestepNumPotentialPlatformErosionCells);
   else
      OutStream << setw(6) << SPACE;

   // Output the this-timestep % of sea cells with actual shore platform erosion ================================================
   OutStream << std::fixed << setprecision(1);
   OutStream << setw(6) << 100 * static_cast<double>(m_ulThisTimestepNumActualPlatformErosionCells) / static_cast<double>(m_ulThisTimestepNumSeaCells);

   // Output per-timestep actual shore platform erosion in mm (average for all sea cells)
   double dThisTimestepActualPlatformErosion = m_dThisTimestepActualPlatformErosionFine + m_dThisTimestepActualPlatformErosionSand + m_dThisTimestepActualPlatformErosionCoarse;
   OutStream << setw(6) << 1000 * dThisTimestepActualPlatformErosion / static_cast<double>(m_ulThisTimestepNumSeaCells);

   // Output per-timestep actual shore platform erosion in mm (average for all cells with actual shore platform erosion)
   OutStream << std::fixed << setprecision(0);
   if (m_ulThisTimestepNumActualPlatformErosionCells > 0)
      OutStream << setw(5) << 1000 * dThisTimestepActualPlatformErosion / static_cast<double>(m_ulThisTimestepNumActualPlatformErosionCells);
   else
      OutStream << setw(5) << SPACE;

   // Output per-timestep actual shore platform erosion in mm (average for all sea cells)
   OutStream << std::fixed << setprecision(1);
   OutStream << setw(4) << 1000 * m_dThisTimestepActualPlatformErosionFine / static_cast<double>(m_ulThisTimestepNumSeaCells);
   OutStream << setw(4) << 1000 * m_dThisTimestepActualPlatformErosionSand / static_cast<double>(m_ulThisTimestepNumSeaCells);
   OutStream << setw(4) << 1000 * m_dThisTimestepActualPlatformErosionCoarse / static_cast<double>(m_ulThisTimestepNumSeaCells);

   // Output the this-timestep % of sea cells with potential beach erosion ======================================================
   OutStream << setw(7) << 100 * static_cast<double>(m_ulThisTimestepNumPotentialBeachErosionCells) / static_cast<double>(m_ulThisTimestepNumSeaCells);

   // Output per-timestep potential beach erosion in mm (average for all sea cells)
   OutStream << std::fixed << setprecision(0);
   assert(m_ulThisTimestepNumSeaCells > 0);
   double dTmp = 1000 * m_dThisTimestepPotentialBeachErosion / static_cast<double>(m_ulThisTimestepNumSeaCells);
   if (dTmp  > 99999)
   {
      OutStream << setw(6) << std::scientific << setprecision(0) << dTmp;
      OutStream << std::fixed;
   }
   else
      OutStream << setw(6) << dTmp;

   // Output per-timestep potential beach erosion in mm (average for all cells with potential beach erosion)
   OutStream << std::fixed << setprecision(0);
   if (m_ulThisTimestepNumPotentialBeachErosionCells > 0)
   {
      dTmp = 1000 * m_dThisTimestepPotentialBeachErosion / static_cast<double>(m_ulThisTimestepNumPotentialBeachErosionCells);
      if (dTmp > 99999)
      {
         OutStream << setw(6) << std::scientific << setprecision(0) << dTmp;
         OutStream << std::fixed;
      }
      else
         OutStream << setw(6) << dTmp;
   }
   else
      OutStream << setw(6) << SPACE;

   // This-timestep % of sea cells with actual beach erosion ====================================================================
   OutStream << std::fixed << setprecision(1);
   OutStream << setw(7) << 100 * static_cast<double>(m_ulThisTimestepNumActualBeachErosionCells) / static_cast<double>(m_ulThisTimestepNumSeaCells);

   // Output per-timestep actual beach erosion in mm (average for all sea cells)
   double dThisTimestepActualBeachErosion = m_dThisTimestepActualBeachErosionFine + m_dThisTimestepActualBeachErosionSand + m_dThisTimestepActualBeachErosionCoarse;
   OutStream << setw(6) << 1000 * dThisTimestepActualBeachErosion / static_cast<double>(m_ulThisTimestepNumSeaCells);

   // Per-timestep actual beach erosion in mm (average for all cells with actual beach erosion)
   OutStream << std::fixed << setprecision(0);
   if (m_ulThisTimestepNumActualBeachErosionCells > 0)
      OutStream << setw(6) << 1000 * dThisTimestepActualBeachErosion / static_cast<double>(m_ulThisTimestepNumActualBeachErosionCells);
   else
      OutStream << setw(6) << SPACE;

   // Per-timestep actual beach erosion in mm (average for all sea cells)
   OutStream << std::fixed << setprecision(1);
   OutStream << setw(4) << 1000 * m_dThisTimestepActualBeachErosionFine / static_cast<double>(m_ulThisTimestepNumSeaCells);
   OutStream << setw(4) << 1000 * m_dThisTimestepActualBeachErosionSand / static_cast<double>(m_ulThisTimestepNumSeaCells);
   OutStream << setw(4) << 1000 * m_dThisTimestepActualBeachErosionCoarse / static_cast<double>(m_ulThisTimestepNumSeaCells);

   // Output the this-timestep % of sea cells with beach deposition =============================================================
   OutStream << setw(7) << 100 * static_cast<double>(m_ulThisTimestepNumBeachDepositionCells) / static_cast<double>(m_ulThisTimestepNumSeaCells);

   // Per-timestep beach deposition in mm (average for all sea cells)
   double dThisTimestepBeachDeposition = m_dThisTimestepBeachDepositionSand + m_dThisTimestepBeachDepositionCoarse;
   OutStream << setw(6) << 1000 * dThisTimestepBeachDeposition / static_cast<double>(m_ulThisTimestepNumSeaCells);

   // Per-timestep beach deposition in mm (average for all cells with beach deposition)
   OutStream << std::fixed << setprecision(0);
   if (m_ulThisTimestepNumBeachDepositionCells > 0)
      OutStream << setw(8) << 1000 * dThisTimestepBeachDeposition / static_cast<double>(m_ulThisTimestepNumBeachDepositionCells);
   else
      OutStream << setw(8) << SPACE;

   // Per-timestep beach deposition in mm (average for all sea cells)
   OutStream << std::fixed << setprecision(1);
   OutStream << setw(4) << 1000 * m_dThisTimestepBeachDepositionSand / static_cast<double>(m_ulThisTimestepNumSeaCells);
   OutStream << setw(4) << 1000 * m_dThisTimestepBeachDepositionCoarse / static_cast<double>(m_ulThisTimestepNumSeaCells);

   // Per-timestep cliff collapse erosion in mm (average for all coast cells) ===========================================================
   OutStream << setw(4) << 1000 * m_dThisTimestepCliffCollapseErosionFine / static_cast<double>(m_ulThisTimestepNumCoastCells);
   OutStream << setw(4) << 1000 * m_dThisTimestepCliffCollapseErosionSand / static_cast<double>(m_ulThisTimestepNumCoastCells);
   OutStream << setw(4) << 1000 * m_dThisTimestepCliffCollapseErosionCoarse / static_cast<double>(m_ulThisTimestepNumCoastCells);

   // Per-timestep cliff collapse deposition in mm (average for all sea cells) ==================================================
   OutStream << setw(4) << 1000 * m_dThisTimestepCliffDepositionSand / static_cast<double>(m_ulThisTimestepNumSeaCells);
   OutStream << setw(4) << 1000 * m_dThisTimestepCliffDepositionCoarse / static_cast<double>(m_ulThisTimestepNumSeaCells);

   // Output per-timestep fine sediment going to suspension, in mm (average for all sea cells) ==================================
   OutStream << setw(6) << 1000 * m_dThisTimestepFineSedimentToSuspension / static_cast<double>(m_ulThisTimestepNumSeaCells);

   OutStream << " ";

   // Finally, set 'markers' for events that have occurred this timestep
   if (m_bSaveGISThisTimestep)
      OutStream << " GIS" << m_nGISSave;

   OutStream << endl;

   // Did a text file write error occur?
   if (OutStream.fail())
      return false;

   return true;
}

/*===============================================================================================================================

 Write the results for this timestep to the time series CSV files

===============================================================================================================================*/
bool CSimulation::bWriteTSFiles(void)
{
   // Sea area
   if (m_bSeaAreaTS)
   {
      // Output in external CRS units
      SeaAreaTSStream << m_dSimElapsed << "\t,\t" << m_dExtCRSGridArea * static_cast<double>(m_ulThisTimestepNumSeaCells) / static_cast<double>(m_ulNumCells) << endl;

      // Did a time series file write error occur?
      if (SeaAreaTSStream.fail())
         return false;
   }

   // Still water level
   if (m_bStillWaterLevelTS)
   {
      // Output as is (m)
      StillWaterLevelTSStream << m_dSimElapsed << "\t,\t" << m_dThisTimestepSWL << endl;

      // Did a time series file write error occur?
      if (StillWaterLevelTSStream.fail())
         return false;
   }

   // Actual platform erosion (fine, sand, and coarse)
   if (m_bActualPlatformErosionTS)
   {
      // Output as is (m depth equivalent)
      ErosionTSStream << m_dSimElapsed  << "\t,\t" << m_dThisTimestepActualPlatformErosionFine << ",\t" << m_dThisTimestepActualPlatformErosionSand << ",\t" << m_dThisTimestepActualPlatformErosionCoarse << endl;

      // Did a time series file write error occur?
      if (ErosionTSStream.fail())
         return false;
   }

   // Cliff collapse erosion (fine, sand, and coarse)
   if (m_bCliffCollapseErosionTS)
   {
      // Output as is (m depth equivalent)
      CliffCollapseErosionTSStream << m_dSimElapsed << "\t,\t" << m_dThisTimestepCliffCollapseErosionFine << ",\t" << m_dThisTimestepCliffCollapseErosionSand << ",\t" << m_dThisTimestepCliffCollapseErosionCoarse << endl;

      // Did a time series file write error occur?
      if (CliffCollapseErosionTSStream.fail())
         return false;
   }

   // Cliff collapse deposition (sand and coarse)
   if (m_bCliffCollapseDepositionTS)
   {
      // Output as is (m depth equivalent)
      CliffCollapseDepositionTSStream << m_dSimElapsed << "\t,\t" << m_dThisTimestepCliffDepositionSand << ",\t" << m_dThisTimestepCliffDepositionCoarse << endl;

      // Did a time series file write error occur?
      if (CliffCollapseDepositionTSStream.fail())
         return false;
   }

   // Cliff collapse net
   if (m_bCliffCollapseNetTS)
   {
      // Output as is (m depth equivalent)
      CliffCollapseNetTSStream << noshowpos << m_dSimElapsed << "\t,\t" << showpos << -m_dThisTimestepCliffErosionFine + (m_dThisTimestepCliffDepositionSand - m_dThisTimestepCliffTalusSandErosion) + (m_dThisTimestepCliffDepositionCoarse - m_dThisTimestepCliffTalusCoarseErosion) << endl;

      // Did a time series file write error occur?
      if (CliffCollapseNetTSStream.fail())
         return false;
   }

   // Beach erosion (fine, sand, and coarse)
   if (m_bBeachErosionTS)
   {
      // Output as is (m depth equivalent)
      BeachErosionTSStream << m_dSimElapsed << "\t,\t" << m_dThisTimestepActualBeachErosionFine << ",\t" << m_dThisTimestepActualBeachErosionSand << ",\t" << m_dThisTimestepActualBeachErosionCoarse << endl;

      // Did a time series file write error occur?
      if (BeachErosionTSStream.fail())
         return false;
   }

   // Beach deposition (sand and coarse)
   if (m_bBeachDepositionTS)
   {
      // Output as is (m depth equivalent)
      BeachDepositionTSStream << m_dSimElapsed << "\t,\t" << m_dThisTimestepBeachDepositionSand << ",\t" << m_dThisTimestepBeachDepositionCoarse << endl;

      // Did a time series file write error occur?
      if (BeachDepositionTSStream.fail())
         return false;
   }

   // Net change in beach sediment
   if (m_bBeachSedimentChangeNetTS)
   {
      // Output as is (m depth equivalent)
      BeachSedimentChangeNetTSStream << noshowpos << m_dSimElapsed << "\t,\t" << showpos << -m_dThisTimestepActualBeachErosionFine + (m_dThisTimestepBeachDepositionSand - m_dThisTimestepActualBeachErosionSand) + (m_dThisTimestepBeachDepositionCoarse - m_dThisTimestepActualBeachErosionCoarse) << endl;

      // Did a time series file write error occur?
      if (BeachSedimentChangeNetTSStream.fail())
         return false;
   }

   if (m_bSuspSedTS)
   {
      // Output as is (m depth equivalent)
      SedLoadTSStream << m_dSimElapsed << "\t,\t" << m_dThisTimestepFineSedimentToSuspension << endl;

      // Did a time series file write error occur?
      if (SedLoadTSStream.fail())
         return false;
   }

   return true;
}

/*==============================================================================================================================

 Output the erosion potential look-up values, for checking purposes

==============================================================================================================================*/
void CSimulation::WriteLookUpData(void)
{
   // Open the output file
   string strLookUpFile = m_strOutPath;
   strLookUpFile.append(EROSION_POTENTIAL_LOOKUP_FILE);
   ofstream LookUpOutStream;
   LookUpOutStream.open(strLookUpFile.c_str(), ios::out | ios::trunc);

   if (LookUpOutStream)
   {
      // File opened OK, so output the values
      LookUpOutStream << "DepthOverDB, \tErosionPotential" << endl;
      double dDepthOverDB = 0.0;
      while (dDepthOverDB <= m_dDepthOverDBMax)
      {
	 double dErosionPotential = interpolate( m_VdDepthOverDB, m_VdErosionPotential, dDepthOverDB, false );
         LookUpOutStream << dDepthOverDB << ",\t" << dErosionPotential << endl;
	 dDepthOverDB += DEPTH_OVER_DB_INCREMENT;
      }
      LookUpOutStream << endl;

      // And close the file
      LookUpOutStream.close();
   }
}


/*==============================================================================================================================

 Save a coastline-normal profile

==============================================================================================================================*/
int CSimulation::nSaveProfile(int const nProfile, int const nCoast, int const nProfSize, vector<double> const* pdVDistXY, vector<double> const* pdVZ, vector<double> const* pdVDepthOverDB, vector<double> const* pdVErosionPotentialFunc, vector<double> const* pdVSlope, vector<double> const* pdVRecessionXY, vector<double> const* pdVChangeElevZ, vector<CGeom2DIPoint>* const pPtVGridProfile, vector<double> const* pdVScapeXY)
{
   // TODO make this more efficient, also give warnings if no profiles will be output
   for (unsigned int i = 0; i < m_VlProfileTimestep.size(); i++)
   {
      for (unsigned int j = 0; j < m_VnProfileToSave.size(); j++)
      {
         if ((m_ulIter == m_VlProfileTimestep[i]) && (nProfile == m_VnProfileToSave[j]))
         {
            if (! bWriteProfileData(nCoast, nProfile, nProfSize, pdVDistXY, pdVZ, pdVDepthOverDB, pdVErosionPotentialFunc, pdVSlope, pdVRecessionXY, pdVChangeElevZ, pPtVGridProfile, pdVScapeXY))
               return RTN_ERR_PROFILEWRITE;
         }
      }
   }

   return RTN_OK;
}


/*==============================================================================================================================

 Writes values for a single profile, for checking purposes

==============================================================================================================================*/
bool CSimulation::bWriteProfileData(int const nCoast, int const nProfile, int const nProfSize, vector<double> const* pdVDistXY, vector<double> const* pdVZ, vector<double> const* pdVDepthOverDB, vector<double> const* pdVErosionPotentialFunc, vector<double> const* pdVSlope, vector<double> const* pdVRecessionXY, vector<double> const* pdVChangeElevZ, vector<CGeom2DIPoint>* const pPtVGridProfile, vector<double> const* pdVScapeXY) const
{
   string strFName = m_strOutPath;
   stringstream ststrTmp;

   strFName.append("profile_");
   ststrTmp << FillToWidth('0', 3) << nProfile;
   strFName.append(ststrTmp.str());

   strFName.append("_timestep_");
   ststrTmp.clear();
   ststrTmp.str(string());
   ststrTmp << FillToWidth('0', 4) << m_ulIter;
   strFName.append(ststrTmp.str());

   strFName.append(".csv");

   ofstream OutProfStream;
   OutProfStream.open(strFName.c_str(), ios::out | ios::trunc);
   if (! OutProfStream)
   {
      // Error, cannot open file
      cerr << ERR << "cannot open " << strFName << " for output" << endl;
      return false;
   }

   OutProfStream << "\"Dist\", \"X\", \"Y\", \"Z (before erosion)\", \"Depth/DB\", \"Erosion Potential\", \"Slope\", \"Recession XY\", \"Change Elev Z\", \"Grid X\",  \"Grid Y\",  \"Weight\",  \"For profile " << nProfile << " from coastline " << nCoast << " at timestep " << m_ulIter << "\"" << endl;
   for (int i = 0; i < nProfSize; i++)
   {
      double dX = dGridCentroidXToExtCRSX(pPtVGridProfile->at(i).nGetX());
      double dY = dGridCentroidYToExtCRSY(pPtVGridProfile->at(i).nGetY());

      OutProfStream << pdVDistXY->at(i) << ",\t" << dX << ",\t" << dY << ",\t" << pdVZ->at(i) << ",\t" << pdVDepthOverDB->at(i) << ",\t" << pdVErosionPotentialFunc->at(i) << ",\t" << pdVSlope->at(i) << ",\t" << pdVRecessionXY->at(i) << ",\t" << pdVChangeElevZ->at(i) << ",\t" <<  pPtVGridProfile->at(i).nGetX() <<  ",\t" << pPtVGridProfile->at(i).nGetY() <<  ", \t" << pdVScapeXY->at(i) << endl;
   }

   OutProfStream.close();

   return true;
}


/*==============================================================================================================================

 Save a coastline-normal parallel profile

==============================================================================================================================*/
int CSimulation::nSaveParProfile(int const nProfile, int const nCoast, int const nParProfSize, int const nDirection, int const nDistFromProfile, vector<double> const* pdVDistXY, vector<double> const* pdVZ, vector<double> const* pdVDepthOverDB, vector<double> const* pdVErosionPotentialFunc, vector<double> const* pdVSlope, vector<double> const* pdVRecessionXY, vector<double> const* pdVChangeElevZ, vector<CGeom2DIPoint>* const pPtVGridProfile, vector<double> const* pdVScapeXY)
{
   // TODO make this more efficient, also give warnings if no profiles will be output
   for (unsigned int i = 0; i < m_VlProfileTimestep.size(); i++)
   {
      for (unsigned int j = 0; j < m_VnProfileToSave.size(); j++)
      {
         if ((m_ulIter == m_VlProfileTimestep[i]) && (nProfile == m_VnProfileToSave[j]))
         {
            if (! bWriteParProfileData(nCoast, nProfile, nParProfSize, nDirection, nDistFromProfile, pdVDistXY, pdVZ, pdVDepthOverDB, pdVErosionPotentialFunc, pdVSlope, pdVRecessionXY, pdVChangeElevZ, pPtVGridProfile, pdVScapeXY))
               return RTN_ERR_PROFILEWRITE;
         }
      }
   }

   return RTN_OK;
}


/*==============================================================================================================================

 Writes values for a single parallel profile, for checking purposes

==============================================================================================================================*/
bool CSimulation::bWriteParProfileData(int const nCoast, int const nProfile, int const nProfSize, int const nDirection, int const nDistFromProfile, vector<double> const* pdVDistXY, vector<double> const* pdVZ, vector<double> const* pdVDepthOverDB, vector<double> const* pdVErosionPotentialFunc, vector<double> const* pdVSlope, vector<double> const* pdVRecessionXY, vector<double> const* pdVChangeElevZ, vector<CGeom2DIPoint>* const pPtVGridProfile, vector<double> const* pdVScapeXY) const
{
   string strFName = m_strOutPath;
   stringstream ststrTmp;

   strFName.append("profile_");
   ststrTmp << FillToWidth('0', 3) << nProfile;
   strFName.append(ststrTmp.str());

   strFName.append("_parallel_");
   ststrTmp.clear();
   ststrTmp.str(string());
   ststrTmp << FillToWidth('0', 3) << nDistFromProfile;
   strFName.append(ststrTmp.str());

   strFName.append((nDirection == 0 ? "_F" : "_B"));

   strFName.append("_timestep_");
   ststrTmp.clear();
   ststrTmp.str(string());
   ststrTmp << FillToWidth('0', 4) << m_ulIter;
   strFName.append(ststrTmp.str());

   strFName.append(".csv");

   ofstream OutProfStream;
   OutProfStream.open(strFName.c_str(), ios::out | ios::trunc);
   if (! OutProfStream)
   {
      // Error, cannot open file
      cerr << ERR << "cannot open " << strFName << " for output" << endl;
      return false;
   }

   OutProfStream << "\"Dist\", \"X\", \"Y\", \"Z (before erosion)\", \"Depth/DB\", \"Erosion Potential\", \"Slope\", \"Recession XY\", \"Change Elev Z\", \"Grid X\",  \"Grid Y\",  \"Weight\",  \"For profile " << nProfile << " from coastline " << nCoast << " at timestep " << m_ulIter << "\"" << endl;
   for (int i = 0; i < nProfSize; i++)
   {
      double dX = dGridCentroidXToExtCRSX(pPtVGridProfile->at(i).nGetX());
      double dY = dGridCentroidYToExtCRSY(pPtVGridProfile->at(i).nGetY());

      OutProfStream << pdVDistXY->at(i) << ",\t" << dX << ",\t" << dY << ",\t" << pdVZ->at(i) << ",\t" << pdVDepthOverDB->at(i) << ",\t" << pdVErosionPotentialFunc->at(i) << ",\t" << pdVSlope->at(i) << ",\t" << pdVRecessionXY->at(i) << ",\t" << pdVChangeElevZ->at(i) << ",\t" <<  pPtVGridProfile->at(i).nGetX() <<  ",\t" << pPtVGridProfile->at(i).nGetY() <<  ", \t" << pdVScapeXY->at(i) << endl;
   }

   OutProfStream.close();

   return true;
}


/*==============================================================================================================================

 Writes end-of-run information to Out, Log and time-series files

==============================================================================================================================*/
int CSimulation::nWriteEndRunDetails(void)
{
   // Final write to time series CSV files
   if (! bWriteTSFiles())
      return (RTN_ERR_TIMESERIES_FILE_WRITE);

   // Save the values from the RasterGrid array into raster GIS files
   if (! bSaveAllRasterGISFiles())
      return (RTN_ERR_RASTER_FILE_WRITE);

   // Save the vector GIS files
   if (! bSaveAllVectorGISFiles())
      return (RTN_ERR_VECTOR_FILE_WRITE);

   OutStream << " GIS" << m_nGISSave << endl;

   // Print out run totals etc.
   OutStream << PERITERHEAD1 << endl;
   OutStream << PERITERHEAD2 << endl;
   OutStream << PERITERHEAD3 << endl;
   OutStream << PERITERHEAD4 << endl;
   OutStream << PERITERHEAD5 << endl;

   OutStream << std::fixed << setprecision(2);
   OutStream << endl << endl;

   // Write out hydrology grand totals etc.
   OutStream << ENDHYDROLOGYHEAD << endl;
   OutStream << "Minimum still water level = " << m_dMinSWL << endl;
   OutStream << "Maximum still water level = " << m_dMaxSWL << endl;
   OutStream << endl;

   // Now write out sediment movement grand totals etc.
   OutStream << ENDSEDIMENTHEAD << endl << endl;

   OutStream << "PLATFORM EROSION" << endl;
   OutStream << "Total potential platform erosion = " << m_ldGTotPotentialPlatformErosion * m_dCellArea << " m^3" << endl << endl;

   OutStream << "Total actual platform erosion = " << (m_ldGTotFineActualPlatformErosion + m_ldGTotSandActualPlatformErosion + m_ldGTotCoarseActualPlatformErosion) * m_dCellArea << " m^3" << endl;
   OutStream << "Total fine actual platform erosion = " << m_ldGTotFineActualPlatformErosion * m_dCellArea << " m^3" << endl;
   OutStream << "Total sand actual platform erosion = " << m_ldGTotSandActualPlatformErosion * m_dCellArea << " m^3" << endl;
   OutStream << "Total coarse actual platform erosion = " << m_ldGTotCoarseActualPlatformErosion * m_dCellArea << " m^3" << endl;
   OutStream << endl;

   OutStream << "CLIFF COLLAPSE" << endl;
   OutStream << "Total cliff collapse = " << (m_ldGTotCliffCollapseFine + m_ldGTotCliffCollapseSand + m_ldGTotCliffCollapseCoarse) * m_dCellArea << " m^3" << endl;
   OutStream << "Total fine cliff collapse = " << m_ldGTotCliffCollapseFine * m_dCellArea << " m^3" << endl;
   OutStream << "Total sand cliff collapse = " << m_ldGTotCliffCollapseSand * m_dCellArea << " m^3" << endl;
   OutStream << "Total coarse cliff collapse = " << m_ldGTotCliffCollapseCoarse * m_dCellArea << " m^3" << endl;
   OutStream << endl;

   OutStream << "DEPOSITION OF CLIFF COLLAPSE TALUS" << endl;
   OutStream << "Total cliff collapse talus deposition = " << (m_ldGTotCliffTalusSandDeposition + m_ldGTotCliffTalusCoarseDeposition) * m_dCellArea << " m^3" << endl;
   OutStream << "Total sand cliff collapse talus deposition = " << m_ldGTotCliffTalusSandDeposition * m_dCellArea << " m^3" << endl;
   OutStream << "Total coarse cliff collapse talus deposition = " << m_ldGTotCliffTalusCoarseDeposition * m_dCellArea << " m^3" << endl;
   OutStream << endl;

   OutStream << "EROSION OF CLIFF COLLAPSE TALUS" << endl;
   OutStream << "Total erosion of cliff collapse talus = " << (m_ldGTotCliffTalusFineErosion + m_ldGTotCliffTalusSandErosion + m_ldGTotCliffTalusCoarseErosion) * m_dCellArea << " m^3" << endl;
   OutStream << "Total fine erosion of cliff collapse talus = " << m_ldGTotCliffTalusFineErosion * m_dCellArea << " m^3" << endl;
   OutStream << "Total sand erosion of cliff collapse talus = " << m_ldGTotCliffTalusSandErosion * m_dCellArea << " m^3" << endl;
   OutStream << "Total coarse erosion of cliff collapse talus = " << m_ldGTotCliffTalusCoarseErosion * m_dCellArea << " m^3" << endl;
   OutStream << endl;

   OutStream << "BEACH EROSION" << endl;
   OutStream << "Total potential beach erosion = " << m_ldGTotPotentialBeachErosion * m_dCellArea << " m^3" << endl << endl;
   OutStream << "Total actual beach erosion = " << (m_ldGTotActualFineBeachErosion + m_ldGTotActualSandBeachErosion + m_ldGTotActualCoarseBeachErosion) * m_dCellArea << " m^3" << endl;
   OutStream << "Total actual fine beach erosion = " << m_ldGTotActualFineBeachErosion * m_dCellArea << " m^3" << endl;
   OutStream << "Total actual sand beach erosion = " << m_ldGTotActualSandBeachErosion * m_dCellArea << " m^3" << endl;
   OutStream << "Total actual coarse beach erosion = " << m_ldGTotActualCoarseBeachErosion * m_dCellArea << " m^3" << endl;
   OutStream << endl;

   OutStream << "BEACH DEPOSITION" << endl;
   OutStream << "Total beach deposition = " << (m_ldGTotSandBeachDeposition + m_ldGTotCoarseBeachDeposition) * m_dCellArea << " m^3" << endl;
   OutStream << "Total sand beach deposition = " << m_ldGTotSandBeachDeposition * m_dCellArea << " m^3" << endl;
   OutStream << "Total coarse beach deposition = " << m_ldGTotCoarseBeachDeposition * m_dCellArea << " m^3" << endl;
   OutStream << endl;

   OutStream << "SUSPENDED SEDIMENT" << endl;
   OutStream << "Suspended fine sediment = " << m_ldGTotSuspendedSediment * m_dCellArea << " m^3" << endl;
   OutStream << endl;

   OutStream << "LOST FROM GRID" << endl;
   OutStream << "Total potential sediment lost from grid due to beach erosion = " << m_ldGTotPotentialSedLostBeachErosion * m_dCellArea << " m^3" << endl << endl;

   OutStream << "Total actual fine sediment lost from grid due to beach erosion = " << m_ldGTotActualFineSedLostBeachErosion * m_dCellArea << " m^3" << endl;
   OutStream << "Total actual sand sediment lost from grid due to beach erosion = " << m_ldGTotActualSandSedLostBeachErosion * m_dCellArea << " m^3" << endl;
   OutStream << "Total actual coarse sediment lost from grid due to beach erosion = " << m_ldGTotActualCoarseSedLostBeachErosion * m_dCellArea << " m^3" << endl << endl;

   OutStream << "Total sand sediment lost from grid due to cliff collapse = " << m_ldGTotSandSedLostCliffCollapse * m_dCellArea << " m^3" << endl;
   OutStream << "Total coarse sediment lost from grid due to cliff collapse = " << m_ldGTotCoarseSedLostCliffCollapse * m_dCellArea << " m^3" << endl;
   OutStream << endl;

   OutStream << "ERRORS" << endl;
   OutStream << "Erosion error = " << m_ldGTotMassBalanceErosionError << " m^3" << endl;
   OutStream << "Deposition error = " << m_ldGTotMassBalanceDepositionError << " m^3" << endl;
   OutStream << endl;

   OutStream << "ALL-PROCESS TOTALS" << endl;
   long double ldActualTotalEroded = m_ldGTotFineActualPlatformErosion + m_ldGTotSandActualPlatformErosion + m_ldGTotCoarseActualPlatformErosion + m_ldGTotCliffCollapseFine + m_ldGTotCliffCollapseSand + m_ldGTotCliffCollapseCoarse + m_ldGTotCliffTalusFineErosion + m_ldGTotCliffTalusSandErosion + m_ldGTotCliffTalusCoarseErosion + m_ldGTotActualFineBeachErosion + m_ldGTotActualSandBeachErosion + m_ldGTotActualCoarseBeachErosion;
   OutStream << "Total sediment eroded (all processes) = " << ldActualTotalEroded * m_dCellArea << " m^3" << endl;

   long double ldTotalDepositedAndSuspension = m_ldGTotCliffTalusSandDeposition + m_ldGTotCliffTalusCoarseDeposition + m_ldGTotSandBeachDeposition + m_ldGTotCoarseBeachDeposition + m_ldGTotSuspendedSediment;
   OutStream << "Total sediment deposited and in suspension (all processes) = " << ldTotalDepositedAndSuspension * m_dCellArea << " m^3" << endl;

   long double ldTotalLost = m_ldGTotActualFineSedLostBeachErosion + m_ldGTotActualSandSedLostBeachErosion + m_ldGTotActualCoarseSedLostBeachErosion + m_ldGTotSandSedLostCliffCollapse + m_ldGTotCoarseSedLostCliffCollapse;
   OutStream << "Total sediment lost from grid (all processes) = " << ldTotalLost * m_dCellArea << " m^3" << endl;
   OutStream << "                                              = " << ldTotalLost * m_dCellArea / m_dSimDuration << " m^3/hour" << endl;
   OutStream << std::fixed << setprecision(4);
   OutStream << "                                              = " << ldTotalLost * m_dCellArea / (m_dSimDuration * 3600) << " m^3/sec" << endl << endl;
   OutStream << std::fixed << setprecision(2);

   OutStream << "NOTE grid edge option is ";
   if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_CLOSED)
      OutStream << "CLOSED, therefore values above are for fine sediment only";
   else if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_OPEN)
      OutStream << "OPEN, therefore values above are for all sediment size classes";
   else if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_RECIRCULATE)
      OutStream << "RE-CIRCULATING, therefore values above are for all sediment size classes";
   OutStream << endl;

   long double
      ldLHS = ldActualTotalEroded + ldTotalLost,
      ldRHS = ldTotalDepositedAndSuspension ;
   OutStream << "Eroded + lost = " << ldLHS * m_dCellArea << " m^3" << endl;
   OutStream << "Deposited + suspension = " << ldRHS * m_dCellArea << " m^3" << endl;

   OutStream << endl << endl;

   // Finally calculate performance details
   OutStream << PERFORMHEAD << endl;

   // Get the time that the run ended
   m_tSysEndTime = time(nullptr);

   OutStream << "Run ended at " << put_time(localtime(&m_tSysEndTime), "%T on %A %d %B %Y") << endl;
   OutStream << "Time simulated: " << strDispSimTime(m_dSimDuration) << endl << endl;

   // Write to log file
   LogStream << "END OF RUN ================================================================================================" << endl << endl;

   LogStream << "ERRORS" << endl;
   LogStream << "Erosion error = " << m_ldGTotMassBalanceErosionError << " m^3" << endl;
   LogStream << "Deposition error = " << m_ldGTotMassBalanceDepositionError << " m^3" << endl;
   LogStream << endl;

   LogStream << "ALL-PROCESS TOTALS" << endl;
   LogStream << "Total sediment eroded (all processes) = " << ldActualTotalEroded * m_dCellArea << " m^3" << endl;

   LogStream << "Total sediment deposited and in suspension (all processes) = " << ldTotalDepositedAndSuspension * m_dCellArea << " m^3" << endl;

   LogStream << "Total sediment lost from grid (all processes) = " << ldTotalLost * m_dCellArea << " m^3" << endl;
   LogStream << "                                              = " << ldTotalLost * m_dCellArea / m_dSimDuration << " m^3/hour" << endl;
   LogStream << std::fixed << setprecision(4);
   LogStream << "                                              = " << ldTotalLost * m_dCellArea / (m_dSimDuration * 3600) << " m^3/sec" << endl << endl;
   LogStream << std::fixed << setprecision(2);

   LogStream << "NOTE grid edge option is ";
   if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_CLOSED)
      LogStream << "CLOSED, therefore values above are for fine sediment only";
   else if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_OPEN)
      LogStream << "OPEN, therefore values above are for all sediment size classes";
   else if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_RECIRCULATE)
      LogStream << "RE-CIRCULATING, therefore values above are for all sediment size classes";
   LogStream << endl;

   LogStream << "Eroded + lost = " << ldLHS * m_dCellArea << " m^3" << endl;
   LogStream << "Deposited + suspension = " << ldRHS * m_dCellArea << " m^3" << endl << endl;

   // Output averages for on-profile and between-profile potential shore platform erosion, ideally these are roughly equal
   LogStream << std::fixed;
   LogStream << endl;
   LogStream << "On-profile average potential shore platform erosion = " << (m_ulTotPotentialPlatformErosionOnProfiles > 0 ? m_dTotPotentialPlatformErosionOnProfiles / static_cast<double>(m_ulTotPotentialPlatformErosionOnProfiles) : 0) << " mm (n = " << m_ulTotPotentialPlatformErosionOnProfiles << ")" << endl;
   LogStream << "Between-profile average potential shore platform erosion = " << (m_ulTotPotentialPlatformErosionBetweenProfiles > 0 ? m_dTotPotentialPlatformErosionBetweenProfiles / static_cast<double>(m_ulTotPotentialPlatformErosionBetweenProfiles) : 0) << " mm (n = " << m_ulTotPotentialPlatformErosionBetweenProfiles << ")" << endl;
   LogStream << endl;

#if ! defined RANDCHECK
   // Calculate length of run, write in file (note that m_dSimDuration is in hours)
   CalcTime(m_dSimDuration * 3600);
#endif

   // Calculate statistics re. memory usage etc.
   CalcProcessStats();
   OutStream << endl << "END OF RUN" << endl;
   LogStream << endl << "END OF RUN" << endl;

   // Need to flush these here (if we don't, the buffer may not get written)
   LogStream.flush();
   OutStream.flush();

   return RTN_OK;
}


