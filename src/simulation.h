/*!
 *
 * \class CSimulation
 * \brief This class runs CoastalME simulations
 * \details TODO This is a more detailed description of the CSimulation class
 * \author David Favis-Mortlock
 * \author Andres Payo

 * \date 2017
 * \copyright GNU General Public License
 *
 * \file simulation.h
 * \brief Contains CSimulation definitions
 *
 */

#ifndef SIMULATION_H
#define SIMULATION_H
/*===============================================================================================================================

 This file is part of CoastalME, the Coastal Modelling Environment.

 CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

===============================================================================================================================*/
#include <fstream>
using std::ofstream;

#include <string>
using std::string;

#include <utility>
using std::pair;

#include <gdal_priv.h>

#include "line.h"
#include "i_line.h"


int const
   NRNG    = 2,
   SAVEMAX = 1000;

class CGeomRasterGrid;               // Forward declarations
class CRWCoast;
class CGeomProfile;
class CGeomCoastPolygon;
class CRWCliff;

class CSimulation
{
private:
   bool
      m_bBasementElevSave,
      m_bSedimentTopSurfSave,
      m_bTopSurfSave,
      m_bSliceSave,
      m_bSeaDepthSave,
      m_bAvgSeaDepthSave,
      m_bWaveHeightSave,
      m_bAvgWaveHeightSave,
      m_bWaveAngleSave,
      m_bAvgWaveAngleSave,
      m_bWaveEnergySinceCollapseSave,
      m_bMeanWaveEnergySave,
      m_bBreakingWaveHeightSave,
      m_bBeachProtectionSave,
      m_bPotentialPlatformErosionSave,
      m_bActualPlatformErosionSave,
      m_bTotalPotentialPlatformErosionSave,
      m_bTotalActualPlatformErosionSave,
      m_bPotentialBeachErosionSave,
      m_bActualBeachErosionSave,
      m_bTotalPotentialBeachErosionSave,
      m_bTotalActualBeachErosionSave,
      m_bLandformSave,
      m_bInterventionClassSave,
      m_bInterventionHeightSave,
      m_bSuspSedSave,
      m_bAvgSuspSedSave,
      m_bFineUnconsSedSave,
      m_bSandUnconsSedSave,
      m_bCoarseUnconsSedSave,
      m_bFineConsSedSave,
      m_bSandConsSedSave,
      m_bCoarseConsSedSave,
      m_bRasterCoastlineSave,
      m_bRasterNormalSave,
      m_bDistWeightSave,
      m_bActiveZoneSave,
      m_bCliffCollapseSave,
      m_bTotCliffCollapseSave,
      m_bCliffCollapseDepositionSave,
      m_bTotCliffCollapseDepositionSave,
      m_bRasterPolygonSave,
      m_bPotentialPlatformErosionMaskSave,
      m_bSeaMaskSave,
      m_bBeachMaskSave,
      m_bShadowZoneCodesSave,
      m_bSaveRegular,
      m_bCoastSave,
      m_bNormalsSave,
      m_bInvalidNormalsSave,
      m_bCoastCurvatureSave,
      m_bPolygonNodeSave,
      m_bPolygonBoundarySave,
      m_bCliffNotchSave,
      m_bShadowZoneLineSave,
      m_bSeaAreaTS,
      m_bStillWaterLevelTS,
      m_bActualPlatformErosionTS,
      m_bDepositionTS,
      m_bPotentialSedLostFromGridTS,
      m_bSuspSedTS,
      m_bSaveGISThisTimestep,
      m_bOutputProfileData,
      m_bOutputParallelProfileData,
      m_bOutputLookUpData,
      m_bOmitSearchNorthEdge,
      m_bOmitSearchSouthEdge,
      m_bOmitSearchWestEdge,
      m_bOmitSearchEastEdge,
      m_bErodeShorePlatformAlternateDirection,
      m_bDoCoastPlatformErosion,
      m_bDoCliffCollapse,
      m_bGDALCanCreate,
      m_bGDALCanWriteFloat,
      m_bGDALCanWriteInt32,
      m_bScaleRasterOutput,
      m_bWorldFile;

   char** m_papszGDALRasterOptions;
   char** m_papszGDALVectorOptions;

   int
      m_nXGridMax,
      m_nYGridMax,
      m_nLayers,
      m_nCoastSmooth,
      m_nCoastSmoothWindow,
      m_nSavGolCoastPoly,
      m_nProfileSmoothWindow,
      m_nCoastNormalAvgSpacing,              // In cells
      m_nCoastCurvatureInterval,             // A length, measured in coastline points
      m_nNaturalCapeNormals,
      m_nGISSave,
      m_nUSave,
      m_nThisSave,
      m_nCoastMax,
      m_nCoastMin,
      m_nNThisTimestepCliffCollapse,
      m_nNTotCliffCollapse,
      m_nCliffDepositionPlanviewWidth,
      m_nGlobalPolygonID,                    // There are m_nGlobalPolygonID + 1 polygons at any time (all coasts)
      m_nUnconsSedimentHandlingAtGridEdges,
      m_nBeachErosionDepositionEquation,
      m_nMissingValue;

   GDALDataType
      m_GDALWriteIntDataType,
      m_GDALWriteFloatDataType;

   long
      m_lGDALMaxCanWrite,
      m_lGDALMinCanWrite;

   unsigned long
      m_ulTimestep,
      m_ulTotTimestep,
      m_ulRandSeed[NRNG],
      m_ulNumCells,
      m_ulThisTimestepNumSeaCells,
      m_ulThisTimestepNumCoastCells,
      m_ulThisTimestepNumPotentialPlatformErosionCells,
      m_ulThisTimestepNumActualPlatformErosionCells,
      m_ulThisTimestepNumPotentialBeachErosionCells,
      m_ulThisTimestepNumActualBeachErosionCells,
      m_ulThisTimestepNumBeachDepositionCells,
      m_ulTotPotentialPlatformErosionOnProfiles,
      m_ulTotPotentialPlatformErosionBetweenProfiles;

   double
      m_dDurationUnitsMult,
      m_dNorthWestXExtCRS,
      m_dNorthWestYExtCRS,
      m_dSouthEastXExtCRS,
      m_dSouthEastYExtCRS,
      m_dExtCRSGridArea,
      m_dCellSide,                     // Length of cell side (in external CRS units)
      m_dCellArea,                     // Area of cell  (in external CRS units)
      m_dCellDiagonal,                 // Length of cell's diagonal (in external CRS units)
      m_dInvCellSide,                  // Inverse of m_dCellSide
      m_dInvCellDiagonal,              // Inverse of m_dCellDiagonal
      m_dSimDuration,                  // Duration of simulation, in hours
      m_dTimeStep,
      m_dSimElapsed,                   // Time simulated so far, in hours
      m_dRSaveTime,
      m_dRSaveInterval,
      m_dUSaveTime[SAVEMAX],
      m_dClkLast,                      // Last value returned by clock()
      m_dCPUClock,                     // Total elapsed CPU time
      m_dGeoTransform[6],
      m_nWavePropagationModel,
      m_dSeaWaterDensity,
      m_dOrigSWL,
      m_dFinalSWL,
      m_dDeltaSWLPerTimestep,
      m_dThisTimestepSWL,
      m_dMinSWL,
      m_dMaxSWL,
      m_dBreakingWaveHeight,
      m_dWavePeriod,
      m_dC_0,                          // Deep water wave speed (m/s)
      m_dL_0,                          // Deep water wave length (m)
      m_dWaveDepthRatioForWaveCalcs,
      m_dDeepWaterWaveHeight,
      m_dDeepWaterWaveOrientation,
      m_dR,
      m_dD50Fine,
      m_dD50Sand,
      m_dD50Coarse,
      m_dBeachSedimentDensity,
      m_dBeachSedimentPorosity,
      m_dFineErodibility,
      m_dSandErodibility,
      m_dCoarseErodibility,
      m_dFineErodibilityNormalized,
      m_dSandErodibilityNormalized,
      m_dCoarseErodibilityNormalized,
      m_dKLS,
      m_dKamphuis,
      m_dG,
      m_dInmersedToBulkVolumetric,
      m_dDepthOfClosure,
      m_dCoastNormalAvgSpacing,        // In m
      m_dCoastNormalLength,
      m_dThisTimestepTotSeaDepth,
      m_dThisTimestepPotentialPlatformErosion,
      m_dThisTimestepActualFinePlatformErosion,
      m_dThisTimestepActualSandPlatformErosion,
      m_dThisTimestepActualCoarsePlatformErosion,
      m_dThisTimestepPotentialBeachErosion,
      m_dThisTimestepActualFineBeachErosion,
      m_dThisTimestepActualSandBeachErosion,
      m_dThisTimestepActualCoarseBeachErosion,
      m_dThisTimestepSandBeachDeposition,
      m_dThisTimestepCoarseBeachDeposition,
      m_dThisTimestepFineSedimentToSuspension,
      m_dThisTimestepPotentialSedLostBeachErosion,
      m_dThisTimestepActualFineSedLostBeachErosion,
      m_dThisTimestepActualSandSedLostBeachErosion,
      m_dThisTimestepActualCoarseSedLostBeachErosion,
      m_dThisTimestepEstimatedActualFineBeachErosion,
      m_dThisTimestepEstimatedActualSandBeachErosion,
      m_dThisTimestepEstimatedActualCoarseBeachErosion,
      m_dThisTimestepCliffTalusFineErosion,
      m_dThisTimestepCliffTalusSandErosion,
      m_dThisTimestepCliffTalusCoarseErosion,
      m_dThisTimestepSandSedLostCliffCollapse,
      m_dThisTimestepCoarseSedLostCliffCollapse,
      m_dThisTimestepMassBalanceErosionError,
      m_dThisTimestepMassBalanceDepositionError,
      m_dDepthOverDBMax,                                    // Used in erosion potential look-up function
      m_dTotPotErosionOnProfiles,
      m_dTotPotErosionBetweenProfiles,
      m_dProfileMaxSlope,
      m_dSimpleSmoothWeight,
      m_dBeachSmoothingVertTolerance,
      m_dCliffErodibility,
      m_dNotchOverhangAtCollapse,
      m_dNotchBaseBelowSWL,
      m_dCliffDepositionA,
      m_dCliffDepositionPlanviewLength,
      m_dCliffDepositionHeightFrac,
      m_dThisTimestepCliffCollapseFine,
      m_dThisTimestepCliffCollapseSand,
      m_dThisTimestepCliffCollapseCoarse,
      m_dThisTimestepCliffTalusSandDeposition,
      m_dThisTimestepCliffTalusCoarseDeposition,
      m_dCoastNormalRandSpaceFact,
      m_dDeanProfileStartAboveSWL,
      m_dMissingValue;

   // These grand totals are all long doubles, the aim is to minimize rounding errors when many very small numbers are added to a single much larger number, see e.g. http://www.ddj.com/cpp/184403224
   long double
      m_ldGTotPotentialPlatformErosion,
      m_ldGTotFineActualPlatformErosion,
      m_ldGTotSandActualPlatformErosion,
      m_ldGTotCoarseActualPlatformErosion,
      m_ldGTotPotentialSedLostBeachErosion,
      m_ldGTotActualFineSedLostBeachErosion,
      m_ldGTotActualSandSedLostBeachErosion,
      m_ldGTotActualCoarseSedLostBeachErosion,
      m_ldGTotSandSedLostCliffCollapse,
      m_ldGTotCoarseSedLostCliffCollapse,
      m_ldGTotCliffCollapseFine,
      m_ldGTotCliffCollapseSand,
      m_ldGTotCliffCollapseCoarse,
      m_ldGTotCliffTalusSandDeposition,
      m_ldGTotCliffTalusCoarseDeposition,
      m_ldGTotCliffTalusFineErosion,
      m_ldGTotCliffTalusSandErosion,
      m_ldGTotCliffTalusCoarseErosion,
      m_ldGTotPotentialBeachErosion,
      m_ldGTotActualFineBeachErosion,
      m_ldGTotActualSandBeachErosion,
      m_ldGTotActualCoarseBeachErosion,
      m_ldGTotSandBeachDeposition,
      m_ldGTotCoarseBeachDeposition,
      m_ldGTotSuspendedSediment,
      m_ldGTotMassBalanceErosionError,
      m_ldGTotMassBalanceDepositionError;

   string
      m_strCMEDir,
      m_strCMEIni,
      m_strMailAddress,
      m_strDataPathName,
      m_strRasterGISOutFormat,
      m_strVectorGISOutFormat,
      m_strInitialBasementDEMFile,
      m_strInitialLandformFile,
      m_strInterventionClassFile,
      m_strInterventionHeightFile,
      m_strInitialSuspSedimentFile,
      m_strInitialCoastlineFile,
      m_strShapeFunctionFile,
//       m_strTideDataFile,
      m_strLogFile,
      m_strOutPath,
      m_strOutFile,
      m_strPalFile,
      m_strGDALBasementDEMDriverCode,           // Basement DEM (raster)
      m_strGDALBasementDEMDriverDesc,
      m_strGDALBasementDEMProjection,
      m_strGDALBasementDEMDataType,
      m_strGDALLDriverCode,                     // Initial Landform Class (raster)
      m_strGDALLDriverDesc,
      m_strGDALLProjection,
      m_strGDALLDataType,
      m_strGDALICDriverCode,                    // Initial Intervention Class (raster)
      m_strGDALICDriverDesc,
      m_strGDALICProjection,
      m_strGDALICDataType,
      m_strGDALIHDriverCode,                    // Initial Intervention Class (raster)
      m_strGDALIHDriverDesc,
      m_strGDALIHProjection,
      m_strGDALIHDataType,
      m_strGDALIWDriverCode,                    // Initial Water Depth (raster)
      m_strGDALIWDriverDesc,
      m_strGDALIWProjection,
      m_strGDALIWDataType,
      m_strGDALISSDriverCode,                   // Initial Suspended Sediment (raster)
      m_strGDALISSDriverDesc,
      m_strGDALISSProjection,
      m_strGDALISSDataType,
      m_strOGRICDriverCode,                     // Initial Coastline (vector)
      m_strOGRICGeometry,
      m_strOGRICDataType,
      m_strOGRICDataValue,
      m_strGDALRasterOutputDriverLongname,
      m_strGDALRasterOutputDriverExtension,
      m_strOGRVectorOutputExtension,
      m_strRunName,
      m_strDurationUnits;

   struct RandState
   {
      unsigned long s1, s2, s3;
   } m_ulRState[NRNG];

   time_t
      m_tSysStartTime,
      m_tSysEndTime;

   ofstream
      OutStream,
      SeaAreaTSStream,
      StillWaterLevelTSStream,
      ErosionTSStream,
      DepositionTSStream,
      SedLostTSStream,
      SedLoadTSStream;

   vector<bool>
      m_bConsChangedThisTimestep,
      m_bUnconsChangedThisTimestep;

   vector<int>
      m_VnProfileToSave,
      m_VnSavGolIndexCoast;            // Savitzky-Golay shift index for the coastline vector(s)

   vector<unsigned long>
      m_VulProfileTimestep;

   vector<double>
      m_VdSliceElev,
      m_VdErosionPotential,            // For erosion potential lookup
      m_VdSavGolFCRWCoast,               // Savitzky-Golay filter coefficients for the coastline vector(s)
      m_VdSavGolFCGeomProfile;             // Savitzky-Golay filter coefficients for the profile vectors
//       m_VdTideData;                    // Tide data: one record per timestep, is the change (m) from still water level for that timestep

   vector<string>
      m_VstrInitialFineUnconsSedimentFile,
      m_VstrInitialSandUnconsSedimentFile,
      m_VstrInitialCoarseUnconsSedimentFile,
      m_VstrInitialFineConsSedimentFile,
      m_VstrInitialSandConsSedimentFile,
      m_VstrInitialCoarseConsSedimentFile,
      m_VstrGDALIUFDriverCode,          // Initial Fine Unconsolidated Sediment (raster)
      m_VstrGDALIUFDriverDesc,
      m_VstrGDALIUFProjection,
      m_VstrGDALIUFDataType,
      m_VstrGDALIUSDriverCode,          // Initial Sand Unconsolidated Sediment (raster)
      m_VstrGDALIUSDriverDesc,
      m_VstrGDALIUSProjection,
      m_VstrGDALIUSDataType,
      m_VstrGDALIUCDriverCode,          // Initial Coarse Unconsolidated Sediment (raster)
      m_VstrGDALIUCDriverDesc,
      m_VstrGDALIUCProjection,
      m_VstrGDALIUCDataType,
      m_VstrGDALICFDriverCode,          // Initial Fine Consolidated Sediment (raster)
      m_VstrGDALICFDriverDesc,
      m_VstrGDALICFProjection,
      m_VstrGDALICFDataType,
      m_VstrGDALICSDriverCode,          // Initial Sand Consolidated Sediment (raster)
      m_VstrGDALICSDriverDesc,
      m_VstrGDALICSProjection,
      m_VstrGDALICSDataType,
      m_VstrGDALICCDriverCode,          // Initial Coarse Consolidated Sediment (raster)
      m_VstrGDALICCDriverDesc,
      m_VstrGDALICCProjection,
      m_VstrGDALICCDataType;

   // The raster grid object
   CGeomRasterGrid* m_pRasterGrid;

   // The coastline objects
   vector<CRWCoast> m_VCoast;

   // Pointers to coast polygon objects
   vector<CGeomCoastPolygon*> m_pVCoastPolygon;

private:
   // Input and output routines
   static int nHandleCommandLineParams(int, char* []);
   bool bReadIni(void);
   bool bReadRunData(void);
   bool bOpenLogFile(void);
   bool bSetUpTSFiles(void);
   void WriteStartRunDetails(void);
   bool bWritePerTimestepResults(void);
   bool bWriteTSFiles(void);
   int nWriteEndRunDetails(void);
   int nReadShapeFunction(void);
//    int nReadTideData(void);
   int nSaveProfile(int const, int const, int const, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<CGeom2DIPoint>* const, vector<double> const*);
   bool bWriteProfileData(int const, int const, int const, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<CGeom2DIPoint>* const, vector<double> const*);
   int nSaveParProfile(int const, int const, int const, int const, int const, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<CGeom2DIPoint>* const, vector<double> const*);
   bool bWriteParProfileData(int const, int const, int const, int const, int const, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<CGeom2DIPoint>* const, vector<double> const*);
   void WriteLookUpData(void);
   bool bWriteRasterGridWaveData2CSV(string&);

   // GIS input and output stuff
   int nReadBasementDEMData(void);
   int nReadRasterGISData(int const, int const);
//    int nReadVectorGISData(int const);        // NO LONGER USED BUT MAY BE USEFUL SOMEDAY
   bool bWriteRasterGISFloat(int const, string const*, int const = 0);
   bool bWriteRasterGISInt(int const, string const*, double const = 0);
   bool bWriteVectorGIS(int const, string const*);
   void GetRasterOutputMinMax(int const, double&, double&, int const, double const);
   void SetRasterFileCreationDefaults(void);
   int  InterpolateWavePropertiesToCellsGdal(string);
   bool  bReadWaveAttributesXYZ(string);

   // Initialization
   bool bCreateErosionPotentialLookUp(vector<double>*, vector<double>*, vector<double>*);

   // Top-level simulation routines
   static int nUpdateIntervention(void);
   int nCalcExternalForcing(void);
   int nInitGridAndCalcStillWaterLevel(void);
   int nLocateSeaAndCoasts(void);
   static int nLocateAllEstuaries(void);
   int nAssignAllCoastalLandforms(void);
   int nAssignNonCoastlineLandforms(void);
   int nDoAllPropagateWaves(void);
   int nDoAllShorePlatFormErosion(void);
   int nDoAllWaveEnergyToCoastLandforms(void);
   int nDoCliffCollapse(CRWCliff*, double const, double&, double&, double&);
   int nDoCliffCollapseDeposition(CRWCliff*, double const, double const, double const);
   int nUpdateGrid(void);

   // Lower-level simulation routines
   void FindAllSeaCells(void);
   void FloodFillSea(int const, int const);
   int nTraceCoastLine(int const, int const, int const, int const);
   int nTraceAllCoasts(void);
   void DoCoastCurvature(int const, int const);
   int nCreateAllNormalProfilesAndCheckForIntersection(void);
   int nCreateAllNormalProfiles(void);
   void CreateNaturalCapeNormals(int const, int&, int const, vector<bool>*, vector<pair<int, double> > const*);
   void CreateRestOfNormals(int const, int&, int const, double const, vector<bool>*, vector<pair<int, double> > const*);
   void CreateInterventionProfiles(int const, int&, int const);
   int nCreateNormalProfile(int const, int const, int&);
   void CreateGridEdgeProfile(bool const, int const, int&);
   int nPutAllProfilesOntoGrid(void);
   int nModifyAllIntersectingProfiles(void);
   static bool bCheckForIntersection(CGeomProfile* const, CGeomProfile* const, int&, int&, double&, double&, double&, double&);
   void MergeProfilesAtFinalLineSegments(int const, int const, int const, int const, int const, double const, double const, double const, double const);
   void TruncateOneProfileRetainOtherProfile(int const, int const, int const, double const, double const, int const, int const, bool const);
   int nInsertPointIntoProfilesIfNeededThenUpdate(int const, int const, double const, double const, int const, int const, int const, bool const);
   void TruncateProfileAndAppendNew(int const, int const, int const, vector<CGeom2DPoint> const*, vector<vector<pair<int, int> > > const*);
   void RasterizeProfile(int const, int const, vector<CGeom2DIPoint>*, vector<bool>*, bool&, bool&, bool&, bool&, bool&);
   int nRasterizeCliffCollapseProfile(vector<CGeom2DPoint> const*, vector<CGeom2DIPoint>*) const;
   int nCalcPotentialPlatformErosionOnProfile(int const, int const);
   int nCalcPotentialPlatformErosionBetweenProfiles(int const, int const, int const);
   void ConstructParallelProfile(int const, int const, int const, int const, int const, vector<CGeom2DIPoint>* const, vector<CGeom2DIPoint>*, vector<CGeom2DPoint>*);
   double dCalcBeachProtectionFactor(int const, int const, double const);
   void FillInBeachProtectionHoles(void);
   void FillPotentialPlatformErosionHoles(void);
   void DoActualShorePlatformErosionOnCell(int const, int const);
   double dLookUpErosionPotential(double const) const;
   static CGeom2DPoint PtChooseEndPoint(int const, CGeom2DPoint const*, CGeom2DPoint const*, double const, double const, double const, double const);
   int nGetCoastNormalEndPoint(int const, int const, int const, CGeom2DPoint const*, double const, CGeom2DPoint*);
   int nLandformToGrid(int const, int const);
   void CalcWavePropertiesOnProfile(int const, int const, int const); 
   int GetThisProfilePointsElevationVectors(int const,int const, int const&, vector<double>&, vector<double>&); // used to run CShore
   int CreateCShoreinfile(double const, double const, double const, double const , double const, double const, vector<double>&, vector<double>&); // used to run CShore
   int dLookUpCShoreouputs(string const, unsigned const, unsigned const, vector<double> const , vector<double>& vdMyInterpolatedValues); // used to run CShore
   double dCalcWaveAngleToCoastNormal(double const, int const);
   void CalcCoastTangents(int const);
   void InterpolateWavePropertiesToCoastline(int const, int const, int const);
   void InterpolateWavePropertiesToCells(int const, int const, int const);
   static double dCalcCurvature(int const, CGeom2DPoint const*, CGeom2DPoint const*, CGeom2DPoint const*);
   void CalcD50AndFillWaveCalcHoles(void);
   int nDoAllShadowZones(void);
   int nFindAllShadowZones(void);   
   int nFloodFillShadowZone(int const, int const, CGeom2DIPoint const*, int const, CGeom2DIPoint const*);
   int nSweepShadowZone(int const, int const, CGeom2DIPoint const*, int const, CGeom2DIPoint const*, int&);
   int nSweepDownDriftFromShadowZone(int const, int const, CGeom2DIPoint const*, int const, CGeom2DIPoint const*, int const);
   int nCreateAllPolygons(void);
   void RasterizePolygonJoiningLine(CGeom2DPoint const*, CGeom2DPoint const*);
   static bool bIsWithinPolygon(CGeom2DPoint const*, vector<CGeom2DPoint> const*);
   static CGeom2DPoint PtFindPointInPolygon(vector<CGeom2DPoint> const*, int const);
   void MarkPolygonCells(void);
   void DoPolygonSharedBoundaries(void);
   void DoAllPotentialBeachErosion(void);
   int nDoAllActualBeachErosionAndDeposition(void);
   int nEstimateActualBeachErosionOnPolygon(int const, int const, double const);
   void EstimateActualBeachErosionOnCell(int const, int const, int const, double const, double&, double&, double&);
   void ErodeBeachConstrained(int const, int const, int const, double const, double&, double&, double&);
   int nRouteActualBeachErosionToAdjacentPolygons(int const, int const);
   int nDoWithinPolygonBeachRedistribution(int const, int const);
   int nDoBeachErosionOnCells(int const, int const, double const);
   int nDoBeachDepositionOnCells(int const, int const, double const);

   // GIS utility routines
   bool bCheckRasterGISOutputFormat(void);
   bool bCheckVectorGISOutputFormat(void);
   bool bSaveAllRasterGISFiles(void);
   bool bSaveAllVectorGISFiles(void);
   bool bIsWithinGrid(int const, int const) const;
   bool bIsWithinGrid(CGeom2DIPoint const*) const;
   double dGridCentroidXToExtCRSX(int const) const;
   double dGridCentroidYToExtCRSY(int const) const;
   double dGridXToExtCRSX(double const) const;
   double dGridYToExtCRSY(double const) const;
   double dExtCRSXToGridCentroidX(double const) const;
   double dExtCRSYToGridCentroidY(double const) const;
   double dExtCRSXToGridX(double const) const;
   double dExtCRSYToGridY(double const) const;
   int nExtCRSXToGridX(double const) const;
   int nExtCRSYToGridY(double const) const;
   CGeom2DIPoint PtiExtCRSToGrid(CGeom2DPoint const*) const;
   CGeom2DPoint PtGridCentroidToExt(CGeom2DIPoint const*) const;
   static double dGetDistanceBetween(CGeom2DPoint const*, CGeom2DPoint const*);
   static double dGetDistanceBetween(CGeom2DIPoint const*, CGeom2DIPoint const*);
   static double dTriangleAreax2(CGeom2DPoint const*, CGeom2DPoint const*, CGeom2DPoint const*);
   void KeepWithinGrid(int const, int const, int&, int&) const;
   void KeepWithinGrid(CGeom2DIPoint const*, CGeom2DIPoint*) const;
   static double dKeepWithin360(double const);
//    vector<CGeom2DPoint> VGetPerpendicular(CGeom2DPoint const*, CGeom2DPoint const*, double const, int const);
   static CGeom2DPoint PtGetPerpendicular(CGeom2DPoint const*, CGeom2DPoint const*, double const, int const);
   static CGeom2DIPoint PtiGetPerpendicular(CGeom2DIPoint const*, CGeom2DIPoint const*, double const, int const);
   static CGeom2DPoint PtAverage(CGeom2DPoint const*, CGeom2DPoint const*);
   static CGeom2DPoint PtAverage(vector<CGeom2DPoint>*);
   static CGeom2DIPoint PtiAverage(CGeom2DIPoint const*, CGeom2DIPoint const*);
   static CGeom2DIPoint PtiWeightedAverage(CGeom2DIPoint const*, CGeom2DIPoint const*, double const);
   static double dAngleSubtended(CGeom2DIPoint const*, CGeom2DIPoint const*, CGeom2DIPoint const*);

   // Utility routines
   static void AnnounceStart(void);
   void AnnounceLicence(void);
   void AnnounceReadBasementDEM(void) const;
   static void AnnounceAddLayers(void);
   static void AnnounceReadRasterFiles(void);
//    static void AnnounceReadVectorFiles(void);
   void AnnounceReadLGIS(void) const;
   void AnnounceReadICGIS(void) const;
   void AnnounceReadIHGIS(void) const;
   static void AnnounceInitializing(void);
   void AnnounceReadInitialSuspSedGIS(void) const;
   void AnnounceReadInitialFineUnconsSedGIS(int const) const;
   void AnnounceReadInitialSandUnconsSedGIS(int const) const;
   void AnnounceReadInitialCoarseUnconsSedGIS(int const) const;
   void AnnounceReadInitialFineConsSedGIS(int const) const;
   void AnnounceReadInitialSandConsSedGIS(int const) const;
   void AnnounceReadInitialCoarseConsSedGIS(int const) const;
//    void AnnounceReadTideData(void) const;
   static void AnnounceReadShapeFunctionFile(void);
   static void AnnounceAllocateMemory(void);
   static void AnnounceIsRunning(void);
   static void AnnounceSimEnd(void);
   void StartClock(void);
   bool bFindExeDir(char* pcArg);
   bool bTimeToQuit(void);
   static int nDoTimeUnits(string const*);
   int nDoSimulationTimeMultiplier(string const*);
   static double dGetTimeMultiplier(string const*);
   void UpdateGrandTotals(void);
   static string strGetBuild(void);
   static string strGetComputerName(void);
   void DoCPUClockReset(void);
   void CalcTime(double const);
   static string strDispTime(double const, bool const, bool const);
   static string strDispSimTime(double const);
   void AnnounceProgress(void);
   static string strGetErrorText(int const);
   string strListRasterFiles(void) const;
   string strListVectorFiles(void) const;
   string strListTSFiles(void) const;
   void CalcProcessStats(void);
   void CalcSavitzkyGolayCoeffs(void);
   CGeomLine LSmoothCoastSavitzkyGolay(CGeomLine*, int const, int const) const;
   CGeomLine LSmoothCoastRunningMean(CGeomLine*, int const, int const) const;
   vector<double> dVSmoothProfileSlope(vector<double>*);
//    vector<double> dVCalCGeomProfileSlope(vector<CGeom2DPoint>*, vector<double>*);
   vector<double> dVSmoothProfileSavitzkyGolay(vector<double>*, vector<double>*);
   vector<double> dVSmoothProfileRunningMean(vector<double>*);
   void CalcSavitzkyGolay(double[], int const, int const, int const, int const, int const);

   // Random number stuff
   static unsigned long ulGetTausworthe(unsigned long const, unsigned long const, unsigned long const, unsigned long const, unsigned long const);
   void InitRand0(unsigned long const);
   void InitRand1(unsigned long const);
   unsigned long ulGetRand0(void);
   unsigned long ulGetRand1(void);
   static unsigned long ulGetLCG(unsigned long const);            // Used by all generators
   double dGetRand0d1(void);
//    int nGetRand0To(int const);
   int nGetRand1To(int const);
//    double dGetRand0GaussPos(double const, double const);
   double dGetRand0Gaussian(void);
//    double dGetCGaussianPDF(double const);
   void Rand1Shuffle(int*, int);
#ifdef RANDCHECK
   void CheckRand(void) const;
#endif

public:
   ofstream LogStream;

   CSimulation(void);
   ~CSimulation(void);

   //! Returns this timestep's still water level
   double dGetThisTimestepSWL(void) const;

   //! Returns the vertical tolerance for beach cells to be included in smoothing
   double dGetBeachSmoothingVertTolerance(void) const;

   //! Returns the cell size
//    double dGetCellSide(void) const;

   //! Returns the size of the grid in the X direction
   int nGetGridXMax(void) const;

   //! Returns the size of the grid in the Y direction
   int nGetGridYMax(void) const;

   //! Returns the global d50 value for fine sediment
   double dGetD50Fine(void) const;

   //! Returns the global d50 value for sand sediment
   double dGetD50Sand(void) const;

   //! Returns the global d50 value for coarse sediment
   double dGetD50Coarse(void) const;

   //! Runs the simulation
   int nDoSimulation(int, char*[]);

   //! Carries out end-of-simulation tidying (error messages etc.)
   void DoSimulationEnd(int const);
};
#endif // SIMULATION_H
