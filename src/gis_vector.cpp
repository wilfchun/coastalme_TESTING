/*!
 *
 * \file gis_vector.cpp
 * \brief These functions use GDAL to read and write vector GIS files in several formats. This version will build with GDAL version 2
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
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::ios;

#include <sstream>
using std::stringstream;

#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include "cme.h"
#include "simulation.h"
#include "coast.h"
#include "cliff.h"


/*==============================================================================================================================

 Reads vector GIS datafiles NOT NOW USED BUT MAY BE SOMEDAY

===============================================================================================================================*/
// int CSimulation::nReadVectorGISData(int const nDataItem)
// {
//    int
//       nMaxLayer = 0,
//       nNeedFieldType = 0,
//       nNeedGeometry = 0;
//    string
//       strGISFile,
//       strDriverCode,
//       strGeometry,
//       strDataType,
//       strDataValue;
//
//    // Set up file name and constraints
//    switch (nDataItem)
//    {
//       case (COAST_VEC):
//          strGISFile = m_strInitialCoastlineFile;
//          nMaxLayer = COAST_VEC_MAX_LAYER;
//          nNeedFieldType = COAST_VEC_FIELD_DATA_TYPE,
//          nNeedGeometry = COAST_VEC_GEOMETRY;
//          break;
//
//       // TODO Others
//    }
//
//    // Open the OGR datasource
//    OGRDataSource* pOGRDataSource = NULL;
//    pOGRDataSource = OGRSFDriverRegistrar::Open(strGISFile.c_str(), FALSE);
//    if (pOGRDataSource == NULL)
//    {
//       // Can't open file (note will already have sent GDAL error message to stdout)
//       cerr << ERR << "cannot open " << strGISFile << " for input: " << CPLGetLastErrorMsg() << endl;
//       return RTN_ERR_VECTOR_FILE_READ;
//    }
//
//    // Opened OK, so get dataset information
//    strDriverCode = pOGRDataSource->GetDriver()->GetName();
//
//    // Find out number of layers, and compare with the required number
//    int nLayer = pOGRDataSource->GetLayerCount();
//    if (nLayer > nMaxLayer)
//       LogStream << WARN << "need " << nMaxLayer << (nMaxLayer > 1 ? "layers" : "layer") << " in " << strGISFile << ", " << nLayer << " found. Only the first " << nMaxLayer << (nMaxLayer > 1 ? "layers" : "layer") << " will be read." << endl;
//
//    for (int n = 0; n < nMaxLayer; n++)
//    {
//       // Open this layer
//       OGRLayer* pOGRLayer;
//       pOGRLayer = pOGRDataSource->GetLayer(n);
//
//       // Get features from the layer
//       OGRFeature* pOGRFeature;
//
//       // Make sure we are at the beginning of the layer, then iterate through all features in the layer
//       pOGRLayer->ResetReading();
//       while ((pOGRFeature = pOGRLayer->GetNextFeature()) != NULL)
//       {
//          OGRFeatureDefn* pOGRFeatureDefn = pOGRLayer->GetLayerDefn();
//          for (int nField = 0; nField < pOGRFeatureDefn->GetFieldCount(); nField++)
//          {
//             OGRFieldDefn* pOGRFieldDefn = pOGRFeatureDefn->GetFieldDefn(nField);
//
//             int nFieldType = pOGRFieldDefn->GetType();
//             int nThisFieldType = 0;
//             switch (nFieldType)
//             {
//                case OFTInteger:
//                   nThisFieldType = VEC_FIELD_DATA_INT;
//                   strDataType = "integer";
//                   break;
//
//                case OFTReal:
//                   nThisFieldType = VEC_FIELD_DATA_REAL;
//                   strDataType = "real";
//                   break;
//
//                case OFTString:
//                   nThisFieldType = VEC_FIELD_DATA_STRING;
//                   strDataType = "string";
//                   break;
//
//                default:
//                   nThisFieldType = VEC_FIELD_DATA_OTHER;
//                   strDataType = "other";
//                   break;
//             }
//
//             // Check whether we have the expected field data type
//             if (nNeedFieldType != VEC_FIELD_DATA_ANY)
//             {
//                if (nThisFieldType != nNeedFieldType)
//                {
//                   // Error: we have not got the expected field type
//                   string strNeedType;
//                   switch (nNeedFieldType)
//                   {
//                      case VEC_FIELD_DATA_INT:
//                         strNeedType = "integer";
//                         break;
//
//                      case VEC_FIELD_DATA_REAL:
//                         strNeedType = "real";
//                         break;
//
//                      case VEC_FIELD_DATA_STRING:
//                         strNeedType = "string";
//                         break;
//
//                      case VEC_FIELD_DATA_OTHER:
//                         strNeedType = "other";
//                         break;
//                   }
//
//                   cerr << ERR << strDataType << " field data found in " << strGISFile << ", but " << strNeedType << " field data needed." << endl;
//                   return RTN_ERR_VECTOR_FILE_READ;
//                }
//             }
//
//             // OK we have the desired field data type, so get the value
//             // TODO WILL WE EVER ACTUALLY USE THIS VALUE? May as well just get it as a string for later display
//             strDataValue = pOGRFeature->GetFieldAsString(nField);
// /*          switch (nFieldType)
//             {
//                case OFTInteger:
//                   pOGRFeature->GetFieldAsInteger(nField);
//                   break;
//
//                case OFTReal:
//                   pOGRFeature->GetFieldAsDouble(nField);
//                   break;
//
//                case OFTString:
//                   pOGRFeature->GetFieldAsString(nField);
//                   break;
//
//                default:
//                   pOGRFeature->GetFieldAsString(nField);
//                   break;
//             } */
//          }
//
//          // Now get the geometry
//          OGRGeometry* pOGRGeometry;
//          pOGRGeometry = pOGRFeature->GetGeometryRef();
//          if (pOGRGeometry == NULL)
//          {
//             cerr << ERR << " null geometry in " << strGISFile << "." << endl;
//             return RTN_ERR_VECTOR_FILE_READ;
//          }
//
//          // And then get the geometry type
//          int nGeometry = wkbFlatten(pOGRGeometry->getGeometryType());
//          int nThisGeometry = 0;
//          switch (nGeometry)
//          {
//             case wkbPoint:
//                nThisGeometry = VEC_GEOMETRY_POINT;
//                strGeometry = "point";
//                break;
//
//             case wkbLineString:
//                nThisGeometry = VEC_GEOMETRY_LINE;
//                strGeometry = "line";
//                break;
//
//             case wkbPolygon:
//                nThisGeometry = VEC_GEOMETRY_POLYGON;
//                strGeometry = "polygon";
//                break;
//
//             default:
//                // NOTE may need wkbMultiLineString or similar for channel network
//                nThisGeometry = VEC_GEOMETRY_OTHER;
//                strGeometry = "other";
//                break;
//          }
//
//          // Have we got the expected geometry?
//          if (nThisGeometry != nNeedGeometry)
//          {
//             // Error, we do not have the desired geometry
//             string strNeedGeometry;
//             switch (nNeedGeometry)
//             {
//                case VEC_FIELD_DATA_INT:
//                   strNeedGeometry = "integer";
//                   break;
//
//                case VEC_FIELD_DATA_REAL:
//                   strNeedGeometry = "real";
//                   break;
//
//                case VEC_FIELD_DATA_STRING:
//                   strNeedGeometry = "string";
//                   break;
//
//                case VEC_FIELD_DATA_OTHER:
//                   strNeedGeometry = "other";
//                   break;
//             }
//
//             cerr << strGeometry << " data found in " << strGISFile << ", but " << strNeedGeometry << " data needed." << endl;
//             return RTN_ERR_VECTOR_FILE_READ;
//          }
//
//          // The geometry is OK, so (at last) process the data
//          int nPoints = 0;
//          OGRPoint* pOGRPoint;
//          OGRLineString* pOGRLineString;
//          switch (nGeometry)
//          {
//             case wkbPoint:
//                // Point data
//                pOGRPoint = (OGRPoint *) pOGRGeometry;
//                // TODO write the rest of this
//                cout << "Point: x = " << pOGRPoint->getX() << ", y = " << pOGRPoint->getY() << endl;
//                break;
//
//             case wkbLineString:
//                // Line data
//                pOGRLineString = (OGRLineString *) pOGRGeometry;
//                nPoints = pOGRLineString->getNumPoints();
//                for (int i = 0; i < nPoints; i++)
//                {
//                   switch (nDataItem)
//                   {
//                      case (COAST_VEC):
//                         // Add this point to the coastline. Note that we are assuming only one coastline object here
//                         m_VCoast[0].AppendToCoastline(pOGRLineString->getX(i), pOGRLineString->getY(i));
//                         break;
//
//                      // TODO others
//                   }
//                }
//                break;
//
//             case wkbPolygon:
//                // Polygon data
//                // TODO write this
//                break;
//
//             default:
//                // NOTE will need wkbMultiLineString or similar for channel network
//                // TODO write this
//                break;
//          }
//
//          // Get rid of the Feature object
//          OGRFeature::DestroyFeature(pOGRFeature);
//
//          // Pass on some info to show in the text output
//          switch (nDataItem)
//          {
//             case (COAST_VEC):
//                m_strOGRICDriverCode = strDriverCode;
//                m_strOGRICDataType   = strDataType;
//                m_strOGRICDataValue  = strDataValue;
//                m_strOGRICGeometry   = strGeometry;
//                break;
//
//             // TODO Others
//
//          }
//       }
//    }
//
//    // Clean up: get rid of the data source object
//    OGRDataSource::DestroyDataSource(pOGRDataSource);
//
//    return RTN_OK;
// }


/*==============================================================================================================================

 Writes vector GIS files using OGR

===============================================================================================================================*/
bool CSimulation::bWriteVectorGIS(int const nDataItem, string const* strPlotTitle)
{
   // Begin constructing the file name for this save
   string strFilePathName(m_strOutPath);

   switch (nDataItem)
   {
      case (PLOT_COAST):
      {
         strFilePathName.append(VECTOR_COAST_NAME);
         break;
      }

      case (PLOT_NORMALS):
      {
         strFilePathName.append(VECTOR_NORMALS_NAME);
         break;
      }

      case (PLOT_INVALID_NORMALS):
      {
         strFilePathName.append(VECTOR_INVALID_NORMALS_NAME);
         break;
      }

      case (PLOT_COAST_CURVATURE):
      {
         strFilePathName.append(VECTOR_COAST_CURVATURE_NAME);
         break;
      }

      case (PLOT_WAVE_AND_HEIGHT):
      {
         strFilePathName.append(VECTOR_WAVE_ANGLE_NAME);
         break;
      }

      case (PLOT_AVG_WAVE_AND_HEIGHT):
      {
         strFilePathName.append(VECTOR_AVG_WAVE_ANGLE_NAME);
         break;
      }

      case (PLOT_WAVE_ENERGY_SINCE_COLLAPSE):
      {
         strFilePathName.append(VECTOR_WAVE_ENERGY_SINCE_COLLAPSE_NAME);
         break;
      }

      case (PLOT_MEAN_WAVE_ENERGY):
      {
         strFilePathName.append(VECTOR_MEAN_WAVE_ENERGY_NAME);
         break;
      }

      case (PLOT_BREAKING_WAVE_HEIGHT):
      {
         strFilePathName.append(VECTOR_BREAKING_WAVE_HEIGHT_NAME);
         break;
      }

      case (PLOT_POLYGON_NODES):
      {
         strFilePathName.append(VECTOR_POLYGON_NODES_NAME);
         break;
      }

      case (PLOT_POLYGON_BOUNDARY):
      {
         strFilePathName.append(VECTOR_POLYGON_BOUNDARY_NAME);
         break;
      }

      case (PLOT_CLIFF_NOTCH_SIZE):
      {
         strFilePathName.append(VECTOR_CLIFF_NOTCH_SIZE_NAME);
         break;
      }
      
      case (PLOT_SHADOW_ZONE_BOUNDARY):
      {
         strFilePathName.append(VECTOR_SHADOW_ZONE_LINE_NAME);
         break;
      }
   }

   // Append the 'save number' to the filename, and prepend zeros to the save number
   strFilePathName.append("_");
   stringstream ststrTmp;
   ststrTmp << FillToWidth('0', MAX_SAVE_DIGITS) << m_nGISSave;
   strFilePathName.append(ststrTmp.str());

   // Make a copy of the filename without any extension
   string strFilePathNameNoExt = strFilePathName;

   // If desired, append an extension
   if (! m_strOGRVectorOutputExtension.empty())
      strFilePathName.append(m_strOGRVectorOutputExtension);

   // Set up the vector driver
   GDALDriver* pGDALDriver = GetGDALDriverManager()->GetDriverByName(m_strVectorGISOutFormat.c_str());
   if (pGDALDriver == NULL)
   {
      cerr << ERR << "vector GIS output driver " << m_strVectorGISOutFormat << CPLGetLastErrorMsg() << endl;
      return false;
   }

   // Now create the dataset
   GDALDataset* pGDALDataSet = NULL;
   pGDALDataSet = pGDALDriver->Create(strFilePathName.c_str(), 0, 0, 0, GDT_Unknown, m_papszGDALVectorOptions);
   if (pGDALDataSet == NULL)
   {
      cerr << ERR << "cannot create " << m_strVectorGISOutFormat << " named " << strFilePathName << "\n" << CPLGetLastErrorMsg() << endl;
      return false;
   }

   // Create the output layer
   OGRLayer* pOGRLayer = NULL;
   OGRSpatialReference* pOGRSpatialRef = NULL;     // TODO add spatial reference
   OGRwkbGeometryType eGType = wkbUnknown;
   string strType = "unknown";

   pOGRLayer = pGDALDataSet->CreateLayer(strFilePathNameNoExt.c_str(), pOGRSpatialRef, eGType, m_papszGDALVectorOptions);
   if (pOGRLayer == NULL)
   {
      cerr << ERR << "cannot create '" << strType << "' layer in " << strFilePathName << "\n" << CPLGetLastErrorMsg() << endl;
      return false;
   }

   switch (nDataItem)
   {
      case (PLOT_COAST):
      {
         eGType = wkbLineString;
         strType = "line";

         // The layer has been created, so create an integer-numbered value (the number of the coast object) for the multi-line
         string strFieldValue1 = "Coast";
         OGRFieldDefn OGRField1(strFieldValue1.c_str(), OFTInteger);
         if (pOGRLayer->CreateField(&OGRField1) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 1 '" << strFieldValue1 << "' in " << strFilePathName << "\n" << CPLGetLastErrorMsg() << endl;
            return false;
         }

         // OK, now do features
         OGRLineString OGRls;

         for (int i = 0; i < static_cast<int>(m_VCoast.size()); i++)
         {
            // Create a feature object, one per coast
            OGRFeature *pOGRFeature = NULL;
            pOGRFeature = OGRFeature::CreateFeature(pOGRLayer->GetLayerDefn());

            // Set the feature's attribute (the coast number)
            pOGRFeature->SetField(strFieldValue1.c_str(), i);

            // Now attach a geometry to the feature object
            for (int j = 0; j < m_VCoast[i].pLGetCoastline()->nGetSize(); j++)
               //  In external CRS
               OGRls.addPoint(m_VCoast[i].pPtGetVectorCoastlinePoint(j)->dGetX(), m_VCoast[i].pPtGetVectorCoastlinePoint(j)->dGetY());

            pOGRFeature->SetGeometry(&OGRls);

            // Create the feature in the output layer
            if (pOGRLayer->CreateFeature(pOGRFeature) != OGRERR_NONE)
            {
               cerr << ERR << "cannot create  " << strType << " feature " << strPlotTitle << " for coast " << i << " in " << strFilePathName << "\n" << CPLGetLastErrorMsg() << endl;
               return false;
            }

            // Tidy up: empty the line string and get rid of the feature object
            OGRls.empty();
            OGRFeature::DestroyFeature(pOGRFeature);
         }

         break;
      }

      case (PLOT_NORMALS):
      case (PLOT_INVALID_NORMALS):
      {
         eGType = wkbLineString;
         strType = "line";

         // The layer has been created, so create an integer-numbered value (the number of the normal) associated with the line
         string strFieldValue1 = "Normal";
         OGRFieldDefn OGRField1(strFieldValue1.c_str(), OFTInteger);
         if (pOGRLayer->CreateField(&OGRField1) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 1 '" << strFieldValue1 << "' in " << strFilePathName << "\n" << CPLGetLastErrorMsg() << endl;
            return false;
         }

         // Also create other integer-numbered values for the category codes of the coastline-normalprofile
         string
            strFieldValue2 = "StartCoast",
            strFieldValue3 = "EndCoast",
            strFieldValue4 = "HitLand",
            strFieldValue5 = "HitCoast",
            strFieldValue6 = "HitNormal";
         OGRFieldDefn
            OGRField2(strFieldValue2.c_str(), OFTInteger),
            OGRField3(strFieldValue3.c_str(), OFTInteger),
            OGRField4(strFieldValue4.c_str(), OFTInteger),
            OGRField5(strFieldValue5.c_str(), OFTInteger),
            OGRField6(strFieldValue6.c_str(), OFTInteger);
         if (pOGRLayer->CreateField(&OGRField2) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 2 '" << strFieldValue2 << "' in " << strFilePathName << "\n" << CPLGetLastErrorMsg() << endl;
            return false;
         }
         if (pOGRLayer->CreateField(&OGRField3) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 3 '" << strFieldValue3 << "' in " << strFilePathName << "\n" << CPLGetLastErrorMsg() << endl;
            return false;
         }
         if (pOGRLayer->CreateField(&OGRField4) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 4 '" << strFieldValue4 << "' in " << strFilePathName << "\n" << CPLGetLastErrorMsg() << endl;
            return false;
         }
         if (pOGRLayer->CreateField(&OGRField5) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 5 '" << strFieldValue5 << "' in " << strFilePathName << "\n" << CPLGetLastErrorMsg() << endl;
            return false;
         }
         if (pOGRLayer->CreateField(&OGRField6) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 6 '" << strFieldValue6 << "' in " << strFilePathName << "\n" << CPLGetLastErrorMsg() << endl;
            return false;
         }

         // OK, now create features
         OGRLineString OGRls;

         for (int i = 0; i < static_cast<int>(m_VCoast.size()); i++)
         {
            for (int j = 0; j < m_VCoast[i].nGetNumProfiles(); j++)
            {
               CGeomProfile* pProfile = m_VCoast[i].pGetProfile(j);

               if (((nDataItem == PLOT_NORMALS) && (pProfile->bOKIncStartAndEndOfCoast())) || ((nDataItem == PLOT_INVALID_NORMALS) && (! pProfile->bOKIncStartAndEndOfCoast())))
               {
                  // Create a feature object, one per profile
                  OGRFeature *pOGRFeature = NULL;
                  pOGRFeature = OGRFeature::CreateFeature(pOGRLayer->GetLayerDefn());

                  // Set the feature's attributes
                  pOGRFeature->SetField(strFieldValue1.c_str(), j);
                  pOGRFeature->SetField(strFieldValue2.c_str(), 0);
                  pOGRFeature->SetField(strFieldValue3.c_str(), 0);
                  pOGRFeature->SetField(strFieldValue4.c_str(), 0);
                  pOGRFeature->SetField(strFieldValue5.c_str(), 0);
                  pOGRFeature->SetField(strFieldValue6.c_str(), 0);
                  if (pProfile->bStartOfCoast())
                     pOGRFeature->SetField(strFieldValue2.c_str(), 1);
                  if (pProfile->bEndOfCoast())
                     pOGRFeature->SetField(strFieldValue3.c_str(), 1);
                  if (pProfile->bHitLand())
                     pOGRFeature->SetField(strFieldValue4.c_str(), 1);
                  if (pProfile->bHitCoast())
                     pOGRFeature->SetField(strFieldValue5.c_str(), 1);
                  if (pProfile->bHitAnotherProfile())
                     pOGRFeature->SetField(strFieldValue6.c_str(), 1);

                  // Now attach a geometry to the feature object
                  for (int k = 0; k < pProfile->nGetProfileSize(); k++)
                     OGRls.addPoint(pProfile->pPtGetPointInProfile(k)->dGetX(), pProfile->pPtGetPointInProfile(k)->dGetY());

                  pOGRFeature->SetGeometry(&OGRls);
                  OGRls.empty();

                  // Create the feature in the output layer
                  if (pOGRLayer->CreateFeature(pOGRFeature) != OGRERR_NONE)
                  {
                     cerr << ERR << "cannot create  " << strType << " feature " << strPlotTitle << " for coast " << i << " and profile " << j << " in " << strFilePathName << "\n" << CPLGetLastErrorMsg() << endl;
                     return false;
                  }

                  // Tidy up: get rid of the feature object
                  OGRFeature::DestroyFeature(pOGRFeature);
               }
            }
         }

         break;
      }

      case (PLOT_COAST_CURVATURE):
      case (PLOT_WAVE_ENERGY_SINCE_COLLAPSE):
      case (PLOT_MEAN_WAVE_ENERGY):
      case (PLOT_BREAKING_WAVE_HEIGHT):
      case (PLOT_POLYGON_NODES):
      case (PLOT_CLIFF_NOTCH_SIZE):
      {
         eGType = wkbPoint;
         strType = "point";

         // The layer has been created, so create a real-numbered value associated with each point
         string strFieldValue1;
         if (nDataItem == PLOT_COAST_CURVATURE)
            strFieldValue1 = "Curve";
         else if ((nDataItem == PLOT_WAVE_ENERGY_SINCE_COLLAPSE) || (nDataItem == PLOT_MEAN_WAVE_ENERGY))
            strFieldValue1 = "Energy";
         else if (nDataItem == PLOT_BREAKING_WAVE_HEIGHT)
            strFieldValue1 = "Height";
         else if (nDataItem == PLOT_POLYGON_NODES)
            strFieldValue1 = "Node";
         else if (nDataItem == PLOT_CLIFF_NOTCH_SIZE)
            strFieldValue1 = "Notch";

         OGRFieldDefn OGRField1(strFieldValue1.c_str(), OFTReal);
         if (pOGRLayer->CreateField(&OGRField1) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 1 '" << strFieldValue1 << "' in " << strFilePathName << "\n" << CPLGetLastErrorMsg() << endl;
            return false;
         }

         // OK, now create features
         OGRLineString OGRls;
         OGRMultiLineString OGRmls;
         OGRPoint OGRPt;

         for (int i = 0; i < static_cast<int>(m_VCoast.size()); i++)
         {
            for (int j = 0; j < m_VCoast[i].pLGetCoastline()->nGetSize(); j++)
            {
               // Create a feature object, one per coastline point
               OGRFeature *pOGRFeature = NULL;
               pOGRFeature = OGRFeature::CreateFeature(pOGRLayer->GetLayerDefn());

               // Set the feature's geometry (in external CRS)
               OGRPt.setX(m_VCoast[i].pPtGetVectorCoastlinePoint(j)->dGetX());
               OGRPt.setY(m_VCoast[i].pPtGetVectorCoastlinePoint(j)->dGetY());
               pOGRFeature->SetGeometry(&OGRPt);

               if (nDataItem == PLOT_COAST_CURVATURE)
               {
                  double dCurvature = m_VCoast[i].dGetDetailedCurvature(j);
                  if (dCurvature == DBL_NODATA)
                        continue;

                  // Set the feature's attribute
                  pOGRFeature->SetField(strFieldValue1.c_str(), dCurvature);
               }
               else if (nDataItem == PLOT_WAVE_ENERGY_SINCE_COLLAPSE)
               {
                  // Set the feature's attribute
                  pOGRFeature->SetField(strFieldValue1.c_str(), m_VCoast[i].dGetWaveEnergy(j));
               }
               else if (nDataItem == PLOT_MEAN_WAVE_ENERGY)
               {
                  // Set the feature's attribute
                  double dEnergy = m_VCoast[i].pGetCoastLandform(j)->dGetTotAccumWaveEnergy();
                  dEnergy *= 24;
                  dEnergy /= m_dSimElapsed;     // Is in energy units per day

                  pOGRFeature->SetField(strFieldValue1.c_str(), dEnergy);
               }
               else if (nDataItem == PLOT_BREAKING_WAVE_HEIGHT)
               {
                  // Set the feature's attribute
                  double dHeight = m_VCoast[i].dGetBreakingWaveHeight(j);
                  pOGRFeature->SetField(strFieldValue1.c_str(), dHeight);
               }
               else if (nDataItem == PLOT_POLYGON_NODES)
               {
                  int nNode = m_VCoast[i].nGetPolygonNode(j);
                  if (nNode == INT_NODATA)
                     continue;

                  // Set the feature's attribute
                  pOGRFeature->SetField(strFieldValue1.c_str(), nNode);
               }
               else if (nDataItem == PLOT_CLIFF_NOTCH_SIZE)
               {
                  CACoastLandform* pCoastLandform = m_VCoast[i].pGetCoastLandform(j);
                  int nCategory = pCoastLandform->nGetLandFormCategory();
                  double dNotchOverhang = DBL_NODATA;

                  if (nCategory == LF_CAT_CLIFF)
                  {
                     CRWCliff* pCliff = reinterpret_cast<CRWCliff*>(pCoastLandform);

                     // Get attribute values from the cliff object
                     dNotchOverhang = pCliff->dGetNotchOverhang();
                  }

                  // Set the feature's attribute
                  pOGRFeature->SetField(strFieldValue1.c_str(), dNotchOverhang);
               }

               // Create the feature in the output layer
               if (pOGRLayer->CreateFeature(pOGRFeature) != OGRERR_NONE)
               {
                  cerr << ERR << "cannot create " << strType << " feature " << strPlotTitle << " for coast " << i << " point " << j << " in " << strFilePathName << "\n" << CPLGetLastErrorMsg() << endl;
                  return false;
               }

               // Get rid of the feature object
               OGRFeature::DestroyFeature(pOGRFeature);
            }
         }

         break;
      }

      case (PLOT_WAVE_AND_HEIGHT):
      {
         eGType = wkbPoint;
         strType = "point";

         // The layer has been created, so create real-numbered values associated with each point
         string
            strFieldValue1 = "Angle",
            strFieldValue2 = "Height";

         // Create the first field
         OGRFieldDefn OGRField1(strFieldValue1.c_str(), OFTReal);
         if (pOGRLayer->CreateField(&OGRField1) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 1 '" << strFieldValue1 << "' in " << strFilePathName << "\n" << CPLGetLastErrorMsg() << endl;
            return false;
         }

         // Create the second field
         OGRFieldDefn OGRField2(strFieldValue2.c_str(), OFTReal);
         if (pOGRLayer->CreateField(&OGRField2) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 2 '" << strFieldValue2 << "' in " << strFilePathName << "\n" << CPLGetLastErrorMsg() << endl;
            return false;
         }

         // OK, now create features
         OGRLineString OGRls;
         OGRMultiLineString OGRmls;
         OGRPoint OGRPt;

         for (int nX = 0; nX < m_nXGridMax; nX++)
         {
            for (int nY = 0; nY < m_nYGridMax; nY++)
            {
               // Only output a value if the cell is a sea cell which is not in the active zone (wave height and angle values are meaningless if in the active zone)
               if ((m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea()) && (! m_pRasterGrid->m_Cell[nX][nY].bIsInActiveZone()))
               {
                  // Create a feature object, one per sea cell
                  OGRFeature *pOGRFeature = NULL;
                  pOGRFeature = OGRFeature::CreateFeature(pOGRLayer->GetLayerDefn());

                  // Set the feature's geometry (in external CRS)
                  OGRPt.setX(dGridCentroidXToExtCRSX(nX));
                  OGRPt.setY(dGridCentroidYToExtCRSY(nY));
                  pOGRFeature->SetGeometry(&OGRPt);

                  double
                     dOrientation = m_pRasterGrid->m_Cell[nX][nY].dGetWaveOrientation(),
                     dHeight = m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight();

                  if ((dHeight == DBL_NODATA) || (dOrientation == DBL_NODATA))
                     continue;

                  // Set the feature's attributes
                  pOGRFeature->SetField(strFieldValue1.c_str(), dOrientation);
                  pOGRFeature->SetField(strFieldValue2.c_str(), dHeight);

                  // Create the feature in the output layer
                  if (pOGRLayer->CreateFeature(pOGRFeature) != OGRERR_NONE)
                  {
                     cerr << ERR << "cannot create " << strType << " feature " << strPlotTitle << " for cell [" << nX << "][" << nY << "] in " << strFilePathName << "\n" << CPLGetLastErrorMsg() << endl;
                     return false;
                  }

                  // Get rid of the feature object
                  OGRFeature::DestroyFeature(pOGRFeature);
               }
            }
         }
      break;
      }

      case (PLOT_AVG_WAVE_AND_HEIGHT):
      {
         eGType = wkbPoint;
         strType = "point";

         // The layer has been created, so create real-numbered values associated with each point
         string
            strFieldValue1 = "Angle",
            strFieldValue2 = "Height";

         // Create the first field
         OGRFieldDefn OGRField1(strFieldValue1.c_str(), OFTReal);
         if (pOGRLayer->CreateField(&OGRField1) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 1 '" << strFieldValue1 << "' in " << strFilePathName << "\n" << CPLGetLastErrorMsg() << endl;
            return false;
         }

         // Create the second field
         OGRFieldDefn OGRField2(strFieldValue2.c_str(), OFTReal);
         if (pOGRLayer->CreateField(&OGRField2) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 2 '" << strFieldValue2 << "' in " << strFilePathName << "\n" << CPLGetLastErrorMsg() << endl;
            return false;
         }

         // OK, now create features
         OGRLineString OGRls;
         OGRMultiLineString OGRmls;
         OGRPoint OGRPt;

         for (int nX = 0; nX < m_nXGridMax; nX++)
         {
            for (int nY = 0; nY < m_nYGridMax; nY++)
            {
               // Create a feature object, one per sea cell
               OGRFeature *pOGRFeature = NULL;
               pOGRFeature = OGRFeature::CreateFeature(pOGRLayer->GetLayerDefn());

               // Set the feature's geometry (in external CRS)
               OGRPt.setX(dGridCentroidXToExtCRSX(nX));
               OGRPt.setY(dGridCentroidYToExtCRSY(nY));
               pOGRFeature->SetGeometry(&OGRPt);

               double
                  dOrientation = m_pRasterGrid->m_Cell[nX][nY].dGetTotWaveOrientation() / m_ulTimestep,
                  dHeight = m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight() / m_ulTimestep;

               if ((dHeight == DBL_NODATA) || (dOrientation == DBL_NODATA))
                  continue;

               // Set the feature's attributes
               pOGRFeature->SetField(strFieldValue1.c_str(), dOrientation);
               pOGRFeature->SetField(strFieldValue2.c_str(), dHeight);

               // Create the feature in the output layer
               if (pOGRLayer->CreateFeature(pOGRFeature) != OGRERR_NONE)
               {
                  cerr << ERR << "cannot create " << strType << " feature " << strPlotTitle << " for cell [" << nX << "][" << nY << "] in " << strFilePathName << "\n" << CPLGetLastErrorMsg() << endl;
                  return false;
               }

               // Get rid of the feature object
               OGRFeature::DestroyFeature(pOGRFeature);
            }
         }
      break;
      }

      case (PLOT_POLYGON_BOUNDARY):
      {
         eGType = wkbPolygon;
         strType = "polygon";

         // The layer has been created, so create two integer-numbered values (the number of the polygon object, and the number of the coast point which is the polygon's node) for the polygon
         string
            strFieldValue1 = "Polygon",
            strFieldValue2 = "CoastNode",
            strFieldValue3 = "TotSedChng",
            strFieldValue4 = "FinSedChng",
            strFieldValue5 = "SndSedChng",
            strFieldValue6 = "CrsSedChng";

         OGRFieldDefn OGRField1(strFieldValue1.c_str(), OFTInteger);
         if (pOGRLayer->CreateField(&OGRField1) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 1 '" << strFieldValue1 << "' in " << strFilePathName << "\n" << CPLGetLastErrorMsg() << endl;
            return false;
         }

         OGRFieldDefn OGRField2(strFieldValue2.c_str(), OFTInteger);
         if (pOGRLayer->CreateField(&OGRField2) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 2 '" << strFieldValue2 << "' in " << strFilePathName << "\n" << CPLGetLastErrorMsg() << endl;
            return false;
         }

         OGRFieldDefn OGRField3(strFieldValue3.c_str(), OFTReal);
         if (pOGRLayer->CreateField(&OGRField3) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 3 '" << strFieldValue3 << "' in " << strFilePathName << "\n" << CPLGetLastErrorMsg() << endl;
            return false;
         }

         OGRFieldDefn OGRField4(strFieldValue4.c_str(), OFTReal);
         if (pOGRLayer->CreateField(&OGRField4) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 4 '" << strFieldValue4 << "' in " << strFilePathName << "\n" << CPLGetLastErrorMsg() << endl;
            return false;
         }

         OGRFieldDefn OGRField5(strFieldValue5.c_str(), OFTReal);
         if (pOGRLayer->CreateField(&OGRField5) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 5 '" << strFieldValue5 << "' in " << strFilePathName << "\n" << CPLGetLastErrorMsg() << endl;
            return false;
         }

         OGRFieldDefn OGRField6(strFieldValue6.c_str(), OFTReal);
         if (pOGRLayer->CreateField(&OGRField6) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 6 '" << strFieldValue6 << "' in " << strFilePathName << "\n" << CPLGetLastErrorMsg() << endl;
            return false;
         }

         // OK, now do features
         OGRLineString OGRls;

         for (int i = 0; i < static_cast<int>(m_VCoast.size()); i++)
         {
            for (int j = 0; j < m_VCoast[i].nGetNumPolygons(); j++)
            {
               // Create a feature object, one per polygon
               OGRFeature *pOGRFeature = NULL;
               pOGRFeature = OGRFeature::CreateFeature(pOGRLayer->GetLayerDefn());

               CGeomCoastPolygon* pPolygon = m_VCoast[i].pGetPolygon(j);

               // Set the feature's attributes
               pOGRFeature->SetField(strFieldValue1.c_str(), j);
               pOGRFeature->SetField(strFieldValue2.c_str(), pPolygon->nGetNodeCoastPoint());
               pOGRFeature->SetField(strFieldValue3.c_str(), pPolygon->dGetDeltaActualTotalSediment());
               pOGRFeature->SetField(strFieldValue4.c_str(), pPolygon->dGetDeltaActualUnconsFine());
               pOGRFeature->SetField(strFieldValue5.c_str(), pPolygon->dGetDeltaActualUnconsSand());
               pOGRFeature->SetField(strFieldValue6.c_str(), pPolygon->dGetDeltaActualUnconsCoarse());

               // Now attach a geometry to the feature object
               for (int n = 0; n < pPolygon->nGetBoundarySize(); n++)
                  //  In external CRS
                  OGRls.addPoint(pPolygon->pPtGetBoundaryPoint(n)->dGetX(), pPolygon->pPtGetBoundaryPoint(n)->dGetY());

               pOGRFeature->SetGeometry(&OGRls);

               // Create the feature in the output layer
               if (pOGRLayer->CreateFeature(pOGRFeature) != OGRERR_NONE)
               {
                  cerr << ERR << "cannot create " << strType << " feature " << strPlotTitle << " for coast " << i << " polygon " << j << " in " << strFilePathName << "\n" << CPLGetLastErrorMsg() << endl;
                  return false;
               }

               // Tidy up: empty the line string and get rid of the feature object
               OGRls.empty();
               OGRFeature::DestroyFeature(pOGRFeature);
            }
         }

         break;
      }
      
      case (PLOT_SHADOW_ZONE_BOUNDARY):
      {
         eGType = wkbLineString;
         strType = "line";
         
         // Create an integer-numbered value (the number of the shadow zone line object) for the multi-line
         string strFieldValue1 = "ShadowLine";
         OGRFieldDefn OGRField1(strFieldValue1.c_str(), OFTInteger);
         if (pOGRLayer->CreateField(&OGRField1) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 1 '" << strFieldValue1 << "' in " << strFilePathName << "\n" << CPLGetLastErrorMsg() << endl;
            return false;
         }
         
         // OK, now do features
         OGRLineString OGRls;
         
         for (int i = 0; i < static_cast<int>(m_VCoast.size()); i++)
         {
            for (int j = 0; j < m_VCoast[i].nGetNumShadowZoneBoundaries(); j++)
            {
               // Create a feature object, one per coast
               OGRFeature *pOGRFeature = NULL;
               pOGRFeature = OGRFeature::CreateFeature(pOGRLayer->GetLayerDefn());
               
               // Set the feature's attribute (the shadow zone line number)
               pOGRFeature->SetField(strFieldValue1.c_str(), j);
               
               // Now attach a geometry to the feature object
               CGeomLine LShadow = *m_VCoast[i].pGetShadowZoneBoundary(j);
               OGRls.addPoint(LShadow.dGetXAt(0), LShadow.dGetYAt(0));
               OGRls.addPoint(LShadow.dGetXAt(1), LShadow.dGetYAt(1));
               
               pOGRFeature->SetGeometry(&OGRls);
               
               // Create the feature in the output layer
               if (pOGRLayer->CreateFeature(pOGRFeature) != OGRERR_NONE)
               {
                  cerr << ERR << "cannot create  " << strType << " feature " << strPlotTitle << j << " for coast " << i << " in " << strFilePathName << "\n" << CPLGetLastErrorMsg() << endl;
                  return false;
               }
               
               // Tidy up: empty the line string and get rid of the feature object
               OGRls.empty();
               OGRFeature::DestroyFeature(pOGRFeature);
            }
         }
         
         break;         
      }
   }

   // Get rid of the dataset object
   GDALClose(pGDALDataSet);

   return true;
}
