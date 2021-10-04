/*!
 *
 * \mainpage
   <b>CoastalME</b> (Coastal Modelling Environment) simulates the long-term behaviour of a coast. This initial version is a prototype which considers simple soft cliff cross-shore effects only\n

   See <a href="https://github.com/coastalme/coastalme_TESTING" target="_blank">https://github.com/coastalme/coastalme_TESTING</a> for the latest version of the source code, and <a href="https://github.com/coastalme/coastalme" target="_blank">https://github.com/coastalme/coastalme</a> for the stable version.\n

 * \section intro_sec Introduction
 * <b>TODO</b> Say more about CoastalME here\n
   \n
   From Shingle Street\n
   To Orford Ness\n
   The waves maraud,\n
   The winds oppress,\n
   The earth can’t help\n
   But acquiesce\n
   For this is east\n
   And east means loss,\n
   A lessening shore, receding ground,\n
   Three feet gone last year, four feet this\n
   Where land runs out and nothing’s sound.\n
   Nothing lasts long on Shingle Street.\n
   \n
   By Blake Morrison (2018). See <a href="http://www.randomhouse.co.uk/editions/shingle-street/9780701188771" target="_blank">http://www.randomhouse.co.uk/editions/shingle-street/9780701188771</a>\n

 * \section install_sec Installation

 * \subsection step1 Step 1: Opening the box

 * \subsection step2 Step 2: Running CoastalME

 * \subsection step3 Step 3: Building datasets
 *
 * \file cme.h
 * \brief This file contains global definitions for CoastalME
 *
 */

#ifndef CME_H
#define CME_H
/*===============================================================================================================================

 This file is part of CoastalME, the Coastal Modelling Environment.

 CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

===============================================================================================================================*/
#include <climits>

#include <sstream>
using std::ostream;
using std::ostringstream;

#include "simulation.h"

//===================================================== platform-specific stuff =================================================
#ifdef _WIN32
#define access _access
#define F_OK 0 // Test for file existence
#endif

#ifdef _MSC_VER
// MS Visual C++, byte order is IEEE little-endian, 32-bit
#ifdef _DEBUG
#include <crtdbg.h> // useful
#endif

// clock_t is a signed long: see <time.h>
long const CLOCK_T_MIN = LONG_MIN;
double const CLOCK_T_RANGE = static_cast<double>(LONG_MAX) - static_cast<double>(CLOCK_T_MIN);
#ifdef _M_ALPHA
string const PLATFORM = "Alpha/MS Visual C++";
#elif defined _M_IX86
string const PLATFORM = "Intel x86/MS Visual C++";
#elif defined _M_MPPC
string const PLATFORM = "Power PC/MS Visual C++";
#elif defined _M_MRX000
string const PLATFORM = "MIPS/MS Visual C++";
#else
string const PLATFORM = "Other/MS Visual C++";
#endif
#endif

#ifdef __GNUG__
// GNU C++
#ifndef CPU
#error CPU not defined
#else
#ifdef x86
// Intel x86, byte order is little-endian, 32-bit
string const PLATFORM = "Intel x86/GNU C++";
// clock_t is an unsigned long: see <time.h>
unsigned long const CLOCK_T_MIN = 0;
double const CLOCK_T_RANGE = static_cast<double>(ULONG_MAX);

#elif defined rs6000
// IBM RS-6000, byte order is big-endian, 32-bit
string const PLATFORM = "IBM RS-6000/GNU C++";
// clock_t is a signed long: see <time.h> NEED TO CHECK
long const CLOCK_T_MIN = LONG_MIN;
double const CLOCK_T_RANGE = static_cast<double>(LONG_MAX) - static_cast<double>(CLOCK_T_MIN);
#elif defined ultrasparc
// Sun UltraSparc, byte order is big-endian, 32-bit
string const PLATFORM = "Sun UltraSPARC/GNU C++";
// clock_t is a signed long: see <time.h>
long const CLOCK_T_MIN = LONG_MIN;
double const CLOCK_T_RANGE = static_cast<double>(LONG_MAX) - static_cast<double>(CLOCK_T_MIN);
#else
// Something else, assume 32-bit
string const PLATFORM = "Other/GNU C++";
// clock_t is a signed long: NEED TO CHECK <time.h>
long const CLOCK_T_MIN = LONG_MIN;
double const CLOCK_T_RANGE = static_cast<double>(LONG_MAX) - static_cast<double>(CLOCK_T_MIN);
#endif
#endif
#endif

#ifdef __MINGW32__
// Minimalist GNU for Windows
//   #define __USE_MINGW_ANSI_STDIO 1        // Fix long doubles output problem, see http://stackoverflow.com/questions/7134547/gcc-printf-and-long-double-leads-to-wrong-output-c-type-conversion-messes-u

#define WEXITSTATUS(x) ((x)&0xff)
#endif

#ifdef __HP_aCC
// HP-UX aCC, byte order is big-endian, can be either 32-bit or 64-bit
string const PLATFORM = "HP-UX aC++";
// clock_t is an unsigned long: see <time.h>
unsigned long const CLOCK_T_MIN = 0;
#ifdef __ia64
// However, clock_t is a 32-bit unsigned long and we are using 64-bit unsigned longs here
double const CLOCK_T_RANGE = 4294967295UL; // crude, improve
#else
double const CLOCK_T_RANGE = static_cast<double>(ULONG_MAX);
#endif
#endif

//===================================================== hard-wired constants ====================================================
char const PATH_SEPARATOR = '/'; // Works for Windows too!
char const SPACE = ' ';
char const QUOTE1 = ';';
char const QUOTE2 = '#';
char const COMMA = ',';
char const DASH = '-';
char const COLON = ':';
char const SLASH = '/';
char const TILDE = '~';

// TESTING options
bool const USE_DEEP_WATER_FOR_SHADOW_LINE = true;       // Use deep water wave orientation in determining shadow line orientation?
bool const CREATE_SHADOW_ZONE_IF_HITS_GRID_EDGE = true; // If shadow line tracing hits grid edge, create shadow zone?
bool const SAVE_CSHORE_OUTPUT = true;                   // #ifdef CSHORE_FILE_INOUT || CSHORE_BOTH, append all CShore output files to a whole-run master

// TODO Make this a user input
bool const ACCEPT_SHORT_PROFILES = true;

int const BUF_SIZE = 2048;     // Max length (inc. terminating NULL) of any C-type string
int const MAX_SAVE_DIGITS = 3; // Maximum number of digits for GIS save number
int const CLOCK_CHECK_ITERATION = 5000;
int const SAVGOL_POLYNOMIAL_MAX_ORDER = 6;   // Maximum order of Savitsky-Golay smoothing polynomial
int const COAST_LENGTH_MAX = 10;             // For safety check when tracing coast
int const COAST_LENGTH_MIN_X_PROF_SPACE = 2; // Ignore very short coasts less than this x profile spacing
int const MAX_NUM_SHADOW_ZONES = 10;         // Consider at most this number of shadow zones
int const GRID_MARGIN = 10;                  // Ignore this many along-coast grid-edge points re. shadow zone calcs
int const MIN_PROFILE_SPACING = 20;          // In cells: profile creation does not work well if profiles are too closely spaced
int const CAPE_POINT_MIN_SPACING = 10;       // In cells: for shadow zone stuff, cape points must not be closer than this
int const FLOOD_FILL_START_OFFSET = 2;       // In cells: flood fill starts this distance inside polygon
int const SHADOW_LINE_MIN_SINCE_HIT_SEA = 5;
int const MAX_LEN_SHADOW_LINE_TO_IGNORE = 200;                // In cells: if can't find flood fill start point, continue if short shadow line
int const MIN_PAR_PROFILE_SIZE = 3;                           // In cells: min size for uncons sed parallel profile
int const MAX_NUM_PREV_ORIENTATION_VALUES = 10;               // Max length of deque used in tracing shadow boundary
int const MIN_INLAND_OFFSET_FOR_BEACH_EROSION_ESTIMATION = 5; // Used in estimation of beach erosion

// Log file detail level
int const LOG_FILE_LEAST_DETAIL = 1;
int const LOG_FILE_MIDDLE_DETAIL = 2;
int const LOG_FILE_MOST_DETAIL = 3;

// TODO Make these user inputs
int const NUMBER_OF_INTERVENTION_CAPES = 4;
int const CSHORE_INTERPOLATION_LINEAR = 0;
int const CSHORE_INTERPOLATION_HERMITE_CUBIC = 1;

int const INT_NODATA = -9999;

int const NO_DIRECTION = 0;
int const NORTH = 1;
int const NORTH_EAST = 2;
int const EAST = 3;
int const SOUTH_EAST = 4;
int const SOUTH = 5;
int const SOUTH_WEST = 6;
int const WEST = 7;
int const NORTH_WEST = 8;

int const DIRECTION_DOWNCOAST = 0; // Down-coast, i.e. along the coast so that the index of coastline points INCREASES
int const DIRECTION_UPCOAST = 1;   // Up-coast, i.e. along the coast so that the index of coastline points DECREASES

// Handedness codes, these show which side the sea is on when travelling down-coast (i.e. in the direction in which coastline point numbers INCREASE)
int const NULL_HANDED = -1;
int const RIGHT_HANDED = 0;
int const LEFT_HANDED = 1;

// Sediment texture codes
int const TEXTURE_FINE = 0;
int const TEXTURE_SAND = 1;
int const TEXTURE_COARSE = 2;

// Time unit codes
int const TIME_UNKNOWN = -1;
int const TIME_HOURS = 0;
int const TIME_DAYS = 1;
int const TIME_MONTHS = 2;
int const TIME_YEARS = 3;

// Intervention input and output codes
int const IO_INTERVENTION_NONE = 0;
int const IO_INTERVENTION_STRUCT = 1;
int const IO_INTERVENTION_NON_STRUCT = 2;

// Generic landform code
int const LF_NONE = 0;

// Landform category codes for cells and coast landform objects (see old source for full list, to be used eventually)
int const LF_CAT_HINTERLAND = 1;
int const LF_CAT_SEA = 2;
int const LF_CAT_CLIFF = 3;
int const LF_CAT_DRIFT = 4;
int const LF_CAT_INTERVENTION = 5;

// Landform sub-category codes for cells, LF_CAT_CLIFF
int const LF_SUBCAT_CLIFF_ON_COASTLINE = 6;
int const LF_SUBCAT_CLIFF_INLAND = 7;

// Landform sub-category codes for cells, for LF_CAT_DRIFT
int const LF_SUBCAT_DRIFT_MIXED = 8;
int const LF_SUBCAT_DRIFT_TALUS = 9;
int const LF_SUBCAT_DRIFT_BEACH = 10;
// TODO
int const LF_SUBCAT_DRIFT_DUNES = 11;

// Landform sub-category codes for cells, for LF_CAT_INTERVENTION. See also "Intervention input and output codes"
int const LF_SUBCAT_INTERVENTION_STRUCT = 12;
int const LF_SUBCAT_INTERVENTION_NON_STRUCT = 13;

// GIS raster input codes
int const FINE_CONS_RASTER = 1;
int const SAND_CONS_RASTER = 2;
int const COARSE_CONS_RASTER = 3;
int const FINE_UNCONS_RASTER = 4;
int const SAND_UNCONS_RASTER = 5;
int const COARSE_UNCONS_RASTER = 6;
int const SUSP_SED_RASTER = 7;
int const LANDFORM_RASTER = 8;
int const INTERVENTION_CLASS_RASTER = 9;
int const INTERVENTION_HEIGHT_RASTER = 10;

// GIS vector data type codes
int const VEC_FIELD_DATA_ANY = 0;
int const VEC_FIELD_DATA_INT = 1;
int const VEC_FIELD_DATA_REAL = 2;
int const VEC_FIELD_DATA_STRING = 3;
int const VEC_FIELD_DATA_OTHER = 4;

// GIS vector geometry codes
int const VEC_GEOMETRY_POINT = 1;
int const VEC_GEOMETRY_LINE = 2;
int const VEC_GEOMETRY_POLYGON = 3;
int const VEC_GEOMETRY_OTHER = 4;

// GIS vector input codes and constraints
int const DEEP_WATER_WAVE_STATIONS_VEC = 1;
int const DEEP_WATER_WAVE_STATIONS_MAX_LAYER = 1;
int const DEEP_WATER_WAVE_STATIONS_POINT_GEOMETRY = VEC_GEOMETRY_POINT;
int const SEDIMENT_INPUT_EVENT_LOCATION_VEC = 2;
int const SEDIMENT_INPUT_EVENT_LOCATION_MAX_LAYER = 1;
int const SEDIMENT_INPUT_EVENT_LOCATION_POINT_GEOMETRY = VEC_GEOMETRY_POINT;
int const SEDIMENT_INPUT_EVENT_LOCATION_LINE_GEOMETRY = VEC_GEOMETRY_LINE;

// GIS raster output codes
int const RASTER_PLOT_ACTIVE_ZONE = 1;
int const RASTER_PLOT_ACTUAL_BEACH_EROSION = 2;
int const RASTER_PLOT_ACTUAL_PLATFORM_EROSION = 3;
int const RASTER_PLOT_AVG_SEA_DEPTH = 4;
int const RASTER_PLOT_AVG_SUSPENDED_SEDIMENT = 5;
int const RASTER_PLOT_AVG_WAVE_HEIGHT = 6;
int const RASTER_PLOT_AVG_WAVE_ORIENTATION = 7;
int const RASTER_PLOT_BASEMENT_ELEVATION = 8;
int const RASTER_PLOT_BEACH_DEPOSITION = 9;
int const RASTER_PLOT_BEACH_MASK = 10;
int const RASTER_PLOT_BEACH_PROTECTION = 11;
int const RASTER_PLOT_CLIFF_COLLAPSE = 12;
int const RASTER_PLOT_CLIFF_COLLAPSE_DEPOSIT = 13;
int const RASTER_PLOT_COARSE_CONSOLIDATED_SEDIMENT = 14;
int const RASTER_PLOT_COARSE_UNCONSOLIDATED_SEDIMENT = 15;
int const RASTER_PLOT_COAST = 16;
int const RASTER_PLOT_DEEP_WATER_WAVE_HEIGHT = 17;
int const RASTER_PLOT_DEEP_WATER_WAVE_ORIENTATION = 18;
int const RASTER_PLOT_DEEP_WATER_WAVE_PERIOD = 19;
int const RASTER_PLOT_FINE_CONSOLIDATED_SEDIMENT = 20;
int const RASTER_PLOT_FINE_UNCONSOLIDATED_SEDIMENT = 21;
int const RASTER_PLOT_INTERVENTION_CLASS = 22;
int const RASTER_PLOT_INTERVENTION_HEIGHT = 23;
int const RASTER_PLOT_INUNDATION_MASK = 24;
int const RASTER_PLOT_LANDFORM = 25;
int const RASTER_PLOT_LOCAL_SLOPE_OF_CONSOLIDATED_SEDIMENT = 26;
int const RASTER_PLOT_NORMAL = 27;
int const RASTER_PLOT_OVERALL_TOP_ELEVATION = 28;
int const RASTER_PLOT_POLYGON = 29;
int const RASTER_PLOT_POLYGON_GAIN_OR_LOSS = 30;
int const RASTER_PLOT_POLYGON_UPDRIFT_OR_DOWNDRIFT = 31;
int const RASTER_PLOT_POTENTIAL_BEACH_EROSION = 32;
int const RASTER_PLOT_POTENTIAL_PLATFORM_EROSION = 33;
int const RASTER_PLOT_SAND_CONSOLIDATED_SEDIMENT = 34;
int const RASTER_PLOT_SAND_UNCONSOLIDATED_SEDIMENT = 35;
int const RASTER_PLOT_SEA_DEPTH = 36;
int const RASTER_PLOT_SEDIMENT_TOP_ELEVATION_ELEV = 37;
int const RASTER_PLOT_SHADOW_DOWNDRIFT_ZONE = 38;
int const RASTER_PLOT_SHADOW_ZONE = 39;
int const RASTER_PLOT_SLICE = 40;
int const RASTER_PLOT_SUSPENDED_SEDIMENT = 41;
int const RASTER_PLOT_TOTAL_ACTUAL_BEACH_EROSION = 42;
int const RASTER_PLOT_TOTAL_ACTUAL_PLATFORM_EROSION = 43;
int const RASTER_PLOT_TOTAL_BEACH_DEPOSITION = 44;
int const RASTER_PLOT_TOTAL_CLIFF_COLLAPSE = 45;
int const RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_DEPOSIT = 46;
int const RASTER_PLOT_TOTAL_POTENTIAL_BEACH_EROSION = 47;
int const RASTER_PLOT_TOTAL_POTENTIAL_PLATFORM_EROSION = 48;
int const RASTER_PLOT_WAVE_HEIGHT = 49;
int const RASTER_PLOT_WAVE_ORIENTATION = 50;
int const RASTER_PLOT_POTENTIAL_PLATFORM_EROSION_MASK = 51;
int const RASTER_PLOT_SEDIMENT_INPUT = 52;

// GIS vector output codes
int const VECTOR_PLOT_AVG_WAVE_ANGLE_AND_HEIGHT = 1;
int const VECTOR_PLOT_BREAKING_WAVE_HEIGHT = 2;
int const VECTOR_PLOT_CLIFF_NOTCH_SIZE = 3;
int const VECTOR_PLOT_COAST = 4;
int const VECTOR_PLOT_COAST_CURVATURE = 5;
int const VECTOR_PLOT_DEEP_WATER_WAVE_ANGLE_AND_HEIGHT = 6;
int const VECTOR_PLOT_DOWNDRIFT_BOUNDARY = 7;
int const VECTOR_PLOT_INVALID_NORMALS = 8;
int const VECTOR_PLOT_MEAN_WAVE_ENERGY = 9;
int const VECTOR_PLOT_NORMALS = 10;
int const VECTOR_PLOT_POLYGON_BOUNDARY = 11;
int const VECTOR_PLOT_POLYGON_NODES = 12;
int const VECTOR_PLOT_SHADOW_BOUNDARY = 13;
int const VECTOR_PLOT_WAVE_ANGLE_AND_HEIGHT = 14;
int const VECTOR_PLOT_WAVE_ENERGY_SINCE_COLLAPSE = 15;
int const VECTOR_PLOT_WAVE_SETUP = 16;
int const VECTOR_PLOT_STORM_SURGE = 17;

// Return codes
int const RTN_OK = 0;
int const RTN_HELPONLY = 1;
int const RTN_CHECKONLY = 2;
int const RTN_USERABORT = 3;
int const RTN_ERR_BADPARAM = 4;
int const RTN_ERR_INI = 5;
int const RTN_ERR_CMEDIR = 6;
int const RTN_ERR_RUNDATA = 7;
int const RTN_ERR_SCAPE_SHAPE_FUNCTION_FILE = 8;
int const RTN_ERR_TIDEDATAFILE = 9;
int const RTN_ERR_LOGFILE = 10;
int const RTN_ERR_OUTFILE = 11;
int const RTN_ERR_TSFILE = 12;
int const RTN_ERR_DEMFILE = 13;
int const RTN_ERR_RASTER_FILE_READ = 14;
int const RTN_ERR_VECTOR_FILE_READ = 15;
int const RTN_ERR_MEMALLOC = 16;
int const RTN_ERR_RASTER_GIS_OUT_FORMAT = 17;
int const RTN_ERR_VECTOR_GIS_OUT_FORMAT = 18;
int const RTN_ERR_TEXT_FILE_WRITE = 19;
int const RTN_ERR_RASTER_FILE_WRITE = 20;
int const RTN_ERR_VECTOR_FILE_WRITE = 21;
int const RTN_ERR_TIMESERIES_FILE_WRITE = 22;
int const RTN_ERR_LINETOGRID = 23;
int const RTN_ERR_PROFILESPACING = 24;
int const RTN_ERR_PROFILE_ENDPOINT_IS_OFFGRID = 25;
int const RTN_ERR_PROFILE_ENDPOINT_IS_INLAND = 26;
int const RTN_ERR_NO_SOLUTION_FOR_ENDPOINT = 27;
int const RTN_ERR_PROFILE_END_INSUFFICIENT_DEPTH = 28;
int const RTN_ERR_NO_PROFILES_1 = 29;
int const RTN_ERR_NO_PROFILES_2 = 30;
int const RTN_ERR_NOSEACELLS = 31;
int const RTN_ERR_GRIDTOLINE = 32;
int const RTN_ERR_TRACECOAST = 33;
int const RTN_ERR_NOCOAST = 34;
int const RTN_ERR_PROFILEWRITE = 35;
int const RTN_ERR_TIMEUNITS = 36;
int const RTN_ERR_CLIFFNOTCH = 37;
int const RTN_ERR_CLIFFDEPOSIT = 38;
int const RTN_ERR_BAD_INDEX = 39;
int const RTN_ERR_EDGEOFGRID = 40;
int const RTN_ERR_NO_SEAWARD_END_OF_PROFILE_1 = 41;
int const RTN_ERR_NO_SEAWARD_END_OF_PROFILE_2 = 42;
int const RTN_ERR_NO_SEAWARD_END_OF_PROFILE_3 = 43;
int const RTN_ERR_NO_SEAWARD_END_OF_PROFILE_4 = 44;
int const RTN_ERR_LANDFORM_TO_GRID = 45;
int const RTN_ERR_NO_TOP_LAYER = 46;
int const RTN_ERR_NO_ADJACENT_POLYGON = 47;
int const RTN_ERR_BAD_MULTILINE = 48;
int const RTN_ERR_CANNOT_INSERT_POINT = 49;
int const RTN_ERR_CANNOT_ASSIGN_COASTAL_LANDFORM = 50;
int const RTN_ERR_SHADOW_ZONE_FLOOD_FILL_NOGRID = 51;
int const RTN_ERR_SHADOW_ZONE_FLOOD_START_POINT = 52;
int const RTN_ERR_CSHORE_EMPTY_PROFILE = 53;
int const RTN_ERR_CSHORE_FILE_INPUT = 54;
int const RTN_ERR_READING_CSHORE_FILE_OUTPUT = 55;
int const RTN_ERR_WAVE_INTERPOLATION_LOOKUP = 56;
int const RTN_ERR_GRIDCREATE = 57;
int const RTN_ERR_COAST_CANT_FIND_EDGE_CELL = 58;
int const RTN_ERR_CSHORE_ERROR = 59;
int const RTN_ERR_NO_CELL_UNDER_COASTLINE = 60;
int const RTN_ERR_OPEN_DEEP_WATER_WAVE_DATA = 61;
int const RTN_ERR_READING_DEEP_WATER_WAVE_DATA = 62;
int const RTN_ERR_BOUNDING_BOX = 63;
int const RTN_ERR_READING_SEDIMENT_INPUT_EVENT = 64;
int const RTN_ERR_SEDIMENT_INPUT_EVENT = 65;

// Elevation and 'slice' codes
int const ELEV_IN_BASEMENT = -1;
int const ELEV_ABOVE_SEDIMENT_TOP = -2;
int const NO_NONZERO_THICKNESS_LAYERS = -3;

// Vector smoothing codes
int const SMOOTH_NONE = 0;
int const SMOOTH_RUNNING_MEAN = 1;
int const SMOOTH_SAVITZKY_GOLAY = 2;

// Grid-edge boundary treatment for unconsolidated sediment movement
int const GRID_EDGE_CLOSED = 0;
int const GRID_EDGE_OPEN = 1;
int const GRID_EDGE_RECIRCULATE = 2;

// Model for wave propagation
int const WAVE_MODEL_COVE = 0;
int const WAVE_MODEL_CSHORE = 1;

// Equation for estimating erosion of unconsolidated sediment
int const UNCONS_SEDIMENT_EQUATION_CERC = 0;
int const UNCONS_SEDIMENT_EQUATION_KAMPHUIS = 1;

unsigned long const MASK = 0xfffffffful;
unsigned long const SEDINPUTEVENTERROR = -1;

double const PI = 3.141592653589793238462643;

double const D50_FINE_DEFAULT = 0.0625; // In mm
double const D50_SAND_DEFAULT = 0.42;   // Ditto
double const D50_COARSE_DEFAULT = 19.0; // Ditto

double const BEACH_PROTECTION_HB_RATIO = 0.23; // The beach protection factor is this times breaking depth
double const WALKDEN_HALL_PARAM_1 = 3.25;      // First param in Equation 4 from Walkden & Hall, 2005
double const WALKDEN_HALL_PARAM_2 = 1.50;      // Second param in Equation 4 from Walkden & Hall, 2005

double const DEPTH_OVER_DB_INCREMENT = 0.001;        // Depth Over DB increment for erosion potential look-up function
double const INVERSE_DEPTH_OVER_DB_INCREMENT = 1000; // Inverse of the above
double const DEAN_POWER = 2.0 / 3.0;                 // Dean profile exponent

// TODO Let the user define these CShore input parameters
double const CSHORE_FRICTION_FACTOR = 0.015; // Friction factor for CShore model
double const CSHORE_SURGE_LEVEL = 0.0;       // Not used, but in the future we might include surge in the calculations

double const TOLERANCE = 1e-4;                // For bFPIsEqual, if too small (e.g. 1e-10), get spurious "rounding" errors
double const SEDIMENT_ELEV_TOLERANCE = 1e-10; // Differences in depth-equivalent sediment amount (m) less than this are ignored
double const STRAIGHT_COAST_MAX_DETAILED_CURVATURE = -5;
double const STRAIGHT_COAST_MAX_SMOOTH_CURVATURE = -1;
double const MIN_LENGTH_OF_SHADOW_ZONE_LINE = 10;     // Used in shadow line tracing
double const MAX_LAND_LENGTH_OF_SHADOW_ZONE_LINE = 5; // Used in shadow line tracing

double const DBL_NODATA = -9999;

string const PROGRAM_NAME = "CoastalME 0.9.9 TESTING: 4 October 2021";
string const PROGRAM_NAME_SHORT = "CME";
string const CME_INI = "cme.ini";

string const COPYRIGHT = "(C) 2021 Andres Payo and David Favis-Mortlock";
string const LINE = "-------------------------------------------------------------------------------";
string const DISCLAIMER1 = "This program is distributed in the hope that it will be useful, but WITHOUT ANY";
string const DISCLAIMER2 = "WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A";
string const DISCLAIMER3 = "PARTICULAR PURPOSE. See the GNU General Public License for more details. You";
string const DISCLAIMER4 = "should have received a copy of the GNU General Public License along with this";
string const DISCLAIMER5 = "program; if not, contact the Free Software Foundation, Inc., 675 Mass Ave,";
string const DISCLAIMER6 = "Cambridge, MA 02139, USA.";

string const ABOUT = "simulates the long-term behaviour of a coast. This initial version considers only simple soft cliff cross-shore effects";
string const THANKS = "Many thanks to:\n\tJim W. Hall\n\tMartin D. Hurst\n\tMike J.A. Walkden\n\tIan Townend\n\tMark Dickson\n\tMatthew Ives\n\tRobert J. Nicholls";
string const GDALDRIVERS = "GDAL drivers";

string const USAGE = "Usage: cme [OPTION]...";
string const USAGE1 = "  --gdal             List GDAL drivers";
string const USAGE2 = "  --about            Information about this program";
string const USAGE3 = "  --help             Display this text";
string const USAGE4 = "  --home=DIRECTORY   Specify the location of the .ini file etc.";
string const USAGE5 = "  --datafile=FILE    Specify the location and name of the main datafile";

string const START_NOTICE = "- Started on ";
string const INITIALIZING_NOTICE = "- Initializing";
string const READING_FILE_LOCATIONS = "  - Reading file locations: ";
string const READING_RUN_DATA = "  - Reading run data file: ";
string const READING_BASEMENT = "  - Reading basement DEM: ";
string const READING_RASTER_FILES = "  - Reading raster GIS files";
string const READING_LANDFORM_FILE = "    - Landform class: ";
string const READING_INTERVENTION_CLASS_FILE = "    - Intervention class: ";
string const READING_INTERVENTION_HEIGHT_FILE = "    - Intervention height: ";
string const READING_SUSPENDED_SEDIMENT_FILE = "    - Suspended sediment: ";
string const READING_UNCONS_FINE_SEDIMENT_FILE = "    - Unconsolidated fine sediment (layer ";
string const READING_UNCONS_SAND_SEDIMENT_FILE = "    - Unconsolidated sand sediment (layer ";
string const READING_UNCONS_COARSE_SEDIMENT_FILE = "    - Unconsolidated coarse sediment (layer ";
string const READING_CONS_FINE_SEDIMENT_FILE = "    - Consolidated fine sediment (layer ";
string const READING_CONS_SAND_SEDIMENT_FILE = "    - Consolidated sand sediment (layer ";
string const READING_CONS_COARSE_SEDIMENT_FILE = "    - Consolidated coarse sediment (layer ";
string const READING_VECTOR_FILES = "  - Reading vector GIS files";
string const READING_DEEP_WATER_WAVE_FILE = "    - Deep water wave values: ";
string const READING_SED_INPUT_EVENT_FILE = "    - Sediment input event values: ";
string const READING_SCAPE_SHAPE_FUNCTION_FILE = "  - Reading SCAPE shape function file";
string const READING_TIDE_DATA_FILE = "  - Reading tide data file: ";
string const ALLOCATE_MEMORY = "  - Allocating memory for raster grid";
string const ADD_LAYERS = "  - Adding sediment layers to raster grid";
string const INITIALIZING = "  - Initializing";
string const RUN_NOTICE = "- Running simulation";
string const SIMULATING = "\r  - Simulating ";
string const FINAL_OUTPUT = "- Writing final output";
string const SEND_EMAIL = "  - Sending email to ";
string const RUN_END_NOTICE = "- Run ended at ";
string const PRESS_KEY = "Press any key to continue...";

string const ERROR_NOTICE = " with error code ";
string const EMAIL_ERROR = "Could not send email";

string const SCAPE_DIR = "scape/";
string const SCAPE_SHAPE_FUNCTION_FILE = "ShapeFunction.dat";
string const EROSION_POTENTIAL_LOOKUP_FILE = "ErosionPotential.csv";

string const CSHORE_DIR = "cshore/";
string const CSHORE_INFILE = "infile";

string const ERR = "*** ERROR ";
string const WARN = "WARNING ";
string const NOTE = "      Note ";

string const PERITERHEAD1 =
    "<-----ELAPSED----><-SEA-><----POTENTIAL---><----------ACTUAL-----------><----POTENTIAL----><------------ACTUAL------------><-----------ACTUAL-----------><-SEDIMENT-><--CLIFF COLLAPSE--><SUSP>";

string const PERITERHEAD2 =
    "       TIME        DEPTH  PLATFORM EROSION        PLATFORM EROSION          BEACH EROSION            BEACH EROSION                BEACH DEPOSITION        INPUT EVENT  EROSION DEPOSITION  SED";

string const PERITERHEAD3 =
    "Time  Hours  Years   Avg  % Sea   All  Erod % Sea   All Erod <-sea avg->  % Sea   All  Erod  % Sea   All Erodng <-sea avg->  % Sea   All  Deposit <-sea->            <-coast avg><--sea->";
string const PERITERHEAD4 =
    "Step                       Area   Sea  Area  Area   Sea Area   F   S   C   Area   Sea  Area   Area   Sea   Area   F   S   C   Area   Sea     Area   S   C   F   S   C   F   S   C   S   C    F";
string const PERITERHEAD5 =
    "                                  Avg   Avg         Avg  Avg                            Avg          Avg    Avg                      Avg      Avg";

string const PERITERHEAD =
    "PER-ITERATION RESULTS =========================================================================================================================================================================";
string const ENDHYDROLOGYHEAD =
    "END OF SIMULATION: HYDROLOGY ==================================================================================================================================================================";
string const ENDSEDIMENTHEAD =
    "END OF SIMULATION: SEDIMENT MOVEMENT ==========================================================================================================================================================";
string const PERFORMHEAD =
    "END OF SIMULATION: PERFORMANCE ================================================================================================================================================================";

string const OUTEXT = ".out";
string const LOGEXT = ".log";
string const CSVEXT = ".csv";

string const DEEP_WATER_WAVE_STATION_ID = "id";
// string const   DEEP_WATER_WAVE_STATIONS_HEIGHT                             = "HEIGHT";
// string const   DEEP_WATER_WAVE_STATIONS_ANGLE                              = "ANGLE";
// string const   DEEP_WATER_WAVE_STATIONS_PERIOD                             = "PERIOD";
string const SEDIMENT_INPUT_EVENT_LOCATION_ID = "id";

// GIS raster output user codes
string const RASTER_ALL_OUTPUT_CODE = "all";
string const RASTER_SEDIMENT_TOP_CODE = "sediment_top_elevation";
string const RASTER_SEDIMENT_TOP_NAME = "sediment_top_elevation";
string const RASTER_TOP_CODE = "top_elevation";
string const RASTER_TOP_NAME = "top_elevation";
string const RASTER_BASEMENT_ELEVATION_CODE = "basement_elevation";
string const RASTER_BASEMENT_ELEVATION_NAME = "basement_elevation";
string const RASTER_LOCAL_SLOPE_CODE = "local_cons_sediment_slope";
string const RASTER_LOCAL_SLOPE_NAME = "local_cons_sediment_slope";
string const RASTER_SEA_DEPTH_CODE = "sea_depth";
string const RASTER_SEA_DEPTH_NAME = "sea_depth";
string const RASTER_AVG_SEA_DEPTH_CODE = "avg_sea_depth";
string const RASTER_AVG_SEA_DEPTH_NAME = "avg_sea_depth";
string const RASTER_INUNDATION_MASK_CODE = "inundation_mask";
string const RASTER_INUNDATION_MASK_NAME = "inundation_mask";
string const RASTER_WAVE_HEIGHT_CODE = "wave_height";
string const RASTER_WAVE_HEIGHT_NAME = "wave_height";
string const RASTER_AVG_WAVE_HEIGHT_CODE = "avg_wave_height";
string const RASTER_AVG_WAVE_HEIGHT_NAME = "avg_wave_height";
string const RASTER_WAVE_ORIENTATION_CODE = "wave_orientation";
string const RASTER_WAVE_ORIENTATION_NAME = "wave_orientation";
string const RASTER_WAVE_PERIOD_CODE = "wave_period";
string const RASTER_WAVE_PERIOD_NAME = "wave_period";
string const RASTER_AVG_WAVE_ORIENTATION_CODE = "avg_wave_orientation";
string const RASTER_AVG_WAVE_ORIENTATION_NAME = "avg_wave_orientation";
string const RASTER_BEACH_MASK_CODE = "beach_mask";
string const RASTER_BEACH_MASK_NAME = "beach_mask";
string const RASTER_BEACH_PROTECTION_CODE = "beach_protection";
string const RASTER_BEACH_PROTECTION_NAME = "beach_protection";
string const RASTER_POTENTIAL_PLATFORM_EROSION_MASK_CODE = "potential_platform_erosion_mask";
string const RASTER_POTENTIAL_PLATFORM_EROSION_MASK_NAME = "potential_platform_erosion_mask";
string const RASTER_POTENTIAL_PLATFORM_EROSION_CODE = "potential_platform_erosion";
string const RASTER_POTENTIAL_PLATFORM_EROSION_NAME = "potential_platform_erosion";
string const RASTER_ACTUAL_PLATFORM_EROSION_CODE = "actual_platform_erosion";
string const RASTER_ACTUAL_PLATFORM_EROSION_NAME = "actual_platform_erosion";
string const RASTER_TOTAL_POTENTIAL_PLATFORM_EROSION_CODE = "total_potential_platform_erosion";
string const RASTER_TOTAL_POTENTIAL_PLATFORM_EROSION_NAME = "total_potential_platform_erosion";
string const RASTER_TOTAL_ACTUAL_PLATFORM_EROSION_CODE = "total_actual_platform_erosion";
string const RASTER_TOTAL_ACTUAL_PLATFORM_EROSION_NAME = "total_actual_platform_erosion";
string const RASTER_POTENTIAL_BEACH_EROSION_CODE = "potential_beach_erosion";
string const RASTER_POTENTIAL_BEACH_EROSION_NAME = "potential_beach_erosion";
string const RASTER_ACTUAL_BEACH_EROSION_CODE = "actual_beach_erosion";
string const RASTER_ACTUAL_BEACH_EROSION_NAME = "actual_beach_erosion";
string const RASTER_TOTAL_POTENTIAL_BEACH_EROSION_CODE = "total_potential_beach_erosion";
string const RASTER_TOTAL_POTENTIAL_BEACH_EROSION_NAME = "total_potential_beach_erosion";
string const RASTER_TOTAL_ACTUAL_BEACH_EROSION_CODE = "total_actual_beach_erosion";
string const RASTER_TOTAL_ACTUAL_BEACH_EROSION_NAME = "total_actual_beach_erosion";
string const RASTER_BEACH_DEPOSITION_CODE = "beach_deposition";
string const RASTER_BEACH_DEPOSITION_NAME = "beach_deposition";
string const RASTER_TOTAL_BEACH_DEPOSITION_CODE = "total_beach_deposition";
string const RASTER_TOTAL_BEACH_DEPOSITION_NAME = "total_beach_deposition";
string const RASTER_LANDFORM_CODE = "landform_class";
string const RASTER_LANDFORM_NAME = "landform_class";
string const RASTER_INTERVENTION_CLASS_CODE = "intervention_class";
string const RASTER_INTERVENTION_CLASS_NAME = "intervention_class";
string const RASTER_INTERVENTION_HEIGHT_CODE = "intervention_height";
string const RASTER_INTERVENTION_HEIGHT_NAME = "intervention_height";
string const RASTER_SUSP_SED_CODE = "susp_sed";
string const RASTER_SUSP_SED_NAME = "susp_sed";
string const RASTER_AVG_SUSP_SED_CODE = "avg_susp_sed";
string const RASTER_AVG_SUSP_SED_NAME = "avg_susp_sed";
string const RASTER_FINE_UNCONS_CODE = "uncons_sed_fine";
string const RASTER_FINE_UNCONS_NAME = "uncons_sed_fine";
string const RASTER_SAND_UNCONS_CODE = "uncons_sed_sand";
string const RASTER_SAND_UNCONS_NAME = "uncons_sed_sand";
string const RASTER_COARSE_UNCONS_CODE = "uncons_sed_coarse";
string const RASTER_COARSE_UNCONS_NAME = "uncons_sed_coarse";
string const RASTER_FINE_CONS_CODE = "cons_sed_fine";
string const RASTER_FINE_CONS_NAME = "cons_sed_fine";
string const RASTER_SAND_CONS_CODE = "cons_sed_sand";
string const RASTER_SAND_CONS_NAME = "cons_sed_sand";
string const RASTER_COARSE_CONS_CODE = "cons_sed_coarse";
string const RASTER_COARSE_CONS_NAME = "cons_sed_coarse";
string const RASTER_COAST_CODE = "rcoast";
string const RASTER_COAST_NAME = "rcoast";
string const RASTER_COAST_NORMAL_CODE = "rcoast_normal";
string const RASTER_COAST_NORMAL_NAME = "rcoast_normal";
string const RASTER_ACTIVE_ZONE_CODE = "active_zone";
string const RASTER_ACTIVE_ZONE_NAME = "active_zone";
string const RASTER_CLIFF_COLLAPSE_CODE = "cliff_collapse";
string const RASTER_CLIFF_COLLAPSE_NAME = "cliff_collapse";
string const RASTER_TOTAL_CLIFF_COLLAPSE_CODE = "total_cliff_collapse";
string const RASTER_TOTAL_CLIFF_COLLAPSE_NAME = "total_cliff_collapse";
string const RASTER_CLIFF_COLLAPSE_DEPOSITION_CODE = "cliff_collapse_deposition";
string const RASTER_CLIFF_COLLAPSE_DEPOSITION_NAME = "cliff_collapse_deposition";
string const RASTER_TOTAL_CLIFF_COLLAPSE_DEPOSITION_CODE = "total_cliff_collapse_deposition";
string const RASTER_TOTAL_CLIFF_COLLAPSE_DEPOSITION_NAME = "total_cliff_collapse_deposition";
string const RASTER_POLYGON_CODE = "polygon_raster";
string const RASTER_POLYGON_NAME = "polygon_raster";
string const RASTER_SLICE_CODE = "slice";
string const RASTER_SLICE_NAME = "slice";
string const RASTER_SHADOW_ZONE_CODE = "shadow_zones";
string const RASTER_SHADOW_ZONE_NAME = "shadow_zones";
string const RASTER_SHADOW_DOWNDRIFT_ZONE_CODE = "shadow_downdrift_zones";
string const RASTER_SHADOW_DOWNDRIFT_ZONE_NAME = "shadow_downdrift_zones";
string const RASTER_DEEP_WATER_WAVE_ORIENTATION_CODE = "deep_water_wave_orientation";
string const RASTER_DEEP_WATER_WAVE_ORIENTATION_NAME = "deep_water_wave_orientation";
string const RASTER_DEEP_WATER_WAVE_HEIGHT_CODE = "deep_water_wave_height";
string const RASTER_DEEP_WATER_WAVE_HEIGHT_NAME = "deep_water_wave_height";
string const RASTER_DEEP_WATER_WAVE_PERIOD_CODE = "deep_water_wave_period";
string const RASTER_DEEP_WATER_WAVE_PERIOD_NAME = "deep_water_wave_period";
string const RASTER_POLYGON_UPDRIFT_OR_DOWNDRIFT_CODE = "polygon_updrift_or_downdrift";
string const RASTER_POLYGON_UPDRIFT_OR_DOWNDRIFT_NAME = "polygon_updrift_or_downdrift";
string const RASTER_POLYGON_GAIN_OR_LOSS_CODE = "polygon_gain_or_loss";
string const RASTER_POLYGON_GAIN_OR_LOSS_NAME = "polygon_gain_or_loss";
string const RASTER_SEDIMENT_INPUT_EVENT_CODE = "sediment_input_total";
string const RASTER_SEDIMENT_INPUT_EVENT_NAME = "sediment_input_total";

// GIS raster output titles
string const RASTER_PLOT_ACTIVE_ZONE_TITLE = "Active zone";
string const RASTER_PLOT_ACTUAL_BEACH_EROSION_TITLE = "Actual (constrained) beach erosion depth";
string const RASTER_PLOT_ACTUAL_PLATFORM_EROSION_TITLE = "Actual (constrained) shore platform erosion depth";
string const RASTER_PLOT_AVG_SEA_DEPTH_TITLE = "Average sea depth";
string const RASTER_PLOT_AVG_SUSPENDED_SEDIMENT_TITLE = "Average depth of suspended sediment";
string const RASTER_PLOT_AVG_WAVE_HEIGHT_TITLE = "Average wave height";
string const RASTER_PLOT_AVG_WAVE_ORIENTATION_TITLE = "Average wave orientation";
string const RASTER_PLOT_BASEMENT_ELEVATION_TITLE = "Basement elevation";
string const RASTER_PLOT_BEACH_DEPOSITION_TITLE = "Beach deposition depth";
string const RASTER_PLOT_BEACH_MASK_TITLE = "Beach mask";
string const RASTER_PLOT_BEACH_PROTECTION_TITLE = "Beach protection factor";
string const RASTER_PLOT_CLIFF_COLLAPSE_DEPOSIT_TITLE = "Cliff collapse deposition depth";
string const RASTER_PLOT_CLIFF_COLLAPSE_TITLE = "Cliff collapse depth";
string const RASTER_PLOT_COARSE_CONSOLIDATED_SEDIMENT_TITLE = "Consolidated coarse sediment depth";
string const RASTER_PLOT_COARSE_UNCONSOLIDATED_SEDIMENT_TITLE = "Unconsolidated coarse sediment depth";
string const RASTER_PLOT_COAST_TITLE = "Rasterized coastline";
string const RASTER_PLOT_DEEP_WATER_WAVE_HEIGHT_TITLE = "Deep water wave height";
string const RASTER_PLOT_DEEP_WATER_WAVE_ORIENTATION_TITLE = "Deep water wave orientation";
string const RASTER_PLOT_DEEP_WATER_WAVE_PERIOD_TITLE = "Deep water wave period";
string const RASTER_PLOT_FINE_CONSOLIDATED_SEDIMENT_TITLE = "Consolidated fine sediment depth";
string const RASTER_PLOT_FINE_UNCONSOLIDATED_SEDIMENT_TITLE = "Unconsolidated fine sediment depth";
string const RASTER_PLOT_INTERVENTION_CLASS_TITLE = "Intervention class";
string const RASTER_PLOT_INTERVENTION_HEIGHT_TITLE = "Intervention height";
string const RASTER_PLOT_INUNDATION_MASK_TITLE = "Inundated area mask";
string const RASTER_PLOT_LANDFORM_TITLE = "Landform class";
string const RASTER_PLOT_LOCAL_SLOPE_OF_CONSOLIDATED_SEDIMENT_TITLE = "Local slope of consolidated sediment";
string const RASTER_PLOT_NORMAL_TITLE = "Rasterized normals to coastline";
string const RASTER_PLOT_OVERALL_TOP_ELEVATION_TITLE = "Elevation of sediment top plus intervention, or sea surface";
string const RASTER_PLOT_POLYGON_GAIN_OR_LOSS_TITLE = "Polygon gain or loss of unconsolidated sediment";
string const RASTER_PLOT_POLYGON_TITLE = "Rasterized polygon boundaries";
string const RASTER_PLOT_POLYGON_UPDRIFT_OR_DOWNDRIFT_TITLE = "Polygon updrift or downdrift movement of unconsolidated sediment";
string const RASTER_PLOT_POTENTIAL_BEACH_EROSION_TITLE = "Potential (unconstrained) beach erosion depth";
string const RASTER_PLOT_POTENTIAL_PLATFORM_EROSION_TITLE = "Potential (unconstrained) shore platform erosion depth";
string const RASTER_PLOT_SAND_CONSOLIDATED_SEDIMENT_TITLE = "Consolidated sand sediment depth";
string const RASTER_PLOT_SAND_UNCONSOLIDATED_SEDIMENT_TITLE = "Unconsolidated sand sediment depth";
string const RASTER_PLOT_SEA_DEPTH_TITLE = "Sea depth";
string const RASTER_PLOT_SEDIMENT_TOP_ELEVATION_ELEV_TITLE = "Elevation of sediment top";
string const RASTER_PLOT_SHADOW_DOWNDRIFT_ZONE_TITLE = "Downdrift of wave shadow zones";
string const RASTER_PLOT_SHADOW_ZONE_TITLE = "Wave shadow zones";
string const RASTER_PLOT_SLICE_TITLE = "Slice though layers at elevation = ";
string const RASTER_PLOT_SUSPENDED_SEDIMENT_TITLE = "Suspended sediment depth";
string const RASTER_PLOT_TOTAL_ACTUAL_BEACH_EROSION_TITLE = "Total actual (constrained) beach erosion depth";
string const RASTER_PLOT_TOTAL_ACTUAL_PLATFORM_EROSION_TITLE = "Total actual (constrained) shore platform erosion depth";
string const RASTER_PLOT_TOTAL_BEACH_DEPOSITION_TITLE = "Total beach deposition depth";
string const RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_DEPOSIT_TITLE = "Total of cliff collapse deposition depth";
string const RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_TITLE = "Total of cliff collapse depth";
string const RASTER_PLOT_TOTAL_POTENTIAL_BEACH_EROSION_TITLE = "Total potential (unconstrained) beach erosion depth";
string const RASTER_PLOT_TOTAL_POTENTIAL_PLATFORM_EROSION_TITLE = "Total potential (unconstrained) shore platform erosion depth";
string const RASTER_PLOT_WAVE_HEIGHT_TITLE = "Wave height";
string const RASTER_PLOT_WAVE_ORIENTATION_TITLE = "Wave orientation";
string const RASTER_PLOT_POTENTIAL_PLATFORM_EROSION_MASK_TITLE = "Potential (unconstrained) shore platform erosion binary mask";
string const RASTER_PLOT_SEDIMENT_INPUT_EVENT_TITLE = "Sediment input event(s) since last GIS save";

// GIS vector output user codes
string const VECTOR_ALL_OUTPUT_CODE = "all";
string const VECTOR_COAST_CODE = "coast";
string const VECTOR_COAST_NAME = "coast";
string const VECTOR_NORMALS_CODE = "normals";
string const VECTOR_NORMALS_NAME = "normals";
string const VECTOR_INVALID_NORMALS_CODE = "invalid_normals";
string const VECTOR_INVALID_NORMALS_NAME = "invalid_normals";
// string const   VECTOR_COLLAPSE_NORMALS_CODE                             = "collapse_normals";
// string const   VECTOR_COLLAPSE_NORMALS_NAME                             = "collapse_normals";
string const VECTOR_COAST_CURVATURE_CODE = "coast_curvature";
string const VECTOR_COAST_CURVATURE_NAME = "coast_curvature";
string const VECTOR_WAVE_ANGLE_AND_HEIGHT_CODE = "wave_angle";
string const VECTOR_WAVE_ANGLE_AND_HEIGHT_NAME = "wave_angle";
string const VECTOR_AVG_WAVE_ANGLE_AND_HEIGHT_NAME = "avg_wave_angle";
string const VECTOR_AVG_WAVE_ANGLE_AND_HEIGHT_CODE = "avg_wave_angle";
string const VECTOR_WAVE_ENERGY_SINCE_COLLAPSE_CODE = "wave_energy";
string const VECTOR_WAVE_ENERGY_SINCE_COLLAPSE_NAME = "wave_energy";
string const VECTOR_MEAN_WAVE_ENERGY_CODE = "mean_wave_energy";
string const VECTOR_MEAN_WAVE_ENERGY_NAME = "mean_wave_energy";
string const VECTOR_BREAKING_WAVE_HEIGHT_CODE = "breaking_wave_height";
string const VECTOR_BREAKING_WAVE_HEIGHT_NAME = "breaking_wave_height";
string const VECTOR_POLYGON_NODE_SAVE_CODE = "node";
string const VECTOR_POLYGON_NODES_NAME = "node";
string const VECTOR_POLYGON_BOUNDARY_SAVE_CODE = "polygon";
string const VECTOR_POLYGON_BOUNDARY_NAME = "polygon";
string const VECTOR_CLIFF_NOTCH_SIZE_CODE = "cliff_notch";
string const VECTOR_CLIFF_NOTCH_SIZE_NAME = "cliff_notch";
string const VECTOR_SHADOW_BOUNDARY_CODE = "shadow_boundary";
string const VECTOR_SHADOW_BOUNDARY_NAME = "shadow_boundary";
string const VECTOR_DOWNDRIFT_BOUNDARY_CODE = "downdrift_boundary";
string const VECTOR_DOWNDRIFT_BOUNDARY_NAME = "downdrift_boundary";
string const VECTOR_DEEP_WATER_WAVE_ANGLE_AND_HEIGHT_CODE = "deep_water_wave_angle";
string const VECTOR_DEEP_WATER_WAVE_ANGLE_AND_HEIGHT_NAME = "deep_water_wave_angle";
string const VECTOR_WAVE_SETUP_CODE = "wave_setup";
string const VECTOR_WAVE_SETUP_NAME = "wave_setup";
string const VECTOR_STORM_SURGE_CODE = "storm_surge";
string const VECTOR_STORM_SURGE_NAME = "storm_surge";

// GIS vector output titles
string const VECTOR_PLOT_AVG_WAVE_ANGLE_AND_HEIGHT_TITLE = "Average wave orientation and height";
string const VECTOR_PLOT_BREAKING_WAVE_HEIGHT_TITLE = "Breaking wave height";
string const VECTOR_PLOT_CLIFF_NOTCH_SIZE_TITLE = "Cliff notch incision";
string const VECTOR_PLOT_COAST_CURVATURE_TITLE = "Coastline curvature";
string const VECTOR_PLOT_COAST_TITLE = "Coastline";
string const VECTOR_PLOT_DEEP_WATER_WAVE_ANGLE_AND_HEIGHT_TITLE = "Deep water wave orientation and height";
string const VECTOR_PLOT_DOWNDRIFT_BOUNDARY_TITLE = "Downdrift zone boundary";
string const VECTOR_PLOT_INVALID_NORMALS_TITLE = "INVALID Coastline-normal profiles";
string const VECTOR_PLOT_MEAN_WAVE_ENERGY_TITLE = "Mean wave energy";
string const VECTOR_PLOT_NORMALS_TITLE = "Coastline-normal profiles";
string const VECTOR_PLOT_POLYGON_BOUNDARY_TITLE = "Polygons";
string const VECTOR_PLOT_POLYGON_NODES_TITLE = "Polygon nodes";
string const VECTOR_PLOT_SHADOW_BOUNDARY_TITLE = "Shadow zone boundary";
string const VECTOR_PLOT_WAVE_ANGLE_AND_HEIGHT_TITLE = "Wave orientation and height";
string const VECTOR_PLOT_WAVE_ENERGY_SINCE_COLLAPSE_TITLE = "Wave energy since collapse";
string const VECTOR_PLOT_WAVE_SETUP_TITLE = "Wave setup";
string const VECTOR_PLOT_STORM_SURGE_TITLE = "Storm surge";

// Time series codes
string const TIME_SERIES_SEA_AREA_NAME = "sea_area";
string const TIME_SERIES_SEA_AREA_CODE = "sea_area";

string const TIME_SERIES_STILL_WATER_LEVEL_NAME = "still_water_level";
string const TIME_SERIES_STILL_WATER_LEVEL_CODE = "water_level";

string const TIME_SERIES_PLATFORM_EROSION_NAME = "platform_erosion";
string const TIME_SERIES_PLATFORM_EROSION_CODE = "platform_erosion";

string const TIME_SERIES_CLIFF_COLLAPSE_EROSION_NAME = "cliff_collapse_erosion";
string const TIME_SERIES_CLIFF_COLLAPSE_EROSION_CODE = "cliff_collapse_erosion";

string const TIME_SERIES_CLIFF_COLLAPSE_DEPOSITION_NAME = "cliff_collapse_deposition";
string const TIME_SERIES_CLIFF_COLLAPSE_DEPOSITION_CODE = "cliff_collapse_deposition";

string const TIME_SERIES_CLIFF_COLLAPSE_NET_NAME = "cliff_collapse_net";
string const TIME_SERIES_CLIFF_COLLAPSE_NET_CODE = "cliff_collapse_net";

string const TIME_SERIES_BEACH_EROSION_NAME = "beach_erosion";
string const TIME_SERIES_BEACH_EROSION_CODE = "beach_erosion";

string const TIME_SERIES_BEACH_DEPOSITION_NAME = "beach_deposition";
string const TIME_SERIES_BEACH_DEPOSITION_CODE = "beach_deposition";

string const TIME_SERIES_BEACH_CHANGE_NET_NAME = "beach_change_net";
string const TIME_SERIES_BEACH_CHANGE_NET_CODE = "beach_change_net";

string const TIME_SERIES_SUSPENDED_SEDIMENT_NAME = "suspended_sediment";
string const TIME_SERIES_SUSPENDED_SEDIMENT_CODE = "suspended";

// CShore codes
string const WAVEENERGYFLUX = "wave_energy_flux";
string const WAVEHEIGHTX = "WAVEHEIGHTX.csv";
string const WAVEHEIGHTY = "WAVEHEIGHTY.csv";
string const ACTIVEZONE = "ACTIVEZONE.csv";

//================================================ Globally-available functions =================================================
template <class T>
T tMax(T a, T b)
{
    return ((a > b) ? a : b);
}

template <class T>
T tMax(T a, T b, T c)
{
    T max = (a < b) ? b : a;
    return ((max < c) ? c : max);
}

template <class T>
T tMin(T a, T b)
{
    return ((a < b) ? a : b);
}

template <class T>
T tMin(T a, T b, T c)
{
    return (a < b ? (a < c ? a : c) : (b < c ? b : c));
}

template <class T>
T tAbs(T a)
{
    // From a posting dated 18 Nov 93 by rmartin@rcmcon.com (Robert Martin), archived in cpp_tips
    return ((a < 0) ? -a : a);
}

template <class T>
bool bIsBetween(T a, T b, T c)
{
    // Assumes b > c
    return ((a >= b) && (a <= c));
}

template <typename T>
string strNumToStr(const T &t)
{
    // From http://stackoverflow.com/questions/2125880/convert-float-to-stdstring-in-c
    ostringstream os;
    os << t;
    return os.str();
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS

// Definitions are in utilsglobal.cpp
double dRound(double const);
int nRound(double const);
// bool bIsWhole(double const);
bool bIsStringValidDouble(string &);
bool bIsStringValidInt(string &);
struct FillToWidth
{
    FillToWidth(char f, int w) : chFill(f), nWidth(w)
    {
    }
    char chFill;
    int nWidth;
};
ostream &operator<<(ostream &, const FillToWidth &);

#endif

//================================================= debugging stuff =============================================================
//#define CLOCKCHECK          // Uncomment to check CPU clock rollover settings
//#define RANDCHECK           // Uncomment to check randomness of random number generator

#endif // CME_H
