/*!
 *
 * \file utils_global.cpp
 * \brief Globally-available utility routines
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
#include <cmath>
#include <cfloat>

#include <sstream>
using std::stringstream;

#include <iomanip>
using std::setw;

#include "cme.h"


/*==============================================================================================================================

 Correctly rounds doubles

==============================================================================================================================*/
double dRound(double const d)
{
   // Rounds positive or negative doubles correctly
   return ((d < 0.0) ? ceil(d - 0.5) : floor(d + 0.5));
}


// bool bIsWhole(double d)
// {
//    // From http://answers.yahoo.com/question/index?qid=20110320132617AAMdb7u
//    return (static_cast<int>(d) == d);
// }


/*==============================================================================================================================

 Checks a double to see if it is NaN. From http://www.johndcook.com/blog/IEEE_exceptions_in_cpp/

==============================================================================================================================*/
bool bDoubleIsValid(double const dX)
{
   // This looks like it should always be true, but it is false if dX is a NaN
   return (dX == dX);
}


/*==============================================================================================================================

 Checks a double to see if it is finite. From http://www.johndcook.com/blog/IEEE_exceptions_in_cpp/

==============================================================================================================================*/
// bool bIsFinite(double const dX)
// {
//    return (dX <= DBL_MAX && dX >= -DBL_MAX);
// }


/*==============================================================================================================================

 Operator that inserts a given fill character, to a given width, into an output stream. From http://stackoverflow.com/questions/2839592/equivalent-of-02d-with-stdstringstream

==============================================================================================================================*/
ostream& operator<< (ostream& ostr, const FillToWidth& args)
{
   ostr.fill(args.chFill);
   ostr.width(args.nWidth);

   return ostr;
}

