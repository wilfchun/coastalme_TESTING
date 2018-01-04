/*!
 *
 * \file 2di_point.cpp
 * \brief Geometry class used to represent 2D point objects with integer co-ordinates
 * \details The CGeom2DIPoint geometry class is used to represent 2D points where the x and y co-ordinates can only be integer values, e.g. points for which the x and y co-ordinates are in the raster-grid CRS (co-ordinate reference system)
 * \author David Favis-Mortlock
 * \author Andres Payo
 * \date 2018
 * \copyright GNU General Public License
 *
 */

/*===============================================================================================================================

 This file is part of CoastalME, the Coastal Modelling Environment.

 CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

===============================================================================================================================*/
#include "2di_point.h"


//! Constructor with no parameters (the X and Y co-ordinates of the CGeom2DIPoint object are set to zero)
CGeom2DIPoint::CGeom2DIPoint(void)
:  nX(0),
   nY(0)
{
}

//! Constructor with two integer parameters, for the X and Y co-ordinates of the CGeom2DIPoint object
CGeom2DIPoint::CGeom2DIPoint(int const nNewX, int const nNewY)
:  nX(nNewX),
   nY(nNewY)
{
}

//! CGeom2DIPoint copy constructor
CGeom2DIPoint::CGeom2DIPoint(CGeom2DIPoint const& Pti)
{
   nX = Pti.nX;
   nY = Pti.nY;
}


//! Returns the CGeom2DIPoint object's integer X co-ordinate
int CGeom2DIPoint::nGetX(void) const
{
   return nX;
}

//! Returns the CGeom2DIPoint object's integer Y co-ordinate
int CGeom2DIPoint::nGetY(void) const
{
   return nY;
}

//! Returns a reference to the CGeom2DIPoint object's integer X co-ordinate
int* CGeom2DIPoint::pnGetX(void)
{
   return &nX;
}

//! Returns a reference to the CGeom2DIPoint object's integer Y co-ordinate
int* CGeom2DIPoint::pnGetY(void)
{
   return &nY;
}

//! The integer parameter sets a value for the CGeom2DIPoint object's X co-ordinate
void CGeom2DIPoint::SetX(int const nNewX)
{
   nX = nNewX;
}

//! The integer parameter sets a value for the CGeom2DIPoint object's Y co-ordinate
void CGeom2DIPoint::SetY(int const nNewY)
{
   nY = nNewY;
}

//! The two integer parameters set values for the CGeom2DIPoint object's X and Y co-ordinates
// void CGeom2DIPoint::SetXY(int const nNewX, int const nNewY)
// {
//    nX = nNewX;
//    nY = nNewY;
// }

//! The parameter is a pointer to a CGeom2DIPoint object, this is used to set values for the CGeom2DIPoint object's X and Y co-ordinates
// void CGeom2DIPoint::SetXY(CGeom2DIPoint const* Pti)
// {
//    nX = Pti->nGetX();
//    nY = Pti->nGetY();
// }


//! Adds the first parameter to the CGeom2DIPoint object's X co-ordinate, adds the second parameter to the CGeom2DIPoint object's Y co-ordinate
void CGeom2DIPoint::AddXAddY(int const nXToAdd, int const nYToAdd)
{
   nX += nXToAdd;
   nY += nYToAdd;
}


//! Sets one CGeom2DIPoint object to be the same as another CGeom2DIPoint object
void CGeom2DIPoint::operator= (CGeom2DIPoint* Pti)
{
   nX = Pti->nGetX();
   nY = Pti->nGetY();
}

// //! Returns true if a pointed-to CGeom2DIPoint object is the same as this CGeom2DIPoint object, returns false otherwise
// bool CGeom2DIPoint::operator== (CGeom2DIPoint* Pti) const
// {
//    if ((Pti->nGetX() == nX) && (Pti->nGetY() == nY))
//       return true;
//
//    return false;
// }

//! Returns true if a CGeom2DIPoint object is the same as this CGeom2DIPoint object, returns false otherwise
bool CGeom2DIPoint::operator== (CGeom2DIPoint Pti) const
{
   if ((Pti.nGetX() == nX) && (Pti.nGetY() == nY))
      return true;

   return false;
}
