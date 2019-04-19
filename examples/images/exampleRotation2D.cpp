/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

/**
 * @file images/exampleRotation2D.cpp
 * @ingroup Examples
 * @author Raphael Lenain (\c raphael.lenain@univ-poitiers.fr )
 * Laboratoire XLIM - Axe ASALI - Equipe IG, Poitiers, France
 *
 * @date 2019/04/16
 *
 * An example file named rotation2D.
 *
 * This file is part of the DGtal library.
 */

/**
*  Example of 2D rotation.
   @see @ref moduleDiscreteTransform
   \image html church_rotated.jpg "Result rotation" 
*  \example images/exampleRotation2D.cpp
**/


///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <cmath>
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/ConstImageAdapter.h"
#include "ConfigExamples.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/Common.h"
#include "DGtal/io/readers/PGMReader.h"
#include "DGtal/io/writers/GenericWriter.h"
//! [include]
#include "DGtal/images/DiscreteTransformation2D.h"
//! [include]
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
using namespace functors;
using namespace Z2i;

///////////////////////////////////////////////////////////////////////////////

int main( int , char** )
{
  typedef ImageSelector<Domain, unsigned char >::Type Image;
  //! [def]
  typedef DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain, unsigned char> Image;
  typedef Image::Domain Domain;
  typedef Domain::Point Point;
  //! [def]
  std::string name = "church";
  trace.beginBlock ( "Example Rotation2D of " + name);
    std::string filename = examplesPath + "samples/" + name + ".pgm";
    Image image = DGtal::PGMReader<Image>::importPGM(filename);

    Domain outputDomain(image.domain().lowerBound(), image.domain().upperBound());
    Image outputImage(outputDomain);
    Point origin((image.domain().lowerBound()[0] + image.domain().upperBound()[0]) / 2,
      (image.domain().lowerBound()[1] + image.domain().upperBound()[1]) / 2);
    double angle = M_PI / 8;
    double parameter = 0.0;
    DGtal::functors::DiscreteRotation2D<DGtal::Z2i::Space> rotation(origin, parameter, angle);
    
    for (Domain::ConstIterator iterator = outputDomain.begin(), iteratorEnd = outputDomain.end(); iteratorEnd != iterator; iterator++)
    {
      Point originalPoint = *iterator, rotatedPoint = rotation(originalPoint);
      unsigned char valeurPoint = !image.domain().isInside(rotatedPoint) ? DGtal::Color::Black.getRGB() : image(rotatedPoint);
      outputImage.setValue(originalPoint, valeurPoint);
    }
    outputImage >> name + "_rotated.pgm";
  trace.endBlock();
  return 0;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
