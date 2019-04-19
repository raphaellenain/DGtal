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
 * @file testRotation2D.cpp
 * @ingroup Tests
 * @author Raphael Lenain (\c raphael.lenain@univ-poitiers.fr )
 * Laboratoire XLIM - Axe ASALI - Equipe IG, Poitiers, France
 *
 * @date 2019/04/16
 *
 * Functions for testing class testRotation2D.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <cmath>
#include <DGtal/images/ImageSelector.h>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include "DGtal/images/ConstImageAdapter.h"
#include "DGtal/base/Common.h"
#include "ConfigTest.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/DiscreteTransformation2D.h"
#include "DGtal/io/readers/PGMReader.h"
#include "DGtal/io/writers/GenericWriter.h"

///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
using namespace functors;
using namespace Z2i;

///////////////////////////////////////////////////////////////////////////////
// Functions for testing class testRotation2D.
///////////////////////////////////////////////////////////////////////////////
class testRotation2D
{
  public:
    bool rotationTest (int testNumber, DGtal::Z2i::Point initialPoint, DGtal::Z2i::Point origin, double parameter, double angle, DGtal::Z2i::Point expectedPoint)
    {
      DGtal::Z2i::Point outputPoint;
      outputPoint = DGtal::functors::DiscreteRotation2D<DGtal::Z2i::Space>(origin, parameter, angle)(initialPoint);
      bool out = outputPoint[0] == expectedPoint[0] && outputPoint[1] == expectedPoint[1];
      if (!out)
      {
        trace.emphase() << "Error on test #" << testNumber << endl << " -> Expected : {" << expectedPoint[0] << ";" << expectedPoint[1] << "} - Result : "<< "{" << outputPoint[0] << ";" << outputPoint[1] << "}"<< endl;
      }
      return out;
    }
};

///////////////////////////////////////////////////////////////////////////////
// Standard services - public :

int main( int, char** )
{
  bool res = true;
  testRotation2D rotation2DTest;
  trace.beginBlock ( "Testing Rotation2D" );
    // Test 1
    res &= rotation2DTest.rotationTest(1, DGtal::Z2i::Point(1.0, -1.0), DGtal::Z2i::Point(0.0, 0.0), 0.0, M_PI / 2, DGtal::Z2i::Point(1.0, 1.0));
    // Test 2
    res &= rotation2DTest.rotationTest(2, DGtal::Z2i::Point(-1.0, -1.0), DGtal::Z2i::Point(0.0, 0.0), 0.0, M_PI / 2, DGtal::Z2i::Point(1.0, -1.0));
    // Test 3
    res &= rotation2DTest.rotationTest(3, DGtal::Z2i::Point(-1.0, 0.0), DGtal::Z2i::Point(0.0, 0.0), 0.0, M_PI / 2, DGtal::Z2i::Point(0.0, -1.0));
    // Test 4
    res &= rotation2DTest.rotationTest(4, DGtal::Z2i::Point(0.0, -1.0), DGtal::Z2i::Point(0.0, 0.0), 0.0, M_PI / 2, DGtal::Z2i::Point(1.0, 0.0));
    // Test 5
    res &= rotation2DTest.rotationTest(5, DGtal::Z2i::Point(0.0, -2.0), DGtal::Z2i::Point(0.0, 1.0), -M_PI / 4, M_PI / 4, DGtal::Z2i::Point(2.0, -1.0));
    // Test 6
    res &= rotation2DTest.rotationTest(6, DGtal::Z2i::Point(2.0, -1.0), DGtal::Z2i::Point(0.0, 1.0), -M_PI / 4, -M_PI / 4, DGtal::Z2i::Point(1.0, -2.0));
    // Test 7
    res &= rotation2DTest.rotationTest(7, DGtal::Z2i::Point(0.0, -2.0), DGtal::Z2i::Point(0.0, 1.0), 7 * M_PI / 4, M_PI / 2, DGtal::Z2i::Point(3.0, 1.0));
    // Test 8
    res &= rotation2DTest.rotationTest(8, DGtal::Z2i::Point(3.0, 1.0), DGtal::Z2i::Point(0.0, 1.0), 7 * M_PI / 4, -M_PI / 2, DGtal::Z2i::Point(0.0, -2.0));
    // Test 9
    res &= rotation2DTest.rotationTest(9, DGtal::Z2i::Point(0.0, -2.0), DGtal::Z2i::Point(0.0, 1.0), 3 * M_PI / 4, M_PI / 8, DGtal::Z2i::Point(2.0, -1.0));
    // Test 10
    res &= rotation2DTest.rotationTest(10, DGtal::Z2i::Point(2.0, -1.0), DGtal::Z2i::Point(0.0, 1.0), 3 * M_PI / 4, -M_PI / 8, DGtal::Z2i::Point(1.0, -1.0));
  trace.emphase() << ( res ? "Passed." : "Error." ) << endl;
  trace.endBlock();
  return res ? 0 : 1;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
