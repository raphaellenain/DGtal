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

#pragma once

/**
* @file DiscreteTransformation2D.h
* @author Eric Andres (\c eric.andres@univ-poitiers.fr )
* Laboratoire XLIM - Axe ASALI - Equipe IG, Poitiers, France
*
* @date 2019/04/16
*
* This file is part of the DGtal library.
*/

#if defined(DiscreteTransformation2D_RECURSES)
#error Recursive header files inclusion detected in DiscreteTransformation2D.h
#else // defined(DiscreteTransformation2D_RECURSES)
/** Prevents recursive inclusion of headers. */
#define DiscreteTransformation2D_RECURSES

#if !defined DiscreteTransformation2D_h
/** Prevents repeated inclusion of headers. */
#define DiscreteTransformation2D_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <cmath>
#include <climits>
#include <utility>
#include "DGtal/base/Common.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/kernel/domains/CDomain.h>
#include <DGtal/kernel/CSpace.h>
//////////////////////////////////////////////////////////////////////////////
#ifndef M_PI_4
#define M_PI_4     0.785398163397448309616  // pi/4
#endif
namespace DGtal
{
	namespace functors
	{
		/////////////////////////////////////////////////////////////////////////////
		// Template class DiscreteReflection2D
		/**
		* Description of template functor like class 'DiscreteReflection2D' <p>
		* \brief Aim: implements reflection of a point in the 2D integer space.
		*
		* @tparam TSpace a 2 dimensional space.
		*
		* @see exampleReflection2D.cpp
		*/
		template <typename TSpace>
		class DiscreteReflection2D
		{
			///Checking concepts
			BOOST_CONCEPT_ASSERT((concepts::CSpace<TSpace>));
			BOOST_STATIC_ASSERT((TSpace::dimension == 2));
			BOOST_STATIC_ASSERT((TSpace::Point::dimension == 2));

			// ----------------------- Types ------------------------------
		public:
			typedef typename TSpace::Point Point;

			// ----------------------- Interface --------------------------------------
		public:
			/**
			* Constructor.
			* @param aOrigin  the center of reflection.
			* @param anAngle  the angle given in radians.
			*/
			DiscreteReflection2D(const Point & aOrigin, const double & anAngle)
				:origin(aOrigin), angle(DiscreteReflection2D::modulo2Pi(anAngle))
			{
				this->sinusAng = DiscreteReflection2D::sinus(angle);
				this->cosinusAng = DiscreteReflection2D::cosinus(angle);
				if (this->angle <= M_PI_4)
				{
					this->orthogonalLinesDivider = this->cosinusAng;
					this->computingFunction = &DiscreteReflection2D::computeX;
					this->indexRoundedCoordinate = 1;
					this->indexComputedCoordinate = 0;
				}
				else
				{
					this->orthogonalLinesDivider = this->sinusAng;
					this->computingFunction = &DiscreteReflection2D::computeY;
					this->indexRoundedCoordinate = 0;
					this->indexComputedCoordinate = 1;
				}
			}

			/**
			* Operator
			* @param point  the point to reflect.
			* @return the transformed point.
			*/
			inline
				Point operator()(const Point & point) const
			{
				Point ret;
				Point pt1, pt2;
				//Computing the index of the orthogonal lines
				double orthogonalLinesIndex = floor(0.5 + (this->cosinusAng * (point[0] - this->origin[0]) + this->sinusAng * (point[1] - this->origin[1])) / this->orthogonalLinesDivider);
				pt1 = this->computePoint(orthogonalLinesIndex, ceil);
				pt2 = this->computePoint(orthogonalLinesIndex, floor);
				double diffPt1 = this->computeDiffPoint(pt1);
				double diffPt2 = this->computeDiffPoint(pt2);
				if (-this->orthogonalLinesDivider / 2.0 <= diffPt1 && diffPt1 < this->orthogonalLinesDivider / 2.0)
				{
					// pt1 is the intersection point which exists in this case
					ret[this->indexRoundedCoordinate] = 2.0 * pt1[this->indexRoundedCoordinate] - point[this->indexRoundedCoordinate];
				}
				else if (-this->orthogonalLinesDivider / 2.0 <= diffPt2 && diffPt2 < this->orthogonalLinesDivider / 2.0)
				{
					// pt2 is the intersection point which exists in this case
					ret[this->indexRoundedCoordinate] = 2.0 * pt2[this->indexRoundedCoordinate] - point[this->indexRoundedCoordinate];
				}
				else
				{
					//there is no intersection point
					ret[this->indexRoundedCoordinate] = (pt1[this->indexRoundedCoordinate] + pt2[this->indexRoundedCoordinate]) - point[this->indexRoundedCoordinate];
				}

				ret[this->indexComputedCoordinate] = this->computingFunction(ret[this->indexRoundedCoordinate], orthogonalLinesIndex, this->origin, this->cosinusAng, this->sinusAng);

				return ret;
			}

			// ------------------------- Protected Methods ------------------------------
		protected:
			//=============================================================
			//Compute the difference of a point
			double computeDiffPoint(Point point) const
			{
				return this->sinusAng * (point[0] - this->origin[0]) - this->cosinusAng * (point[1] - this->origin[1]);
			}

			//=============================================================
			//Compute X knowing Y
			static double computeX(double y, double orthogonalLinesIndex, Point origin, double cosinusAng, double sinusAng)
			{
				return ceil(origin[0] + ((2.0 * orthogonalLinesIndex - 1.0) * cosinusAng - 2.0 * sinusAng * (y - origin[1])) / (2.0 * cosinusAng));
			}

			//=============================================================
			//Compute Y knowing X
			static double computeY(double x, double orthogonalLinesIndex, Point origin, double cosinusAng, double sinusAng)
			{
				return ceil(origin[1] + ((2.0 * orthogonalLinesIndex - 1.0) * sinusAng - 2.0 * cosinusAng * (x - origin[0])) / (2.0 * sinusAng));
			}

			//=============================================================
			//Compute a point
			Point computePoint(double orthogonalLinesIndex, double(*roundFunction)(double)) const
			{
				Point retour;
				retour[this->indexRoundedCoordinate] = roundFunction(this->sinusAng * this->cosinusAng * orthogonalLinesIndex + this->origin[this->indexRoundedCoordinate]);
				retour[this->indexComputedCoordinate] = this->computingFunction(retour[this->indexRoundedCoordinate], orthogonalLinesIndex, this->origin, this->cosinusAng, this->sinusAng);

				return retour;
			}

			//=============================================================
			//Compute angle % 2PI
			static double modulo2Pi(double angle)
			{
				double out = angle;
				double added = 2 * M_PI;
				while (out < 0)
				{
					// Adding 2 * PI
					out += added;
				}
				while (out >= added)
				{
					// Substracting 2 * PI
					out -= added;
				}
				return out;
			}

			//=============================================================
			//Compute cosinus(angle) with rounded values of M_PI suggested by the library
			// angle has to be "moduloed" with 2PI
			static double cosinus(double angle)
			{
				double out = 0;
				if (angle == 0 || angle == 2 * M_PI)
				{
					out = 1.0;
				}
				else if (angle == M_PI_2 || angle == 3 * M_PI_2)
				{
					out = 0.0;
				}
				else if (angle == M_PI)
				{
					out = -1.0;
				}
				else
				{
					out = std::cos(angle);
				}
				return out;
			}

			//=============================================================
			//Compute sinus(angle) with rounded values of M_PI suggested by the library
			// angle has to be "moduloed" with 2PI
			static double sinus(double angle)
			{
				double out = 0;
				if (angle == 0 || angle == 2 * M_PI || angle == M_PI)
				{
					out = 0.0;
				}
				else if (angle == M_PI_2)
				{
					out = 1.0;
				}
				else if (angle == 3 * M_PI_2)
				{
					out = -1.0;
				}
				else
				{
					out = std::sin(angle);
				}
				return out;
			}

			// ------------------------- Protected Datas ------------------------------
		protected:
			Point origin;
			double sinusAng;
			double cosinusAng;
			double angle;
			double orthogonalLinesDivider;
			int indexRoundedCoordinate;
			int indexComputedCoordinate;
			double (*computingFunction)(double, double, Point, double, double);
		};


		/////////////////////////////////////////////////////////////////////////////
		// Template class DiscreteRotation2D
		/**
		* Description of template functor like class 'DiscreteRotation2D' <p>
		* \brief Aim: implements rotation of a point in the 2D integer space.
		*
		* @tparam TSpace a 2 dimensional space.
		*
		* @see exampleRotation2D.cpp
		*/
		template <typename TSpace>
		class DiscreteRotation2D
		{
			///Checking concepts
			BOOST_CONCEPT_ASSERT((concepts::CSpace<TSpace>));
			BOOST_STATIC_ASSERT((TSpace::dimension == 2));
			BOOST_STATIC_ASSERT((TSpace::Point::dimension == 2));

			// ----------------------- Types ------------------------------
		public:
			typedef typename TSpace::Point Point;
			// ----------------------- Interface --------------------------------------
		public:
			/**
			* Constructor.
			* @param aOrigin  the center of rotation.
			* @param parametricAngle  the parametric angle given in radians.
			* @param anAngle  the angle given in radians.
			*/
			DiscreteRotation2D(const Point & aOrigin, const double & parametricAngle, const double & anAngle)
			{
				this->firstReflection = new DiscreteReflection2D<TSpace>(aOrigin, parametricAngle);
				this->secondReflection = new DiscreteReflection2D<TSpace>(aOrigin, parametricAngle + (anAngle / 2.0));
			}

			/**
			* Destructor.
			*/
			~DiscreteRotation2D()
			{
				delete this->firstReflection;
				delete this->secondReflection;
			}

			/**
			* Operator
			* @param point  the point to rotate.
			* @return the transformed point.
			*/
			inline
				Point operator()(const Point & point) const
			{
				return (*this->secondReflection)((*this->firstReflection)(point));
			}

			// ------------------------- Protected Datas ------------------------------
		protected:
			DiscreteReflection2D<TSpace> *firstReflection, *secondReflection;
		};

	}// namespace DGtal::functors
}// namespace DGtal

#endif // !defined DiscreteTransformation2D_h

#undef DiscreteTransformation2D_RECURSES
#endif // else defined(DiscreteTransformation2D_RECURSES)

