/*---------------------------------------------------------------------------*\
License
    This file is part of iRheoFoam which is an extension to OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "rheoResiduals.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rheoResiduals::rheoResiduals
(
    const scalar& tol
)
:
    tol_(tol),
    gammaHat0_(0.0),
    UsHat0_(0.0),
    gammaHat_(0.0),
    UsHat_(0.0),
    hasConverged_(false),
    gRes_(GREAT),
    uRes_(GREAT)   
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::rheoResiduals::calculateInit
(
    const areaScalarField& gamma,
    const areaVectorField& Us
)
{
    scalar maxg = max(gamma.internalField()).value();
    scalar maxUs = max(max(mag(Us.internalField())).value(), SMALL);
    reduce(maxg, maxOp<scalar>());
    reduce(maxUs, maxOp<scalar>());

    scalar avg = gAverage(gamma.internalField());
    scalar avUs = mag(gAverage(Us.internalField()));

    gammaHat0_ = avg/maxg;
    UsHat0_ = avUs/maxUs;
    hasConverged_ = false;
}

void Foam::rheoResiduals::calculate
(
    const areaScalarField& gamma,
    const areaVectorField& Us
)
{
    scalar maxg = max(gamma.internalField()).value();
    scalar maxUs = max(max(mag(Us.internalField())).value(), SMALL);
    reduce(maxg, maxOp<scalar>());
    reduce(maxUs, maxOp<scalar>());

    scalar avg = gAverage(gamma.internalField());
    scalar avUs = mag(gAverage(Us.internalField()));

    gammaHat_ = avg/maxg;
    UsHat_ = avUs/maxUs;

    gRes_ = mag(gammaHat_ - gammaHat0_);
    uRes_ = mag(UsHat_ - UsHat0_);
    if (gRes_ < tol_ && uRes_ < tol_)
    {
        hasConverged_ = true;
    }
}

// ************************************************************************* //
