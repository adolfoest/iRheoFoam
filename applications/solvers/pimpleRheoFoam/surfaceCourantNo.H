/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
-------------------------------------------------------------------------------
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

Global
    surfaceCourantNo

Author
    Hrvoje Jasak, Wikki Ltd.

Description
    Calculates and outputs the mean and maximum Courant Numbers for the
    Finite Area method.

\*---------------------------------------------------------------------------*/

scalar CoNums = 0.0;
scalar meanCoNums = 0.0;
scalar velMags = 0.0;

if (aMesh.nInternalEdges())
{
    edgeScalarField SfUfbyDelta
    (
        aMesh.edgeInterpolation::deltaCoeffs()*mag(phis)
    );

    CoNums = max(SfUfbyDelta/aMesh.magLe())
        .value()*runTime.deltaT().value();

    meanCoNums = (sum(SfUfbyDelta)/sum(aMesh.magLe()))
        .value()*runTime.deltaT().value();

    velMags = max(mag(phis)/aMesh.magLe()).value();
}

Info<< "Interface Courant Number mean: " << meanCoNums
    << " max: " << CoNums
    << " velocity magnitude: " << velMags << endl;

CoNum = max(CoNum, CoNums);
// ************************************************************************* //
