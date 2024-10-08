/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenFOAM Foundation
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

Class
    Foam::CorrectPhiRheo

Description
    Minimal modification of CorrectPhi class to include pimpleControlRheo class
    instead of the original pimpleControl class.

    Reference:
    \verbatim
        Esteban A, Hernández J, Tajuelo J, Rubio MÁ, 
        A numerical model for the simulation of complex 
        planar Newtonian interfaces, 
        Applied Mathematical Modelling, 2025, 
        doi 10.1016/j.apm.2024.115653.
    \endverbatim

SourceFiles
    CorrectPhiRheo.C

\*---------------------------------------------------------------------------*/

#ifndef CorrectPhiRheo_H
#define CorrectPhiRheo_H
 
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
    class pimpleControlRheo;
 
    //- If the mesh is moving correct the velocity BCs on the moving walls to
    //  ensure the corrected fluxes and velocity are consistent
    void correctUphiBCsRheo
    (
        volVectorField& U,
        surfaceScalarField& phi
    );
 
    //- If the mesh is moving correct the velocity BCs on the moving walls to
    //  ensure the corrected fluxes and velocity are consistent
    void correctUphiBCsRheo
    (
        const volScalarField& rho,
        volVectorField& U,
        surfaceScalarField& phi
    );
 
    template<class RAUfType, class DivUType>
    void CorrectPhiRheo
    (
        volVectorField& U,
        surfaceScalarField& phi,
        const volScalarField& p,
        const RAUfType& rAUf,
        const DivUType& divU,
        pimpleControlRheo& pimple
    );
 
    template<class RAUfType, class DivRhoUType>
    void CorrectPhiRheo
    (
        volVectorField& U,
        surfaceScalarField& phi,
        const volScalarField& p,
        const volScalarField& rho,
        const volScalarField& psi,
        const RAUfType& rAUf,
        const DivRhoUType& divRhoU,
        pimpleControlRheo& pimple
    );
}
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
#ifdef NoRepository
    #include "CorrectPhiRheo.C"
#endif
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
#endif
 
// ************************************************************************* //