/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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
    Foam::interfaceVelocityFvPatchVectorField

Group
    grpWallBoundaryConditions

Description
    This boundary condition provides a velocity condition for cases with
    a planar interface.

    Reference:
    \verbatim
        Esteban A, Hernández J, Tajuelo J, Rubio MÁ, 
        A numerical model for the simulation of complex 
        planar Newtonian interfaces, 
        Applied Mathematical Modelling, 2025, 
        doi 10.1016/j.apm.2024.115653.
    \endverbatim

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            interfaceVelocity;
        value           uniform (0 0 0);    // Initial value
        Us              Us;
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchVectorField
    Foam::MRFZone

SourceFiles
    interfaceVelocityFvPatchVectorField.C

Author
    Adolfo Esteban, UNED (2023)


\*---------------------------------------------------------------------------*/

#ifndef interfaceVelocityFvPatchVectorField_H
#define interfaceVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "areaFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class interfaceVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class interfaceVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // private data
    word UsName_;
public:

    //- Runtime type information
    TypeName("interfaceVelocity");


    // Constructors

        //- Construct from patch and internal field
        interfaceVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        interfaceVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given interfaceVelocityFvPatchVectorField
        //  onto a new patch
        interfaceVelocityFvPatchVectorField
        (
            const interfaceVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        interfaceVelocityFvPatchVectorField
        (
            const interfaceVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new interfaceVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        interfaceVelocityFvPatchVectorField
        (
            const interfaceVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new interfaceVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
