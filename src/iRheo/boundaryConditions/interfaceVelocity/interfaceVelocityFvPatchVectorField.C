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

\*---------------------------------------------------------------------------*/

#include "interfaceVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceVelocityFvPatchVectorField::
interfaceVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    UsName_("Us")
{}


Foam::interfaceVelocityFvPatchVectorField::
interfaceVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    UsName_(dict.getOrDefault<word>("Us", "Us"))
{}


Foam::interfaceVelocityFvPatchVectorField::
interfaceVelocityFvPatchVectorField
(
    const interfaceVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    UsName_(ptf.UsName_)
{}


Foam::interfaceVelocityFvPatchVectorField::
interfaceVelocityFvPatchVectorField
(
    const interfaceVelocityFvPatchVectorField& ivpvf
)
:
    fixedValueFvPatchVectorField(ivpvf),
    UsName_(ivpvf.UsName_)
{}


Foam::interfaceVelocityFvPatchVectorField::
interfaceVelocityFvPatchVectorField
(
    const interfaceVelocityFvPatchVectorField& ivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ivpvf, iF),
    UsName_(ivpvf.UsName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::interfaceVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    const areaVectorField& Us =
        db().lookupObject<areaVectorField>("Us");
    vectorField usv = Us;

    vectorField::operator=(usv);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::interfaceVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeEntryIfDifferent<word>("Us", "Us", UsName_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        interfaceVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
