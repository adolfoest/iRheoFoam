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

#include "vanderWaalsIsothermModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace isothermModels
{
    defineTypeNameAndDebug(vanderWaalsIsothermModel, 0);

    addToRunTimeSelectionTable
    (
        isothermModel,
        vanderWaalsIsothermModel,
        components
    );
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isothermModels::vanderWaalsIsothermModel::vanderWaalsIsothermModel
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    isothermModel
    (
        typeName,
        mesh,
        dict
    ),
    subDict_(dict.subDict(typeName + "Coeffs")),
    gamma0_(subDict_.get<scalar>("gamma0")),
    gammaInf_(subDict_.get<scalar>("gammaInf")),
    T_(subDict_.get<scalar>("T")),
    Beta_(subDict_.get<scalar>("Beta"))
{
    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::isothermModels::vanderWaalsIsothermModel::~vanderWaalsIsothermModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::areaScalarField>
Foam::isothermModels::vanderWaalsIsothermModel::KPi
(
    const areaScalarField& gamma
) const
{
    return
        tmp<areaScalarField>::New
        (
            IOobject
            (
                "KPi",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            (kB*T_*gamma0_*gamma)/sqr(1. - (gamma0_*gamma)/gammaInf_)
            - Beta_*sqr(gamma0_*gamma)
        );
}

Foam::tmp<Foam::areaScalarField>
Foam::isothermModels::vanderWaalsIsothermModel::Pi
(
    const areaScalarField& gamma,
    const scalar& Ma
) const
{
    return
        tmp<areaScalarField>::New
        (
            IOobject
            (
                "KPi",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            (kB*T_*gamma0_*gamma)/(1. - (gamma0_*gamma)/gammaInf_)
            - Beta_*sqr(gamma0_*gamma)
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
