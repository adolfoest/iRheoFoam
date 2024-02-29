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

#include "idealGasIsothermModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace isothermModels
{
    defineTypeNameAndDebug(idealGasIsothermModel, 0);

    addToRunTimeSelectionTable
    (
        isothermModel,
        idealGasIsothermModel,
        components
    );
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isothermModels::idealGasIsothermModel::idealGasIsothermModel
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
)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::isothermModels::idealGasIsothermModel::~idealGasIsothermModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::areaScalarField>
Foam::isothermModels::idealGasIsothermModel::KPi
(
    const areaScalarField& gamma
) const
{
    return tmp<areaScalarField>::New
    (
        "KPi", gamma
    );
}

Foam::tmp<Foam::areaScalarField>
Foam::isothermModels::idealGasIsothermModel::Pi
(
    const areaScalarField& gamma,
    const scalar& Ma
) const
{
    return tmp<areaScalarField>::New
    (
        "Pi", (Ma*gamma)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
