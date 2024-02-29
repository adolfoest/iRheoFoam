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

#include "isothermModel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::isothermModel> Foam::isothermModel::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    word modelType("idealGas");

    if (!dict.readIfPresent("isothermModel", modelType))
    {
        Warning
            << "Entry '"
            << "isothermModel" << "' not found in dictionary "
            << dict.name() << nl
            << "using default" << nl;
    }

    Info<< "Selecting isothermModel: " << modelType << endl;

    auto* ctorPtr = componentsConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "isothermModel",
            modelType,
            *componentsConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<isothermModel>(ctorPtr(mesh, dict));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
