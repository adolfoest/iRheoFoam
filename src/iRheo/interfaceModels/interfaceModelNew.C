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

#include "interfaceModel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::interfaceModel> Foam::interfaceModel::New
(
    const fvMesh& mesh,
    const word& dictName
)
{
    dictionary dict
    (
        IOdictionary
        (
            IOobject
            (
                dictName,
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        )
    );

    word modelType("BoussinesqScriven");

    if (!dict.readIfPresent("interfaceModel", modelType))
    {
        Warning
            << "Entry '"
            << "interfaceModel" << "' not found in dictionary "
            << dict.name() << nl
            << "using default" << nl;
    }

    Info<< "Selecting interfaceModel: " << modelType << endl;

    auto* ctorPtr = componentsConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "interfaceModel",
            modelType,
            *componentsConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<interfaceModel>(ctorPtr(mesh, dictName));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
