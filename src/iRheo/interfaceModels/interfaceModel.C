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


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(interfaceModel, 0);
    defineRunTimeSelectionTable(interfaceModel, components);
}

const Foam::word Foam::interfaceModel::dictName("transportProperties");

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
void Foam::interfaceModel::calcIs()
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const polyPatch& pp = pbm[patchID_];
    if (isA<emptyPolyPatch>(pp) || pp.empty())
    {
        n_ = modelDict().get<vector>("normal");
    }
    else
    {
        const fvPatch& patch = mesh_.boundary()[patchID_];
        const vectorField nf = patch.nf();
        n_ = nf[0];
    }
    
    Is_ = dimensionedTensor("Is", dimless, tensor::I - sqr(n_));

    /*const areaVectorField nAvf =
    aMesh.faceAreaNormals() / mag(aMesh.faceAreaNormals());
    const vectorField n = nAvf.internalField();
    areaTensorField Is
    (
        IOobject
        (
            "Is",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh,
        dimensionedTensor("Is", dimless, tensor::I)
    );
    tensorField& Ispfr = Is.primitiveFieldRef();
    Ispfr = Ispfr - sqr(n);

    const faBoundaryMesh& aPatches = aMesh.boundary();
    forAll(aPatches, patchI)
    {
        const vectorField& nb = nAvf.boundaryField()[patchI];
        tensorField& Isb = Is.boundaryFieldRef()[patchI];
        Isb = Isb - sqr(nb);
    }
    */
}

Foam::areaTensorField Foam::interfaceModel::Hadamard(const areaTensorField& t) const
{
    const vector e({1,1,1});
    areaTensorField th = t;
    forAll(t, i)
    {
        vector d = th[i].diag();
        tensor de = d*e;
        th[i] -= de.schur(Is_.value());
    }
    const faMesh& mesh = t.mesh();
    const faBoundaryMesh& patches = mesh.boundary();
    forAll(patches, patchI)
    {
        tensorField& tb = th.boundaryFieldRef()[patchI];
        forAll(tb, i)
        {
            vector d = tb[i].diag();
            tensor de = d*e;
            tb[i] -= de.schur(Is_.value());
        }
    }
    return th;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceModel::interfaceModel
(
    const word& type,
    const fvMesh& mesh,
    const word& dictName
)
:
    IOdictionary
    (
        IOobject
        (
            "interfaceModel",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    interfaceModelCoeffs_(mesh.lookupObject<IOdictionary>(dictName)),
    mesh_(mesh),
    Ma_("Ma", dimArea/sqr(dimTime), 0.0, modelDict()),
    Bq_("Bq", dimViscosity, 0.0, modelDict()),
    Bqk_("Bqk", dimViscosity, 0.0, modelDict()),
    Theta_("Theta", dimViscosity, 0.0, modelDict()),
    Tr_("Tr", dimViscosity, 0.0, modelDict()),
    Pes_("Pes", dimTime/dimArea, 0.0, modelDict()),
    Re_("Re", dimless, 1.0, modelDict()),
    n_({0,0,1}),
    Is_(tensor::I),
    patchID_(-1),
    coupledToBulk_(true),
    isotherm_(nullptr)
{
    // get interface patch
    word interfacePatch(modelDict().get<word>("interfacePatch"));
    patchID_ = mesh_.boundaryMesh().findPatchID(interfacePatch);

    if (patchID_ == -1)
    {
        FatalErrorIn("interfaceModel::interfaceModel(...)")
            << "Interface patch not defined."
            << abort(FatalError);
    }
    calcIs();

    isotherm_ = isothermModel::New(mesh_, interfaceModelCoeffs_);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceModel::~interfaceModel()
{}

// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * * * * * * * * //

const Foam::dictionary& Foam::interfaceModel::modelDict() const
{
    return interfaceModelCoeffs_;
}

Foam::dictionary& Foam::interfaceModel::modelDict()
{
    return interfaceModelCoeffs_;
}

void Foam::interfaceModel::coupledToBulk(const bool& ctb)
{
    coupledToBulk_ = ctb;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //
