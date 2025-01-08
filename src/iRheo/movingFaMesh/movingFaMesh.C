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

#include "movingFaMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
void Foam::movingFaMesh::setBfLabels()
{
    labelUList mpf = aMesh_.boundary().edgeFaces()[movingAPatch_];
    const labelList& faceOwner = mesh_.faceOwner();
    bfLabels_.setSize(mpf.size(), -1);

    forAll(mpf, i)
    {
        const label& faceI = mesh_.boundaryMesh()[patchID_].start() + mpf[i];
        const label& cellI = faceOwner[faceI];

        forAll(mesh_.cells()[cellI], j)
        {
            const label& faceJ = mesh_.cells()[cellI][j];
            if (!mesh_.isInternalFace(faceJ))
            {
                const label& pf = mesh_.boundaryMesh().whichPatch(faceJ);
                if (pf == movingPatch_)
                {
                    const label& k = mesh_.boundaryMesh()[pf].whichFace(faceJ);
                    bfLabels_[i] = k;
                }
            }
        }
    }
}
Foam::direction Foam::movingFaMesh::cmpt()
{
    const dynamicMotionSolverFvMesh& dynMsh =
        refCast<const dynamicMotionSolverFvMesh>(mesh_);
    motionSolver& motion =
        const_cast<motionSolver&>(dynMsh.motion());
    
    const dictionary& dict = motion.coeffDict();
    word cmptName(dict.lookup("component"));

    if (cmptName == "x")
    {
        return vector::X;
    }
    else if (cmptName == "y")
    {
        return vector::Y;
    }
    else if (cmptName == "z")
    {
        return vector::Z;
    }
    else
    {
        FatalErrorInFunction
            << "Given component name " << cmptName << " should be x, y or z"
            << exit(FatalError);

        return 0;
    }
    
}

void Foam::movingFaMesh::update()
{
    const dynamicMotionSolverFvMesh& dynMsh =
        refCast<const dynamicMotionSolverFvMesh>(mesh_);
    motionSolver& motion =
        const_cast<motionSolver&>(dynMsh.motion());

    velocityComponentLaplacianFvMotionSolver& vlMotion =
        refCast<velocityComponentLaplacianFvMotionSolver>
        (
            motion
        );
    const pointScalarField& pointMotionU = 
        vlMotion.pointMotionU();

    const polyPatch& pp = mesh_.boundaryMesh()[patchID_];

    vectorField up(pp.localPoints().size(), Zero);
    forAll(pp.meshPoints(), i)
    {
        const label& pi = pp.meshPoints()[i];
        up[i].replace(cmpt_, pointMotionU[pi]);
    }

    const edgeList& aEdges = aMesh_.edges();
    const edgeVectorField& le = aMesh_.Le();
    forAll(aEdges, i)
    {
        const label& p0 = aEdges[i][0];
        const label& p1 = aEdges[i][1];

        const vector& u0 = up[p0];
        const vector& u1 = up[p1];
        vector uav = 0.5*(u0 + u1);
        if (aMesh_.isInternalEdge(i))
        {
            meshPhis()[i] = uav & le[i];
        }
        else
        {
            const label& patchi = aMesh_.boundary().whichPatch(i);
            const label& patchEdgei =
                aMesh_.boundary()[patchi].whichEdge(i);
            meshPhis().boundaryFieldRef()[patchi][patchEdgei] =
                uav & le.boundaryField()[patchi][patchEdgei];
            if (patchi == movingAPatch_)
            {
                const label& patchEdgei =
                    aMesh_.boundary()[patchi].whichEdge(i);
                ue_[patchEdgei] = uav;
            } 
        }
    }
    
}

void Foam::movingFaMesh::init()
{
    movingPatch_  = mesh_.boundaryMesh().findPatchID(movingPatchName_);
    movingAPatch_ = aMesh_.boundary().findPatchID(movingPatchName_);
    cmpt_ = cmpt();
    setBfLabels();

    meshPhisPtr_ = new edgeScalarField
    ( 
        IOobject
        (
            "meshPhis",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh_,
        dimensionedScalar(dimArea/dimTime, 0.0)
    );
    ue_.setSize(aMesh_.boundary()[movingAPatch_].size());
    ue_ = Zero;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingFaMesh::movingFaMesh
(
    const dynamicFvMesh& mesh,
    const faMesh& aMesh,
    const IOdictionary& dict,
    const label& patchID
)
:
    mesh_(mesh),
    aMesh_(aMesh),
    meshPhisPtr_(nullptr),
    movingPatchName_(dict.lookupOrDefault<word>("movingPatch", "movingWire")),
    movingPatch_(-1),
    movingAPatch_(-1),
    patchID_(patchID),
    bfLabels_(),
    ubf_(),
    nbf_(),
    timeIndex_(-1),
    ue_(),
    cmpt_()
{
    init();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::movingFaMesh::makeRelative
(
    edgeScalarField& phis
)
{
    if (timeIndex_ != mesh_.time().timeIndex())
    {
        update();
        timeIndex_ = mesh_.time().timeIndex();
    }
    phis -= meshPhis();
}

void Foam::movingFaMesh::correct
(
    edgeVectorField& Usf,
    const edgeScalarField& phis
)
{
    edgeVectorField ns(aMesh_.Le()/aMesh_.magLe());
    Usf += ns*(phis/aMesh_.magLe() - (ns & Usf));    
}

void Foam::movingFaMesh::correctBC
(
    areaVectorField& Us,
    const volVectorField& U
)
{
    if (timeIndex_ != U.mesh().time().timeIndex())
    {
        //setbfLabels(); // no topology changes
        labelUList mpf = aMesh_.boundary().edgeFaces()[movingAPatch_];
        ubf_.setSize(mpf.size(), Zero);
        nbf_.setSize(mpf.size(), Zero);
        forAll(bfLabels_, i)
        {
            const label& j = bfLabels_[i];
            ubf_[i] = U.boundaryField()[movingPatch_][j];
            nbf_[i] = aMesh_.faceAreaNormals()[mpf[i]];
        }
        ubf_ -= nbf_*(nbf_&ubf_);
        timeIndex_ = U.mesh().time().timeIndex();
    }
    Us.boundaryFieldRef()[movingAPatch_] == ubf_;
}

void Foam::movingFaMesh::correctBC
(
    areaVectorField& Us
)
{
    if (timeIndex_ != mesh_.time().timeIndex())
    {
        update();
        timeIndex_ = mesh_.time().timeIndex();   
    }
    Us.boundaryFieldRef()[movingAPatch_] == ue_;
}

Foam::edgeScalarField& Foam::movingFaMesh::meshPhis()
{
    return *meshPhisPtr_;
}

// ************************************************************************* //
