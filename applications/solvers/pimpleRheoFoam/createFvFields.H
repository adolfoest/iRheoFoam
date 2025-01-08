#include "createRDeltaT.H"

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

volScalarField gammaVol
(
    IOobject
    (
        "gammaVol",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimless, 0)
);

volVectorField UsVol
(
    IOobject
    (
        "UsVol",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimless, Zero)
);


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


#include "createMRF.H"
#include "createFvOptions.H"

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

autoPtr<interfaceModel> interfaceModel
(
    interfaceModel::New(mesh, "transportProperties")
);

// Create volume-to surface mapping object
volSurfaceMapping vsm(aMesh);

movingFaMesh* mfamPtr(nullptr);
if (mesh.dynamic())
{
    Info<< "\nCreate movingFaMesh object" << endl;
    mfamPtr = new movingFaMesh
    (
        mesh,
        aMesh,
        transportProperties,
        interfaceModel->patchID()
    );
    mfamPtr->correctBC(Us);
}

