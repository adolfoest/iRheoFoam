Info<< "Reading field Us" << endl;
areaVectorField Us
(
    IOobject
    (
        "Us",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    aMesh
);

edgeScalarField phis
(
    IOobject
    (
        "phis",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fac::interpolate(Us) & aMesh.Le()
);

Info<< "Reading field gamma" << endl;
areaScalarField gamma
(
    IOobject
    (
        "gamma",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    aMesh
);
gamma.storePrevIter();

areaScalarField divUs
(
    IOobject
    (
        "divUs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fac::div(phis)
);

edgeVectorField Usf
(
    IOobject
    (
        "Usf",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fac::interpolate(Us)
);

