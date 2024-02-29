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

#include "extensionalViscosity.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace interfaceModels
{
    defineTypeNameAndDebug(extensionalViscosity, 0);

    addToRunTimeSelectionTable
    (
        interfaceModel,
        extensionalViscosity,
        components
    );
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceModels::extensionalViscosity::extensionalViscosity
(
    const fvMesh& mesh,
    const word& dictName
)
:
interfaceModel
(
    typeName,
    mesh,
    dictName
)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceModels::extensionalViscosity::~extensionalViscosity()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::tmp<Foam::faVectorMatrix> 
Foam::interfaceModels::extensionalViscosity::divTaus
(
    areaVectorField& Us
) const
{
    tmp<areaTensorField> tDs = Ds(Us);
    areaTensorField& Ds = tDs.ref();

    areaTensorField Dshear = Hadamard(Ds);
    areaTensorField gradUsMod = Hadamard(fac::grad(Us));

    return
    (
      - fac::div((Theta_ - 0.25*Tr_)*tr(Ds)*Is_)
      - fac::div(dimensionedScalar(dimViscosity, 2.0)*Dshear)
      - fam::laplacian(2.*0.25*Tr_, Us)
      + fac::div(2.*0.25*Tr_*gradUsMod)
    );
}

Foam::tmp<Foam::areaVectorField>
Foam::interfaceModels::extensionalViscosity::marangoni
(
    const areaScalarField& gamma
) const
{
    tmp<areaScalarField> tKPi = isotherm_->KPi(gamma);
    areaScalarField& KPi = tKPi.ref();
    return 
    (
        Ma_*KPi*fac::grad(log(gamma))
    );
}

Foam::tmp<Foam::areaVectorField>
Foam::interfaceModels::extensionalViscosity::nTaut
(
    volVectorField& U,
    areaVectorField& Us
) const
{
    tensorField guvb(fvc::grad(U)().boundaryField()[patchID_]);

    // Get the projection of grad(U) tangent to the interface
    guvb -= n_*(n_ & guvb);
    areaTensorField tau(fac::grad(Us)*dimensionedScalar(dimVelocity, 0.0));
    if (coupledToBulk_)
    {
        forAll(tau, i)
        {
            tau[i] = guvb[i] + T(guvb[i]);
        }
        tau.correctBoundaryConditions();
    }
    
    return 
    (
        (n_ & -tau)/Bq_.value()
    );
}

Foam::scalarField
Foam::interfaceModels::extensionalViscosity::snGradTau
(
    const volVectorField& U
) const
{
    volSymmTensorField tau(twoSymm(fvc::grad(U)));
    volVectorField gTau(fvc::grad(((n_ & tau) & n_)));
    vectorField gTauB(gTau.boundaryField()[patchID_]);
    return 
    (
        (n_ & gTauB)
    );
}

Foam::tmp<Foam::areaTensorField>
Foam::interfaceModels::extensionalViscosity::taus
(
    areaVectorField& Us
) const
{
    tmp<areaTensorField> tDs = Ds(Us);
    areaTensorField& Ds = tDs.ref();

    areaTensorField Dshear = Hadamard(Ds);
    areaTensorField gradUs = fac::grad(Us);
    areaTensorField Dext = gradUs - Hadamard(gradUs);

    return
        tmp<areaTensorField>::New
        (
            IOobject
            (
                "taus",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            Bq_.value()*
            (
                (Theta_ - 0.25*Tr_)*tr(Ds)*Is_ 
              + dimensionedScalar(dimViscosity, 2.0)*Dshear 
              + 0.5*Tr_*Dext
            )
        );
}

Foam::tmp<Foam::areaTensorField>
Foam::interfaceModels::extensionalViscosity::Ds
(
    areaVectorField& Us
) const
{
    return
        tmp<areaTensorField>::New
        (
            IOobject
            (
                "Ds",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            0.5*((fac::grad(Us) & Is_) + (Is_ & T(fac::grad(Us))))
        );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
