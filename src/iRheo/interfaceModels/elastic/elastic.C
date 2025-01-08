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

#include "elastic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace interfaceModels
{
    defineTypeNameAndDebug(elastic, 0);

    addToRunTimeSelectionTable
    (
        interfaceModel,
        elastic,
        components
    );
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceModels::elastic::elastic
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

Foam::interfaceModels::elastic::~elastic()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::tmp<Foam::faVectorMatrix> 
Foam::interfaceModels::elastic::divTaus
(
    areaVectorField& Us
) const
{
    tmp<areaTensorField> tDs = Ds(Us);
    areaTensorField& Ds = tDs.ref();
    return
    (
      - fam::laplacian(Bq_, Us)
      - fac::div(Bqk_*tr(Ds)*Is_)
    );
}

Foam::tmp<Foam::areaVectorField>
Foam::interfaceModels::elastic::marangoni
(
    const areaScalarField& gamma
) const
{
    tmp<areaScalarField> tKPi = KPi(gamma);
    areaScalarField& KPi = tKPi.ref();
    return 
    (
        Ma_*KPi*fac::grad(log(gamma))
    );
}

Foam::tmp<Foam::areaVectorField>
Foam::interfaceModels::elastic::nTaut
(
    volVectorField& U,
    areaVectorField& Us
) const
{
    tensorField guvb(fvc::grad(U)().boundaryField()[patchID_]);
    
    // Get the projection of grad(U) tangent to the interface
    guvb -= (guvb & n_)*n_;
    areaTensorField tau(fac::grad(Us)*dimensionedScalar(dimVelocity, 0.0));
    if (coupledToBulk_)
    {
        forAll(tau, i)
        {
            tau[i] = guvb[i] + T(guvb[i]);
        }
    }
    return 
    (
        (n_ & -tau)
    );
}

Foam::scalarField
Foam::interfaceModels::elastic::snGradTau
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
Foam::interfaceModels::elastic::taus
(
    areaVectorField& Us
) const
{
    tmp<areaTensorField> tDs = Ds(Us);
    areaTensorField& Ds = tDs.ref();
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
            ((Bqk_ - Bq_)*(Is_ && Ds)*Is_ + 2.*Bq_*Ds)
        );
}

Foam::tmp<Foam::areaTensorField>
Foam::interfaceModels::elastic::Ds
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
