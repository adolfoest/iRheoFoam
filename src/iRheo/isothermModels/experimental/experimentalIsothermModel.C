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

#include "experimentalIsothermModel.H"
#include "addToRunTimeSelectionTable.H"
#include "interpolateSplineXY.H"
#include "tableReader.H"
#include "zeroGradientFaPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace isothermModels
{
    defineTypeNameAndDebug(experimentalIsothermModel, 0);

    addToRunTimeSelectionTable
    (
        isothermModel,
        experimentalIsothermModel,
        components
    );
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isothermModels::experimentalIsothermModel::experimentalIsothermModel
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
    fileName_(subDict_.get<fileName>("file")),
    Gamma_(),
    Pi_()
{
    readTable();

    h_ = (Gamma_[Gamma_.size() - 1] - Gamma_[0])/(Gamma_.size()*100);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::isothermModels::experimentalIsothermModel::~experimentalIsothermModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::isothermModels::experimentalIsothermModel::readTable()
{
    autoPtr<tableReader<scalar>> pReader = 
        tableReader<scalar>::New(subDict_);
    
    List<Tuple2<scalar, scalar>> data;
    fileName_.expand();

    auto& reader = pReader.ref();
    reader(fileName_, data);

    if (data.empty())
    {
        FatalErrorInFunction
            << "table read from " << fileName_ << " is empty" << nl
            << exit(FatalError);
    }

    // Resize to experimental data 
    Gamma_.resize(data.size());
    Pi_.resize(data.size());

    forAll(data, i)
    {
        const auto& dataTuple = data[i]; 
        Gamma_[i] = dataTuple.first();
        Pi_[i] = dataTuple.second();
    } 
}

void Foam::isothermModels::experimentalIsothermModel::checkGamma
(
    const areaScalarField& gamma
)const
{
    scalar gammaMin = gMin(gamma.primitiveField());
    scalar gammaMax = gMax(gamma.primitiveField());

    if (gammaMin < Gamma_[0] || gammaMax > Gamma_[Gamma_.size() - 1])
    {
        FatalErrorInFunction
            << "gamma field out of bounds "
            << "of the experimental isotherm" 
            << nl
            << exit(FatalError);
    }
}

Foam::scalar Foam::isothermModels::experimentalIsothermModel::partialPiGammaForw
(
    const scalar& gammaI
)const
{
    // Second-order forward differentiation
    scalar fx = interpolateSplineXY
    (
        gammaI, 
        Gamma_, 
        Pi_
    );

    scalar fx1 = interpolateSplineXY
    (
        gammaI+h_, 
        Gamma_, 
        Pi_
    );

    scalar fx2 = interpolateSplineXY
    (
        gammaI+2*h_, 
        Gamma_, 
        Pi_
    );

    return ((-fx2 + 4*fx1 - 3*fx)/(2*h_));
}

Foam::scalar Foam::isothermModels::experimentalIsothermModel::partialPiGammaBack
(
    const scalar& gammaI
)const
{
    // Second-order backward differentiation
    scalar fx = interpolateSplineXY
    (
        gammaI, 
        Gamma_, 
        Pi_
    );

    scalar fx1 = interpolateSplineXY
    (
        gammaI-h_, 
        Gamma_, 
        Pi_
    );

    scalar fx2 = interpolateSplineXY
    (
        gammaI-2*h_, 
        Gamma_, 
        Pi_
    );

    return ((fx2 - 4*fx1 + 3*fx)/(2*h_));
}

Foam::tmp<Foam::areaScalarField>
Foam::isothermModels::experimentalIsothermModel::KPi
(
    const areaScalarField& gamma
) const
{
    checkGamma(gamma);
    auto tKPi = tmp<areaScalarField>::New
    (
        IOobject
        (
            "KPi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        gamma.mesh(),
        dimensionedScalar("KPi", dimless, 0.0), 
        zeroGradientFaPatchScalarField::typeName
    );

    auto& KPi = tKPi.ref();

    forAll(KPi, i)
    {
        const scalar& gammaI = gamma[i]; 
        if (i == 0)
        {
            KPi[i] = partialPiGammaForw(gammaI);
        }
        else
        {
            KPi[i] = partialPiGammaBack(gammaI);
        }
        KPi[i] *= gammaI;
    }
    KPi.correctBoundaryConditions();

    return tKPi;
}

Foam::tmp<Foam::areaScalarField>
Foam::isothermModels::experimentalIsothermModel::Pi
(
    const areaScalarField& gamma,
    const scalar& Ma
) const
{
    checkGamma(gamma);
    auto tPi = tmp<areaScalarField>::New
    (
        IOobject
        (
            "Pi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        gamma.mesh(),
        dimensionedScalar("Pi", dimless, 0.0), 
        zeroGradientFaPatchScalarField::typeName
    );

    auto& Pi = tPi.ref();

    forAll(Pi, i)
    {
        const scalar& gammaI = gamma[i]; 
        scalar piInterp = interpolateSplineXY
        (
            gammaI, 
            Gamma_, 
            Pi_
        );
        Pi[i] = piInterp;
    }
    Pi.correctBoundaryConditions();

    return tPi;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
