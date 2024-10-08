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

Class
    Foam::interfaceModels::BoussinesqScriven

Description
    Boussinesq-Scriven interfacial model

    Reference:
    \verbatim
        Esteban A, Hernández J, Tajuelo J, Rubio MÁ, 
        A numerical model for the simulation of complex 
        planar Newtonian interfaces, 
        Applied Mathematical Modelling, 2025, 
        doi 10.1016/j.apm.2024.115653.
    \endverbatim

SourceFiles
    BoussinesqScriven.C

Author
    Adolfo Esteban, UNED (2023)

\*---------------------------------------------------------------------------*/

#ifndef BoussinesqScriven_H
#define BoussinesqScriven_H

#include "interfaceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace interfaceModels
{

/*---------------------------------------------------------------------------*\
                      Class BoussinesqScriven Declaration
\*---------------------------------------------------------------------------*/

class BoussinesqScriven
:
    public interfaceModel
{
        
    // Private member functions

        //- No copy construct
        BoussinesqScriven(const BoussinesqScriven&) = delete;

        //- No copy assignment
        void operator=(const BoussinesqScriven&) = delete;


protected:

    // Protected data
        

public:

    //- Runtime type information
    TypeName("BoussinesqScriven");


    // Constructors

        //- Construct from surface film model
        BoussinesqScriven
        (
            const fvMesh& mesh,
            const word& dictName = interfaceModel::dictName
        );


    //- Destructor
    virtual ~BoussinesqScriven();


    // Member Functions

        //- 
        virtual tmp<faVectorMatrix> divTaus(areaVectorField& Us) const;

        virtual tmp<areaVectorField> marangoni(const areaScalarField& gamma) const;

        virtual tmp<areaVectorField> nTaut
        (
            volVectorField& U,
            areaVectorField& Us
        ) const;

        virtual scalarField snGradTau
        (
            const volVectorField& U
        ) const;

        virtual tmp<areaTensorField> taus(areaVectorField& Us) const;

        virtual tmp<areaTensorField> Ds(areaVectorField& Us) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interfaceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
