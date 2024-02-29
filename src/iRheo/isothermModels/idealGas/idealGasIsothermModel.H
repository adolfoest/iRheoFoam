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
    Foam::isothermModels::idealGasIsothermModel

Description
    Ideal gas isotherm model. The Marangoni modulus is calculated
    following
    \verbatim
        K_\Pi = k_B T \Gamma
    \verbatim

SourceFiles
    idealGasIsothermModel.C

Author
    Adolfo Esteban, UNED (2023)

\*---------------------------------------------------------------------------*/

#ifndef idealGasIsothermModel_H
#define idealGasIsothermModel_H

#include "isothermModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace isothermModels
{

/*---------------------------------------------------------------------------*\
                      Class idealGasIsothermModel Declaration
\*---------------------------------------------------------------------------*/

class idealGasIsothermModel
:
    public isothermModel
{
        
    // Private member functions

        //- No copy construct
        idealGasIsothermModel(const idealGasIsothermModel&) = delete;

        //- No copy assignment
        void operator=(const idealGasIsothermModel&) = delete;


protected:

    // Protected data
        

public:

    //- Runtime type information
    TypeName("idealGas");


    // Constructors

        //- Construct from surface film model
        idealGasIsothermModel
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~idealGasIsothermModel();


    // Member Functions

        //- 
        virtual tmp<areaScalarField> KPi(const areaScalarField& gamma) const;

        virtual tmp<areaScalarField> Pi
        (
            const areaScalarField& gamma,
            const scalar& Ma = 0.
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace isothermModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
