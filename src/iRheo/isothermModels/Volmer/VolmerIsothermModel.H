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
    Foam::isothermModels::VolmerIsothermModel

Description
    Volmer isotherm model. The Marangoni modulus is calculated
    following
    \f[
        K_\Pi = \frac{k_B T \Gamma}{(1-\Gamma/\Gamma_\infty)^2}
    \f]
Usage
    Minimal example by using \c constant/transportProperties:
    \verbatim
        isothermModel       Volmer;
        VolmerCoeffs
        {
            gamma0      <value>; // [m^-2]
            gammaInf    <value>; // [m^-2]
            T           <value>; // [K]
        }
    \endverbatim

    where the entries mean:
    \table
      Property           | Description          | Type   | Reqd  | Dflt
      isothermModel      | Type name: Volmer    | word   | yes   | -
      gamma0             | gamma0               | scalar | yes   | -
      gammaInf           | gamma infinity       | scalar | yes   | -
      T                  | Absolute temperature | scalar | yes   | -
    \endtable

SourceFiles
    VolmerIsothermModel.C

Author
    Adolfo Esteban, UNED (2023)

\*---------------------------------------------------------------------------*/

#ifndef VolmerIsothermModel_H
#define VolmerIsothermModel_H

#include "isothermModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace isothermModels
{

/*---------------------------------------------------------------------------*\
                      Class VolmerIsothermModel Declaration
\*---------------------------------------------------------------------------*/

class VolmerIsothermModel
:
    public isothermModel
{
        
    dictionary subDict_;

    scalar gamma0_;

    scalar gammaInf_;

    scalar T_;
    
    // Private member functions

        //- No copy construct
        VolmerIsothermModel(const VolmerIsothermModel&) = delete;

        //- No copy assignment
        void operator=(const VolmerIsothermModel&) = delete;


protected:

    // Protected data
        

public:

    //- Runtime type information
    TypeName("Volmer");


    // Constructors

        //- Construct from surface film model
        VolmerIsothermModel
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~VolmerIsothermModel();


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
