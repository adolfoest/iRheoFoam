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
    Foam::isothermModel

Description
    Base class for isotherm models

SourceFiles
    isothermModel.C
    isothermModelNew.C
    
Author
    Adolfo Esteban, UNED (2023)

\*---------------------------------------------------------------------------*/

#ifndef isothermModel_H
#define isothermModel_H

#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "autoPtr.H"
#include "IOdictionary.H"
#include "areaFields.H"
#include "faCFD.H"
#include "fvc.H"
#include "fvm.H"
#include "isothermModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class isothermModel Declaration
\*---------------------------------------------------------------------------*/

class isothermModel
: 
    public IOdictionary
{
protected:
    // Protected data

        //- Const reference to mesh
        const fvMesh& mesh_;

    // Protected Member Functions


    //- No copy construct
    isothermModel(const isothermModel &) = delete;

    //- No copy assignment
    void operator=(const isothermModel &) = delete;

public:

    //- Runtime type information
    TypeName("isothermModel");

    // Boltzmann constant, [N m K^-1]
    static const scalar kB;

    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             isothermModel,
             components,
             (
                const fvMesh& mesh,
                const dictionary& dict
             ),
             (mesh, dict)
         );

    //- Constructor
        //- Construct from components
        isothermModel
        (
            const word& type,
            const fvMesh& mesh,
            const dictionary& dict
        );

    // Selectors

        //- Return a reference to the selected interface model
        static autoPtr<isothermModel> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
    virtual ~isothermModel();


    // Member Functions

        //- Return K_\Pi
        virtual tmp<areaScalarField> KPi(const areaScalarField& gamma) const = 0;

        //- Return \Pi
        virtual tmp<areaScalarField> Pi
        (
            const areaScalarField& gamma,
            const scalar& Ma = 0.
        ) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
