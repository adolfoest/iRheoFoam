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
    Foam::isothermModels::experimentalIsothermModel

Description
    Experimental isotherm model. The Marangoni modulus is calculated
    following
    \f[
        K_\Pi = \Gamma \frac{\partial \Pi}{\Partial \Gamma}
    \f]
    with the \f$ \Gamma \f$ and \f$ \Pi \f$ values obtained experimentally
    stored in a table.
Usage
    Reads the filename from the dictionary and has an optional readerType. 
    Default is to read OpenFOAM format. The only other format is csv 
    (comma separated values). Minimal example by using 
    \c constant/transportProperties:
    \verbatim
        isothermModel       experimental;
        experimentalCoeffs
        {
            readerType          csv;
            file                "<constant>/isotherm.csv";
            hasHeaderLine       true;   // skip first line
            refColumn           0;      // \Gamma is in column 0, 
                                        // must be dimensionless with gamma0
            componentColumns    (1);    // \Pi values starts in column 1, 
                                        // must be dimensionless with \eta U
            separator           ",";
        }
    \endverbatim

SourceFiles
    experimentalIsothermModel.C

Author
    Adolfo Esteban, UNED (2023)

\*---------------------------------------------------------------------------*/

#ifndef experimentalIsothermModel_H
#define experimentalIsothermModel_H

#include "isothermModel.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace isothermModels
{

/*---------------------------------------------------------------------------*\
                      Class experimentalIsothermModel Declaration
\*---------------------------------------------------------------------------*/

class experimentalIsothermModel
:
    public isothermModel
{
        dictionary subDict_;

        // Name of the file to be read containg the table
        fileName fileName_;

        // Values of \Gamma read from the file
        scalarField Gamma_;

        // Values of Pi_ read from the file
        scalarField Pi_;

        // Differentiation interval
        scalar h_;
    
    // Private member functions

        //- No copy construct
        experimentalIsothermModel(const experimentalIsothermModel&) = delete;

        //- No copy assignment
        void operator=(const experimentalIsothermModel&) = delete;


    // Read the \Gamma-\Pi table
        void readTable();

        void checkGamma(const areaScalarField& gamma) const;

        scalar partialPiGammaForw(const scalar& gammaI) const;

        scalar partialPiGammaBack(const scalar& gammaI) const;
        

public:

    //- Runtime type information
    TypeName("experimental");


    // Constructors

        //- Construct from surface film model
        experimentalIsothermModel
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~experimentalIsothermModel();


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
