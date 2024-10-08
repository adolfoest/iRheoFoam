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
    Foam::rheoResiduals

Description
    Class to obtain the residuals when solving the surfactant transport
    and the interfacial momentum balance equations.

    Reference:
    \verbatim
        Esteban A, Hernández J, Tajuelo J, Rubio MÁ, 
        A numerical model for the simulation of complex 
        planar Newtonian interfaces, 
        Applied Mathematical Modelling, 2025, 
        doi 10.1016/j.apm.2024.115653.
    \endverbatim

SourceFiles
    rheoResiduals.C

Author
    Adolfo Esteban, UNED (2023)

\*---------------------------------------------------------------------------*/

#ifndef rheoResiduals_H
#define rheoResiduals_H

#include "fvCFD.H"
#include "faCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class rheoResiduals Declaration
\*---------------------------------------------------------------------------*/

class rheoResiduals
{
    // Private Data
        scalar tol_;
        scalar gammaHat0_;
        scalar UsHat0_;
        scalar gammaHat_;
        scalar UsHat_;
        bool hasConverged_;
        scalar gRes_;
        scalar uRes_;
    
    // Private Member Functions

        //- No copy construct
        rheoResiduals(const rheoResiduals&) = delete;

        //- No copy assignment
        void operator=(const rheoResiduals&) = delete;

public:

    //- Constructor
    rheoResiduals
    (
        const scalar& tol
    );

    //- Destructor
    virtual ~rheoResiduals() = default;


    // Member Functions

        void calculateInit
        (
            const areaScalarField& gamma,
            const areaVectorField& Us
        );

        void calculate
        (
            const areaScalarField& gamma,
            const areaVectorField& Us
        );

        const bool& hasConverged() const
        {
            return hasConverged_;
        }

        const scalar& gamma() const
        {
            return gRes_;
        }
        const scalar& Us() const
        {
            return uRes_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
