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
    Foam::interfaceModel

Description
    Base class for interface models

    Reference:
    \verbatim
        Esteban A, Hernández J, Tajuelo J, Rubio MÁ, 
        A numerical model for the simulation of complex 
        planar Newtonian interfaces, 
        Applied Mathematical Modelling, 2025, 
        doi 10.1016/j.apm.2024.115653.
    \endverbatim

SourceFiles
    interfaceModel.C
    interfaceModelNew.C
    
Author
    Adolfo Esteban, UNED (2023)

\*---------------------------------------------------------------------------*/

#ifndef interfaceModel_H
#define interfaceModel_H

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
                      Class interfaceModel Declaration
\*---------------------------------------------------------------------------*/

class interfaceModel
:
    public IOdictionary
{

private:
    //-
    dictionary  interfaceModelCoeffs_;

protected:
    // Protected data

        //- Const reference to mesh
        const fvMesh& mesh_;

        //- Marangoni number
        dimensionedScalar Ma_;

        //- Boussinesq_1
        dimensionedScalar Bq_;

        //- Boussinesq_2
        dimensionedScalar Bqk_;

        //- \Theta = Bq_2/Bq_1
        dimensionedScalar Theta_;

        //- Trouton
        dimensionedScalar Tr_;

        //- Peclet
        dimensionedScalar Pes_;

        //- Reynolds
        dimensionedScalar Re_;

        //- Interface unit normal vector (uniform, planar interface)
        vector n_;

        //- Interface projection tensor (uniform, planar interface)
        dimensioned<tensor> Is_;
    
        //- Interface patch label  
        label patchID_;

        //- Interface is coupled to bulk
        bool coupledToBulk_;

        autoPtr<isothermModel> isotherm_;

    // Protected Member Functions

        //- Calculates Is
        void calcIs();

        //- Returns \vec t - (\mathrm{diag}(\vec t) \otimes \vec e )\circ \vec I_s
        areaTensorField Hadamard(const areaTensorField& t) const;

    //- No copy construct
    interfaceModel(const interfaceModel &) = delete;

    //- No copy assignment
    void operator=(const interfaceModel &) = delete;

public:
    //- Access to model dictionary 
    dictionary& modelDict();

    //- Const access to model dictionary
    const dictionary& modelDict() const;

public:

    //- Runtime type information
    TypeName("interfaceModel");

    //- Default name of the interfaceModel dictionary
    static const word dictName;

    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             interfaceModel,
             components,
             (
                const fvMesh& mesh,
                const word& dictName
             ),
             (mesh, dictName)
         );

    //- Constructor
        //- Construct from components
        interfaceModel
        (
            const word& type,
            const fvMesh& mesh,
            const word& dictName
        );

    // Selectors

        //- Return a reference to the selected interface model
        static autoPtr<interfaceModel> New
        (
            const fvMesh& mesh,
            const word& dictName
        );

    //- Destructor
    virtual ~interfaceModel();


    // Member Functions

        //- Return \div \cdot \vec{tau}_s
        virtual tmp<faVectorMatrix> divTaus(areaVectorField& Us) const = 0;

        //- Return Ma term 
        virtual tmp<areaVectorField> marangoni(const areaScalarField& gamma) const = 0;

        //- Return \vec{n} \cdot \vec{tau} \cdot \vec{t}
        virtual tmp<areaVectorField> nTaut
        (
            volVectorField& U,
            areaVectorField& Us
        ) const = 0;

        //- Return \nabla(\vec{n} \cdot \vec{tau} \cdot \vec{n}) \cdot \vec{n}
        virtual scalarField snGradTau
        (
            const volVectorField& U
        ) const = 0;

        //- Return \vec{tau}_s
        virtual tmp<areaTensorField> taus(areaVectorField& Us) const = 0;

        //- Return \vec{D}_s
        virtual tmp<areaTensorField> Ds(areaVectorField& Us) const = 0;

        //- Acces function to interface patch label
        const label& patchID() const
        {
            return patchID_;
        }

        const dimensionedScalar& Ma() const
        {
            return Ma_;
        }

        const dimensionedScalar& Bq() const
        {
            return Bq_;
        }

        const dimensionedScalar& Bqk() const
        {
            return Bqk_;
        }

        const dimensionedScalar& Theta() const
        {
            return Theta_;
        }

        const dimensionedScalar& Tr() const
        {
            return Tr_;
        }

        const dimensionedScalar& Pes() const
        {
            return Pes_;
        }

        const dimensionedScalar& Re() const
        {
            return Re_;
        }

        const dimensionedScalar invRe() const
        {
            return dimensionedScalar
            (
                "1/Re",
                dimViscosity, 
                1.0/Re_.value()
            );
        }

        const dimensionedTensor& Is() const
        {
            return Is_;
        }

        void coupledToBulk(const bool& ctb);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
