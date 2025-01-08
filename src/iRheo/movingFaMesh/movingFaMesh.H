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
    Foam::movingFaMesh

Description
    Class to obtain the finite area mesh movement based on the finite
    volume mesh movement for its corresponding moving patches.
    For the moment, works with velocityComponentLaplacian only.

    Reference:
    \verbatim
        Esteban A, Hernández J, Tajuelo J, Rubio MÁ, 
        A numerical model for the simulation of complex 
        planar Newtonian interfaces, 
        Applied Mathematical Modelling, 2025, 
        doi 10.1016/j.apm.2024.115653.
    \endverbatim

SourceFiles
    movingFaMesh.C

Author
    Adolfo Esteban, UNED (2023)

\*---------------------------------------------------------------------------*/

#ifndef movingFaMesh_H
#define movingFaMesh_H

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "faCFD.H"
#include "regIOobject.H"

#include "dynamicMotionSolverFvMesh.H"
#include "motionSolver.H"
#include "velocityComponentLaplacianFvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class movingFaMesh Declaration
\*---------------------------------------------------------------------------*/

class movingFaMesh
{
    // Private Data

        const dynamicFvMesh& mesh_;

        const faMesh& aMesh_;

        // meshPhi for aMesh
        mutable edgeScalarField* meshPhisPtr_;

        // Name of moving patch in mesh
        word movingPatchName_;

        // Label of moving patch in mesh
        label movingPatch_;

        // Corresponding label in aMesh assuming 
        // same name as movingPatch
        label movingAPatch_;

        // Patch label of mesh in which 
        // aMesh is constructed
        const label& patchID_;

        // List of movingPatch faces owning the 
        // corresponding aMesh edge
        labelList bfLabels_;

        // Velocity of the moving patch
        vectorField ubf_;

        // Normals of the moving patch
        vectorField nbf_;

        label timeIndex_;

        vectorField ue_;

        direction cmpt_;
    
    // Private Member Functions

        //- No copy construct
        movingFaMesh(const movingFaMesh&) = delete;

        //- No copy assignment
        void operator=(const movingFaMesh&) = delete;

        void setBfLabels();

        // Get the velocity of the moving patch points
        void update();

        direction cmpt();

        void init();


public:

    //- Constructor
    movingFaMesh
    (
        const dynamicFvMesh& mesh,
        const faMesh& aMesh,
        const IOdictionary& dict,
        const label& patchID
    );

    //- Destructor
    virtual ~movingFaMesh() = default;


    // Member Functions

        void makeRelative
        (
            edgeScalarField& phis
        );

        void correct
        (
            edgeVectorField& Usf,
            const edgeScalarField& phis
        );

        // Set the velocity of the mesh moving patch to
        // the corresponding aMesh moving patch
        void correctBC
        (
            areaVectorField& Us,
            const volVectorField& U
        );

        void correctBC
        (
            areaVectorField& Us
        );

        const edgeScalarField& meshPhis() const
        {
            return *meshPhisPtr_;
        }

        edgeScalarField& meshPhis();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
