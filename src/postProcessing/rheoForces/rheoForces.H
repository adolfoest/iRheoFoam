/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
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
    Foam::functionObjects::rheoForces

Group
    grpForcesFunctionObjects

Description
    Calculates the forces by integrating the pressure and
    the stress tensor divergence rheoForces over a given list of patches.

    Results are written to multiple files as a function of time in the
    postProcessing/\<functionObjectName\> directory:
    - force.dat          : rheoForces

    Reference:
    \verbatim
        Esteban A, Hernández J, Tajuelo J, Rubio MÁ, 
        A numerical model for the simulation of complex 
        planar Newtonian interfaces, 
        Applied Mathematical Modelling, 2025, 
        doi 10.1016/j.apm.2024.115653.
    \endverbatim

Usage
    Example of function object specification:
    \verbatim
    rheoForces
    {
        type        rheoForces;
        libs        (rheoForces);
        ...
        log         yes;
        writeFields yes;
        patches     (walls);
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: rheoForces       | yes         |
        log          | Write force data to standard output | no | no
        writeFields  | Write the force and moment fields | no | no
        patches      | Patches included in the rheoForces calculation | yes |
        p            | Pressure field name     | no          | p
        U            | Velocity field name     | no          | U
    \endtable

Note
  - Derived from forces.H

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile
    Foam::functionObjects::timeControl

SourceFiles
    rheoForces.C

Author
    Adolfo Esteban, UNED (2023)

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_rheoForces_H
#define functionObjects_rheoForces_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"
#include "HashSet.H"
#include "Tuple2.H"
#include "OFstream.H"
#include "Switch.H"
#include "point.H"
#include "vectorField.H"
#include "faCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class rheoForces Declaration
\*---------------------------------------------------------------------------*/

class rheoForces
:
    public fvMeshFunctionObject,
    public writeFile
{

protected:

    // Protected data

        //- Pressure and viscous force
        List<vector> force_;

        // File streams

            //- Forces
            autoPtr<OFstream> forceFilePtr_;

        // Read from dictionary

            //- Patches to integrate rheoForces over
            labelHashSet patchSet_;

            wordRes areaPatches_;

            //- Name of pressure field
            word pName_;

            //- Name of velocity field
            word UName_;

            //- The name of the surface velocity field
            word UsName_;

            word gammaName_;

            //- Write fields flag
            bool writeFields_;

            //- Initialised flag
            bool initialised_;

            bool calcBulkForce_;

            bool isCIT_;

            scalar x0_;

            scalar L0_;

            scalar Ma_;

            scalar Theta_;

            scalar Tr_;

            tensor Is_;

            scalar Bq_;

            bool isExtensional_;


    // Protected Member Functions

        //- Create the output files
        void createFiles();

        //- Write header for integrated data
        void writeIntegratedHeader(const word& header, Ostream& os) const;

        //- Initialise the fields
        void initialise();

        //- Reset the fields prior to accumulation of forces
        void resetFields();

        //- Return \vect{tau}_s
        tmp<areaTensorField> taus() const;

        //- Return \Pi
        tmp<areaScalarField> Pi() const;

        void getInterfaceParam();

        //- Add patch contributions to force field
        void addToFields
        (
            const label patchi,
            const vectorField& fBulk
        );

        //- Add patch contributions to force field
        void addToFields
        (
            const label patchi,
            const vectorField& fMa,
            const vectorField& fs
        );

        //- Helper function to write integrated rheoForces and moments
        void writeIntegratedForce
        (
            const string& descriptor,
            const vector& fm0,
            const vector& fm1,
            const vector& fm2,
            autoPtr<OFstream>& osPtr
        ) const;

        //- Write force data
        void writeForces();

        //- No copy construct
        rheoForces(const rheoForces&) = delete;

        //- No copy assignment
        void operator=(const rheoForces&) = delete;


public:

    //- Runtime type information
    TypeName("rheoForces");


    // Constructors

        //- Construct from Time and dictionary
        rheoForces
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool readFields = true
        );

        //- Construct from objectRegistry and dictionary
        rheoForces
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const bool readFields = true
        );


    //- Destructor
    virtual ~rheoForces() = default;


    // Member Functions

        //- Read the rheoForces data
        virtual bool read(const dictionary&);

        //- Calculate the rheoForces and moments
        virtual void calcForces();

        //- Return the total force
        virtual vector forceEff() const;

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the rheoForces
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
