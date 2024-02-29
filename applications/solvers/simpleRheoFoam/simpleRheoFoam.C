/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Application
    simpleRheoFoam

Description
    Steady-state solver for incompressible flow of a Newtonian fluid
    with a compressible planar interface on a static mesh.

Notes
    - Derived from simpleFoam.C

Author
    Adolfo Esteban, UNED (2023)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"
#include "fvOptions.H"

#include "faCFD.H" 
#include "loopControl.H" 
#include "interfaceModel.H" 
#include "rheoResiduals.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Steady-state solver for incompressible flow of a Newtonian fluid"
        " with a compressible planar interface on a static mesh."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H" 
    #include "createFaMesh.H" 
    #include "createControl.H"
    #include "createFaFields.H" 
    #include "createFvFields.H"
    #include "createTimeControls.H" 
    #include "initContinuityErrs.H"
    #include "readSolutionControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        loopControl looping(runTime, aMesh.solutionDict(),"solution");
        #include "surfaceCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;
        p.storePrevIter(); 

        while (looping.loop())
        {
            residuals.calculateInit(gamma, Us);
            if(coupledToBulk)
            {
                #include "UEqn.H"
                #include "pEqn.H"
            }
             
            #include "interfStressBalance.H"
            #include "gammaEqn.H"
            residuals.calculate(gamma, Us);

            if (residuals.hasConverged())
            {
                break;
            }
        }
        
        #include "smoothFields.H"    
        #include "mapToVolFields.H"

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
