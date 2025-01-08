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
    /*scalar ttb = 0;
    scalar tts = 0;
    scalar ttg = 0;
    label titer = 0;
    label steps = 0;*/
    while (runTime.run())
    {
        loopControl looping(runTime, aMesh.solutionDict(),"solution");
        #include "surfaceCourantNo.H"
        #include "setDeltaT.H"
        #include "checkDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;
        p.storePrevIter(); 
        /*scalar tb = 0;
        scalar ts = 0;
        scalar tg = 0;
        label iter = 0;*/

        while (looping.loop())
        {
            residuals.calculateInit(gamma, Us);
            //scalar t1 = mesh.time().elapsedCpuTime();
            if(coupledToBulk)
            {
                #include "UEqn.H"
                #include "pEqn.H"
            }
            //scalar t2 = mesh.time().elapsedCpuTime();
            
            #include "interfStressBalance.H"
            //scalar t3 = mesh.time().elapsedCpuTime();
            #include "gammaEqn.H"
            //scalar t4 = mesh.time().elapsedCpuTime();
            residuals.calculate(gamma, Us);

            /*tb += t2-t1;
            ts += t3-t2;
            tg += t4-t3;
            iter++;*/

            if (residuals.hasConverged())
            {
                break;
            }
        }
        /*reduce(tb, sumOp<scalar>());
        reduce(ts, sumOp<scalar>());
        reduce(tg, sumOp<scalar>());
        reduce(iter, sumOp<label>());
        Info<<"cpuT: " << runTime.timeName() << tab
            << tb << tab << ts << tab << tg << tab << iter << endl;
        ttb += tb;
        tts += ts;
        ttg += tg;
        titer += iter;*/
        
        #include "smoothFields.H"    
        #include "mapToVolFields.H"

        runTime.write();

        runTime.printExecutionTime(Info);
        //steps++;
    }
    /*Info<<"totalTimes:" << tab
        << ttb << tab << tts << tab << ttg << tab 
        << titer << tab << steps << endl;*/

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
