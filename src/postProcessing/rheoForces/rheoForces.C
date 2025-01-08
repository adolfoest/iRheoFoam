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

\*---------------------------------------------------------------------------*/

#include "rheoForces.H"
#include "fvcGrad.H"
#include "interfaceModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(rheoForces, 0);
    addToRunTimeSelectionTable(functionObject, rheoForces, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::rheoForces::createFiles()
{
    if (writeToFile() && !forceFilePtr_)
    {
        forceFilePtr_ = createFile("force");
        writeIntegratedHeader("Force", forceFilePtr_());
    }

}


void Foam::functionObjects::rheoForces::writeIntegratedHeader
(
    const word& header,
    Ostream& os
) const
{
    writeHeader(os, header);
    writeHeader(os, "");
    writeCommented(os, "Time");
    writeTabbed(os, "(total_x total_y total_z)");
    writeTabbed(os, "(ma_x ma_y ma_z)");
    writeTabbed(os, "(s_x s_y s_z)");
    writeTabbed(os, "(Bulk_x Bulk_y Bulk_z)");

    os  << endl;

}


void Foam::functionObjects::rheoForces::initialise()
{
    if (initialised_)
    {
        return;
    }

    if
    (
        !foundObject<areaVectorField>(UsName_)
      || !foundObject<areaScalarField>(gammaName_)

    )
    {
        FatalErrorInFunction
            << "Could not find Us: " << UsName_ << " or gamma: " << gammaName_
            << " in database"
            << exit(FatalError);
    }

    if (calcBulkForce_)
    {     
        if
        (
            !foundObject<volVectorField>(UName_)
         || !foundObject<volScalarField>(pName_)

        )
        {
            FatalErrorInFunction
                << "Could not find U: " << UName_ << " or p:" << pName_
                << " in database"
                << exit(FatalError);
        }

    }
    
    getInterfaceParam();
    initialised_ = true;
}


void Foam::functionObjects::rheoForces::resetFields()
{
    force_[0] = Zero;
    force_[1] = Zero;
    force_[2] = Zero;

    if (writeFields_)
    {
        volVectorField& force =
            lookupObjectRef<volVectorField>(scopedName("force"));

        force == dimensionedVector(force.dimensions(), Zero);
    }

}


Foam::tmp<Foam::areaTensorField>
Foam::functionObjects::rheoForces::taus() const
{
    areaVectorField Us = lookupObject<areaVectorField>(UsName_);

    if (foundObject<interfaceModel>("interfaceModel"))
    {
        const interfaceModel& im =
            lookupObject<interfaceModel>("interfaceModel");
        return im.taus(Us);
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for tau_s calculation"
            << exit(FatalError);     
        return areaTensorField::null();
    } 
}

Foam::tmp<Foam::areaScalarField>
Foam::functionObjects::rheoForces::Pi() const
{
    areaScalarField gamma = lookupObject<areaScalarField>(gammaName_);

    if (foundObject<interfaceModel>("interfaceModel"))
    {
        const interfaceModel& im =
            lookupObject<interfaceModel>("interfaceModel");
        return im.Pi(gamma);
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for Pi calculation"
            << exit(FatalError);     
        return areaScalarField::null();
    } 
}

void Foam::functionObjects::rheoForces::getInterfaceParam()
{
    if (foundObject<interfaceModel>("interfaceModel"))
    {
        const interfaceModel& im =
            lookupObject<interfaceModel>("interfaceModel");
        Ma_ = im.Ma().value();
        Theta_ = im.Theta().value();
        Tr_ = im.Tr().value();
        Is_ = im.Is().value(); 
        Bq_ = im.Bq().value();
    }
}

void Foam::functionObjects::rheoForces::addToFields
(
    const label patchi,
    const vectorField& fBulk
)
{
    force_[2] += sum(fBulk);
    if (!writeFields_)
    {
        return;
    }

    auto& force = lookupObjectRef<volVectorField>(scopedName("force"));
    vectorField& pf = force.boundaryFieldRef()[patchi];
    pf += fBulk; 
}

void Foam::functionObjects::rheoForces::addToFields
(
    const label patchi,
    const vectorField& fMa,
    const vectorField& fs
)
{
    force_[0] += sum(fMa);
    force_[1] += sum(fs);
    if (!writeFields_)
    {
        return;
    }

    if(!calcBulkForce_)
    {
        auto& force = lookupObjectRef<volVectorField>(scopedName("force"));
        vectorField& pf = force.boundaryFieldRef()[patchi];
        pf += fMa + fs;
    }

}

void Foam::functionObjects::rheoForces::writeIntegratedForce
(
    const string& descriptor,
    const vector& fm0,
    const vector& fm1,
    const vector& fm2,
    autoPtr<OFstream>& osPtr
) const
{
    vector total = fm0 + fm1 + fm2;

    if (writeToFile())
    {
        Ostream& os = osPtr();

        writeCurrentTime(os);

        os  << tab << total
            << tab << fm0
            << tab << fm1
            << tab << fm2;

        os  << endl;
    }
}


void Foam::functionObjects::rheoForces::writeForces()
{
    //Log << type() << " " << name() << " write:" << nl
    //    << force_[0] << nl
    //    << force_[1] << nl
    //    << force_[2] << nl;

    writeIntegratedForce
    (
        "forces",
        force_[0],
        force_[1],
        force_[2],
        forceFilePtr_
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::rheoForces::rheoForces(
    const word &name,
    const Time &runTime,
    const dictionary &dict,
    bool readFields)
    : fvMeshFunctionObject(name, runTime, dict),
      writeFile(mesh_, name),
      force_(3),
      forceFilePtr_(),
      patchSet_(),
      areaPatches_(),
      pName_("p"),
      UName_("U"),
      UsName_("Us"),
      gammaName_("gamma"),
      writeFields_(false),
      initialised_(false),
      calcBulkForce_(true),
      isCIT_(false),
      x0_(-1.0),
      L0_(-1.0),
      Ma_(0.0),
      Theta_(0.0),
      Tr_(0.0),
      Is_(tensor::I),
      Bq_(0.0),
      isExtensional_(false)
{
    if (readFields)
    {
        read(dict);
        Log << endl;
    }
}


Foam::functionObjects::rheoForces::rheoForces
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    bool readFields
)
:
    fvMeshFunctionObject(name, obr, dict),
    writeFile(mesh_, name),
    force_(3),
    forceFilePtr_(),
    patchSet_(),
    areaPatches_(),
    pName_("p"),
    UName_("U"),
    UsName_("Us"),
    gammaName_("gamma"),
    writeFields_(false),
    initialised_(false),
    calcBulkForce_(true),
    isCIT_(false),
    x0_(-1.0),
    L0_(-1.0),
    Ma_(0.0),
    Theta_(0.0),
    Tr_(0.0),
    Is_(tensor::I),
    Bq_(0.0),
    isExtensional_(false)
{
    if (readFields)
    {
        read(dict);
        Log << endl;
    }

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::rheoForces::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    initialised_ = false;

    Info<< type() << " " << name() << ":" << nl;

    areaPatches_ = dict.get<wordRes>("areaPatches");

    calcBulkForce_ =
        dict.get<bool>("calcBulkForce");

    if (calcBulkForce_)
    {
        patchSet_ =
            mesh_.boundaryMesh().patchSet
            (
                dict.get<wordRes>("patches")
            );
    }

    // Optional field name entries
    if (dict.readIfPresent<word>("p", pName_))
    {
        Info<< "    p: " << pName_ << endl;
    }
    if (dict.readIfPresent<word>("U", UName_))
    {
        Info<< "    U: " << UName_ << endl;
    }
    if (dict.readIfPresent<word>("Us", UsName_))
    {
        Info << "  Us: " << UsName_ << endl;
    }
    if (dict.readIfPresent<word>("gamma", gammaName_))
    {
        Info << "  gamma: " << gammaName_ << endl;
    }
    if (dict.readIfPresent<scalar>("x0", x0_))
    {
        Info << "  x0: " << x0_ << endl;
    }

    writeFields_ = dict.getOrDefault("writeFields", false);

    if (writeFields_)
    {
        Info<< "    Fields will be written" << endl;

        volVectorField* forcePtr
        (
            new volVectorField
            (
                IOobject
                (
                    scopedName("force"),
                    time_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedVector(dimless, Zero)
            )
        );

        mesh_.objectRegistry::store(forcePtr);
    }

    if (foundObject<interfaceModel>("interfaceModel"))
    {
        const interfaceModel& im =
            lookupObject<interfaceModel>("interfaceModel");
        
        const dictionary& imDict = im.modelDict();
        const word& isothermName = imDict.get<word>("interfaceModel");
        if (isothermName == "extensionalViscosity")
        {
            isExtensional_ = true;
        }
    }

    return true;
}


void Foam::functionObjects::rheoForces::calcForces()
{
    initialise();
    resetFields();
    
    const HashTable<const faMesh*> meshes =
        obr().lookupClass<faMesh>();
    wordList available = obr().sortedNames<faMesh>();
    if (!available.size())
    {
        FatalErrorInFunction
            << "No valid number of area meshes"
            << exit(FatalError);
    }
    const word& areaName = available.first();
    const faMesh& aMesh = *meshes[areaName];

    const edgeVectorField::Boundary& Leb = aMesh.Le().boundaryField();

    // CIT: calculate rod perimeter at the interface
    if (x0_ > SMALL && L0_ < -SMALL)
    {
        L0_ = 0.;
        for (const word& apName : areaPatches_)
        {
            const label& patchi = aMesh.boundary().findPatchID(apName);
            scalarField sA(mag(Leb[patchi]));
            L0_ += sum(x0_ * sA);
        }  
        reduce(L0_, sumOp<scalar>());
        Info<< type() << " " << name() << ":" << nl
            << "    L0 = " << L0_ << endl;
    }   

    if (calcBulkForce_)
    {
        // Bulk
        const volVectorField& U = lookupObject<volVectorField>(UName_);
        const volScalarField& p = lookupObject<volScalarField>(pName_);
        const surfaceVectorField::Boundary& Sfb = mesh_.Sf().boundaryField();
        const volTensorField tau(fvc::grad(U) + T(fvc::grad(U)));
        const volTensorField::Boundary& taub = tau.boundaryField();
        const volScalarField::Boundary& pb = p.boundaryField();

        for (const label& patchi : patchSet_)
        {
            vectorField fBulk = 
                -Sfb[patchi] & (-pb[patchi]*tensor::I + taub[patchi]);
                
            // Lateral surface of the disc
            /*forAll(fBulk, i)
            {
            	if (mag(Sfb[patchi][i][1])/mag(Sfb[patchi][i]) < 1e-4)
            	{
            		fBulk[i] *= 0.0;
            	}
            }*/
            // CIT: Make dS_f dimensionless with L_0^2
            if (x0_ > SMALL)
            {
                scalarField SfbNew 
                    = mag(Sfb[patchi])*sqr(x0_/L0_);
                fBulk = fBulk*SfbNew/mag(Sfb[patchi]);
            }
            addToFields(patchi, fBulk);
        }
    }
    // Interface    
    tmp<areaTensorField> tTaus = taus();
    const areaTensorField::Boundary& tausb = tTaus().boundaryField();

    tmp<areaScalarField> tPi = Pi();
    const areaScalarField::Boundary& Pib = tPi().boundaryField();

    //const areaScalarField& gamma = lookupObject<areaScalarField>(gammaName_);
    //const areaScalarField::Boundary& gsb = gamma.boundaryField();

    for (const word& apName : areaPatches_)
    {
        const label& patchi = aMesh.boundary().findPatchID(apName);
        scalarField dl(mag(Leb[patchi]));
        vectorField n(-Leb[patchi]/dl);

        // Marangoni
        vectorField fMa = (-Pib[patchi] * Is_) & n;

        // Surface stress tensor
        vectorField fs = tausb[patchi] & n;

        // CIT: make dl dimensionless with L_0
        if (x0_ > SMALL)
        {
            dl = dl*x0_/L0_;
        }
        
        fs *= dl;
        fMa *= dl;

        if (isExtensional_)
        {
            fs /= Bq_;
        }
        addToFields(patchi, fMa, fs);
    }      
}

Foam::vector Foam::functionObjects::rheoForces::forceEff() const
{
    return force_[0] + force_[1] + force_[2];
}


bool Foam::functionObjects::rheoForces::execute()
{
    calcForces();

    reduce(force_[0], sumOp<vector>());
    reduce(force_[1], sumOp<vector>());
    reduce(force_[2], sumOp<vector>());

    if (Pstream::master())
    {
        createFiles();

        writeForces();

        Log << endl;
    }

    // Write state/results information
    setResult("F_Ma", mag(force_[0]));
    setResult("F_s", mag(force_[1]));
    setResult("F_Bulk", mag(force_[2]));

    return true;
}


bool Foam::functionObjects::rheoForces::write()
{
    if (writeFields_)
    {
        lookupObject<volVectorField>(scopedName("force")).write();
    }

    return true;
}


// ************************************************************************* //
