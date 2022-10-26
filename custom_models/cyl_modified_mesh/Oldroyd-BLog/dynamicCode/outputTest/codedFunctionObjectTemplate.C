/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "codedFunctionObjectTemplate.H"
#include "fvCFD.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(outputTestFunctionObject, 0);

addRemovableToRunTimeSelectionTable
(
    functionObject,
    outputTestFunctionObject,
    dictionary
);


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = ba4dbdd3ab8bce4f927d847044fc82b1a9aedbbf
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void outputTest_ba4dbdd3ab8bce4f927d847044fc82b1a9aedbbf(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const fvMesh& outputTestFunctionObject::mesh() const
{
    return refCast<const fvMesh>(obr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

outputTestFunctionObject::outputTestFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::regionFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

outputTestFunctionObject::~outputTestFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool outputTestFunctionObject::read(const dictionary& dict)
{
    if (false)
    {
        Info<<"read outputTest sha1: ba4dbdd3ab8bce4f927d847044fc82b1a9aedbbf\n";
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool outputTestFunctionObject::execute()
{
    if (false)
    {
        Info<<"execute outputTest sha1: ba4dbdd3ab8bce4f927d847044fc82b1a9aedbbf\n";
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool outputTestFunctionObject::write()
{
    if (false)
    {
        Info<<"write outputTest sha1: ba4dbdd3ab8bce4f927d847044fc82b1a9aedbbf\n";
    }

//{{{ begin code
    #line 67 "/home/clove/rheoTool/of90/custom_models/cyl_modified_mesh/Oldroyd-BLog/system/controlDict/functions/outputTest"
// Lookup/create variable 
	     
	   const volVectorField& U = mesh().lookupObject<volVectorField>("U");
           const volSymmTensorField& tau = mesh().lookupObject<volSymmTensorField>("tau");
           const volScalarField& p = mesh().lookupObject<volScalarField>("p");
           const dictionary& constDict = mesh().lookupObject<IOdictionary>("constitutiveProperties");
           dimensionedScalar rho_(constDict.subDict("parameters").lookup("rho"));
           dimensionedScalar etaS_(constDict.subDict("parameters").lookup("etaS"));
           dimensionedScalar etaP_(constDict.subDict("parameters").lookup("etaP"));
  
           label cyl = mesh().boundaryMesh().findPatchID("cylinder");
           scalarList list;
    
          // Compute cd
 
           volTensorField L(fvc::grad(U));

	   volSymmTensorField F(tau + symm( L + L.T() ) * etaS_ - p * symmTensor::I * rho_);

           vector Fpatch = gSum( ( -mesh().boundaryMesh()[cyl].faceAreas() ) & F.boundaryField()[cyl] )/(etaS_ + etaP_).value();
  
           list.append(mesh().time().value());  // Time (col 0)  
           list.append(Fpatch.x());             // Cd   (col 1)  

          // Write data

           string comsh;           
           string filename("test.txt");
	   std::stringstream doub2str; doub2str.precision(12);

           comsh = "./writeData " + filename;
           forAll(list, id)
            {
              doub2str.str(std::string());
              doub2str << list[id]; 
              comsh += " " + doub2str.str();
            }
           
           if (Pstream::master())
            {
	      system(comsh);
            }
//}}} end code

    return true;
}


bool outputTestFunctionObject::end()
{
    if (false)
    {
        Info<<"end outputTest sha1: ba4dbdd3ab8bce4f927d847044fc82b1a9aedbbf\n";
    }

//{{{ begin code
    
//}}} end code

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

