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

defineTypeNameAndDebug(outputForceFunctionObject, 0);

addRemovableToRunTimeSelectionTable
(
    functionObject,
    outputForceFunctionObject,
    dictionary
);


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = 7306e0537ad00ff8242d5996a9e3fef2fa3dbafc
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void outputForce_7306e0537ad00ff8242d5996a9e3fef2fa3dbafc(bool load)
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

const fvMesh& outputForceFunctionObject::mesh() const
{
    return refCast<const fvMesh>(obr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

outputForceFunctionObject::outputForceFunctionObject
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

outputForceFunctionObject::~outputForceFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool outputForceFunctionObject::read(const dictionary& dict)
{
    if (false)
    {
        Info<<"read outputForce sha1: 7306e0537ad00ff8242d5996a9e3fef2fa3dbafc\n";
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool outputForceFunctionObject::execute()
{
    if (false)
    {
        Info<<"execute outputForce sha1: 7306e0537ad00ff8242d5996a9e3fef2fa3dbafc\n";
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool outputForceFunctionObject::write()
{
    if (false)
    {
        Info<<"write outputForce sha1: 7306e0537ad00ff8242d5996a9e3fef2fa3dbafc\n";
    }

//{{{ begin code
    #line 65 "/home/clove/rheoTool/of90/custom_models/sphere_con/Oldroyd-BLog/system/controlDict/functions/outputForce"
//const volScalarField& p = mesh().lookupObject<volScalarField>("p");
      label sphere = mesh().boundaryMesh().findPatchID("sphere");
      //const volVectorField& U = mesh().lookupObject<volVectorField>("U");
      const volSymmTensorField& tau = mesh().lookupObject<volSymmTensorField>("tau");
      //const volScalarField& p = mesh().lookupObject<volScalarField>("p");
      //const dictionary& constDict = mesh().lookupObject<IOdictionary>("constitutiveProperties");

      //dimensionedScalar rho_(constDict.subDict("parameters").lookup("rho"));
      //dimensionedScalar etaS_(constDict.subDict("parameters").lookup("etaS"));
      //dimensionedScalar etaP_(constDict.subDict("parameters").lookup("etaP"));
  
      //label cyl = mesh().boundaryMesh().findPatchID("cylinder");  //this appears to be searching a subset of the mesh() along the boundaries, locating the patches of interest
                                                                  //by the patch ID "cylinder" which is probably the patches located along the cylinder boundary condition.
                                                                  //the created label is later used to index through the boundaryField and obtain values relevant to the calculations
      scalarList list;
    
      // Compute cd (drag coefficient) this is the original computation
 
      //volTensorField L(fvc::grad(U));  //using fvc to access the grad function to get the gradient of velocity vector field U
	    //volSymmTensorField F(tau); //symm returns the symmetric part of a second rank tensor. so this is shear stress added to 
      //vector Fpatch = (F.boundaryField()[cyl]);                                                                                  //the symmetric velocity gradient multiplied by a viscosity coefficient minus pressure*density*Identity matrix
                                                                                        //thus this appears to be a force balance constructing a force tensor, F
      //vector Fpatch = gSum( ( -mesh().boundaryMesh()[cyl].faceAreas() ) & F.boundaryField()[cyl] )/(etaS_ + etaP_).value(); summation of force times area over etaS and etaP to obtain drag coefficient
                                                                                                                        // note use of boundaryMesh()[cyl] only looks at values along cylinder interface
      
      //compute force applied to cylinder by pressure (neglects shear stress and viscous forces)
      //vectorField area_normal = mesh().Sf().boundaryField()[cyl]; /// mesh().magSf().boundaryField()[cyl];
      //vector Fpatch = F.boundaryField()[sphere]; //F appears to have units of force/area, multiplying by area normal should give force vector
      //vector Fpatch = gSum((area_normal)*(p));
      //gSum((F.boundaryField()[cyl]).component(symmTensor::XX))
      //volScalarField x2(F.component(symmTensor::XX)); //take the XX component of the symmTensor
      
      
      //tau::XY;
      //olScalarField tau_xy(tau.component(symmTensor::XY));
      //vectorField area_normal = mesh().Sf().boundaryField()[sphere];
      volSymmTensorField F(tau);
      vectorField area_normal = mesh().Sf().boundaryField()[sphere];
      vector Fpatch = gSum(-area_normal & F.boundaryField()[sphere]);
      //float tau_sum = gSum(tau_xy_sphere(tau_xy.boundaryField()[sphere]));
      //int len_tau = sizeof(tau_xy.boundaryField()[sphere]);
      //float sum_force_xy = 0;
      //for (int i=0;i<len_tau;i++){
        //sum_force_xy+=(tau_xy.boundaryField()[sphere][i]*area_normal*-1); //sum all shear stress acting along surface of the sphere
      //}
      //float tau_sum=gSum(tau_xy_sphere);
      //vector Fpatch = F.boundaryField()[sphere];
      //vector tautype = tau[sphere];
      list.append(mesh().time().value());  // Time (col 0)  
      list.append(Fpatch.x());             // xx component of force   (col 1)  this appears to be using .x() to reference the x component of the vector Fpatch, Fpatch being declared as a vector object

          // Write data: this function is expected to access the writeData file which is responsible for constructing all simulation result files.
          //without that, no data is saved from simulation run. 

      string comsh;           
      string filename("sphere_force.txt");
	    std::stringstream doub2str; doub2str.precision(12);

      comsh = "./writeData " + filename; //set up to write data to cyl_force.txt
      forAll(list, id)
      {
        doub2str.str(std::string());
        doub2str << list[id]; 
        comsh += " " + doub2str.str(); //write time signature for controlDict function and force 
        //comsh += "  " + list[id];
      }
           
      if (Pstream::master())
      {
	      system(comsh);
      }
//}}} end code

    return true;
}


bool outputForceFunctionObject::end()
{
    if (false)
    {
        Info<<"end outputForce sha1: 7306e0537ad00ff8242d5996a9e3fef2fa3dbafc\n";
    }

//{{{ begin code
    
//}}} end code

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

