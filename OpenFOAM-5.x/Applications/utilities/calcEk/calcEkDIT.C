/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

Application
    calcEk

Description
    Utility to calculate turbulent energy spectra for Decaying Isotropic 
    Turbulence (DIT) test.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentTransportModel.H"

#include "Kmesh.H"
#include "fft.H"
#include "calcEk.H"
#include "graph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    
    #include "createFields.H"

    #include "readTransportProperties.H"
    #include "readTurbulenceProperties.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating Ek...\n" << endl;

    calcEk(U, K).write
    (
        runTime.path()/"Ek"/runTime.timeName(),
        "Ek",
        runTime.graphFormat()
    );
    
    runTime.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
