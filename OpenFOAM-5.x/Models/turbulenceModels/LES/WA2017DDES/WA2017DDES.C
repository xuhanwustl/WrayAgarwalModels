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

\*---------------------------------------------------------------------------*/

#include "WA2017DDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> WA2017DDES<BasicTurbulenceModel>::rd
(
    const volScalarField& nur,
    const volScalarField& magGradU
) const
{
    return min
    (
        nur
       /(
           max
           (
               magGradU,
               dimensionedScalar("SMALL", magGradU.dimensions(), SMALL)
           )*sqr(this->kappa_*this->y_)
       ),
       scalar(10)
    );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WA2017DDES<BasicTurbulenceModel>::fd
(
    const volScalarField& magGradU
) const
{
    return 1 - tanh(pow3(Cd1_*rd(this->nuEff(), magGradU)));
}

template<class BasicTurbulenceModel>
void WA2017DDES<BasicTurbulenceModel>::calc_fdes()
{
    // Calculate RANS length scale lrans
    const volScalarField lrans = max
                                (
                                    sqrt(Rnu_/S_),
                                    dimensionedScalar("SMALL", dimLength, SMALL)
                                );
                                
    // Calculate the DDES length scale lddes
    fd_ = fd(mag(fvc::grad(this->U_)));
    const volScalarField lddes = max
                                 (
                                    lrans - fd_*max(dimensionedScalar("ZERO", dimLength, 0), 
                                                    lrans - CDES_*this->delta()),
                                    dimensionedScalar("SMALL", dimLength, SMALL)
                                 );

    fdes_ = lrans / lddes;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
WA2017DDES<BasicTurbulenceModel>::WA2017DDES
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    WA2017DES<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    ),
    
    Cd1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cd1",
            this->coeffDict_,
            4.0
        )
    ),

    fd_
    (
        IOobject
        (
            "fd",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), 0)
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool WA2017DDES<BasicTurbulenceModel>::read()
{
    if (WA2017DES<BasicTurbulenceModel>::read())
    {
        Cd1_.readIfPresent(this->coeffDict());
        
        return true;
    }
    else
    {
        return false;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
