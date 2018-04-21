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

#include "WA2017IDDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> WA2017IDDES<BasicTurbulenceModel>::alpha() const
{
    return max
    (
        0.25 - this->y_ / static_cast<const volScalarField&>(IDDESDelta_.hmax()),
        scalar(-5)
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> WA2017IDDES<BasicTurbulenceModel>::ft
(
    const volScalarField& magGradU
) const
{
    return tanh(pow3(sqr(ct_) * this->rd(this->nut_, magGradU)));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> WA2017IDDES<BasicTurbulenceModel>::fl
(
    const volScalarField& magGradU
) const
{
    return tanh(pow(sqr(cl_)*this->rd(this->nu(), magGradU), 10));
}

template<class BasicTurbulenceModel>
void WA2017IDDES<BasicTurbulenceModel>::calc_fdes()
{
    // Calculate RANS length scale lrans
    const volScalarField lrans = max
                                (
                                    sqrt(Rnu_/S_),
                                    dimensionedScalar("SMALL", dimLength, SMALL)
                                );
                                
    // Calculate the IDDES length scale liddes

    const volScalarField alpha(this->alpha());
    
    // e^(alpha^2)
    const volScalarField expTerm(exp(sqr(alpha)));

    const volScalarField magGradU(mag(fvc::grad(this->U_)));

    // fe1
    tmp<volScalarField> fHill = 2*(pos(alpha)*pow(expTerm, -11.09) + 
                                neg(alpha)*pow(expTerm, -9.0));
    
    // fb
    tmp<volScalarField> fStep = min(2*pow(expTerm, -9.0), scalar(1));
    
    fd_ = this->fd(magGradU);

    fdtilda_ = max(1 - fd_, fStep);

    // fe2
    tmp<volScalarField> fAmp = 1 - max(ft(magGradU), fl(magGradU));
  
    // fe
    fe_ = max(fHill - 1, scalar(0))*fAmp;
    
    const volScalarField liddes = max
                                  (
                                      fdtilda_*(1 + fe_)*lrans +
                                      (1 - fdtilda_)*CDES_*this->delta(),
                                      dimensionedScalar("SMALL", dimLength, SMALL)
                                  );

    fdes_ = lrans / liddes;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
WA2017IDDES<BasicTurbulenceModel>::WA2017IDDES
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
    WA2017DDES<BasicTurbulenceModel>
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

    fwStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "fwStar",
            this->coeffDict_,
            0.424
        )
    ),

    cl_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cl",
            this->coeffDict_,
            3.55
        )
    ),

    ct_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ct",
            this->coeffDict_,
            1.63
        )
    ),
    
    fdtilda_
    (
        IOobject
        (
            "fdtilda",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), 0)
    ),

    fe_
    (
        IOobject
        (
            "fe",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), 0)
    ),
    
    IDDESDelta_(refCast<IDDESDelta>(this->delta_()))
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool WA2017IDDES<BasicTurbulenceModel>::read()
{
    if (WA2017DDES<BasicTurbulenceModel>::read())
    {
        fwStar_.readIfPresent(this->coeffDict());
        cl_.readIfPresent(this->coeffDict());
        ct_.readIfPresent(this->coeffDict());
        
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
