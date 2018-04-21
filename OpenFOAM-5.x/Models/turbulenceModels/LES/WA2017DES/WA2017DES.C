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

#include "WA2017DES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void WA2017DES<BasicTurbulenceModel>::calc_f1()
{    
    tmp<volScalarField> eta = y_*sqrt(Rnu_*S_)/(20.0*this->nu());
    f1_ = tanh(pow((1.0 + 20.0*eta) / 
                   (1.0 + sqr(max(y_*sqrt(Rnu_*S_), 1.5*Rnu_) / 
                              (20.0*this->nu())))
                   ,
                   4.0));  
    f1_ = min(f1_,0.9);
    bound(f1_,SMALL);
}

template<class BasicTurbulenceModel>
void WA2017DES<BasicTurbulenceModel>::calc_fdes()
{
    fdes_ = max
            (
                sqrt(Rnu_) / (sqrt(S_) * CDES_*this->delta()),
                scalar(1)
            );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WA2017DES<BasicTurbulenceModel>::blendFactor() const
{
    return scalar(1.0) - scalar(1.0) / fdes_;
}

template<class BasicTurbulenceModel>
void WA2017DES<BasicTurbulenceModel>::calcBlendFactors()
{
    blendfactor_ = blendFactor();
    
    UBlendingFactor_ = fvc::interpolate(blendfactor_);
    RnuBlendingFactor_ = UBlendingFactor_;
    pBlendingFactor_ = UBlendingFactor_;
    KBlendingFactor_ = UBlendingFactor_;
    eBlendingFactor_ = UBlendingFactor_;
    hBlendingFactor_ = UBlendingFactor_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
WA2017DES<BasicTurbulenceModel>::WA2017DES
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
    WrayAgarwalLESModel<BasicTurbulenceModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),
    
    CDES_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CDES",
            this->coeffDict_,
            0.41
        )
    ),

    y_(wallDist::New(this->mesh_).y()),

    fdes_
    (
        IOobject
        (
            "fdes",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0.0", dimensionSet(0, 0, 0, 0, 0), 0.0)
    ),

    outDelta_
    (
        IOobject
        (
            "outDelta",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0.0", dimensionSet(0, 1, 0, 0, 0), 0.0)
    ),

    blendfactor_
    (
        IOobject
        (
            "blendfactor",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), 0)
    ),

    UBlendingFactor_
    (
        IOobject
        (
            "UBlendingFactor",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(blendfactor_)
    ),

    RnuBlendingFactor_
    (
        IOobject
        (
            "RnuBlendingFactor",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(blendfactor_)
    ),

    pBlendingFactor_
    (
        IOobject
        (
            "pBlendingFactor",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(blendfactor_)
    ),

    KBlendingFactor_
    (
        IOobject
        (
            "KBlendingFactor",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(blendfactor_)
    ),

    eBlendingFactor_
    (
        IOobject
        (
            "eBlendingFactor",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(blendfactor_)
    ),

    hBlendingFactor_
    (
        IOobject
        (
            "hBlendingFactor",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(blendfactor_)
    )
{
    outDelta_ = this->delta();
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool WA2017DES<BasicTurbulenceModel>::read()
{
    if (WrayAgarwalLESModel<BasicTurbulenceModel>::read())
    {
        CDES_.readIfPresent(this->coeffDict());
        
        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicTurbulenceModel>
void WA2017DES<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }
    
    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;

    LESeddyViscosity<BasicTurbulenceModel>::correct();
    
    // Calculate Strain Rate Magnitude S_
    volScalarField S2(2.0*magSqr(symm(fvc::grad(U))));
    S_ = sqrt(S2);
    bound(S_, dimensionedScalar("0", S_.dimensions(), SMALL));
    bound(S2, dimensionedScalar("0", S2.dimensions(), SMALL));

    // Calculate switch function f1_
    calc_f1();

    // Calculate hybrid switch function fdes_
    calc_fdes();
    const volScalarField fdes2 = sqr(fdes_);
    
    // Blend Scheme
    calcBlendFactors();
    
    // Define and solve Rnu Equation
    tmp<fvScalarMatrix> RnuEqn
    (
        fvm::ddt(alpha, rho, Rnu_)
      + fvm::div(alphaRhoPhi, Rnu_)
      - fvm::laplacian(alpha*rho*this->DRnuEff(f1_), Rnu_)
     ==
        alpha*rho*this->C1(f1_)*fvm::Sp(S_, Rnu_)
      + alpha*rho*f1_*C2kw_*fvm::Sp((fvc::grad(Rnu_)&fvc::grad(S_))/S_/fdes2, Rnu_)
      - alpha*rho*(1.0-f1_)*C2ke_*fvm::Sp(Rnu_*magSqr(fvc::grad(S_))/S2/fdes2, Rnu_)
    );

    RnuEqn.ref().relax();
    solve(RnuEqn);
    bound(Rnu_, dimensionedScalar("0", Rnu_.dimensions(), 0.0));
    Rnu_.correctBoundaryConditions();

    this->correctNut();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> WA2017DES<BasicTurbulenceModel>::LESRegion() const
{
    tmp<volScalarField> tLESRegion
    (
        new volScalarField
        (
            IOobject
            (
                "DES::LESRegion",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            neg(scalar(1) - fdes_)
        )
    );

    return tLESRegion;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
