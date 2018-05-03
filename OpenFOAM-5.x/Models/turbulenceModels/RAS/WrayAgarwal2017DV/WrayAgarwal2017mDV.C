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

#include "WrayAgarwal2017mDV.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
WrayAgarwal2017mDV<BasicTurbulenceModel>::WrayAgarwal2017mDV
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
    WrayAgarwal2017m<BasicTurbulenceModel>
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
    
    gamma_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma",
            this->coeffDict_,
            1.4
        )
    ),

    Rsp_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Rsp",
            this->coeffDict_,
            286.9
        )
    ),

    Cr1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cr1",
            this->coeffDict_,
            0.01
        )
    ),

    Cr2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cr2",
            this->coeffDict_,
            0.01
        )
    ),

    Crho1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Crho1",
            this->coeffDict_,
            1.4
        )
    ),

    Crho2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Crho2",
            this->coeffDict_,
            2.0
        )
    ),

    sigmaRho_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaRho",
            this->coeffDict_,
            1.0
        )
    ),

    T_
    (
        IOobject
        (
            "T",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    // Density variance fields
    rhoPrime2_
    (
        IOobject
        (
            "rhoPrime2",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    fcomp_
    (
        IOobject
        (
            "fcomp",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0.0", dimensionSet(0, 2, -2, 0, 0), 0.0)
    ),
    
    gradUSum_
    (
        IOobject
        (
            "gradUSum",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0.0", dimensionSet(0, 0, -1, 0, 0), 0.0)
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool WrayAgarwal2017mDV<BasicTurbulenceModel>::read()
{
    if (WrayAgarwal2017m<BasicTurbulenceModel>::read())
    {
        gamma_.readIfPresent(this->coeffDict());
        Rsp_.readIfPresent(this->coeffDict());
        Cr1_.readIfPresent(this->coeffDict());
        Cr2_.readIfPresent(this->coeffDict());
        Crho1_.readIfPresent(this->coeffDict());
        Crho2_.readIfPresent(this->coeffDict());
        sigmaRho_.readIfPresent(this->coeffDict());
        
        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwal2017mDV<BasicTurbulenceModel>::DRhoPrime2Eff(volScalarField Switch) const
{
    return tmp<volScalarField>
    (
        new volScalarField("DRhoPrime2Eff", sigmaRho_*this->nut_ + this->nu())
    );
}

template<class BasicTurbulenceModel>
void WrayAgarwal2017mDV<BasicTurbulenceModel>::correct()
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

    eddyViscosity<RASModel<BasicTurbulenceModel> >::correct();

    // Calculate Strain Rate Magnitude S_
    volScalarField S2(2.0*magSqr(symm(fvc::grad(U))));
    S_ = sqrt(S2);
    bound(S_, dimensionedScalar("0", S_.dimensions(), SMALL));
    bound(S2, dimensionedScalar("0", S2.dimensions(), SMALL));

    // Calculate switch function f1_
    this->calc_f1();
    
    // Calculate the speed of sound
    volScalarField a = sqrt(gamma_ * T_ * Rsp_ * 
                            dimensionedScalar("RspUnit", 
                                              dimensionSet(0, 2, -2, -1, 0), 
                                              1));

    // Calculate the duj/dxj
    const volTensorField tgradU = fvc::grad(this->U_);    
    gradUSum_ = tr(tgradU);

    // Define and solve rhoPrime2 Equation
    tmp<fvScalarMatrix> rhoPrime2Eqn
    (
        fvm::ddt(alpha, rhoPrime2_)
      + fvm::div(alphaRhoPhi/fvc::interpolate(rho), rhoPrime2_)
      - fvm::laplacian(alpha*DRhoPrime2Eff(f1_), rhoPrime2_)
     ==
        2.0*Crho1_*alpha*this->nut_*magSqr(fvc::grad(rho))
      - alpha*fvm::Sp(2.0*gradUSum_ + Crho2_*S_, rhoPrime2_)
    );
    
    rhoPrime2Eqn.ref().relax();
    solve(rhoPrime2Eqn);
    bound(rhoPrime2_, dimensionedScalar("0", rhoPrime2_.dimensions(), 0.0));
    rhoPrime2_.correctBoundaryConditions();
    
    // Calculate fcomp
    fcomp_ = rhoPrime2_*sqr(a/rho)*(Cr1_*gradUSum_/S_ - Cr2_);

    // Define and solve Rnu Equation
    tmp<fvScalarMatrix> RnuEqn
    (
        fvm::ddt(alpha, rho, Rnu_)
      + fvm::div(alphaRhoPhi, Rnu_)
      - fvm::laplacian(alpha*rho*this->DRnuEff(f1_), Rnu_)
     ==
        alpha*rho*this->C1(f1_)*fvm::Sp(S_, Rnu_)
      + alpha*rho*f1_*C2kw_*fvm::Sp((fvc::grad(Rnu_)&fvc::grad(S_))/S_, Rnu_)
      - alpha*rho*(1.0-f1_)*min(C2ke_*Rnu_*Rnu_*magSqr(fvc::grad(S_))/S2,
                                Cm_*magSqr(fvc::grad(Rnu_)))
    );

    RnuEqn.ref().relax();
    solve(RnuEqn);
    bound(Rnu_, dimensionedScalar("0", Rnu_.dimensions(), 0.0));
    Rnu_.correctBoundaryConditions();

    this->correctNut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
