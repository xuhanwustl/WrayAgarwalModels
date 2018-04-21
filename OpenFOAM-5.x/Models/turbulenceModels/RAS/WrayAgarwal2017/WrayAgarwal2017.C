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

#include "WrayAgarwal2017.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void WrayAgarwal2017<BasicTurbulenceModel>::calc_f1()
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
WrayAgarwal2017<BasicTurbulenceModel>::WrayAgarwal2017
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
    WrayAgarwalRASModel<BasicTurbulenceModel>
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

    y_(wallDist::New(this->mesh_).y())
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void WrayAgarwal2017<BasicTurbulenceModel>::correct()
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
    calc_f1();

    // Define and solve Rnu Equation
    tmp<fvScalarMatrix> RnuEqn
    (
        fvm::ddt(alpha, rho, Rnu_)
      + fvm::div(alphaRhoPhi, Rnu_)
      - fvm::laplacian(alpha*rho*this->DRnuEff(f1_), Rnu_)
     ==
        alpha*rho*this->C1(f1_)*fvm::Sp(S_, Rnu_)
      + alpha*rho*f1_*C2kw_*fvm::Sp((fvc::grad(Rnu_)&fvc::grad(S_))/S_, Rnu_)
      - alpha*rho*(1.0-f1_)*C2ke_*fvm::Sp(Rnu_*magSqr(fvc::grad(S_))/S2, Rnu_)
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
