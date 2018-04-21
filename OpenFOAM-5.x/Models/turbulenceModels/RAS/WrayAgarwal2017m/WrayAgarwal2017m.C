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

#include "WrayAgarwal2017m.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
WrayAgarwal2017m<BasicTurbulenceModel>::WrayAgarwal2017m
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
    WrayAgarwal2017<BasicTurbulenceModel>
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

    Cm_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cm",
            this->coeffDict_,
            8.0
        )
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool WrayAgarwal2017m<BasicTurbulenceModel>::read()
{
    if (WrayAgarwal2017<BasicTurbulenceModel>::read())
    {        
        Cm_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicTurbulenceModel>
void WrayAgarwal2017m<BasicTurbulenceModel>::correct()
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
