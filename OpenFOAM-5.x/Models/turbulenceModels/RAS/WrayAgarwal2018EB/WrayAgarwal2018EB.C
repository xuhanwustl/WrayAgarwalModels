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

#include "WrayAgarwal2018EB.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwal2018EB<BasicTurbulenceModel>::LR2() const
{
    const volScalarField Clnu = (4.0+sqrt(this->chi()))*this->nu();

    return max(C3kw_*Rnu_, Clnu) / (S_ + Clnu/sqr(dimensionedScalar("Lref", 
                                                                    dimensionSet(0, 1, 0, 0, 0), 
                                                                    Lref_.value())
                                                 )
                                   );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
WrayAgarwal2018EB<BasicTurbulenceModel>::WrayAgarwal2018EB
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
    WrayAgarwal2018<BasicTurbulenceModel>
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
    
    C3kw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3kw",
            this->coeffDict_,
            0.171
        )
    ),
    
    Lref_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Lref",
            this->coeffDict_,
            0.0
        )
    ),
    
    PR_
    (
        IOobject
        (
            "PR",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    // Set the coefficients to fit the WrayAgarwal2018EBEB model 
    this->modifyCoeff(this->Cw_, 5.97);
    this->modifyCoeff(this->C1ke_, 0.094);
    this->modifyCoeff(this->C1kw_, 0.2);
    this->modifyCoeff(this->C2ke_, 1.24);
    this->modifyCoeff(this->C2kw_, 2.63);
    
    // Make sure the reference length scale is valid
    if (Lref_.value() <= SMALL) 
    {
        FatalErrorIn(type)
            << "Please assign a postive value to coefficient Lref."
            << exit(FatalError);
    }
    
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool WrayAgarwal2018EB<BasicTurbulenceModel>::read()
{
    if (WrayAgarwal2018<BasicTurbulenceModel>::read())
    {
        C3kw_.readIfPresent(this->coeffDict());
        Lref_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicTurbulenceModel>
void WrayAgarwal2018EB<BasicTurbulenceModel>::correct()
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

    // Calculate Vorticity Magnitude W_
    volScalarField W2(2.0*magSqr(skew(fvc::grad(U))));
    W_ = sqrt(W2);
    bound(W_, dimensionedScalar("0", W_.dimensions(), SMALL));
    bound(W2, dimensionedScalar("0", W2.dimensions(), SMALL));
    
    // Calculate switch function f1_
    this->calc_f1();

    // Define and solve PR Equation
    tmp<fvScalarMatrix> PR_Eqn
    (
      - LR2()*fvm::laplacian(PR_)
      + fvm::Sp(1.0, PR_)
     ==
        S_*Rnu_
    );

    PR_Eqn.ref().relax();
    solve(PR_Eqn);
    PR_.correctBoundaryConditions();

    // Define and solve Rnu Equation
    tmp<fvScalarMatrix> RnuEqn
    (
        fvm::ddt(alpha, rho, Rnu_)
      + fvm::div(alphaRhoPhi, Rnu_)
      - fvm::laplacian(alpha*rho*this->DRnuEff(f1_), Rnu_)
     ==
        alpha*rho*(this->C1(f1_)-1.0)*fvm::Sp(S_, Rnu_)
      + alpha*rho*f1_*C2kw_*fvm::Sp((fvc::grad(Rnu_)&fvc::grad(S_))/S_, Rnu_)
      - alpha*rho*(1.0-f1_)*min(C2ke_*Rnu_*Rnu_*magSqr(fvc::grad(S_))/S2,
                                Cm_*magSqr(fvc::grad(Rnu_)))
      + alpha*rho*PR_
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
