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

#include "WrayAgarwalBase.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class TurbulenceModel, class BasicTurbulenceModel>
void WrayAgarwalBase<TurbulenceModel, BasicTurbulenceModel>::modifyCoeff
(
    dimensionedScalar& coeff, 
    scalar value
)
{
    if (this->coeffDict_.remove(coeff.name())) 
    {
        coeff = dimensioned<scalar>::lookupOrAddToDict
                (
                    coeff.name(),
                    this->coeffDict_,
                    value
                );
        
        // Change the C2ke_ accordingly
        if (coeff.name() == "C1ke") 
        {
              this->coeffDict_.remove(C2ke_.name());
              C2ke_ = dimensioned<scalar>::lookupOrAddToDict
                      (
                          "C2ke",
                          this->coeffDict_,
                          C1ke_.value()/sqr(kappa_.value())+sigmake_.value()
                      );
        }

        // Change the C2kw_ accordingly
        if (coeff.name() == "C1kw") 
        {
              this->coeffDict_.remove(C2kw_.name());
              C2kw_ = dimensioned<scalar>::lookupOrAddToDict
                      (
                          "C2kw",
                          this->coeffDict_,
                          C1kw_.value()/sqr(kappa_.value())+sigmakw_.value()
                      );
        }
    }
    else
    {
        FatalError
            << "Trying to assign an unknown coefficient " << coeff.name() << "."
            << exit(FatalError);
    }
}

template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalBase<TurbulenceModel, BasicTurbulenceModel>::blend
(
    const volScalarField& Switch,
    const dimensionedScalar& psi1,
    const dimensionedScalar& psi2
) const
{
    return Switch*(psi1-psi2) + psi2;
}

template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalBase<TurbulenceModel, BasicTurbulenceModel>::sigmaR
(
    const volScalarField& Switch
) const
{
    return blend(Switch, sigmakw_, sigmake_);
}

template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalBase<TurbulenceModel, BasicTurbulenceModel>::C1
(
    const volScalarField& Switch
) const
{
    return blend(Switch, C1kw_, C1ke_);
}

template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalBase<TurbulenceModel, BasicTurbulenceModel>::chi() const
{
    return Rnu_/this->nu();
}

template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalBase<TurbulenceModel, BasicTurbulenceModel>::fmu
(
    const volScalarField& chi
) const
{
    const volScalarField chi3(pow3(chi));
    return chi3/(chi3 + pow3(Cw_));
}

template<class TurbulenceModel, class BasicTurbulenceModel>
void WrayAgarwalBase<TurbulenceModel, BasicTurbulenceModel>::correctNut()
{
    this->nut_ = Rnu_*fmu(this->chi());
    this->nut_.correctBoundaryConditions();

    BasicTurbulenceModel::correctNut();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TurbulenceModel, class BasicTurbulenceModel>
WrayAgarwalBase<TurbulenceModel, BasicTurbulenceModel>::WrayAgarwalBase
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    TurbulenceModel
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

    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),
    
    Cw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw",
            this->coeffDict_,
            8.54
        )
    ),

    C1ke_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1ke",
            this->coeffDict_,
            0.1127
        )
    ),

    C1kw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1kw",
            this->coeffDict_,
            0.0829
        )
    ),

    sigmake_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmake",
            this->coeffDict_,
            1.0
        )
    ),

    sigmakw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmakw",
            this->coeffDict_,
            0.72
        )
    ),

    C2ke_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2ke",
            this->coeffDict_,
            C1ke_.value()/sqr(kappa_.value())+sigmake_.value()
        )
    ),

    C2kw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2kw",
            this->coeffDict_,
            C1kw_.value()/sqr(kappa_.value())+sigmakw_.value()
        )
    ),

    Rnu_
    (
        IOobject
        (
            "Rnu",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    f1_
    (
        IOobject
        (
            "f1",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0.0", dimensionSet(0, 0, 0, 0, 0), 0.0)
    ),

    S_
    (
        IOobject
        (
            "S",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0.0", dimensionSet(0, 0, -1, 0, 0), 0.0)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TurbulenceModel, class BasicTurbulenceModel>
bool WrayAgarwalBase<TurbulenceModel, BasicTurbulenceModel>::read()
{
    if (TurbulenceModel::read())
    {        
        kappa_.readIfPresent(this->coeffDict());
        Cw_.readIfPresent(this->coeffDict());
        C1ke_.readIfPresent(this->coeffDict());
        C1kw_.readIfPresent(this->coeffDict());
        sigmake_.readIfPresent(this->coeffDict());
        sigmakw_.readIfPresent(this->coeffDict());
        C2ke_ = C1ke_ / sqr(kappa_) + sigmake_;
        C2kw_ = C1kw_ / sqr(kappa_) + sigmakw_;

        return true;
    }
    else
    {
        return false;
    }
}


template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalBase<TurbulenceModel, BasicTurbulenceModel>::DRnuEff(volScalarField Switch) const
{
    return tmp<volScalarField>
    (
        new volScalarField("DRnuEff", Rnu_*sigmaR(Switch) + this->nu())
    );
}


template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalBase<TurbulenceModel, BasicTurbulenceModel>::k() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "k",
                this->runTime_.timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensionedScalar("0", dimensionSet(0, 2, -2, 0, 0), 0)
        )
    );
}


template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalBase<TurbulenceModel, BasicTurbulenceModel>::epsilon() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "epsilon",
                this->runTime_.timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensionedScalar("0", dimensionSet(0, 2, -3, 0, 0), 0)
        )
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
