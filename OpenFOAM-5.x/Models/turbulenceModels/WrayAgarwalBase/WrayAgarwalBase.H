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

Class
    Foam::RASorLESModels::WrayAgarwalBase

Group
    grpRASTurbulence

Description
    Abstract base class for Wray-Agarwal Turbulence Model Family
    Includes two alias template classes:
        WrayAgarwalRASModel for RANS models:
            WrayAgarwal2017
            WrayAgarwal2017m
            WrayAgarwal2018
            WrayAgarwal2018EB
        WrayAgarwalLESModel for LES models:
            WA2017DES
            WA2017DESDIT
            WA2017DDES
            WA2017IDDES

SourceFiles
    WrayAgarwalBase.C

\*---------------------------------------------------------------------------*/

#ifndef WrayAgarwalBase_H
#define WrayAgarwalBase_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class WrayAgarwalBase Declaration
\*---------------------------------------------------------------------------*/

template<class TurbulenceModel, class BasicTurbulenceModel>
class WrayAgarwalBase
:
    public TurbulenceModel
{
    // Private Member Functions


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar kappa_;
            dimensionedScalar Cw_;
            dimensionedScalar C1ke_;
            dimensionedScalar C1kw_;
            dimensionedScalar sigmake_;
            dimensionedScalar sigmakw_;
            dimensionedScalar C2ke_;
            dimensionedScalar C2kw_;

        // Fields

            volScalarField Rnu_;
            volScalarField f1_;
            volScalarField S_;


    // Protected Member Functions
    
        virtual void modifyCoeff
        (
            dimensionedScalar& coeff, 
            scalar value
        );
        
        virtual void calc_f1() = 0;

        tmp<volScalarField> blend
        (
            const volScalarField& Switch,
            const dimensionedScalar& psi1,
            const dimensionedScalar& psi2
        ) const;

        tmp<volScalarField> sigmaR(const volScalarField& Switch) const;

        tmp<volScalarField> C1(const volScalarField& Switch) const;

        tmp<volScalarField> chi() const;

        tmp<volScalarField> fmu(const volScalarField& chi) const;
        
        virtual void correctNut();

public:

    using typename BasicTurbulenceModel::alphaField;
    using typename BasicTurbulenceModel::rhoField;
    using typename BasicTurbulenceModel::transportModel;


    // Constructors

        //- Construct from components
        WrayAgarwalBase
        (
            const word& type,
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName
        );

        // Disallow default bitwise copy construct and assignment
        WrayAgarwalBase(const WrayAgarwalBase&) = delete;
        WrayAgarwalBase& operator=(const WrayAgarwalBase&) = delete;
     

    //- Destructor
    virtual ~WrayAgarwalBase()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for Rnu
        tmp<volScalarField> DRnuEff(volScalarField Switch) const;

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const;

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "WrayAgarwalBase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
