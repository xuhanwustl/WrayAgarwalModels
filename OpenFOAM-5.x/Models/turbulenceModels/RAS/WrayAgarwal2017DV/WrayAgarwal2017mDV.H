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
    Foam::RASModels::WrayAgarwal2017mDV

Group
    grpRASTurbulence

Description
    Wray-Agarwal one equation Turbulence Model
    Enhanced version of the WA-2017 model on NASA Turbulence Modeling 
    Resource (TMR) website
    
    Coupled a trasport equation of density fluctuation to account for 
    extra compressibility effects on turbulence for hypersonic flow.
    
    Reference:
    \verbatim
        https://turbmodels.larc.nasa.gov/wray_agarwal.html
        Han, X., and Agarwal, R. K., “A Density Variance Correction Based on 
        Wray-Agarwal Turbulence Model for Simulation of Hypersonic Flow,” 
        AIAA Aviation Meeting, Atlanta, GA, June 2018.
    \endverbatim
    
    The default model coefficients are
        WrayAgarwal2017mDVCoeffs
        {
            kappa       0.41
            Cw          8.54
            C1ke        0.1127
            C1kw        0.0829
            sigmake     1.0
            sigmakw     0.72
            C2ke        1.6704  (C1ke/sqr(kappa) + sigmake)
            C2kw        1.2132  (C1kw/sqr(kappa) + sigmakw)
            Cm          8.0
            gamma       1.4
            Rsp         286.9
            Cr1         0.01
            Cr2         0.01
            Crho1       1.4
            Crho2       2.0
            sigmaRho    1.0
        }
    
SourceFiles
    WrayAgarwal2017mDV.C

\*---------------------------------------------------------------------------*/

#ifndef WrayAgarwal2017mDV_H
#define WrayAgarwal2017mDV_H

#include "WrayAgarwal2017m.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                      Class WrayAgarwal2017mDV Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class WrayAgarwal2017mDV
:
    public WrayAgarwal2017m<BasicTurbulenceModel>
{
    // Private Member Functions


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar gamma_;
            dimensionedScalar Rsp_;
            dimensionedScalar Cr1_;
            dimensionedScalar Cr2_;
            dimensionedScalar Crho1_;
            dimensionedScalar Crho2_;
            dimensionedScalar sigmaRho_;
            
        // Fields
            volScalarField T_;
            volScalarField rhoPrime2_;
            volScalarField fcomp_;
            volScalarField gradUSum_;


public:

    using typename BasicTurbulenceModel::transportModel;
    using typename BasicTurbulenceModel::rhoField;
    using typename BasicTurbulenceModel::alphaField;

    using WrayAgarwal2017m<BasicTurbulenceModel>::Cm_;
    using WrayAgarwalRASModel<BasicTurbulenceModel>::C2kw_;
    using WrayAgarwalRASModel<BasicTurbulenceModel>::C2ke_;
    using WrayAgarwalRASModel<BasicTurbulenceModel>::Rnu_;
    using WrayAgarwalRASModel<BasicTurbulenceModel>::f1_;
    using WrayAgarwalRASModel<BasicTurbulenceModel>::S_;
    
    
    //- Runtime type information
    TypeName("WrayAgarwal2017mDV");


    // Constructors

        //- Construct from components
        WrayAgarwal2017mDV
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );

        // Disallow default bitwise copy construct and assignment
        WrayAgarwal2017mDV(const WrayAgarwal2017mDV&) = delete;
        WrayAgarwal2017mDV& operator=(const WrayAgarwal2017mDV&) = delete;


    //- Destructor
    virtual ~WrayAgarwal2017mDV()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for rhoPrime2_
        tmp<volScalarField> DRhoPrime2Eff(volScalarField Switch) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "WrayAgarwal2017mDV.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
