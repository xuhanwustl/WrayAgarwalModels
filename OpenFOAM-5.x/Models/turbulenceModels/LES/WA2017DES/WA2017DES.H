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
    Foam::LESModels::WA2017DES

Group
    grpDESTurbulence

Description
    Detached Eddy Simulation (DES) version of the Wray-Agarwal model
    The RANS part is the WA-2017 version of the Wray-Agarwal model on the 
    NASA Turbulence Modeling Resource (TMR) website

    Field blendfactor_ is used to interpolate different schemes between LES 
    and RANS region
    
    Reference:
    \verbatim
        https://turbmodels.larc.nasa.gov/wray_agarwal.html
        Han, X., Wray, T. J., Agarwal, R. K., "Application of a New DES Model 
        Based on Wray-Agarwal Turbulence Model for Simulation of Wall-Bounded 
        Flows with Separation," AIAA Paper 2017-3966, June 2017.
    \endverbatim
    
    The default model coefficients are
        WA2017DESCoeffs
        {
            kappa       0.41
            Cw          8.54
            C1ke        0.1127
            C1kw        0.0829
            sigmake     1.0
            sigmakw     0.72
            C2ke        1.6704  (C1ke/sqr(kappa) + sigmake)
            C2kw        1.2132  (C1kw/sqr(kappa) + sigmakw)
            CDES        0.41
        } 

SourceFiles
    WA2017DES.C

\*---------------------------------------------------------------------------*/

#ifndef WA2017DES_H
#define WA2017DES_H

#include "WrayAgarwalLESModel.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                        Class WA2017DES Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class WA2017DES
:
    public WrayAgarwalLESModel<BasicTurbulenceModel>
{
    // Private Member Functions


protected:

    // Protected data

        // Model constants
            
            dimensionedScalar CDES_;

        // Fields

            // Wall distance
            const volScalarField& y_;
            
            // LES fields
            volScalarField fdes_;
            volScalarField outDelta_;
            volScalarField blendfactor_;
            surfaceScalarField UBlendingFactor_;
            surfaceScalarField RnuBlendingFactor_;
            surfaceScalarField pBlendingFactor_;
            surfaceScalarField KBlendingFactor_;
            surfaceScalarField eBlendingFactor_;
            surfaceScalarField hBlendingFactor_;
            

    // Protected Member Functions
        
        virtual void calc_f1();
        
        virtual void calc_fdes();
        
        virtual tmp<volScalarField> blendFactor() const;
        
        virtual void calcBlendFactors();


public:

    using typename BasicTurbulenceModel::alphaField;
    using typename BasicTurbulenceModel::rhoField;
    using typename BasicTurbulenceModel::transportModel;
    
    using WrayAgarwalLESModel<BasicTurbulenceModel>::C2kw_;
    using WrayAgarwalLESModel<BasicTurbulenceModel>::C2ke_;
    using WrayAgarwalLESModel<BasicTurbulenceModel>::Rnu_;
    using WrayAgarwalLESModel<BasicTurbulenceModel>::f1_;
    using WrayAgarwalLESModel<BasicTurbulenceModel>::S_;


    // Runtime type information
    TypeName("WA2017DES");


    // Constructors

        // Construct from components
        WA2017DES
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
        WA2017DES(const WA2017DES&) = delete;
        WA2017DES& operator=(const WA2017DES&) = delete;

    // Destructor
    virtual ~WA2017DES()
    {}


    // Member Functions

        // Re-read model coefficients if they have changed
        virtual bool read();

        // Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        // Return the LES field indicator
        virtual tmp<volScalarField> LESRegion() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "WA2017DES.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
