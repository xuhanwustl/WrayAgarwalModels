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
    Foam::LESModels::WA2017IDDES

Group
    grpDESTurbulence

Description
    Improved Delayed-Detached Eddy Simulation (IDDES) version of the WA2017DES model
    The RANS part is the WA-2017 version of the Wray-Agarwal model on the 
    NASA Turbulence Modeling Resource (TMR) website

    Reference:
    \verbatim
        https://turbmodels.larc.nasa.gov/wray_agarwal.html
        Han, X., and Agarwal, R. K., “Application of a New IDDES Model Based 
        on Wray-Agarwal Turbulence Model for Simulation of Wall-Bounded Flows 
        with Separation,” AIAA Scitech Meeting, Kissimmee, FL, January 2018.
    \endverbatim
    
    The default model coefficients are
        WA2017IDDESCoeffs
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
            Cd1         4.0
            fwStar      0.424
            cl          3.55
            ct          1.63
        } 
        
        delta               maxDeltaxyz;
        maxDeltaxyzCoeffs
        {
            deltaCoeff      2;
        }


SourceFiles
    WA2017IDDES.C

\*---------------------------------------------------------------------------*/

#ifndef WA2017IDDES_H
#define WA2017IDDES_H

#include "WA2017DES.H"
#include "IDDESDelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                        Class WA2017IDDES Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class WA2017IDDES
:
    public WA2017DDES<BasicTurbulenceModel>
{
    // Private Member Functions
        

protected:

    // Protected data

        // Model constants
            
            dimensionedScalar fwStar_;
            dimensionedScalar cl_;
            dimensionedScalar ct_;

        // Fields
            
            volScalarField fdtilda_;
            volScalarField fe_;
            const IDDESDelta& IDDESDelta_;
            

    // Protected Member Functions

        tmp<volScalarField> alpha() const;
        
        tmp<volScalarField> ft(const volScalarField& magGradU) const;
        
        tmp<volScalarField> fl(const volScalarField& magGradU) const;
        
        virtual void calc_fdes();


public:

    using typename BasicTurbulenceModel::alphaField;
    using typename BasicTurbulenceModel::rhoField;
    using typename BasicTurbulenceModel::transportModel;
    
    using WrayAgarwalLESModel<BasicTurbulenceModel>::Rnu_;
    using WrayAgarwalLESModel<BasicTurbulenceModel>::S_;
    
    using WA2017DES<BasicTurbulenceModel>::CDES_;
    using WA2017DES<BasicTurbulenceModel>::fdes_;
    
    using WA2017DDES<BasicTurbulenceModel>::fd_;


    // Runtime type information
    TypeName("WA2017IDDES");


    // Constructors

        // Construct from components
        WA2017IDDES
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
        WA2017IDDES(const WA2017IDDES&) = delete;
        WA2017IDDES& operator=(const WA2017IDDES&) = delete;


    // Destructor
    virtual ~WA2017IDDES()
    {}


    // Member Functions

        // Re-read model coefficients if they have changed
        virtual bool read();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "WA2017IDDES.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
