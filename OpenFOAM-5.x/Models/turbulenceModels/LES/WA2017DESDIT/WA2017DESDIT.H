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
    Foam::LESModels::WA2017DESDIT

Group
    grpDESTurbulence

Description
    Wray-Agarwal model based DES type turbulence model 
    The RANS part is the WA-2017 version of the Wray-Agarwal model on the 
    NASA Turbulence Modeling Resource (TMR) website

    This model is the same as the WA2017DES model except the LES mode is 
    enforced to perform the Decaying Isotropic Turbulence (DIT) test
    
    Reference:
    \verbatim
        https://turbmodels.larc.nasa.gov/wray_agarwal.html
        Han, X., Wray, T. J., Agarwal, R. K., "Application of a New DES Model 
        Based on Wray-Agarwal Turbulence Model for Simulation of Wall-Bounded 
        Flows with Separation," AIAA Paper 2017-3966, June 2017.
    \endverbatim
    
    The default model coefficients are
        WA2017DESDITCoeffs
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
    WA2017DESDIT.C

\*---------------------------------------------------------------------------*/

#ifndef WA2017DESDIT_H
#define WA2017DESDIT_H

#include "WA2017DES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                        Class WA2017DESDIT Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class WA2017DESDIT
:
    public WA2017DES<BasicTurbulenceModel>
{
    // Private Member Functions


protected:

    // Protected data
            

    // Protected Member Functions
        
        virtual void calc_fdes();


public:

    using typename BasicTurbulenceModel::alphaField;
    using typename BasicTurbulenceModel::rhoField;
    using typename BasicTurbulenceModel::transportModel;

    using WA2017DES<BasicTurbulenceModel>::CDES_;
    using WA2017DES<BasicTurbulenceModel>::fdes_;
    using WA2017DES<BasicTurbulenceModel>::Rnu_;
    using WA2017DES<BasicTurbulenceModel>::S_;
    

    //- Runtime type information
    TypeName("WA2017DESDIT");


    // Constructors

        // Construct from components
        WA2017DESDIT
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
        WA2017DESDIT(const WA2017DESDIT&) = delete;
        WA2017DESDIT& operator=(const WA2017DESDIT&) = delete;


    // Destructor
    virtual ~WA2017DESDIT()
    {}


    // Member Functions
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "WA2017DESDIT.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
