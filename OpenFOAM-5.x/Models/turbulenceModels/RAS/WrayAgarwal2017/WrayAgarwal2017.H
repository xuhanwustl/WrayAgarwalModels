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
    Foam::RASModels::WrayAgarwal2017

Group
    grpRASTurbulence

Description
    Wray-Agarwal one equation Turbulence Model 
    WA-2017 version on NASA Turbulence Modeling Resource (TMR) website
    
    Reference:
    \verbatim
        https://turbmodels.larc.nasa.gov/wray_agarwal.html
        Han, X., Wray, T. J., Agarwal, R. K., "Application of a New DES Model 
        Based on Wray-Agarwal Turbulence Model for Simulation of Wall-Bounded 
        Flows with Separation," AIAA Paper 2017-3966, June 2017. 
    \endverbatim
    
    The default model coefficients are
        WrayAgarwal2017Coeffs
        {
            kappa       0.41
            Cw          8.54
            C1ke        0.1127
            C1kw        0.0829
            sigmake     1.0
            sigmakw     0.72
            C2ke        1.6704  (C1ke/sqr(kappa) + sigmake)
            C2kw        1.2132  (C1kw/sqr(kappa) + sigmakw)
        }   
    
SourceFiles
    WrayAgarwal2017.C

\*---------------------------------------------------------------------------*/

#ifndef WrayAgarwal2017_H
#define WrayAgarwal2017_H

#include "WrayAgarwalRASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                      Class WrayAgarwal2017 Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class WrayAgarwal2017
:
    public WrayAgarwalRASModel<BasicTurbulenceModel>
{
    // Private Member Functions
    

protected:

    // Protected data

        // Model coefficients

        // Fields

            //- Wall distance
            const volScalarField& y_;


    // Protected Member Functions
    
        virtual void calc_f1();

public:

    using typename BasicTurbulenceModel::transportModel;
    using typename BasicTurbulenceModel::rhoField;
    using typename BasicTurbulenceModel::alphaField;
    
    using WrayAgarwalRASModel<BasicTurbulenceModel>::C2kw_;
    using WrayAgarwalRASModel<BasicTurbulenceModel>::C2ke_;
    using WrayAgarwalRASModel<BasicTurbulenceModel>::Rnu_;
    using WrayAgarwalRASModel<BasicTurbulenceModel>::f1_;
    using WrayAgarwalRASModel<BasicTurbulenceModel>::S_;
    
    
    //- Runtime type information
    TypeName("WrayAgarwal2017");


    // Constructors

        //- Construct from components
        WrayAgarwal2017
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
        WrayAgarwal2017(const WrayAgarwal2017&) = delete;
        WrayAgarwal2017& operator=(const WrayAgarwal2017&) = delete;


    //- Destructor
    virtual ~WrayAgarwal2017()
    {}


    // Member Functions

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "WrayAgarwal2017.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
