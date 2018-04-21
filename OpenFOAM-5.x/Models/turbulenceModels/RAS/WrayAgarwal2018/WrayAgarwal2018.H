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
    Foam::RASModels::WrayAgarwal2018

Group
    grpRASTurbulence

Description
    Wray-Agarwal one equation Turbulence Model
    A wall distance free version of the WA-2017 model on NASA Turbulence 
    Modeling Resource (TMR) website
    
    Re-defined the f1 function in the WA-2017 model without using the wall 
    distance
    
    Added a bound on the second destruction term in the Rnu transport 
    equation to fix the significant eddy viscosity drop in the zero strain 
    rate region (e.g. channel center)
    
    Reference:
    \verbatim
        https://turbmodels.larc.nasa.gov/wray_agarwal.html
        Han, X., Rahman, M. M., and Agarwal, R. K., “Development and 
        Application of a Wall Distance Free Wray-Agarwal Turbulence Model,”
        AIAA Scitech Meeting, Kissimmee, FL, January 2018.
    \endverbatim

    The default model coefficients are
        WrayAgarwal2018Coeffs
        {
            kappa       0.41
            Cw          8.54
            C1ke        0.1284
            C1kw        0.0829
            sigmake     1.0
            sigmakw     0.72
            C2ke        1.7638  (C1ke/sqr(kappa) + sigmake)
            C2kw        1.2132  (C1kw/sqr(kappa) + sigmakw)
            Cmu         0.09
            Cm          8.0
        }
    
SourceFiles
    WrayAgarwal2018.C

\*---------------------------------------------------------------------------*/

#ifndef WrayAgarwal2018_H
#define WrayAgarwal2018_H

#include "WrayAgarwalRASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                      Class WrayAgarwal2018 Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class WrayAgarwal2018
:
    public WrayAgarwalRASModel<BasicTurbulenceModel>
{
    // Private Member Functions


protected:

    // Protected data

        // Model coefficients
        
            dimensionedScalar Cmu_;
            dimensionedScalar Cm_;

        // Fields

            volScalarField W_;


    // Protected Member Functions
    
        tmp<volScalarField> WDF_R() const;

        tmp<volScalarField> WDF_omega() const;

        tmp<volScalarField> WDF_k
        (
            const volScalarField& omega 
        ) const;

        tmp<volScalarField> arg1() const;
        
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
    TypeName("WrayAgarwal2018");


    // Constructors

        //- Construct from components
        WrayAgarwal2018
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
        WrayAgarwal2018(const WrayAgarwal2018&) = delete;
        WrayAgarwal2018& operator=(const WrayAgarwal2018&) = delete;


    //- Destructor
    virtual ~WrayAgarwal2018()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "WrayAgarwal2018.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
