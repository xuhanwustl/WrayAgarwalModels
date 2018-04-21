#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"
#include "makeTurbulenceModel.H"

namespace Foam
{

    typedef ThermalDiffusivity<CompressibleTurbulenceModel<fluidThermo> >
        fluidThermoCompressibleTurbulenceModel;

    typedef RASModel<EddyDiffusivity<fluidThermoCompressibleTurbulenceModel> >
        RASfluidThermoCompressibleTurbulenceModel;
        
    typedef LESModel<EddyDiffusivity<fluidThermoCompressibleTurbulenceModel> >
        LESfluidThermoCompressibleTurbulenceModel;

}

///////////////////     RANS models      ///////////////////

#include "WrayAgarwal2017.H"
makeTemplatedTurbulenceModel
(
    fluidThermoCompressibleTurbulenceModel,
    RAS,
    WrayAgarwal2017
);

#include "WrayAgarwal2017m.H"
makeTemplatedTurbulenceModel
(
    fluidThermoCompressibleTurbulenceModel,
    RAS,
    WrayAgarwal2017m
);

#include "WrayAgarwal2018.H"
makeTemplatedTurbulenceModel
(
    fluidThermoCompressibleTurbulenceModel,
    RAS,
    WrayAgarwal2018
);

#include "WrayAgarwal2018EB.H"
makeTemplatedTurbulenceModel
(
    fluidThermoCompressibleTurbulenceModel,
    RAS,
    WrayAgarwal2018EB
);

///////////////////     LES models      ///////////////////

#include "WA2017DES.H"
makeTemplatedTurbulenceModel
(
    fluidThermoCompressibleTurbulenceModel,
    LES,
    WA2017DES
);

#include "WA2017DDES.H"
makeTemplatedTurbulenceModel
(
    fluidThermoCompressibleTurbulenceModel,
    LES,
    WA2017DDES
);

#include "WA2017IDDES.H"
makeTemplatedTurbulenceModel
(
    fluidThermoCompressibleTurbulenceModel,
    LES,
    WA2017IDDES
);
