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
/*
#include "WD2017.H"
makeTemplatedTurbulenceModel
(
    fluidThermoCompressibleTurbulenceModel,
    RAS,
    WD2017
);
*/
///////////////////     LES models      ///////////////////
/*
#include "WA2017IDDES.H"
makeTemplatedTurbulenceModel
(
    fluidThermoCompressibleTurbulenceModel,
    LES,
    WA2017IDDES
);
*/
