#include "IncompressibleTurbulenceModel.H"
#include "incompressible/transportModel/transportModel.H"
#include "addToRunTimeSelectionTable.H"
#include "makeTurbulenceModel.H"

#include "RASModel.H"
#include "LESModel.H"

namespace Foam
{
    typedef IncompressibleTurbulenceModel<transportModel> 
        transportModelIncompressibleTurbulenceModel;
        
    typedef RASModel<transportModelIncompressibleTurbulenceModel>
        RAStransportModelIncompressibleTurbulenceModel;
        
    typedef LESModel<transportModelIncompressibleTurbulenceModel>
        LEStransportModelIncompressibleTurbulenceModel;
}

///////////////////     RANS models      ///////////////////

#include "WrayAgarwal2017.H"
makeTemplatedTurbulenceModel
(
    transportModelIncompressibleTurbulenceModel,
    RAS,
    WrayAgarwal2017
);

#include "WrayAgarwal2017m.H"
makeTemplatedTurbulenceModel
(
    transportModelIncompressibleTurbulenceModel,
    RAS,
    WrayAgarwal2017m
);

#include "WrayAgarwal2018.H"
makeTemplatedTurbulenceModel
(
    transportModelIncompressibleTurbulenceModel,
    RAS,
    WrayAgarwal2018
);

#include "WrayAgarwal2018EB.H"
makeTemplatedTurbulenceModel
(
    transportModelIncompressibleTurbulenceModel,
    RAS,
    WrayAgarwal2018EB
);

///////////////////     LES models      ///////////////////

#include "WA2017DES.H"
makeTemplatedTurbulenceModel
(
    transportModelIncompressibleTurbulenceModel,
    LES,
    WA2017DES
);

#include "WA2017DESDIT.H"
makeTemplatedTurbulenceModel
(
    transportModelIncompressibleTurbulenceModel,
    LES,
    WA2017DESDIT
);

#include "WA2017DDES.H"
makeTemplatedTurbulenceModel
(
    transportModelIncompressibleTurbulenceModel,
    LES,
    WA2017DDES
);

#include "WA2017IDDES.H"
makeTemplatedTurbulenceModel
(
    transportModelIncompressibleTurbulenceModel,
    LES,
    WA2017IDDES
);

