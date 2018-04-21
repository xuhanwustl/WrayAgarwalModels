/* UDF file of Wray-Agarwal Turbulence Model for ANSYS fluent

Description:
    Wray-Agarwal one equation Turbulence Model
    WA-2017 version on NASA Turbulence Modeling Resource (TMR) website

    The default model coefficients are
        kappa       0.41
        Cw          8.54
        C1ke        0.1127
        C1kw        0.0829
        sigmake     1.0
        sigmakw     0.72
        C2ke        1.6704  (C1ke/sqr(kappa) + sigmake)
        C2kw        1.2132  (C1kw/sqr(kappa) + sigmakw)

Reference:
    https://turbmodels.larc.nasa.gov/wray_agarwal.html
    Han, X., Wray, T. J., Agarwal, R. K., "Application of a New DES Model 
    Based on Wray-Agarwal Turbulence Model for Simulation of Wall-Bounded 
    Flows with Separation," AIAA Paper 2017-3966, June 2017. 
    
==========================================================================*/

#include "udf.h"
#include "mem.h"
#include "math.h"

// Coefficients
#define Kappa       0.41
#define Cw          8.54
#define C1ke        0.1127
#define C1kw        0.0829
#define Sigmake     1.0
#define Sigmakw     0.72
#define C2ke        (C1ke/Kappa/Kappa+Sigmake)
#define C2kw        (C1kw/Kappa/Kappa+Sigmakw)
#define WA_SMALL    1e-15

// Fields
enum
{
    R,
    SRM,
    d,
    nu,
    f1
};

#define mut_            C_MU_T(c,t)
#define mu_             C_MU_L(c,t)
#define rho_            C_R(c,t)

#define R_              C_UDSI(c,t,R)
#define S_              C_UDSI(c,t,SRM)
#define GradR_          C_UDSI_G(c,t,R)
#define GradS_          C_UDSI_G(c,t,SRM)

#define d_              C_WALL_DIST(c,t)
#define nu_             C_UDMI(c,t,nu)
#define f1_             C_UDMI(c,t,f1)

// Healper functions
real chi(cell_t c, Thread *t)
{
    return R_/nu_;
}

real fmu(cell_t c, Thread *t)
{
    real chi3 = pow(chi(c,t),3.0);
    return chi3 / (chi3 + pow(Cw,3.0));
}

real calculate_f1(cell_t c, Thread *t)
{
    real sqrtRS = sqrt(R_*S_);
    real arg1 = (1.0+d_*sqrtRS/nu_)
                /
                (1.0+SQR(MAX(d_*sqrtRS,1.5*R_)/20.0/nu_));

    return MIN(tanh(pow(arg1, 4.0)), 0.9);
}

real blend(cell_t c, Thread *t, real Switch, real psi1, real psi2)
{
    return Switch*(psi1-psi2) + psi2;
}

real C1(cell_t c, Thread *t)
{
    return blend(c, t, f1_, C1kw, C1ke);
}

real sigmaR(cell_t c, Thread *t)
{
    return blend(c, t, f1_, Sigmakw, Sigmake);
}

// Initialization functions
DEFINE_ON_DEMAND(WA_setnames) 
{
    Set_User_Scalar_Name(R, "R");
    Set_User_Scalar_Name(SRM, "S");

    Set_User_Memory_Name(d, "d");
    Set_User_Memory_Name(nu, "nu");
    Set_User_Memory_Name(f1, "f1");
}

// Iteration functions
DEFINE_ADJUST(WA_adjust, d)
{
    Thread *t;
    cell_t c;
    
    int n = 0;

    thread_loop_c(t, d)
    {
        begin_c_loop(c, t)
        {
            nu_ = mu_/rho_;

            // Calculate strain rate magnitude S_
            S_ = C_STRAIN_RATE_MAG(c,t);

            // Bound R_ and S_
            R_ = MAX(R_, 0.0);
            S_ = MAX(S_, WA_SMALL);

            // Calculate and bound the switch function f1
            f1_ = calculate_f1(c,t);
        }
        end_c_loop(c, t)
    }

    // Calculate derivatives
    // Using recontruct gradient in R-Equation
    for (n=0; n<2; ++n)
    {
        MD_Alloc_Storage_Vars(d, SV_UDSI_RG(n), SV_UDSI_G(n), SV_NULL);
        Scalar_Reconstruction(d, SV_UDS_I(n), -1, SV_UDSI_RG(n), NULL);
        Scalar_Derivatives(d, SV_UDS_I(n), -1, SV_UDSI_G(n), SV_UDSI_RG(n), NULL);
    }
}

// Turbulent viscosity
DEFINE_TURBULENT_VISCOSITY(WA_mut, c, t)
{
    return rho_*fmu(c,t)*R_;
}

// Transport terms
DEFINE_DIFFUSIVITY(WA_diffusivity, c, t, eqn)
{
    return rho_*sigmaR(c, t)*R_ + mu_;
}

DEFINE_SOURCE(WA_source_prod, c, t, dS, eqn)
{
    real source;

    source = rho_ * C1(c, t) * R_ * S_;
    dS[eqn] = rho_ * C1(c, t) * S_;

    return source;
}

DEFINE_SOURCE(WA_source_dest1, c, t, dS, eqn)
{
    real source;

    source = rho_*(
        f1_*C2kw*R_/S_*NV_DOT(GradR_,GradS_)
        );

    dS[eqn] = rho_*(
        f1_*C2kw/S_*NV_DOT(GradR_,GradS_)
        );

    return source;
}

DEFINE_SOURCE(WA_source_dest2, c, t, dS, eqn)
{
    real source;
    real S2 = MAX(SQR(S_), WA_SMALL);

    source = rho_*(
        -(1.0-f1_)*C2ke*SQR(R_)/S2*NV_MAG2(GradS_)
        );

    dS[eqn] = rho_*(
        -(1.0-f1_)*C2ke*R_/S2*NV_MAG2(GradS_)
        );

    return source;
}
