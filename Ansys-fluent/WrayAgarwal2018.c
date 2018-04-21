/* UDF file of Wray-Agarwal Turbulence Model for ANSYS fluent

Description
    Wray-Agarwal one equation Turbulence Model
    A wall distance free version of the WA-2017 model on NASA Turbulence 
    Modeling Resource (TMR) website
    
    Re-defined the f1 function in the WA-2017 model without using the wall 
    distance
    
    Added a bound on the second destruction term in the Rnu transport 
    equation to fix the significant eddy viscosity drop in the zero strain 
    rate region (e.g. channel center)

    The default model coefficients are
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

Reference:
    https://turbmodels.larc.nasa.gov/wray_agarwal.html
    Han, X., Rahman, M. M., and Agarwal, R. K., “Development and 
    Application of a Wall Distance Free Wray-Agarwal Turbulence Model,”
    AIAA Scitech Meeting, Kissimmee, FL, January 2018.
    
==========================================================================*/

#include "udf.h"
#include "mem.h"
#include "math.h"

// Coefficients
#define Kappa       0.41
#define Cw          8.54
#define C1ke        0.1284
#define C1kw        0.0829
#define Sigmake     1.0
#define Sigmakw     0.72
#define C2ke        (C1ke/Kappa/Kappa+Sigmake)
#define C2kw        (C1kw/Kappa/Kappa+Sigmakw)
#define Cm          8.0
#define Cmu         0.09
#define WA_SMALL    1e-15

// Fields
enum
{
    R,
    SRM,
    W,
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

#define W_              C_UDMI(c,t,W)
#define d_              C_WALL_DIST(c,t)
#define nu_             C_UDMI(c,t,nu)
#define f1_             C_UDMI(c,t,f1)

#define WA_SMALL        1e-15

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

void calculate_SW(cell_t c, Thread *t)
{
    int i, j;
    real GradU[ND_ND][ND_ND];
    real Symm[ND_ND][ND_ND];
    real Skew[ND_ND][ND_ND];
    real S2, W2;

    // Initialize velocity gradient tensor
    #if RP_3D
    GradU[0][0] = C_DUDX(c,t); GradU[0][1] = C_DUDY(c,t); GradU[0][2] = C_DUDZ(c,t);
    GradU[1][0] = C_DVDX(c,t); GradU[1][1] = C_DVDY(c,t); GradU[1][2] = C_DVDZ(c,t);
    GradU[2][0] = C_DWDX(c,t); GradU[2][1] = C_DWDY(c,t); GradU[2][2] = C_DWDZ(c,t);
    #else
    GradU[0][0] = C_DUDX(c,t); GradU[0][1] = C_DUDY(c,t);
    GradU[1][0] = C_DVDX(c,t); GradU[1][1] = C_DVDY(c,t);                  
    #endif

    // Calculate symmetrical tensor and antisymmetrical tensor 
    for (i=0; i<ND_ND; ++i)
    {
        for(j=0; j<ND_ND; ++j)
        {
            Symm[i][j] = 1.0/2.0*(GradU[i][j]+GradU[j][i]);
            Skew[i][j] = 1.0/2.0*(GradU[i][j]-GradU[j][i]);
        }
    }

    // Calculate strain rate magnitude and vorticity magnitude
    S2 = 0;
    W2 = 0;
    for (i=0; i<ND_ND; ++i)
    {
        for(j=0; j<ND_ND; ++j)
        {
            S2 += Symm[i][j]*Symm[i][j];
            W2 += Skew[i][j]*Skew[i][j];
        }
    }
    S_ = sqrt(2.0*S2);
    W_ = sqrt(2.0*W2);
};

real calculate_f1(cell_t c, Thread *t)
{
    // Wall Distance Free
    real WDF_R = W_/S_;
    real WDF_omega = S_/sqrt(Cmu);
    real WDF_k = mut_/rho_*WDF_omega;
    
    real eta = S_*MAX(1.0, WDF_R);
    
    real arg1 = (nu_+R_)/2.0 * SQR(eta) / MAX(Cmu*WDF_k*WDF_omega, WA_SMALL);
    
    return tanh(pow(arg1, 4.0));
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

    Set_User_Memory_Name(W, "W");
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

            // Calculate strain rate magnitude S_ and vorticity magnitude W_
            calculate_SW(c, t);

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
        -(1.0-f1_)*MIN(C2ke*SQR(R_)/S2*NV_MAG2(GradS_),
                     Cm*NV_MAG2(GradR_))
        );

    // The second destruction term is treated as explicit
    dS[eqn] = 0.0;

    return source;
}
