/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Application
    twoSpeciesIcoFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluid with an additional Darcy term
    and hindered transport of 2 platelet species.
    Eqns:
    u_t + div(u grad u) = - grad p + nu laplacian(u) - nu alpha(theta_B) u
    Pm_t = - div( W(theta_T) (uPm - Dp grad(Pm)) ) + Rm(Pm,Pb)
    Pb_t = -Rm(Pm,Pb)
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    // Diffusion rate for virtual substance eta
    dimensionedScalar pseudoDt ("pseudoDt", dimTime, 1.0);
    dimensionedScalar Dn ( "Dn", dimViscosity, 0.0 );
    Dn = pow(L_eta, 2.0) / pseudoDt / 4.0;
    
    // Calculate g0 for the binding affinity function g(eta)
    dimensionedScalar g0
    (
        "g0",
        dimless,
        0.0
    );
    g0 = ( pow(eta_ast,3.0) + pow(1.0 - eta_t,3.0) ) / pow(1.0-eta_t, 3.0);

    // Calculate Theta_T, Theta_B based on initial conditions
    Theta_T = (Pb + Pm)/Pmax;
    Theta_B = Pb/Pmax;
    volScalarField Theta_T0 (Theta_T);
    volScalarField Theta_B0 (Theta_B);              
    
    // Initialize dt0 for linear extrapolation approx of Theta's
    scalar dt0 = runTime.controlDict().lookupOrDefault("deltaT",1.0);


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
	
	    const bool adjustTimeStep =
	        runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

	    scalar maxCo =
	        runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);

	    scalar maxDeltaT =
	        runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);	
        
        #include "CourantNo.H"
        #include "setDeltaT.H"

       // Calculate approximate Theta's for current time step
       volScalarField apprxTheta_T 
       (   
            Theta_T0 + (runTime.deltaTValue()+dt0) / dt0
            * (Theta_T - Theta_T0)
       );       
           
       volScalarField apprxTheta_B
       (   
            Theta_B0 + (runTime.deltaTValue()+dt0) / dt0
            * (Theta_B - Theta_B0)
       );       
 
        // Calculate the frictional resistance (Darcy Term)
        volScalarField alpha
        (
            nu*C_CK*sqr(0.6*apprxTheta_B) / (pow(1.-0.6*apprxTheta_B,3.))
        );

        // Momentum predictor
        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi,U)
          - fvm::laplacian(nu,U)
          + fvm::Sp(alpha,U)
        );

        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }

        // --- PISO loop
        while (piso.correct())
        {
            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            adjustPhi(phiHbyA, U, p);

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAU);

            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                // Pressure corrector

                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);

                pEqn.solve();

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }

        // Solve the Platelet Equations
        #include "plateletEqns.H"

        // Update the platelet fractions
        volScalarField Theta_T0 (Theta_T);
        volScalarField Theta_B0 (Theta_B);
        Theta_T = (Pm + Pb)/Pmax;
        Theta_B = Pb/Pmax;
        dt0 = runTime.deltaTValue();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
