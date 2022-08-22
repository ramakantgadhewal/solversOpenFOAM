{
    // Virtual substance eta
    eta = max(0*Theta_B, min(1.+0*Theta_B,Theta_B));
    solve(
        fvm::Sp(1.0/pseudoDt,eta)
        - fvc::Sp(1.0/pseudoDt,eta)
        -fvm::laplacian(Dn,eta)
    );
    
    // Calculate the binding affinity function g(eta)
    g_eta = 
        g0 * 
	    max
	    (
	       0*eta, 
	       min
	       (
	          1. + 0*eta, 
	          pow(eta-eta_t,3.0) / 
              ( pow(eta_ast,3.0) + pow(eta - eta_t,3.0) )
	       )
	    );

    // Interpolate Theta_T to the faces
    surfaceScalarField Theta_Tf
    (
	    fvc::interpolate(apprxTheta_T,"interpolate(Theta_T)")
    );

    // Calculate the hindered transport function on the faces
    surfaceScalarField W
    ( 
	    max
	    (
	        0*Theta_Tf,
	        min
	        ( 
	           1. + 0*Theta_Tf, 
	           tanh( constant::mathematical::pi * (1. - Theta_Tf))
	        )
	    )
    );

    // Advective Flux Field
    surfaceScalarField phiPlts( phi*W );

    // Diffusive Flux Field
    surfaceScalarField DpW( Dp*W );

    // Build the Matrix
    fvScalarMatrix PmEqn
    (
        fvm::ddt(Pm)
        + fvm::div(phiPlts,Pm,"div(phiPlts,Pm)")
        - fvm::laplacian(DpW,Pm)
        + fvm::Sp(Kadh*max(0*Pb,min(Pmax-Pb,Pmax+0*Pb)),Pm) 
        + fvm::Sp(KcohPmax*g_eta,Pm)
    );

    PmEqn.relax();
    PmEqn.solve();
    Pm.correctBoundaryConditions();

    fvScalarMatrix PbEqn
    (
        fvm::ddt(Pb)
        ==
        Kadh*max(0*Pb,min(Pmax-Pb,Pmax+0*Pb))*Pm
        + KcohPmax*g_eta*Pm
    );

    PbEqn.solve();
    Pb.correctBoundaryConditions();
}

