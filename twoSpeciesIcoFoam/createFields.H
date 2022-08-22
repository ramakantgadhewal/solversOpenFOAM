Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Dynamic Viscosity
dimensionedScalar nu
(
    "nu",
    dimViscosity,
    transportProperties.lookup("nu")
);

// Carman-Kozeny Consant
dimensionedScalar C_CK
(
    "C_CK",
    dimless/dimArea,
    transportProperties.lookup("C_CK")
);

// Maximum Packing Limit of Platelets
dimensionedScalar Pmax
(
    "Pmax",
    dimless/dimVolume,
    transportProperties.lookup("Pmax")
);

// Diameter of a Platelet
dimensionedScalar P_diam
(
    "P_diam",
    dimLength,
    transportProperties.lookup("P_diam")
);

// Platelet Diffusion Constant
dimensionedScalar Dp
(
    "Dp",
    dimViscosity,
    transportProperties.lookup("Dp")
);

// Length of eta diffusion L_eta
dimensionedScalar L_eta
(
    "L_eta",
    dimLength,
    transportProperties.lookup("L_eta")
);

// Virtual substance Eta threshold
dimensionedScalar eta_t
(
    "eta_t",
    dimless,
    transportProperties.lookup("eta_t")
);

// Virtual substance Eta_ast
dimensionedScalar eta_ast
(
    "eta_ast",
    dimless,
    transportProperties.lookup("eta_ast")
);

// Cohesion binding rate
dimensionedScalar KcohPmax
(
    "KcohPmax",
    dimless/dimTime,
    transportProperties.lookup("KcohPmax")
);

// Pressure Field
Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Velocity Field
Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Mobile Platelet Field
Info<< "Reading field Pm\n" << endl;
volScalarField Pm
(
    IOobject
    (
        "Pm",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Bound Platelet Field
Info<< "Reading field Pb\n" << endl;
volScalarField Pb
(
    IOobject
    (
        "Pb",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Total Platelet Fraction
Info<< "Reading field Theta_T\n" << endl;
volScalarField Theta_T
(
    IOobject
    (
        "Theta_T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Bound Platelet Fraction
Info<< "Reading field Theta_B\n" << endl;
volScalarField Theta_B
(
    IOobject
    (
        "Theta_B",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Virtual substance secreted from Bound Platelets
Info<< "Reading field eta\n" << endl;
volScalarField eta
(
    IOobject
    (
        "eta",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Virtual substance secreted from Bound Platelets
Info<< "Reading field g_eta\n" << endl;
volScalarField g_eta
(
    IOobject
    (
        "g_eta",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// k_adh Field
Info<< "Reading field Kadh\n" << endl;
volScalarField Kadh
(
    IOobject
    (
        "Kadh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());
