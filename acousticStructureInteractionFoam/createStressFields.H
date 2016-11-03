    Info<< "Reading field U\n" << endl;
    volVectorField Usolid
    (
        IOobject
        (
            "Usolid",
            runTime.timeName(),
            stressMesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        stressMesh
    );

autoPtr<volScalarField> Tptr(NULL);

if (thermalStress)
{
    Info<< "Reading field T\n" << endl;
    Tptr.reset
    (
        new volScalarField
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
}


Info<< "Calculating stress field sigmaU\n" << endl;
volSymmTensorField sigmaU
(
    IOobject
    (
        "sigmaU",
        runTime.timeName(),
        stressMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mu*twoSymm(fvc::grad(Usolid)) + lambda*(I*tr(fvc::grad(Usolid)))
);

Info<< "Calculating explicit part of div(sigma) divSigmaExp\n" << endl;
volVectorField divSigmaExp
(
    IOobject
    (
        "divSigmaExp",
        runTime.timeName(),
        stressMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(sigmaU)
);

if (compactNormalStress)
{
    divSigmaExp -= fvc::laplacian(2*mu + lambda, Usolid, "laplacian(DUsolid,Usolid)");
}
else
{
    divSigmaExp -= fvc::div((2*mu + lambda)*fvc::grad(Usolid), "div(sigmaU)");
}
