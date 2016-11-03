{
	// Euler explicit
//	phi.boundaryField()[acousticPatchID].setDisplacement();  
	solve(fvm::d2dt2(phi) == c02*fvc::laplacian(phi));
	phi.boundaryField().updateCoeffs();
	p = fvc::ddt(phi) * rhoAcoustic;
//#include "writeAverages.H"

}
