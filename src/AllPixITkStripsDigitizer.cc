/**
 *  Author:
 *    Dino Tahirovic <dino.tahirovic@cern.ch>
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <benoit@lal.in2p3.fr>
 */


#include "AllPixDebug.h"
#include "AllPixITkStripsDigitizer.hh"

#include <array>
#include <climits>
#include "CLHEP/Random/RandGauss.h"

#include "AllPixTrackerHit.hh"
#include "G4DigiManager.hh"
#include "G4VDigitizerModule.hh"
#include "AllPixGeoDsc.hh"

#include "TMath.h"
#include "TF1.h"
#include "TString.h"

const double AllPixITkStripsDigitizer::kPairEnergy = 3.6*eV;
const double AllPixITkStripsDigitizer::fC = CLHEP::coulomb*1e-15;
const double AllPixITkStripsDigitizer::SiPermittivity = 11.8*8.854187817e-12/m;
const double AllPixITkStripsDigitizer::k_B = 1.38064852e-23*m*m*kg/s/s;

const double AllPixITkStripsDigitizer::kMinEnergy = 0.1*eV;
const double AllPixITkStripsDigitizer::kMaxEnergy = 14*TeV;

const double AllPixITkStripsDigitizer::kIterations = 1000;
const double AllPixITkStripsDigitizer::kMaxError = 1e-4 ;


extern DebugLevel debug;


AllPixITkStripsDigitizer::AllPixITkStripsDigitizer(G4String modName, G4String hitsColName, G4String digitColName)
: AllPixDigitizerInterface (modName) ,
  m_primaryVertex(0),
  m_doFast(false),
  m_minTime(-25.0*ns),
  m_maxTime(50.0*ns)
{

	// Registration of digits collection name
	collectionName.push_back(digitColName);
	m_hitsColName.push_back(hitsColName);

	// Detector description handle
	// provided by the interface
	m_geometry = GetDetectorGeoDscPtr();
	m_detectorWidth = m_geometry->GetSensorZ();

	// threshold from gear file, in units of charge
	SetThreshold(m_geometry->GetThreshold());

	m_depletionVoltage = m_geometry->GetDepletionVoltage();
	m_biasVoltage = m_geometry->GetBiasVoltage();
	m_Neff = 5.3e12 * 1/TMath::Power(cm, 3); // 1/cm3 FIXME
	m_depletionDepth = TMath::Sqrt((2*SiPermittivity*m_biasVoltage)/(e_SI * TMath::Abs(m_Neff)));
	//INCLUDED FLUENCE FROM m_GEOMETRY
	m_fluence = m_geometry->GetFluence();

	if (m_depletionDepth > m_detectorWidth)
		m_depletionDepth = m_detectorWidth;
	SetTemperatureK(m_geometry->GetSensorTemperature());

	m_diffusion[Electron] = 36.62*cm2/s; // Electron mobility (cm2/s)
	m_diffusion[Hole]=18*cm2/s;// Hole mobility (cm2/s

	m_pitchX = m_geometry->GetPixelX();
	m_nStrips = m_geometry->GetNPixelsX();
	m_nRows = m_geometry->GetNPixelsY();
	m_pitchY = m_geometry->GetPixelY();

	////////////////////////////////////
	// Numerical integration accuracy //
	////////////////////////////////////

	m_tStepL =  0.001*ns;
	m_tStepU=	0.5*ns;
	m_stepSize = m_tStepU/10.;

	//m_precision = 50; // Precision of charge propagation
	m_precision = 10;
	if (m_precision == 1 and debug==ERROR) {
		G4cout << " [AllPixITkStripDigitizer] ERROR precision set to 1. Too low for real simulation. Exiting." << G4endl;
		exit(-1);
	}

	if (debug>ERROR) {

		G4cout << " [AllPixITkStripDigitizer] Bias voltage " << m_biasVoltage << "V" << endl;
		G4cout << " [AllPixITkStripDigitizer] Depletion voltage " << m_depletionVoltage << "V" << endl;
		G4cout << " [AllPixITkStripDigitizer] threshold " << m_inputParameter.thl/eV << " eV " << TString::Format(" %d pairs\n", getThresholdPairs());
		G4cout << " [AllPixITkStripDigitizer] Sensor temperature " << m_temperature << " K " << G4endl;
		G4cout << " [AllPixITkStripDigitizer] Detector w = " << m_detectorWidth/mm << " mm" << G4endl;
		G4cout << " [AllPixITkStripDigitizer] Depletion d = " << m_depletionDepth/mm << " mm" << G4endl;
		G4cout << " [AllPixITkStripDigitizer] Detector T = "<< m_temperature << " K " << G4endl;
		G4cout << " [AllPixITkStripDigitizer] No. of strips in X = "<< m_nStrips<< G4endl;
		G4cout << " [AllPixITkStripDigitizer] No. of strips in Y = "<< m_nRows<< G4endl;
		G4cout << " [AllPixITkStripDigitizer] pitchX [um] = "<< m_pitchX/um<< G4endl;
		G4cout << " [AllPixITkStripDigitizer] pitchY [um] = "<< m_pitchY/um<< G4endl;
		G4cout << " [AllPixITkStripDigitizer] fluence [neq/cm2] = "<<m_fluence*1/cm2<<" (neq/cm2)"<< G4endl;
			}





	if (debug == DEBUG)
		m_testFile = fopen("output/testfile.dat","w"); // FIXME try and catch
}

AllPixITkStripsDigitizer::~AllPixITkStripsDigitizer(){
	if (debug == DEBUG)
		fclose(m_testFile); // FIXME try and catch

}

G4double AllPixITkStripsDigitizer::getMobility(double x, double y, double z, ElectricField _E, CarrierType carrier){

	// Initialize variables so they have the right scope
	G4double saturationV = 0;
	G4double criticalE = 0;
	G4double beta = 0;

	//These parameterizations come from C. Jacoboni et al., Solid‐State Electronics 20 (1977) 77‐89. (see also https://cds.cern.ch/record/684187/files/indet-2001-004.pdf).
	if(carrier == Electron){
		saturationV = 1.53e9 * TMath::Power(m_temperature, -0.87) * cm / s;
		criticalE = 1.01 * pow(m_temperature, 1.55) / cm;
		beta = 2.57E-2 * pow(m_temperature,0.66);
	}
	else if(carrier == Hole){
		saturationV = 1.62e8*pow(m_temperature,-0.52) * cm / s;
		criticalE = 1.24E-7*pow(m_temperature,1.68) / cm;
		beta = 0.46*pow(m_temperature,0.17);
	}
	else {
		cout << " Carrier undefinied!" << endl;
		exit(-1);
	}

	//ElectricField _E = getEField(x, y, z);
	G4double mobility = (saturationV/criticalE) / pow(1+pow((_E.Mag()/criticalE), beta), (1/beta));
	/*if (debug>INFO) {
		if (carrier == Electron)
			G4cout << " [AllPixITkStripDigitizer] mu_e [cm2/Vs] = " <<
			mobility/ cm2 * s << G4endl;
		else
			G4cout << " [AllPixITkStripDigitizer] mu_h [cm2/Vs] = " <<
			mobility /cm2 * s << G4endl;
	}*/

	if (carrier == Electron) mobility *= -1;

	return mobility;
}

G4double AllPixITkStripsDigitizer::GetDriftVelocity(double x, double y, double z, CarrierType carrier){

	ElectricField electricField = getEField(x, y, z);
	double Mobility = getMobility(x, y, z, electricField, carrier);
	G4double driftVelocity = Mobility*electricField.Mag();
	if(carrier == Hole) driftVelocity = -1*driftVelocity; // Drifts in opposite direction
	//if (debug > ERROR) G4cout << " [AllPixITkStripDigitizer] v = " << driftVelocity/cm*s << " cm/s " << G4endl;
	return driftVelocity;
}

G4double AllPixITkStripsDigitizer::GetMeanFreePath(G4double driftVelocity, G4bool isHoleBit){

	G4double meanFreePath = 0;
	if(!isHoleBit) meanFreePath = driftVelocity*m_trappingTimeElectrons; // mm
	if(isHoleBit) meanFreePath = driftVelocity*m_trappingTimeHoles; // mm
	return meanFreePath;
}


ElectricField AllPixITkStripsDigitizer::getEField(double x, double y, double z){

	if (z<-kMaxError or z>(m_detectorWidth+kMaxError) ) {

		if (debug==DEBUG) {
			G4cout << TString::Format(" [ITkStripDigitizer::EField1D] E = 0 V/cm") << G4endl;
			G4cout << " [ITkStripDigitizer::getEField1D] Carrier out of sensitive strip: " << z/um << G4endl;
		}

		return ElectricField(0, 0, 0);

	}
	
	//if(readoutType==ELECTRON)
	//{

	////////////////////////////////////////////
	//// Electric Field from TCAD MAP//////////
	//////////////////////////////////////////
	TFile* efile = new ("/afs/cern.ch/user/m/minano/work/public/allpix/tcad_map/EField-map-strips-310um-f0e15-400V.root");
	//For all maps, 0 is at the collecting electrode and L is the far side.
       

	double electricFieldX = 0.0;
	double electricFieldY = 0.0;
	double electricFieldZ=(m_biasVoltage+m_depletionVoltage)/m_detectorWidth+
			2*m_depletionVoltage/pow(m_detectorWidth,2)*z;

	if(z>m_depletionDepth){
		electricFieldZ=0;
		//cout << "bip" << endl;
	}

	//double electricField = std::sqrt(electricFieldX*electricFieldX + electricFieldY*electricFieldY + electricFieldZ*electricFieldZ);
	//if (debug==DEBUG) G4cout << TString::Format(" [ITkStripDigitizer::EField1D] E = %g V/cm", electricFieldZ*cm)<< G4endl;
	return ElectricField(electricFieldX, electricFieldY, electricFieldZ);
}

G4double AllPixITkStripsDigitizer::getDriftTime(G4double z, CarrierType carrier, bool doFast)
{
	// drift in a uniform electric field from hit position to surface
	double driftTime = 0;

	if (doFast){
		ElectricField _E = getEField(0,0,z);
		if (_E.Mag() != 0) {

			if (debug==DEBUG) G4cout << " [ComputeDriftTime] Fast ballistic calculation of drift time." << endl;
			// tz is transformation vector to c.s. I quadrant
			//double tz = z + m_detectorWidth/2.;
			// w is the traversed path by electron
			//double w = m_detectorWidth - tz;
			double w = m_detectorWidth - z;
			if (debug == DEBUG) G4cout<< TString::Format(" [ComputeDriftTime] Traversed path by electron %f um\n", w/um);

			if (w > kMaxError) 
				driftTime = w / (GetDriftVelocity(0, 0, z, carrier));
			else
				driftTime = 0.0*ns;
		}

		else
			driftTime = 1000*ns;
	}

	else
	{

		if (debug==DEBUG) G4cout << " [ComputeDriftTime] Numeric integration." << endl;

		array<double,4> finalVector = getDriftVector(0,0,z);
		driftTime = finalVector.at(3);

	}
	if (debug==DEBUG) G4cout << TString::Format(" [ComputeDriftTime] t_drift: %g [ns]\n", driftTime/ns);

	return driftTime;
}

array<double,4>  AllPixITkStripsDigitizer::getDriftVector(double x, double y, double z)
{
	G4double driftTime=0;
	G4double dt=m_stepSize;
	G4double xtemp=x;
	G4double ytemp=y;
	G4double ztemp=z;
	//cout << "z : " << z/um << endl;
	//vector<G4double> step;
	array<double,4> step;
	int counter=0;

	//vector<G4double> vx,vy,vz,vt;

	if (debug>INFO) cout << " [AllPixITkStripDigitizer::getDriftVector]" << endl;

	int iter =0;
	while( ztemp>kMaxError and ztemp<(m_detectorWidth-kMaxError) and iter<kIterations){

		driftTime+=dt;

		/*if(B_Field!=0.0){

			if(readoutType==ELECTRON){
				step=RKF5IntegrationBFieldElectrons(xtemp,ytemp,ztemp,dt);
			}
			else{
				step=RKF5IntegrationBFieldHoles(xtemp,ytemp,ztemp,dt);
			}
		}
		else {*/

		step=RKF5IntegrationElectrons(xtemp, ytemp, ztemp, dt);


		xtemp+=step.at(0);
		ytemp+=step.at(1);
		ztemp+=step.at(2);

		// Check the boundaries
		//if( getENorm(ztemp)==0 ) ztemp = m_detectorWidth/2.0;
		if (ztemp < 0) ztemp = 0.0;
		if (ztemp > m_detectorWidth) ztemp = m_detectorWidth;

		// Adaptable step size
		dt=setDt(dt, step.at(3));


		counter++;

		iter++;
		// for debugging purposes - dots can be plotted
		if (debug == DEBUG)
			fprintf(m_testFile, "ixyz %d %f %f %f\n",iter, xtemp/um, ytemp/um, ztemp/um); // FIXME try and catch
	}
	if (debug == DEBUG) {
		cout << " [AllPixITkStripDigitizer::getDriftVector] No. of iterations : " << iter << endl;
		cout << " [AllPixITkStripDigitizer::getDriftVector] After drift x : " << xtemp/um << " y : " << ytemp/um << " z : " << ztemp/um << " t : " << driftTime/ns << endl;
	}

	m_stepSize = dt;

	array<double,4> output;
	output[0]=xtemp;
	output[1]=ytemp;
	output[2]=ztemp;
	output[3]=driftTime;

	return output;
}

G4double  AllPixITkStripsDigitizer::setDt(G4double dt, G4double ErreurMoy)
{

	double Dt=dt;

	if (debug == DEBUG) cout << " [AllPixITkStripDigitizer::setDt] Starting dt : " << Dt/ns << endl;

	//if(isnan(ErreurMoy)){Dt=tup;}
	if (ErreurMoy != ErreurMoy)
		Dt = m_tStepU;
	else if (fabs(ErreurMoy) > kMaxError) { Dt*=0.9;}
	else if (fabs(ErreurMoy) < kMaxError) { Dt*=1.1;};


	if(Dt<m_tStepL) Dt = m_tStepL;
	if(Dt>m_tStepU) Dt = m_tStepU;

	if(debug==DEBUG) cout << " [AllPixITkStripDigitizer::setDt] Changed to dt : " << Dt/ns << endl;
	return Dt;

}

G4double AllPixITkStripsDigitizer::getDiffusionLRMS(G4double tDrift)
{
	//if(readoutType==ELECTRON){
	if ((tDrift < 0.01*ns) or (tDrift > 100*ns)) return 0.0; //!FIXME bounds hardcoded
	//m_diffusion[Electron] = k_B*m_temperature/e_SI*GetMobility(Electron);
	if (debug == DEBUG) G4cout << TString::Format(" [AllPixITkStripDigitizer::getDiffusionLRMS] diffusion const D_e :%f cm2/s\n", m_diffusion[Electron]/cm/cm*s);
	double _diffusion = TMath::Sqrt(2.*m_diffusion[Electron]*tDrift);
	//}

	//return TMath::Sqrt(2.*Default_Hole_D*tDrift);
	if (debug ==DEBUG) G4cout<<" [AllPixITkStripDigitizer::getDiffusionLRMS] Diffusion length [um]" << _diffusion/um << G4endl;
	return _diffusion;

}



G4double AllPixITkStripsDigitizer::IntegrateGaussian(G4double xhit, G4double yhit, G4double Sigma,
		G4double x1, G4double x2, G4double y1, G4double y2, G4double Energy )
{
	//G4cout << TString::Format("Integration borns xhit=%f yhit=%f x1=%f x2=%f y1=%f y2=%f",xhit/1e3,yhit/1e3,x1/1e3,x2/1e3,y1/1e3,y2/1e3) << endl;
	//G4cout << TString::Format("TMath::Erf test %f %f %f %f",TMath::Erf((x1 - xhit)/(Sqrt(2.)*Sigma)),TMath::Erf((x2 - xhit)/(TMath::Sqrt(2.)*Sigma)),TMath::Erf((y1 - yhit)/(TMath::Sqrt(2.)*Sigma)),TMath::Erf((y2 - yhit)/(TMath::Sqrt(2.)*Sigma)))<< endl;

	//G4cout << " [AllPixITkStripDigitizer::IntegrateGaussian]" << G4endl <<
	//	TString::Format("                             xhit %f yhit %f L %f x1 %f x2 %f y1 %f y2 %f", xhit, yhit, Sigma, x1, x2, y1, y2)<< G4endl;

	// Boundary checks
	if (Energy < kMinEnergy) return 0;
	Energy /= eV; /* convert to eV */

	if (Sigma < kMaxError) return 0;

	double energybis= (-TMath::Erf((x1 - xhit)/(TMath::Sqrt(2.)*Sigma)) + TMath::Erf((x2 - xhit)/(TMath::Sqrt(2.)*Sigma)));
	energybis *=      (-TMath::Erf((y1 - yhit)/(TMath::Sqrt(2.)*Sigma)) + TMath::Erf((y2 - yhit)/(TMath::Sqrt(2.)*Sigma)));
	energybis *= Energy/4.0;

	if (debug==DEBUG) {
		G4cout << " [AllPixITkStripsDigitizer::IntegrateGaussian] Energy : " << energybis/eV << " [eV] Pairs : " << floor(energybis/kPairEnergy) << endl;
	}
	/* convert back to system units */
	return energybis*eV;

}

void AllPixITkStripsDigitizer::Digitize(){

	// create the digits collection
	m_digitsCollection = new AllPixITkStripsDigitsCollection("AllPixITkStripsDigitizer", collectionName[0] );

	// get the digiManager
	G4DigiManager * digiMan = G4DigiManager::GetDMpointer();

	// BoxSD_0_HitsCollection
	G4int hcID = digiMan->GetHitsCollectionID(m_hitsColName[0]);

	AllPixTrackerHitsCollection * hitsCollection = 0;
	hitsCollection = (AllPixTrackerHitsCollection*)(digiMan->GetHitsCollection(hcID));

	// temporary data structure
	map<pair<G4int, G4int>, G4double > stripContent; //! ((x,y), depositedEnergy)
	pair<G4int, G4int> currentStrip;

	G4int nEntries = hitsCollection->entries();

	// Loop over hits, deposition of energy along particle path

	float depositedEnergy = 0;

	for(G4int iHit  = 0 ; iHit < nEntries ; iHit++) {
		AllPixTrackerHit* hit = (*hitsCollection)[iHit];
		// Perform range checking FIXME

		// Locate the strip where the simulated hit is
		currentStrip.first  = hit->GetPixelNbX();
		currentStrip.second = hit->GetPixelNbY();
		double hitEnergy = hit->GetEdep();

		if ( (hitEnergy < kMinEnergy) or (hitEnergy > kMaxEnergy) )
			hitEnergy = 0.0;
		//pixelsContent[currentStrip] += hitEnergy;

		// Calculate drift time from the hit position to the electrode
		// Local frame of reference as if strip defines Ist quadrant
		const double sensorHalfX = m_pitchX * m_nStrips /2.;
		const double sensorHalfY = m_pitchY * m_nRows/2.;

		const double xLocal = (hit->GetPosInLocalReferenceFrame()).x() + sensorHalfX;
		const double yLocal = (hit->GetPosInLocalReferenceFrame()).y() + sensorHalfY;

		const double pixelX = m_pitchX * TMath::FloorNint(xLocal/m_pitchX); // pixel frame of reference
		const double pixelY = m_pitchY * TMath::FloorNint(yLocal/m_pitchY);

		const G4double hitPixelX =  xLocal - pixelX;
		const G4double hitPixelY =  yLocal - pixelY;
		const G4double hitPixelZ = (hit->GetPosInLocalReferenceFrame()).z() + m_detectorWidth/2.;

		if (debug>=INFO) {
			G4ThreeVector position = hit->GetPos();
			G4cout << " [ITkStripDigitizer::Digitize] Hit = "<< " X Y Z "<< position.x() << " "<< position.y() << " "<< position.z()<< G4endl;
			G4cout << " [ITkStripDigitizer::Digitize] Strip = "<< currentStrip.first<< " Row = "<< currentStrip.second<< " Energy [eV]= "<< hit->GetEdep()/eV << G4endl;
			G4cout << " [ITkStripDigitizer::Digitize] Local x [um] = "<< hitPixelX/um << G4endl;
			G4cout << " [ITkStripDigitizer::Digitize] Local y [um] = "<< hitPixelY/um << G4endl;
			G4cout << " [ITkStripDigitizer::Digitize] Local z [um] = "<< hitPixelZ/um << G4endl;
			G4cout << " [ITkStripDigitizer::Digitize] Process name = "<< hit->GetProcessName() << G4endl;
		}

		if (not(iHit%10)) hit->Draw();

		// Should the energy be spread?
		//const double sigma = 3.0;
		G4double eHitTotal = hitEnergy;
		// underdepletion
		eHitTotal *= m_depletionDepth/m_detectorWidth;

		pair<G4int, G4int> extraStrip; /* In case there is charge sharing. */

		// Propagate electrons to electrode in bunches instead of individual pairs
		for(G4int nQ  = 0 ; nQ < m_precision ; nQ++) {

			double eHit = double(eHitTotal)/m_precision; // Should be 2*precision, if including holes

			if (debug==DEBUG) {
				G4cout << " [ITkStripDigitizer::Digitize] Bunch " << nQ <<" dE [eV]   = " << eHit/eV << G4endl;
			}
			if (eHit/eV < kMinEnergy) eHit=0;

			const double driftTime = getDriftTime(hitPixelZ, Electron, m_doFast);


			// charge share integration boundaries +-2*DiffusionSigma
			const double diffusionWidth = 3 * getDiffusionLRMS(driftTime);

			// Is there charge sharing - diffusion
			if (hitPixelX + diffusionWidth > m_pitchX) {
				// Main strip
				double stripEnergy = IntegrateGaussian(hitPixelX/um, hitPixelY/um, diffusionWidth/um,
						(0)/um, (m_pitchX)/um,
						(0)/um, (m_pitchY)/um, eHit);
				stripContent[currentStrip] += stripEnergy;
				depositedEnergy += stripEnergy;
				if (debug>=INFO) G4cout << TString::Format(" [ITkStripDigitizer::Digitize] Main strip %d Bunch %d Deposited energy: %f eV\n", currentStrip.first, nQ, stripContent[currentStrip]/eV);


				/* There is charge sharing into the right strip */
				extraStrip.first = currentStrip.first + 1;
				extraStrip.second = currentStrip.second;

				if(extraStrip.first >= 0 && extraStrip.second>=0 && extraStrip.first < m_nStrips && extraStrip.second < m_nRows)
				{
				double sharedEnergy = IntegrateGaussian(hitPixelX/um, hitPixelY/um, diffusionWidth/um,
							(m_pitchX)/um, (2*m_pitchX)/um,
							(0)/um, (m_pitchY)/um, eHit);
				stripContent[extraStrip] += sharedEnergy;
				depositedEnergy += sharedEnergy;

				if (debug>=INFO) G4cout << TString::Format(" [ITkStripDigitizer::Digitize] Shared into strip %d Bunch %d Shared energy: %f eV\n", extraStrip.first, nQ, stripContent[extraStrip]/eV);
			}
				else {
					//stripContent[extraStrip] = 0; // I think this causes to have strips with negative index
					G4cout<<" [] Charge loss at the edge of sensor."<<G4endl;
				}
			}

			else if (hitPixelX - diffusionWidth < 0.0) {

				// Main strip
				double stripEnergy = IntegrateGaussian(hitPixelX/um, hitPixelY/um, diffusionWidth/um,
						(0)/um, (m_pitchX)/um,
						(0)/um, (m_pitchY)/um, eHit);
				stripContent[currentStrip] += stripEnergy;
				depositedEnergy += stripEnergy;
				if (debug>=INFO) G4cout << TString::Format(" [ITkStripDigitizer::Digitize] Main strip %d Bunch %d Deposited energy: %f eV\n", currentStrip.first, nQ, stripContent[currentStrip]/eV);


				// There is charge sharing into the left strip
				// Charge sharing between strips, only in X

				extraStrip.first = currentStrip.first - 1;
				extraStrip.second = currentStrip.second;

				if(extraStrip.first >= 0 && extraStrip.second>=0 && extraStrip.first < m_nStrips && extraStrip.second < m_nRows)
				{
					double sharedEnergy = IntegrateGaussian(hitPixelX/um, hitPixelY/um, diffusionWidth/um,
							(- m_pitchX)/um, (0)/um,
							(0)/um, (m_pitchY)/um, eHit);
					stripContent[extraStrip] += sharedEnergy;
					depositedEnergy += sharedEnergy;

					if (debug>=INFO) G4cout << TString::Format(" [ITkStripDigitizer::Digitize] Share into strip %d Bunch %d Shared energy: %f eV\n", extraStrip.first, nQ, stripContent[extraStrip]/eV);
				}
				else {
					//stripContent[extraStrip] = 0; // Strip with negative index?
					G4cout<<" [] Charge loss at the edge of sensor."<<G4endl;
				}

	   } // Charge sharing
	else {
		// No charge sharing, add this bunch to the current strip
		stripContent[currentStrip] += eHit;
		depositedEnergy += eHit;

		if (debug>=INFO) G4cout << TString::Format(" [AllPixITkStripDigitizer::Digitize] No charge sharing. Bunch %d. Adding %f eV to strip %d row %d\n",
				nQ, eHit/eV, currentStrip.first,currentStrip.second);
	}
} // nQ

} // Hits in one strip
if (debug>=INFO) G4cout << " [AllPixITkStripDigitizer::Digitize] Total deposited energy keV " << depositedEnergy/keV << G4endl;

//G4cout << TString::Format("Collected energy in pixel : %f eV", );

// Now create digits, one per strips
map<pair<G4int, G4int>, G4double >::iterator iCount = stripContent.begin();
// Loop over added strips
for( ; iCount != stripContent.end() ; iCount++)
{
	if (debug==DEBUG) G4cout << " [AllPixITkStripDigitizer::Digitize] iCount " << (*iCount).first.first << G4endl;

	double stripEnergy = (*iCount).second;
	double depositedCharge = stripEnergy/kPairEnergy * CLHEP::e_SI * coulomb;
	if (debug >= INFO) G4cout << " Deposited charge [fC]" << depositedCharge/fC <<G4endl;

	const double tau = 21.0 / 3.0 *ns; //ns CR-RC3
	//const double tau = 22.0 / 2.523 * ns; // CRRC2

	bool crossed = false;
	for (int timeSlice=m_minTime; timeSlice<m_maxTime; timeSlice++) {
		double t = timeSlice/tau; // that's what I found in Athena SCT code
		double Y = 0;
		// CR-RC3 response, as in SCT
		if (t>0) Y = depositedCharge / 27.0 * exp(3) * t*t*t * exp(-t);
		// CR-RC2 response to step
		//if (t>0) Y = depositedCharge / 6.0 * t*t*t * exp(-t);

		if (crossed == false and Y >= getThresholdCharge() and timeSlice < 50) crossed = true; // FIXME hardcoded timebin

		if (m_testFile) fprintf(m_testFile, "FE %d %f\n", timeSlice, Y/fC);
		if (debug == DEBUG) {
			if (Y>getThresholdCharge()) printf("Above threshold = %f fC\n", getThresholdCharge()/fC);
			else printf ("Below threshold = %f fC\n", getThresholdCharge()/fC);
		}
	}

	// Create one digit per pixel
	AllPixITkStripsDigit * digit = new AllPixITkStripsDigit;
	digit->SetPixelIDX((*iCount).first.first);
	digit->SetPixelIDY((*iCount).first.second);
	digit->SetPixelCounts(crossed ? 1 : 0); // Binary readout, no ToT
	digit->SetPixelEnergyDep(stripEnergy);
	digit->SetPrimaryVertex(m_primaryVertex->GetPosition());

	m_digitsCollection->insert(digit);

		if (debug>=INFO) {
			G4cout << " [ITkStripsDigitizer::Digitize] Total energy in strip " << digit->GetPixelIDX()
									<< " = " << digit->GetPixelEnergyDep()/keV << " keV, counts = " <<digit->GetPixelCounts() << endl;
			//G4cout << " [ITkStripsDigitizer::Digitize] Collection size " << collectionSize << G4endl;
		}
}

G4int dc_entries = m_digitsCollection->entries();
if(dc_entries > 0){
	if (debug>ERROR) G4cout << "--------> Digits Collection : " << collectionName[0]
	                                                             << "(" << m_hitsColName[0] << ")"
	                                                             << " contains " << dc_entries
	                                                             << " digits" << G4endl;
}
StoreDigiCollection(m_digitsCollection);

}

array<double,4>  AllPixITkStripsDigitizer::RKF5IntegrationElectrons(G4double x, G4double y, G4double z,G4double dt)
{
	// This function transport using Euler integration, for field (Ex,Ey,Ez),
	// considered constant over time dt. The movement equation are those
	// of charges in semi-conductors, sx= mu*E*dt;;
	double k1x,k2x,k3x,k4x,k5x,k6x;
	double k1y,k2y,k3y,k4y,k5y,k6y;
	double k1z,k2z,k3z,k4z,k5z,k6z;
	double dx,dy,dz;

	ElectricField E = getEField(x,y,z);

	k1x=-getMobility(x, y, z, E, Electron)*E.X()*dt;
	k1y=-getMobility(x, y, z, E, Electron)*E.Y()*dt;
	k1z=-getMobility(x, y, z, E, Electron)*E.Z()*dt;
	if (debug == DEBUG) cout << " [AllPixITkStripDigitizer::RKF5IntegrationElectrons] k1x : "<<k1x<<", k1y : " << k1y << ", k1z : " << k1z << endl;

	E = getEField(x+k1x/4,y+k1y/4,z+k1z/4);

	k2x=-getMobility(x+k1x/4,y+k1y/4,z+k1z/4, E, Electron)*E.X()*dt;
	k2y=-getMobility(x+k1x/4,y+k1y/4,z+k1z/4, E, Electron)*E.Y()*dt;
	k2z=-getMobility(x+k1x/4,y+k1y/4,z+k1z/4, E, Electron)*E.Z()*dt;

	E = getEField(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z);

	k3x=-getMobility(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z, E, Electron)*E.X()*dt;
	k3y=-getMobility(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z, E, Electron)*E.Y()*dt;
	k3z=-getMobility(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z, E, Electron)*E.Z()*dt;

	E = getEField(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z);

	k4x=-getMobility(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z, E, Electron)*E.X()*dt;
	k4y=-getMobility(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z, E, Electron)*E.Y()*dt;
	k4z=-getMobility(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z, E, Electron)*E.Z()*dt;

	E = getEField(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z);

	k5x=-getMobility(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z, E, Electron)*E.X()*dt;
	k5y=-getMobility(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z, E, Electron)*E.Y()*dt;
	k5z=-getMobility(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z, E, Electron)*E.Z()*dt;

	E = getEField(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
			y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
			z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z);

	k6x=-getMobility(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
			y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
			z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z, E, Electron)*E.X()*dt;
	k6y=-getMobility(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
			y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
			z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z, E, Electron)*E.Y()*dt;
	k6z=-getMobility(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
			y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
			z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z, E, Electron)*E.Z()*dt;

	dx=((16./135)*k1x+(6656./12825)*k3x+(28561./56430)*k4x-(9./50)*k5x+(2./55)*k6x);
	dy=((16./135)*k1y+(6656./12825)*k3y+(28561./56430)*k4y-(9./50)*k5y+(2./55)*k6y);
	dz=((16./135)*k1z+(6656./12825)*k3z+(28561./56430)*k4z-(9./50)*k5z+(2./55)*k6z);

	double Ex,Ey,Ez,Erreur;
	Ex=((1./360)*k1x-(128./4275)*k3x-(2197./75240)*k4x-(1./50)*k5x+(2./55)*k6x);
	Ey=((1./360)*k1y-(128./4275)*k3y-(2197./75240)*k4y-(1./50)*k5y+(2./55)*k6y);
	Ez=((1./360)*k1z-(128./4275)*k3z-(2197./75240)*k4z-(1./50)*k5z+(2./55)*k6z);
	Erreur=sqrt(Ex*Ex+Ey*Ey+Ez*Ez);

	if (debug == DEBUG) cout << " [AllPixITkStripDigitizer::RKF5IntegrationElectrons] k6x : "<<k6x<<", k6y : " << k6y << ", k6z : " << k6z << endl;

	array<double,4> newpoint;
	newpoint.at(0)=(dx);
	newpoint.at(1)=(dy);
	newpoint.at(2)=(dz);
	newpoint.at(3)=(Erreur);

	return newpoint;

}
