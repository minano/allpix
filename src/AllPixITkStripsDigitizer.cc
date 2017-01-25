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

const double AllPixITkStripsDigitizer::kMinEnergy = 1*eV;
const double AllPixITkStripsDigitizer::kMaxEnergy = 14*TeV;

const double kIterations = 1000; // maximal number of iter for num solving of Euler eq.


extern DebugLevel debug;

AllPixITkStripsDigitizer::AllPixITkStripsDigitizer(G4String modName, G4String hitsColName, G4String digitColName) 
: AllPixDigitizerInterface (modName) ,
  m_digitsCollection(0),
  m_primaryVertex(0),
  m_doFast(false),
  m_minTime(-25.0*ns),
  m_maxTime(150.0*ns)
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
	m_maxError = 1e-4 ;
	m_tStepL =  0.001*ns;
	m_tStepU=	0.5*ns;
	m_stepSize = m_tStepU/10.;

	m_precision = 10; // Precision of charge propagation

	if (debug>ERROR) {
		G4cout << " [AllPixITkStripDigitizer] Bias voltage " << m_biasVoltage << "V" << G4endl;
		G4cout << " [AllPixITkStripDigitizer] Depletion voltage " << m_depletionVoltage << "V" << G4endl;
		G4cout << " [AllPixITkStripDigitizer] threshold " << m_inputParameter.thl/eV << " eV " << TString::Format(" %d pairs\n", getThresholdPairs());
		G4cout << " [AllPixITkStripDigitizer] Sensor temperature " << m_temperature << " K " << G4endl;
		G4cout << " [AllPixITkStripDigitizer] Detector w = " << m_detectorWidth/mm << " mm" << G4endl;
		G4cout << " [AllPixITkStripDigitizer] Depletion d = " << m_depletionDepth/mm << " mm" << G4endl;
		G4cout << " [AllPixITkStripDigitizer] Detector T = "<< m_temperature << " K " << G4endl;
		G4cout << " [AllPixITkStripDigitizer] No. of strips in X = "<< m_nStrips<< G4endl;
		G4cout << " [AllPixITkStripDigitizer] No. of strips in Y = "<< m_nRows<< G4endl;
		G4cout << " [AllPixITkStripDigitizer] pitchX = "<< m_pitchX<< G4endl;
		G4cout << " [AllPixITkStripDigitizer] pitchY = "<< m_pitchY<< G4endl;
	}

	if (debug == DEBUG)
		m_testFile = fopen("output/testfile.dat","w");
}

AllPixITkStripsDigitizer::~AllPixITkStripsDigitizer(){
	if (debug>=INFO) G4cout << " [AllPixITkStripsDigitizer] Destructor." << G4endl;
	if (debug == DEBUG)
		fclose(m_testFile);

}

G4double AllPixITkStripsDigitizer::getMobility(double x, double y, double z, CarrierType carrier){

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
		G4cout << " Carrier undefinied!" << G4endl;
		exit(-1);
	}

	ElectricField _E = getEField1D(x, y, z);
	G4double mobility = (saturationV/criticalE) / pow(1+pow((_E/criticalE), beta), (1/beta));
	if (debug>INFO) {
		if (carrier == Electron)
			G4cout << " [AllPixITkStripDigitizer] mu_e [cm2/Vs] = " <<
			mobility/ cm2 * s << G4endl;
		else
			G4cout << " [AllPixITkStripDigitizer] mu_h [cm2/Vs] = " <<
			mobility /cm2 * s << G4endl;
	}

	if (carrier == Electron) mobility *= -1;

	return mobility;
}

G4double AllPixITkStripsDigitizer::GetDriftVelocity(double x, double y, double z, CarrierType carrier){

	ElectricField electricField = getEField1D(x, y, z);
	double Mobility = getMobility(x, y, z, carrier);
	G4double driftVelocity = Mobility*electricField;
	if(carrier == Hole) driftVelocity = -1*driftVelocity; // Drifts in opposite direction
	if (debug > ERROR) G4cout << " [AllPixITkStripDigitizer] v = " << driftVelocity/cm*s << " cm/s " << G4endl;
	return driftVelocity;
}

G4double AllPixITkStripsDigitizer::GetMeanFreePath(G4double driftVelocity, G4bool isHoleBit){

	G4double meanFreePath = 0;
	if(!isHoleBit) meanFreePath = driftVelocity*m_trappingTimeElectrons; // mm
	if(isHoleBit) meanFreePath = driftVelocity*m_trappingTimeHoles; // mm
	return meanFreePath;
}


double AllPixITkStripsDigitizer::getEField1D(double x, double y, double z){

	if (z<-m_detectorWidth/2. or z>m_detectorWidth/2.) {
		if (debug==DEBUG) {
			G4cout << TString::Format(" [ITkStripDigitizer::EField1D] E = 0 V/cm") << G4endl;
		G4cout << " [ITkStripDigitizer::getEField1D] Carrier out of sensitive strip." << G4endl;
		}
				return 0.0;
	}
	//if(readoutType==ELECTRON)
	//{
	//double electricFieldX = 0.0;
	//double electricFieldY = 0.0;
	double electricFieldZ=(m_biasVoltage+m_depletionVoltage)/m_detectorWidth+
			2*m_depletionVoltage/pow(m_detectorWidth,2)*z;

	if(z>m_depletionDepth){
		electricFieldZ=0;
		//G4cout << "bip" << G4endl;
	}

	//double electricField = std::sqrt(electricFieldX*electricFieldX + electricFieldY*electricFieldY + electricFieldZ*electricFieldZ);
	if (debug==DEBUG) G4cout << TString::Format(" [ITkStripDigitizer::EField1D] E = %g V/cm", electricFieldZ*cm)
	<< G4endl;
	return electricFieldZ;
}

/*G4double AllPixTimepix3Digitizer::GetElectricFieldNorm(G4double , G4double , G4double ){

	return TMath::Sqrt(electricFieldX*electricFieldX +electricFieldY*electricFieldY +electricFieldZ*electricFieldZ );
}*/

double AllPixITkStripsDigitizer::getENorm(double z){
	double Ez = getEField1D(0,0,z);
	return TMath::Sqrt(Ez*Ez);
}

G4double AllPixITkStripsDigitizer::getDriftTime(G4double z, CarrierType carrier, bool doFast)
{
	// drift in a uniform electric field from hit position to surface
	double drift = 0;
	if (doFast){
		ElectricField _E = getEField1D(0,0,z);

		if (_E != 0) {
			if (debug==DEBUG) G4cout << " [ComputeDriftTime] Fast ballistic calculaton of drift time." << G4endl;
			// tz is transformation vector to c.s. I quadrant
			double tz = z + m_detectorWidth/2.;
			// w is the traversed path by electron
			double w = m_detectorWidth - tz;
			if (debug == DEBUG) G4cout<< TString::Format(" [ComputeDriftTime] Traversed path by electron %f um\n", w/um);

			if (w > 1E-6) // FIXME hardcoded precision
				drift = w / (GetDriftVelocity(0, 0, z, carrier));
			else
				drift = 0.0*ns;
		}
		else
			drift = 1000*ns;
	}
	else
	{
		if (debug==DEBUG) G4cout << " [ComputeDriftTime] Numeric integration." << G4endl;
		array<double,4> finalVector = getDriftVector(0,0,z);
		drift = finalVector.at(3);

	}
	if (debug==DEBUG) G4cout << TString::Format(" [ComputeDriftTime] t_drift: %g [ns]\n", drift/ns);

	return drift;
}

array<double,4>  AllPixITkStripsDigitizer::getDriftVector(double x, double y, double z)
{
	G4double driftTime=0;
	G4double dt=m_stepSize;
	G4double xtemp=x;
	G4double ytemp=y;
	G4double ztemp=z;
	//G4cout << "z : " << z/um << G4endl;
	//vector<G4double> step;
	array<double,4> step;
	int counter=0;

	//vector<G4double> vx,vy,vz,vt;

	if (debug>INFO) G4cout << " [AllPixITkStripDigitizer::getDriftVector]" << G4endl;

	int iter =0;
	while( ztemp>-m_detectorWidth/2.0 && ztemp<m_detectorWidth/2.0 && iter<kIterations){

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
		if (ztemp < -m_detectorWidth/2.0) ztemp = -m_detectorWidth/2.0;
		if (ztemp > m_detectorWidth/2.0) ztemp = m_detectorWidth/2.0;

		// Adaptable step size
		dt=setDt(dt, step.at(3));


		counter++;

		iter++;
		// for debugging purposes - dots can be plotted
		if (debug == DEBUG)
			fprintf(m_testFile, "ixyz %d %f %f %f\n",iter, xtemp/um, ytemp/um, ztemp/um);
	}
	if (debug == DEBUG) {
		G4cout << " [AllPixITkStripDigitizer::getDriftVector] No. of iterations : " << iter << G4endl;
		G4cout << " [AllPixITkStripDigitizer::getDriftVector] After drift x : " << xtemp/um << " y : " << ytemp/um << " z : " << ztemp/um << " t : " << driftTime/ns << G4endl;
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

	if (debug == DEBUG) G4cout << " [AllPixITkStripDigitizer::setDt] Starting dt : " << Dt/ns << G4endl;

	//if(isnan(ErreurMoy)){Dt=tup;}
	if (ErreurMoy != ErreurMoy)
		Dt = m_tStepU;
	else if (fabs(ErreurMoy) > m_maxError) { Dt*=0.9;}
	else if (fabs(ErreurMoy) < m_maxError) { Dt*=1.1;};


	if(Dt<m_tStepL) Dt = m_tStepL;
	if(Dt>m_tStepU) Dt = m_tStepU;

	if(debug==DEBUG) G4cout << " [AllPixITkStripDigitizer::setDt] Changed to dt : " << Dt/ns << G4endl;
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



G4double AllPixITkStripsDigitizer::IntegrateGaussian(G4double xhit, G4double yhit, G4double Sigma, G4double x1, G4double x2, G4double y1, G4double y2, G4double Energy )
{
	//G4cout << TString::Format("Integration borns xhit=%f yhit=%f x1=%f x2=%f y1=%f y2=%f",xhit/1e3,yhit/1e3,x1/1e3,x2/1e3,y1/1e3,y2/1e3) << G4endl;
	//G4cout << TString::Format("TMath::Erf test %f %f %f %f",TMath::Erf((x1 - xhit)/(Sqrt(2.)*Sigma)),TMath::Erf((x2 - xhit)/(TMath::Sqrt(2.)*Sigma)),TMath::Erf((y1 - yhit)/(TMath::Sqrt(2.)*Sigma)),TMath::Erf((y2 - yhit)/(TMath::Sqrt(2.)*Sigma)))<< G4endl;

	//G4cout << " [AllPixITkStripDigitizer::IntegrateGaussian]" << G4endl <<
	//	TString::Format("                             xhit %f yhit %f L %f x1 %f x2 %f y1 %f y2 %f", xhit, yhit, Sigma, x1, x2, y1, y2)<< G4endl;

	double energybis= (Energy*(-TMath::Erf((x1 - xhit)/(TMath::Sqrt(2.)*Sigma)) + TMath::Erf((x2 - xhit)/(TMath::Sqrt(2.)*Sigma)))
			*(-TMath::Erf((y1 - yhit)/(TMath::Sqrt(2.)*Sigma)) + TMath::Erf((y2 - yhit)/(TMath::Sqrt(2.0)*Sigma))))/4.;

	if (debug> INFO) {
		G4cout << " Energy : " << energybis/eV << " [eV] Pairs : " << floor(energybis/kPairEnergy) << G4endl;
	}
	return energybis;

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
	map<pair<G4int, G4int>, G4double > pixelsContent; //! ((x,y), depositedEnergy)
	pair<G4int, G4int> currentStrip;

	G4int nEntries = hitsCollection->entries();

	// Loop over hits, deposition of energy along particle path

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

		// Hit position in the global system
		if (debug>=INFO) {
			G4ThreeVector position = hit->GetPos();
			G4cout << " [ITkStripDigitizer::Digitize] Hit = "<< " X Y Z "<< position.x() << " "<< position.y() << " "<< position.z()<< G4endl;
			G4cout << " [ITkStripDigitizer::Digitize] Strip = "<< currentStrip.first<< " Row = "<< currentStrip.second << G4endl;
			G4cout << " [ITkStripDigitizer::Digitize] Energy [keV]= "<< hit->GetEdep()/keV << " Process name = "<< hit->GetProcessName() << G4endl;
		}

		// Hit in local c. s. of a strip
		const G4double hitLocalX = (hit->GetPosWithRespectToPixel()).x();
		const G4double hitLocalY = (hit->GetPosWithRespectToPixel()).y();
		const G4double hitLocalZ = (hit->GetPosWithRespectToPixel()).z();

		if (debug>=INFO) {
			G4cout << " [ITkStripDigitizer::Digitize] Local x [um] = "<< hitLocalX/um << G4endl;
			G4cout << " [ITkStripDigitizer::Digitize] Local y [um] = "<< hitLocalY/um << G4endl;
			G4cout << " [ITkStripDigitizer::Digitize] Local z [um] = "<< hitLocalZ/um << G4endl;
		}

		// Should the energy be spread?
		//const double sigma = 3.0;
		G4double eHitTotal = hitEnergy;
		// underdepletion
		eHitTotal *= m_depletionDepth/m_detectorWidth;

		// Propagate electrons to electrode in bunches instead of individual pairs
		for(G4int nQ  = 0 ; nQ < m_precision ; nQ++) {

			double eHit = double(eHitTotal)/m_precision; // Should be 2*precision, if including holes

			if (debug==DEBUG) G4cout << " [ITkStripDigitizer::Digitize] Bunch [eV]   = " << eHit/eV << G4endl;
			if (eHit/eV < kMinEnergy) continue;

			const double driftTime = getDriftTime(hitLocalZ, Electron, m_doFast);


			// charge share integration boundaries +-2*DiffusionSigma
			const double diffusionWidth = getDiffusionLRMS(driftTime);

			// Is there charge sharing - diffusion
			if ((hitLocalX > m_pitchX/2. - 3 * diffusionWidth) or (hitLocalX < - m_pitchX/2. + 3 * diffusionWidth)) {
				// Loop over neighbour strips
				// Charge sharing between strips, only in X
				pair<G4int, G4int> extraStrip = currentStrip;

				for(int i=-1;i<=1;i++){
					extraStrip.first = currentStrip.first + i;
					//for(int j=-nseek;j<=nseek;j++){
						if (debug==DEBUG) G4cout << TString::Format(" [ITkStripDigitizer::Digitize] strip %d row %d \n", extraStrip.first, extraStrip.second);

						//extraStrip.second += j;
						if(extraStrip.first >= 0 && extraStrip.second>=0 && extraStrip.first < m_nStrips && extraStrip.second < m_nRows)
						{
							double sharedEnergy = IntegrateGaussian(hitLocalX/um, hitLocalY/um, diffusionWidth/um,
									(-m_pitchX/2.0 + i*m_pitchX)/um, (-m_pitchX/2.+(i+1)*m_pitchX)/um,
									//(-m_pitchY/2 + j*m_pitchY)/nm, (-m_pitchY/2 + (j+1)*m_pitchY)/nm, eHit);
									(-m_pitchY/2)/um, (m_pitchY/2)/um, eHit);
							pixelsContent[extraStrip] += sharedEnergy;

							if (debug==DEBUG) G4cout << TString::Format(" [ITkStripDigitizer::Digitize] Bunch %d Deposited energy: %f eV\n", nQ, pixelsContent[extraStrip]/eV);
						}
					//}
				} // neighbour stirps
			} // Charge sharing
			else {
				// No charge sharing, add this bunch to current strip
				pixelsContent[currentStrip] += eHit;

				if (debug==DEBUG) G4cout << TString::Format(" [AllPixITkStripDigitizer::Digitize] No charge sharing. Bunch %d. Adding %f eV to strip %d row %d\n",
						nQ, eHit/eV, currentStrip.first,currentStrip.second);
			}
		} // nQ
	} // Hits in one strip


	//G4cout << TString::Format("Collected energy in pixel : %f eV", );

	// Now create digits, one per strips
		map<pair<G4int, G4int>, G4double >::iterator iCount = pixelsContent.begin();
		// Loop over added strips
		for( ; iCount != pixelsContent.end() ; iCount++)
		{
			double depositedEnergy = (*iCount).second;
			double depositedCharge = depositedEnergy/kPairEnergy * CLHEP::e_SI * coulomb;
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

				if (crossed == false and Y >= getThresholdCharge()) crossed = true;

				if (m_testFile) fprintf(m_testFile, "FE %d %f\n", timeSlice, Y/fC);
				if (debug == DEBUG) {
					if (Y>getThresholdCharge()) printf("Above threshold = %f fC\n", getThresholdCharge()/fC);
					else printf ("Below threshold = %f fC\n", getThresholdCharge()/fC);
				}
			}

			if( crossed ) // over threshold !
			{
				// Create one digit per pixel
				AllPixITkStripsDigit * digit = new AllPixITkStripsDigit(
						(*iCount).first.first,
						(*iCount).first.second,
						(*iCount).second,
						m_primaryVertex->GetPosition());

				int collectionSize = m_digitsCollection->insert(digit);

				if (debug>=INFO) {
					G4cout << " [ITkStripsDigitizer::Digitize] Total energy in strip " << digit->GetPixelIDX()
							<< " = " << digit->GetPixelEnergyDep()/MeV << " MeV, counts = " <<digit->GetPixelCounts() << G4endl;
					G4cout << " Collection size " << collectionSize << G4endl;
				}
			}
			else if (debug>=INFO) G4cout << " [ITkStripsDigitizer::Digitize] No digits in strip " << (*iCount).first.first << " in this event." << G4endl;
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

	ElectricField electricFieldZ = getEField1D(x,y,z);

	double electricFieldX = 0;
	double electricFieldY = 0;

	k1x=-getMobility(x, y, z, Electron)*electricFieldX*dt;
	k1y=-getMobility(x, y, z, Electron)*electricFieldY*dt;
	k1z=-getMobility(x, y, z, Electron)*electricFieldZ*dt;
	if (debug == DEBUG) G4cout << " [AllPixITkStripDigitizer::RKF5IntegrationElectrons] k1x : "<<k1x<<", k1y : " << k1y << ", k1z : " << k1z << G4endl;

	electricFieldZ = getEField1D(x+k1x/4,y+k1y/4,z+k1z/4);

	k2x=-getMobility(x+k1x/4,y+k1y/4,z+k1z/4, Electron)*electricFieldX*dt;
	k2y=-getMobility(x+k1x/4,y+k1y/4,z+k1z/4, Electron)*electricFieldY*dt;
	k2z=-getMobility(x+k1x/4,y+k1y/4,z+k1z/4, Electron)*electricFieldZ*dt;

	electricFieldZ = getEField1D(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z);

	k3x=-getMobility(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z, Electron)*electricFieldX*dt;
	k3y=-getMobility(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z, Electron)*electricFieldY*dt;
	k3z=-getMobility(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z, Electron)*electricFieldZ*dt;

	electricFieldZ = getEField1D(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z);

	k4x=-getMobility(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z, Electron)*electricFieldX*dt;
	k4y=-getMobility(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z, Electron)*electricFieldY*dt;
	k4z=-getMobility(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z, Electron)*electricFieldZ*dt;

	electricFieldZ = getEField1D(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z);

	k5x=-getMobility(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z, Electron)*electricFieldX*dt;
	k5y=-getMobility(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z, Electron)*electricFieldY*dt;
	k5z=-getMobility(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z, Electron)*electricFieldZ*dt;

	electricFieldZ = getEField1D(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
			y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
			z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z);

	k6x=-getMobility(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
			y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
			z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z, Electron)*electricFieldX*dt;
	k6y=-getMobility(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
			y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
			z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z, Electron)*electricFieldY*dt;
	k6z=-getMobility(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
			y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
			z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z, Electron)*electricFieldZ*dt;

	dx=((16./135)*k1x+(6656./12825)*k3x+(28561./56430)*k4x-(9./50)*k5x+(2./55)*k6x);
	dy=((16./135)*k1y+(6656./12825)*k3y+(28561./56430)*k4y-(9./50)*k5y+(2./55)*k6y);
	dz=((16./135)*k1z+(6656./12825)*k3z+(28561./56430)*k4z-(9./50)*k5z+(2./55)*k6z);

	double Ex,Ey,Ez,Erreur;
	Ex=((1./360)*k1x-(128./4275)*k3x-(2197./75240)*k4x-(1./50)*k5x+(2./55)*k6x);
	Ey=((1./360)*k1y-(128./4275)*k3y-(2197./75240)*k4y-(1./50)*k5y+(2./55)*k6y);
	Ez=((1./360)*k1z-(128./4275)*k3z-(2197./75240)*k4z-(1./50)*k5z+(2./55)*k6z);
	Erreur=sqrt(Ex*Ex+Ey*Ey+Ez*Ez);

	if (debug == DEBUG) G4cout << " [AllPixITkStripDigitizer::RKF5IntegrationElectrons] k6x : "<<k6x<<", k6y : " << k6y << ", k6z : " << k6z << G4endl;

	array<double,4> newpoint;
	newpoint.at(0)=(dx);
	newpoint.at(1)=(dy);
	newpoint.at(2)=(dz);
	newpoint.at(3)=(Erreur);

	return newpoint;

}
