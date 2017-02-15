/**
 * Author:
 *    Dino Tahirovic <dino.tahirovic@cern.ch>
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <Mathieu.Benoit@cern.ch>
 */

#ifndef AllPixITkStripsDigitizer_h
#define AllPixITkStripsDigitizer_h 1


// allpix Interface
#include "AllPixDigitizerInterface.hh"
// digits for this digitizer
#include "AllPixITkStripsDigit.hh"
#include "G4PrimaryVertex.hh"

#include <map>
#include <vector>
#include <array>

using namespace std;

/**
 *  Digitizer AllPixITkStrip implementation
 */
class AllPixITkStripsDigitizer : public  AllPixDigitizerInterface {

public:

  //! Constructors
  AllPixITkStripsDigitizer(G4String, G4String, G4String);
  virtual ~AllPixITkStripsDigitizer();

  //! Constants
  static const double kPairEnergy; // in eV
  static const double fC;
  static const double SiPermittivity;
  static const double k_B;

  static const double kMinEnergy;
  static const double kMaxEnergy;

  // maximal number of iterations for num solving of Euler eq.
  static const double kIterations;


  // Enums
  enum CarrierType {Electron=1, Hole=0};
  typedef double ElectricField;

  //! Methods

  /** Calculate electric field in 1-dim
   *
   * @param x : x position in ??? co-ord system - where is the centre?
   * @param y : y -II-
   * @param z : position where e-h pair is created
   */
  ElectricField getEField1D(double x, double y, double z);

  double getENorm(double Ez);

  /**
   *
   * @param electricField
   * @param Temperature
   * @param isHoleBit
   * @return
   */
  G4double getMobility(double x, double y, double z, CarrierType carrier);

  /** Compute uniform electric field in Z direction
   *
   * @param carrier type electron or hole
   * @param z : position where e-h pair is created
   * @return drift time from the current position of the e/h
   */
  G4double getDriftTime(G4double z, CarrierType,  bool doFast);

  /** Integration step (time) for eq. of motion in E field
   *
   * @param dt : current step size
   * @param ErreurMoy : current error
   * @return New step size
   */
  G4double setDt(G4double dt, G4double ErreurMoy);

  array<double,4>  getDriftVector(G4double x, G4double y, G4double z);

  /**
   *
   * @param electricField
   * @param Mobility
   * @param isHoleBit
   * @return : drift velocity
   */
  G4double GetDriftVelocity(double x, double y, double z, CarrierType carrier);

  /**
   *
   * @param driftVelocity
   * @param isHoleBit
   * @return
   */
  G4double GetMeanFreePath(G4double driftVelocity, G4bool isHoleBit);

  /** Width of diffusion distribution sqrt(2Dt)
   *
   * @param tDrift drift time
   * @return Diffusion length in 3d (RMS)
   */
  G4double getDiffusionLRMS(G4double tDrift);

  /**
   * Returns the discrimination threshold in number of e-h pairs.
   * Internally, threshold is stored as energy needed to produce
   * a certain number of e-h pairs
   *
   * @return roughly 1 fC is 6250 pairs, or 22keV
   */
  inline int getThresholdPairs() {return m_inputParameter.thl/kPairEnergy;}

  //! Setters
  void SetPrimaryVertex(G4PrimaryVertex * pv) {m_primaryVertex = pv;};
  void Digitize ();
  void SetDetectorDigitInputs(G4double){};

  /** Set discriminator threshold. Internally threshold is in eV
   *
   * @param _thlCharge : threshold given in charge units, ex. 1.0*fC
   */
  inline void SetThreshold(G4double _thlCharge) {m_inputParameter.thl = _thlCharge/CLHEP::e_SI / coulomb * kPairEnergy;};

  /**
   *
   * @return threshold as energy in internal units
   */
  inline double getThresholdEnergy() { return m_inputParameter.thl;}

  /**
   *
   * @return
   */
  inline double getThresholdCharge() { return m_inputParameter.thl/kPairEnergy * CLHEP::e_SI * coulomb; }

  /** Set tempereature of sensor, stored in K
   *
   * @param _temperatureC Input is temp in Celsius
   */
  inline void SetTemperatureK(G4double _temperatureC){m_temperature = 273.15 + _temperatureC;};

  G4double MyErf(G4double x);
  G4double IntegrateGaussian(G4double xhit,G4double yhit,G4double Sigma, G4double x1, G4double x2, G4double y1, G4double y2, G4double Energy );

  array<double,4> RKF5IntegrationElectrons(G4double x, G4double y, G4double z,G4double dt);



private:

  // digitInput typedef is defined in AllPixDigitizerInterface.hh
  digitInput m_inputParameter; //! Very efficient struct that has only one field; .thl energy in eV

  AllPixITkStripsDigitsCollection * m_digitsCollection; //! Pointer to output collection of digitized hits
  vector<G4String> m_hitsColName; //! Input hits name from ?

  AllPixGeoDsc * m_geometry; // Handler to geometry from database (xml file)

  G4PrimaryVertex * m_primaryVertex; //! information from EventAction

  double m_detectorWidth; //! total detector thickness, mm
  double m_depletionVoltage; //! From tests
  double m_biasVoltage; //! Bias voltage, should be pass from xml gear file, V
  double m_Neff; //! Effective carrier concentration, 1/cm3
  double m_depletionDepth; //! theoretical depth of depleted zone from flat diode approx, mm

  double m_temperature; //! From gear file

  double m_trappingTimeElectrons; //! From gear file
  double m_trappingTimeHoles; //! From gear file
  //pair<CarrierType, int> m_diffusion;
  std::map<CarrierType, double> m_diffusion;

  double m_pitchX; //! Pitch between strips, needed for charge sharing
  int m_nStrips; //! Number of strips
  int m_nRows; //! Number of strip rows (mini 1, LS 2, SS 4)
  int m_pitchY; //! Vertical pitch between rows (2.5 mm?)

  double m_precision; //! Number of bunches for charge transportation

  bool m_doFast; //! 1 : Ballistic, 2 : numerically solve eq. of motion

  double m_tStepL; //! Lower step bound for time for RK5 integration
  double m_tStepU; //! Upper step bound for time for RK5 integration
  double m_maxError; //! target error for RK5 integration
  double m_stepSize; //! keep the initial step size for the next bunch propagation in RK5

  double m_minTime; //! earliest time for FE discriminator, example -25 ns
  double m_maxTime; //! latest time for FE discr., example 75 ns

  FILE* m_testFile; //! output of debugging information for plotting

};
#endif
