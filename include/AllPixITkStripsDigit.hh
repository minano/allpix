/**
 *  Author:
 *    Dino Tahirovic <dino.tahirovic@cern.ch>
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <benoit@lal.in2p3.fr>
 */

#ifndef AllPixITkStripsDigit_h
#define AllPixITkStripsDigit_h 1

#include "AllPixDigitInterface.hh"
#include "G4TDigiCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

/**
 *  Digit AllPixITkStrips class. Implementation of interface AllPixInterfaceDigit class.
 */
class AllPixITkStripsDigit : public AllPixDigitInterface {

public:

	/** Default empty constructor
	 *
	 */
  AllPixITkStripsDigit();

  /**
   *  Default constructor. Sets strip count to 1.
   *
   * @param xPixel Strip number
   * @param yPixel Row number
   * @param energy Deposited energy. should be ToT? FIXME
   * @param vertex ?
   */
  AllPixITkStripsDigit(G4int xPixel, G4int yPixel, G4double energy, G4ThreeVector vertex);

  ~AllPixITkStripsDigit();

  /** Constructor from another digit
   *
   * @param AllPixITkStripDigit Digit to copy content from
   */
  AllPixITkStripsDigit(const AllPixITkStripsDigit&);

  /** Assignment operator =
   *
   */
  const AllPixITkStripsDigit& operator=(const AllPixITkStripsDigit&);

  /** Comparison operator
   *
   */
  int operator==(const AllPixITkStripsDigit&) const;
  
  /** Reserve memory
   *
   */
  inline void* operator new(size_t);

  /** Free memory
   *
   */
  inline void  operator delete(void*);
  
  /** Draw digit - not implemented
   *
   */
  void Draw();

  /** Print details - not implemented
   *
   */
  void Print();

private:
  
  G4int m_pixelIDX; //! Pixel ID in x-direction
  G4int m_pixelIDY;   //! Pixel ID in y-direction
  G4int m_pixelCounts; //! Number of pixel counts in one event/frame
  G4double m_depositedEnergy; //! MC, Corrected MC charge (detector effects included, at Digi step)
  G4ThreeVector m_primaryVertex; //! MC Truth primary vertex from hit
  
public:
  
  inline void SetPixelIDX(G4int pidX)   {m_pixelIDX = pidX;};
  inline void SetPixelIDY(G4int pidY)   {m_pixelIDY = pidY;};
  inline void SetPixelCounts(G4int pc)  {m_pixelCounts = pc;};
  inline void SetPixelEnergyDep(G4double ed)  {m_depositedEnergy = ed;}; // MC // Corrected MC charge (detector effects included, at Digi step)
  inline void SetPrimaryVertex(G4ThreeVector pv)  {m_primaryVertex = pv;}; // MC vertex //
  inline void IncreasePixelCounts(void)     {m_pixelCounts++;};

  inline G4int GetPixelIDX() const  {return m_pixelIDX;};
  inline G4int GetPixelIDY() const  {return m_pixelIDY;};
  inline G4int GetPixelCounts() const {return m_pixelCounts;};
  inline G4double GetPixelEnergyDep() const {return m_depositedEnergy;}; // MC //
  inline G4ThreeVector GetPrimaryVertex() const {return m_primaryVertex;}; // MC //

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4TDigiCollection<AllPixITkStripsDigit> AllPixITkStripsDigitsCollection;

extern G4Allocator<AllPixITkStripsDigit> AllPixITkStripsDigitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* AllPixITkStripsDigit::operator new(size_t)
{
  void* aDigi;
  aDigi = (void*) AllPixITkStripsDigitAllocator.MallocSingle();
  return aDigi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void AllPixITkStripsDigit::operator delete(void* aDigi)
{
  AllPixITkStripsDigitAllocator.FreeSingle((AllPixITkStripsDigit*) aDigi);
}

#endif

