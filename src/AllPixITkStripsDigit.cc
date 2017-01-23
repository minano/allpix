/**
 *  Author:
 *    Dino Tahirovic <dino.tahirovic@cern.ch>
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <Mathieu.Benoit@cern.ch>
 */

#include "AllPixITkStripsDigit.hh"

G4Allocator<AllPixITkStripsDigit> AllPixITkStripsDigitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AllPixITkStripsDigit::AllPixITkStripsDigit() :
	m_pixelIDX(-1),
	m_pixelIDY(-1),
	m_pixelCounts(-1),
	m_depositedEnergy(-1),
	m_primaryVertex(G4ThreeVector())
	{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AllPixITkStripsDigit::AllPixITkStripsDigit(G4int xPixel, G4int yPixel, G4double energy, G4ThreeVector vertex) :
	m_pixelIDX(xPixel),
	m_pixelIDY(yPixel),
	m_pixelCounts(1),
	m_depositedEnergy(energy),
	m_primaryVertex(vertex)
{
};

AllPixITkStripsDigit::~AllPixITkStripsDigit()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AllPixITkStripsDigit::AllPixITkStripsDigit(const AllPixITkStripsDigit& right)
: AllPixDigitInterface()
{
	m_pixelIDX = right.m_pixelIDX;
	m_pixelIDY = right.m_pixelIDY;
	m_pixelCounts = right.m_pixelCounts;
	m_depositedEnergy = right.m_depositedEnergy;
	m_primaryVertex = right.m_primaryVertex;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const AllPixITkStripsDigit& AllPixITkStripsDigit::operator=(const AllPixITkStripsDigit& right)
{

	m_pixelIDX = right.m_pixelIDX;
	m_pixelIDY = right.m_pixelIDY;
	m_pixelCounts = right.m_pixelCounts;
	m_depositedEnergy = right.m_depositedEnergy;
    m_primaryVertex = right.m_primaryVertex;
	return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int AllPixITkStripsDigit::operator==(const AllPixITkStripsDigit& right) const
		{
	return ((m_pixelIDX==right.m_pixelIDX)&&
			(m_pixelIDY==right.m_pixelIDY)&&
			(m_pixelCounts==right.m_pixelCounts) &&
			(m_depositedEnergy == right.m_depositedEnergy) &&
			(m_primaryVertex == right.m_primaryVertex));
		}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void AllPixITkStripsDigit::Draw()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void AllPixITkStripsDigit::Print()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
