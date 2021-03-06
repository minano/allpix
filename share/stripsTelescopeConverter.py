'''
Created on June 17, 2013

Tool to convert allpix telescope simulation files into LCIO EUTelescope format.

@author: <a href="mailto:samir.arfaoui@cern.ch">Samir Arfaoui</a>
'''

from pyLCIO import  UTIL
from ROOT import TVector3, TLorentzVector, TRandom3, TMath, std, TH1D, TCanvas

#Compiled with Cython if available, comment otherwise
#import pyximport; pyximport.install(pyimport=True)

from pyLCIO import EVENT,IMPL,IOIMPL
from time import time
from datetime import date
import sys, math, glob, tarfile
from collections import defaultdict

telescopeTestDict = {}
telescopeTestDict['300'] = 0
telescopeTestDict['301'] = 0
telescopeTestDict['302'] = 0
telescopeTestDict['303'] = 0
telescopeTestDict['304'] = 0
telescopeTestDict['305'] = 0

print sys.argv[3:]

mySensorIDlist = sys.argv[3:]
for mySensorID in mySensorIDlist:
    telescopeTestDict[mySensorID] = 0
    
stringDUTID='100008'
stringAPIXID='200'
#########################################################################################
def parseFrameFile( inputFile ):

    # open file
    #inputFile = open( fileName, 'r' )

    # each header should look like this, with different value of ChipboardID
    '''
    # Start time (string) : Jun 27 11:24:21.000 2013 # Start time : 1372335861.000 # Acq time : 0.000000 # ChipboardID : Chip_300 # DACs : 5 100 255 127 127 0 405 7 130 128 80 62 128 128 # Mpx type : 3 # Timepix clock : 40  # Eventnr 2 # RelaxD 2 devs 4 TPX DAQ = 0x110402 = START_HW STOP_HW MPX_PAR_READ COMPRESS=1
    492     231     1
    '''

    isNewDet = False
    sensorID = None
    telescopeData = defaultdict(list)
    DUTData       = defaultdict(list)
    APIXData      = defaultdict(list)

    # read contents
    while inputFile:
        
        line = inputFile.readline()
        s = line.split()
        n = len( s )
        if n == 0: break # end of file

        # read data for new detector
        if isNewDet and s[0] != '#':
            
            x = int( s[0] )
            y = int( s[1] )
            e = int( s[2] )
            d = int('0')

            if sensorID[0] == '3':
                #print 'Creating collection ', sensorID

                telescopeData[sensorID].append( x )
                telescopeData[sensorID].append( y )
                telescopeData[sensorID].append( e )
                telescopeData[sensorID].append( d )
                
                telescopeTestDict[sensorID] += 1

            if sensorID == stringAPIXID:
                #print 'Creating collection ', sensorID
                #print 20, x, y, e, d
                APIXData[sensorID].append(x)
                APIXData[sensorID].append(y)
                APIXData[sensorID].append(e)
                APIXData[sensorID].append(d)
				
                telescopeTestDict[sensorID] += 1
				
            if sensorID == stringDUTID:
                #print 'Creating collection ', sensorID
				
                DUTData[sensorID].append( x )
                DUTData[sensorID].append( y )
                DUTData[sensorID].append( e )
                DUTData[sensorID].append( d )
                
                telescopeTestDict[sensorID] += 1

        else:

            isNewDet = False            

        # check if new detector
        if s[0] == '#':
            
            isNewDet = True
            chipName = s[22] # format: Chip_XXX
            sensorID = chipName.split('_')[-1]

    return telescopeData, DUTData, APIXData

#########################################################################################
def convertRun( inputTarFile, outputFileName ):
    
    # read file names from given path
    #fileNames = sorted( glob.glob( inputPath+'/mpx*.txt' ) )
    #nEvents = len( fileNames )

    inputFiles = []
    tar = tarfile.open( inputTarFile, 'r:*' )
    count = int(0)
    for member in tar:
        if (count>1):
            inputFiles.append( tar.extractfile(member) )
            if (not (count%100)) : print "tar string: ", member
        count = (count+1)
        

    # get run number from first file (same for the rest) assuming filename is this form: mpx-YYMMDD-HHmmSS-RUN_FRAME.txt
    runNumber = int( inputFiles[0].name.split('-')[-1].split('_')[0] )
    runNumber = int( inputTarFile.split('n')[-1].split('.')[0] )
    #runNumber = 999999

    # define detector name
    detectorName = 'EUTelescope'

    # create a writer and open the output file
    writer = IOIMPL.LCFactory.getInstance().createLCWriter()
    writer.open( outputFileName, EVENT.LCIO.WRITE_NEW )

    # create a run header and add it to the file
    dateString = date.today().timetuple()
    run = IMPL.LCRunHeaderImpl()
    run.setRunNumber( runNumber )
    run.setDetectorName( detectorName )
    run.parameters().setValue  ( 'GeoID'            , 0 )
    run.parameters().setValues ( 'MaxX'             , std.vector(int)(6,1151) ) 
    run.parameters().setValues ( 'MaxY'             , std.vector(int)(6,575) )
    run.parameters().setValues ( 'MinX'             , std.vector(int)(6,0) ) 
    run.parameters().setValues ( 'MinY'             , std.vector(int)(6,0) )
    run.parameters().setValue  ( 'NoOfDetector'     , 6 )
    run.parameters().setValues ( 'AppliedProcessor' , std.vector('string')(1,'') )
    run.parameters().setValue  ( 'DAQHWName'        , 'EUDRB' )
    run.parameters().setValue  ( 'DAQSWName'        , 'EUDAQ' )
    run.parameters().setValue  ( 'DataType'         , 'SimData' )
    run.parameters().setValue  ( 'DateTime'         , '{}-{}-{}'.format(dateString[2], dateString[1], dateString[0]) )
    run.parameters().setValue  ( 'EUDRBDet'         , 'MIMOSA26' )
    run.parameters().setValue  ( 'EUDRBMode'        , 'ZS2' )
    writer.writeRunHeader( run )
  
    MAXEVENTS = 1000000
    NEVENTS   = 0

    # event loop
    for eventFile in inputFiles:

        if NEVENTS == MAXEVENTS: break
        NEVENTS += 1

        # get event number from file name ( i.e. frame ID ) assuming file name format above
        try:
            iEvent = int( eventFile.name.split('_')[-1].split('.')[0] )
        except ValueError:
            print 'Exception after ', NEVENTS, ' events.'
            print eventFile.name
            print 'Continuing to next event.'
            continue
			
        if ( NEVENTS%1000 == 0 ):
            print 'Events processed: %i ...' % NEVENTS

        # create an event and set its parameters
        event = IMPL.LCEventImpl()
        event.setEventNumber( iEvent )
        event.setDetectorName( detectorName )
        event.setRunNumber( runNumber )
        event.setTimeStamp( int( time() * 1000000000. ) )
        event.parameters().setValue( 'EventType', 2 )

        # parse input file
        telescopeData, DUTData, APIXData = parseFrameFile( eventFile )

        # is there DUT data?
        containsDUT = False
        if len( DUTData ) > 0: 
            containsDUT = True
            
        containsAPIX = False
        if len( APIXData ) > 0:
			containsAPIX = True 

        # if first event, create additional setup collection(s)
        if iEvent == 0:
            
            eudrbSetup = IMPL.LCCollectionVec( EVENT.LCIO.LCGENERICOBJECT )
            
            # collection parameters
            eudrbSetup.parameters().setValue( 'DataDescription', 'type:i,mode:i,spare1:i,spare2:i,spare3:i' )
            eudrbSetup.parameters().setValue( 'TypeName', 'Setup Description' )
            
            # create on setup object per Telescope plane
            for sensorID in sorted( telescopeData.iterkeys() ):
                print " create setup ", sensorID                    
                setupObj = IMPL.LCGenericObjectImpl(5,0,0)
                setupObj.setIntVal( 0, 102 )
                setupObj.setIntVal( 1, 101 )     
                eudrbSetup.addElement( setupObj )

            event.addCollection ( eudrbSetup, 'eudrbSetup' )

            # check if there is a DUT
            if containsDUT:
                print " Create DUT setup"
                DUTSetup = IMPL.LCCollectionVec( EVENT.LCIO.LCGENERICOBJECT )
                event.addCollection ( DUTSetup, 'DUTSetup' )
                
            if containsAPIX:
                print " Create APIX setup"				
                APIXSetup = IMPL.LCCollectionVec( EVENT.LCIO.LCGENERICOBJECT )
                event.addCollection ( APIXSetup, 'APIXSetup' )

        # ID encoder info
        encodingString = 'sensorID:7,sparsePixelType:5'

        # Telescope data collection
        trackerDataColl = IMPL.LCCollectionVec( EVENT.LCIO.TRACKERDATA )
        idEncoder_Telescope = UTIL.CellIDEncoder( IMPL.TrackerDataImpl )( encodingString, trackerDataColl )

        # check if there is a DUT
        if (containsDUT):

            # DUT data collection
            DUTDataColl = IMPL.LCCollectionVec( EVENT.LCIO.TRACKERDATA )
            idEncoder_DUT = UTIL.CellIDEncoder( IMPL.TrackerDataImpl )( encodingString, DUTDataColl )

            for i,sensorID in enumerate( sorted( DUTData.iterkeys() ) ):
              if (sensorID == stringDUTID):
             
                 planeData = IMPL.TrackerDataImpl()
            
                 idEncoder_DUT.reset()        
                 #idEncoder_DUT['sensorID'] = int( sensorID ) - 500 + 6 # cannot fit 500 in 5 bits!! FIXME
                 #idEncoder_DUT['sensorID'] = 8+i # cannot fit 500 in 5 bits!! FIXME
                 idEncoder_DUT['sensorID'] = int(sensorID) - 100000
                 idEncoder_DUT['sparsePixelType'] = 2
                 idEncoder_DUT.setCellID( planeData )
            
                 chargeVec = std.vector(float)()
                 for val in DUTData[sensorID]:
                     chargeVec.push_back( val )
                     if val < 0:
                         print 'Negative number in Event %i' % iEvent
                    
                 planeData.setChargeValues( chargeVec )

                 DUTDataColl.addElement( planeData )

            event.addCollection( DUTDataColl, 'zsdata_strip' )
                   
        if (containsAPIX):       
            # FEI4 data collection
             APIXDataColl = IMPL.LCCollectionVec( EVENT.LCIO.TRACKERDATA )
             idEncoder_APIX = UTIL.CellIDEncoder( IMPL.TrackerDataImpl ) ( encodingString, APIXDataColl )
             for i,sensorID in enumerate( sorted( APIXData.iterkeys() ) ):
               if (sensorID == stringAPIXID) :
				
                 planeData = IMPL.TrackerDataImpl()
				
                 idEncoder_APIX.reset()
                 idEncoder_APIX['sensorID']= 20
                 idEncoder_APIX['sparsePixelType'] = 2
                 idEncoder_APIX.setCellID( planeData )
                
                 chargeVec = std.vector(float)()
                 for val in APIXData[sensorID] :
                    chargeVec.push_back( val )
                    if val < 0:
                        print 'Negative number in Event %i' % iEvent
                 planeData.setChargeValues(chargeVec )
                
                 APIXDataColl.addElement(planeData)
            
             event.addCollection( APIXDataColl, 'zsdata_apix' )

        # fill telescope collection
        for sensorID in sorted( telescopeData.iterkeys() ):
            
            planeData = IMPL.TrackerDataImpl()

            idEncoder_Telescope.reset()
            idEncoder_Telescope['sensorID'] = int( sensorID ) - 300 # cannot fit 300 in 5 bits!! FIXME
            idEncoder_Telescope['sparsePixelType'] = 2
            idEncoder_Telescope.setCellID( planeData )
            
            # loop over hits
            chargeVec = std.vector(float)()
            for val in telescopeData[sensorID]:
                chargeVec.push_back( val )

            planeData.setChargeValues( chargeVec )

            trackerDataColl.addElement( planeData )

        event.addCollection( trackerDataColl, 'zsdata_m26' )

        writer.writeEvent( event )
    
    writer.flush()
    writer.close()

#########################################################################################
def usage():
    print 'Converts allpix generated telescope files into LCIO format'
    print 'Usage:\n python %s <inputTarball> <outputFile> <sensorID1> <sensorID2> <sensorIDn>' % ( sys.argv[0] )

#########################################################################################
if __name__ == '__main__':
    if len( sys.argv ) < 3:
        usage()
        sys.exit( 1 )
    convertRun( sys.argv[1], sys.argv[2] )

    for sensorID in sorted( telescopeTestDict.iterkeys() ):
        print sensorID, telescopeTestDict[sensorID]
