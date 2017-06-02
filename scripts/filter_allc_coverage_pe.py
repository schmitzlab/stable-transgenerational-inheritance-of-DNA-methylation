import sys, math, glob, multiprocessing, subprocess, os, bisect, random
from bioFiles import *

# Usage: python filter_allc_coverage_pe.py [-v=min_cov] [-num_proc] <allc_path> <sample1> [sampleN]*
# creates new allc files for all input samples that only includes information
# about positions which have at least minCov reads for each sample
# expects all chromosomes in one allC file

PICKLE=False

def processInputs( allcPath, sampleNamesAr, minCov, numProc ):
	print( 'AllC path:', allcPath )
	print( 'Samples:', ', '.join(sampleNamesAr) )
	print( 'Minimum coverage:', minCov )
	info = '#from_script: filter_allc_coverage_pe.py; min_cov: {:d}; samples_included: {:s}\n'.format( minCov, ','.join(sampleNamesAr) )
	
	# loop through all samples and get sets of positions with minCov
	print( 'Analyzing samples for positions covered with {:d} processes'.format( numProc) )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async(getCovPositionSample, args=(allcPath, sampleName, minCov) ) for sampleName in sampleNamesAr ]
	posDictAr = [p.get() for p in results]
	
	# combine by chromosome
	print( 'Combining covered positions' )
	posDict = combineCovSamples( posDictAr )
	for chrm in sorted(posDict.keys()):
		print( ' {:s}: {:d} positions'.format( chrm, len(posDict[chrm]) ))
	
	# loop through samples, directly read allC files and only write
	# lines included in posDict
	print( 'Writing new allC files' )
	results2 = [ pool.apply_async( writeAllcFile, args=(allcPath, sampleName, minCov, posDict, info) ) for sampleName in sampleNamesAr ]
	sub = [p.get() for p in results2 ]
	print( 'Done' )

def getCovPositionSample( allcPath, sampleName, minCov ):
	mFile = FileAllC_full( os.path.normpath( '{:s}/allc_{:s}.tsv'.format( allcPath, sampleName ) ) )
	mDict = mFile.getAllCDict( mtypes = 'C', isPickle=PICKLE )
	outDict = {}
	# loop through chromosomes
	for chrm in mDict.keys():
		tmpAr = []
		for pos in mDict[chrm].keys():
			# check coverage
			if mDict[chrm][pos][1] >= minCov:
				tmpAr += [ pos ]
		# end for pos
		# add to outDict as a set
		outDict[chrm] = set( tmpAr )
	# end for chrm
	return outDict

def combineCovSamples( posDictAr ):
	# posDictAr is array of dictionaries with chrms
	chrmList = list( posDictAr[0].keys() )
	outDict = {}
	# loop through chrm
	for chrm in chrmList:
		# start with set 0 and intersection remaining
		cSet = posDictAr[0][chrm]
		for i in range(1, len(posDictAr) ):
			newSet = posDictAr[i][chrm]
			cSet = cSet & newSet
		# end for i
		outDict[chrm] = cSet
	return outDict
	
def writeAllcFile( allcPath, sampleName, minCov, posDict, info ):
	inFileStr = os.path.normpath( '{:s}/allc_{:s}.tsv'.format( allcPath, sampleName, minCov ) )
	inFileObj = FileBio( inFileStr )
	inFile = inFileObj.fbOpen()
	outFileStr = inFileObj.fbBasename() + '_cov{:d}.tsv'.format( minCov )
	print( 'Writing output to', os.path.basename( outFileStr ) )
	outFile = open( outFileStr, 'w' )
	outFile.write( info )
	curChrm = None
	curChrmSet = None
	
	for line in inFile:
		if line.startswith( '#' ):
			outFile.write( line )
			continue
		lineAr = line.rstrip().split( '\t' )
		if len( lineAr ) < 7 or lineAr[6].isdigit() == False:
			continue
		if lineAr[0] != curChrm:
			curChrm = lineAr[0]
			curChrmSet = posDict[curChrm]
		pos = int( lineAr[1] )
		if pos in curChrmSet:
			outFile.write( line )
	# end for line
	inFile.close()
	outFile.close()

def parseInputs( argv ):
	minCov = 3
	numProc = 1
	startInd = 0
	
	for i in range(min(2,len(argv)-2)):
		if argv[i].startswith( '-v=' ):
			try:
				minCov = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: minimum coverage must be integer' )
				exit()
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of processors must be integer' )
				exit()
		elif argv[i] in [ '-h', '--help', '-help']:
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			exit()
	# end for
	allcPath = argv[startInd]
	if os.path.isdir( allcPath ) == False:
		print( 'ERROR: {:s} is not a path to a directory for allC files'.format( allcPath ) )
		exit()
	sampleNamesAr = []
	for i in range(startInd+1, len(argv)):
		sampleNamesAr += [ argv[i] ]
	processInputs( allcPath, sampleNamesAr, minCov, numProc )

def printHelp():
	print( 'Usage: python filter_allc_coverage.py [-v=min_cov] <allc_path> <sample1> [sampleN]*' )
	print( 'Required:' )
	print( 'allc_path\tpath to allC files' )
	print( 'sampleN\tname of sample; used to find allC files' )
	print( 'Optional:' )
	print( '-v=min_cov\tmin coverage for positions to include [default 3]' )
	print( '-p=num_proc\tnumber of processors to use [default 1]' )

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
