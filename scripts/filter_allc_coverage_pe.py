import sys, math, glob, multiprocessing, subprocess, os
from bioFiles import *

# Usage: python filter_allc_coverage_pe.py [-v=min_cov] [-num_proc] <allc_path> <sample1> [sampleN]*
# creates new allc files for all input samples that only includes information
# about positions which have at least minCov reads for each sample
# expects all chromosomes in one allC file

PICKLE=False
MINCOV=3
NUMPROC=1

def processInputs( allcPath, sampleNamesAr, minCov, numProc, isPrint ):
	if isPrint:
		print( 'AllC path:', allcPath )
		print( 'Samples:', ', '.join(sampleNamesAr) )
		print( 'Minimum coverage:', minCov )

	info = '#from_script: filter_allc_coverage_pe.py; min_cov: {:d}; samples_included: {:s}\n'.format( minCov, ','.join(sampleNamesAr) )

	# loop through all samples and get sets of positions with minCov
	if isPrint:
		print( 'Analyzing samples for positions covered with {:d} processes'.format( numProc) )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async(getCovPositionSample, args=(allcPath, sampleName, minCov) ) for sampleName in sampleNamesAr ]
	posDictAr = [p.get() for p in results]

	# combine by chromosome
	if isPrint:
		print( 'Combining covered positions' )
	posDict = combineCovSamples( posDictAr )
	for chrm in sorted(posDict.keys()):
		if isPrint:
			print( ' {:s}: {:d} positions'.format( chrm, len(posDict[chrm]) ))

	# loop through samples, directly read allC files and only write
	# lines included in posDict
	if isPrint:
		print( 'Writing new allC files' )
	results2 = [ pool.apply_async( writeAllcFile, args=(allcPath, sampleName, minCov, posDict, info, isPrint) ) for sampleName in sampleNamesAr ]
	sub = [p.get() for p in results2 ]
	if isPrint:
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

def writeAllcFile( allcPath, sampleName, minCov, posDict, info, isPrint ):
	inFileStr = os.path.normpath( '{:s}/allc_{:s}.tsv'.format( allcPath, sampleName, minCov ) )
	inFileObj = FileBio( inFileStr )
	inFile = inFileObj.fbOpen()
	outFileStr = inFileObj.fbBasename() + '_cov{:d}.tsv'.format( minCov )
	if isPrint:
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
	minCov = MINCOV
	numProc = NUMPROC
	isPrint = True
	startInd = 0

	for i in range(min(4,len(argv)-2)):
		if argv[i].startswith( '-v=' ):
			try:
				minCov = int( argv[i][3:] )

			except ValueError:
				print( 'WARNING: minimum coverage must be integer...using default', MINCOV )
				minCov = MINCOV
			startInd += 1
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
			except ValueError:
				print( 'WARNING: number of processors must be integer...using default', NUMPROC )
				numProc = NUMPROC
			startInd += 1
		elif argv[i] == '-q':
			isPrint = False
			startInd += 1
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
	processInputs( allcPath, sampleNamesAr, minCov, numProc, isPrint )

def printHelp():
	print( 'Usage:\tpython filter_allc_coverage_pe.py [-q] [-h] [-v=min_cov] <allc_path>\n\t<sample_name> [sample_name]*' )
	print()
	print( 'Required:' )
	print( 'allc_path\tpath to allC files' )
	print( 'sampleN\t\tname of sample; used to find allC files' )
	print()
	print( 'Optional:' )
	print( '-q\t\tquiet; do not print progress' )
	print( '-h\t\tprint help and exit' )
	print( '-v=min_cov\tmin coverage for positions to include [default {:d}]'.format( MINCOV) )
	print( '-p=num_proc\tnumber of processors to use [default {:d}]'.format( NUMPROC) )

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
