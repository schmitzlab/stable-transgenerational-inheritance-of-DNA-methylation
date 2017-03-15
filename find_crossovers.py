import sys, math, glob, multiprocessing, subprocess, os, bisect, random
import pandas as pd

# Usage: python find_crossovers.py [-c=prediction_column] [-o=out_id] <input_file>

def processInputs(inFileStr, outID, predCol):
	
	info = '#from_script: find_crossovers.py; in_file: {:s}; pred_col: {:s}'.format( os.path.basename( inFileStr ), predCol )
	
	inTable = pd.read_table(inFileStr, skiprows=1)
	if predCol not in list(inTable):
		print( 'ERROR: input file does not contain predicition column {:s}...exiting'.format( predCol ) )
		exit()
	df = inTable[['bin', 'sample', predCol]]
	
	print( 'Finding breakpoints' )
	dfGroup = df.groupby( 'sample',sort=False )
	dfBreakpoints = dfGroup.apply( processSample, predCol )
	
	print( 'Counting breakpoints' )
	dfGroup2 = dfBreakpoints.groupby(level=0)
	dfCountBP = pd.DataFrame(dfGroup2['count'].agg( lambda x: max(x) - 1 ))
	bpTotal = dfCountBP.sum()
	#dfCountBP.reset_index(inplace=True)
	#dfCountBP.loc['total','count'] = bpTotal
	
	if outID == None:
		outID = os.path.basename( inFileStr ).replace('.tsv','').replace('.gz','').replace('updated_epigenotype', 'breakpoints' )
	outFileStr1 = outID + '_pos.tsv'
	outFileStr2 = outID + '_count.tsv'
	
	print( 'Writing output to {:s}_*'.format( outID ) )
	with open( outFileStr1, 'w' ) as f:
		f.write(info+'\n')
	dfBreakpoints.to_csv( outFileStr1, sep='\t', mode='a' )
	
	with open( outFileStr2, 'w' ) as f:
		f.write(info+'\n')
	dfCountBP.to_csv( outFileStr2, sep='\t', mode='a' )
	print( 'Done' )
	
def processSample( df, predCol ):
	outDf = pd.DataFrame(columns=['count','bin','from','to', 'length'])
	prevPred = 'start'
	count = 0
	length = -1
	bin = -1
	for row in df.itertuples():
		# (0) index (1) bin (2) sample (3) predCol
		bin = row[1]
		newPred = row[3]
		length += 1
		# is a change
		if newPred != prevPred:
			#print( row,prevPred )
			outDf = outDf.append( {'count':count, 'bin':bin, 'from':prevPred, 'to':newPred, 'length':length}, ignore_index=True )
			count += 1
			length = 0
		prevPred = newPred
	# end for
	# add end
	#print( count, bin, prevPred, 'end', length )
	outDf = outDf.append( {'count':count, 'bin':bin, 'from':prevPred, 'to':'end', 'length':length}, ignore_index=True )
	return outDf
	

def parseInputs( argv ):
	outID = None
	predCol = 'vit.prediction'
	startInd = 0
	
	for i in range(min(2,len(argv)-1)):
		if argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-c=' ):
			predCol = argv[i][3:]
			startInd += 1
	# end for i
	inFileStr = argv[startInd]
	processInputs( inFileStr, outID, predCol )


if __name__ == "__main__":
	if len(sys.argv) < 2 :
		print ("python find_crossovers.py [-c=prediction_column] [-o=out_id] <input_file>")
	else:
		parseInputs( sys.argv[1:] )
