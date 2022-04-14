
#!/usr/bin/python
__author__ = "FirstName LastName"
__email__ = "first.last@yale.edu"
__copyright__ = "Copyright 2021"
__license__ = "GPL"
__version__ = "1.0.0"

### Usage: python hw1.py -i <input file> -s <score file>
### Example: python hw1.py -i input.txt -s blosum62.txt
### Note: Smith-Waterman Algorithm

import argparse, re
import numpy as np

### This is one way to read in arguments in Python. 
parser = argparse.ArgumentParser(description='Smith-Waterman Algorithm')
parser.add_argument('-i', '--input', help='input file', required=True)
parser.add_argument('-s', '--score', help='score file', required=True)
parser.add_argument('-o', '--opengap', help='open gap', required=False, default=-2)
parser.add_argument('-e', '--extgap', help='extension gap', required=False, default=-1)

args = parser.parse_args()


### Function reads the input file and splits the lines into the sequneces
def readInputSeq(file):
    with open(file, 'r') as f:
        lines = f.read().splitlines()
    
    seq1 = lines[0]
    seq2 = lines[1]

    f.close()

    return seq1,seq2


### Read scores from the blosum62.txt file into a dictionary so that we have O(1) lookup time. 
### Example Dictionary Entry: mat[('A', 'A')] = 4. 
def readScoreMat(file):
    mat = {}
    with open(file, 'r') as f:
        # read the amnio acid in the y column
        y_acids = re.split('[ \t]+', f.readline().strip())
        for row in f:
            # read values from matrix
            scores = re.split('[ \t]+', row.strip())
            
            #skip values of zero 
            if len(scores) == 0:
                continue 

            x_acid = scores[0] #first value in array it the amino acid

            for i in range(1, len(scores)):
                y_acid = y_acids[i-1]
                score = int(scores[i])

                #place tuple as key and score as value
                mat[(x_acid, y_acid)] = score
    
    f.close() 

    return mat

def tracebackFunction(score_mat, traceback_mat, bestscore, best_idx, seq1, seq2):

	#extract the row/col producing the best score
	bestrow = best_idx[0]
	bestcol = best_idx[1]

	#initialize arrays for the sequences to print
	#we have to do this and then convert to string because strings aren't mutable
	topvals = ["(",")"] #first sequence corresponding to the matrix rows
	matching = [" "," "] #the "|"'s for matching values
	bottomvals = ["(",")"] #second sequence corresponding to the matrix columns

	#vars to hold current col/row, just so i don't overwrite bestrow and bestcol
	row = bestrow
	col = bestcol

	#row and col are currently at the end of the alignment so we can add everything after it to the alignment arrays
	topvals.insert(1,seq2[bestrow-1])
	bottomvals.insert(1,seq1[bestcol-1])
    	#run the traceback -- until you hit a cell with score 0 keep going

	while score_mat[row][col] != 0: 
	
		#the [0] value will be that of the current cell coords so choose the one after it
		#I know this is inefficient, it's just how I wrote the traceback matrix
		[r,c] = traceback_mat[row][col][1]
		#diagonal
		if r == row-1 and c == col-1: 
			
			#extract the cell producing the [r,c] score
			#we're working right->left on the alignment so each new char gets placed immediately after the "("
			topvals.insert(1,seq2[r])
			bottomvals.insert(1,seq1[c])

			#check to see if these values are equivalent, if so make a note in the matching array
			if seq2[r] == seq1[c]: 
				matching.insert(1,"|")	
			else: 
				matching.insert(1," ")
		
		#horizontal gap 
		elif r == row and c < col:
			
			#if i'm at cell [1,5] and the traceback produces cell [1,3] i need to indicate that i have a gap of length 2
			#I tried range(col,c) and (col-1,c-1) worked ... weirdly
			#we have to count backwards because we're starting from the end of the alignment--hence the -1
			for i in range(col-1,c-1,-1): 
				topvals.insert(1,'-')
				bottomvals.insert(1,seq1[i])
				matching.insert(1," ")

		#vertical gap--same idea as horiz gap
		elif r < row and c == col: 
			#again we have to count backwards
			for j in range(row-1,r-1,-1):
				topvals.insert(1,seq2[j])
				bottomvals.insert(1,'-')
				matching.insert(1," ")

		row = r
		col = c


def calcGapPenalty(open_pen, ext_pen, length):
    return open_pen + (ext_pen * (length - 1))

HORZ,VERT,DIAG = range(3)
### Implement your Smith-Waterman Algorithm
def runSW(inputFile, scoreFile, openGap, extGap):
    
    ### STEP 1: Read input file into sequnces
    seq1, seq2 = readInputSeq(inputFile)

    ### STEP 2: Read score file into matrix
    sim_matrix = readScoreMat(scoreFile)
    ### STEP 3: Calculate scores 

    #initalize score matrix and best values
    score_matrix = np.zeros((len(seq2)+1,len(seq1)+1))
    bestscore = 0
    best_idx = [0,0]

    #initalize traceback matrix
    tb_matrix = [[[[i,j]] for j in range(len(seq1)+1)] for i in range(len(seq2)+1)]

    for i in range(1, len(seq2)+1):
        for j in range(1, len(seq1)+1):
            
            # determine the match/mismatch score 
            diag = score_matrix[i-1][j-1] + sim_matrix[(seq2[i-1], seq1[j-1])]

            # determine vertical score
            vert_score = score_matrix[:i, j] + [calcGapPenalty(openGap, extGap, l) for l in range(i, 0, -1)]
            if len(vert_score) != 0:
                max_vert_score = max(score_matrix[:i, j] + [calcGapPenalty(openGap, extGap, l) for l in range(i, 0, -1)])
            

            #determine horizontal gap score

            horz_score = score_matrix[i, :j] + [calcGapPenalty(openGap, extGap, l) for l in range(j, 0, -1)]
            if len(horz_score) != 0:
                max_horz_score = max(score_matrix[i, :j] + [calcGapPenalty(openGap, extGap, l) for l in range(j, 0, -1)])
            
            
            # place max (diag, horizontal gap, vertical gap) into score matrix 
            score = int(max(diag, max_horz_score, max_vert_score, 0))
            score_matrix[i][j] = score

            # update new best score
            if score > bestscore:
                bestscore = score
                best_idx = [i,j]

            
            # STEP 4: add values to traceback matrix so that we can track what the best alignment is
             ## We must determine which (horizontal, vertical, or diag) was the best score 
             ## and use that to add the values to the traceback matrix
             ## for horizontal and vertical we ensure that we add all of the places the max score is
             ## found -- may not be required however
            if score == diag: 
                tb_matrix[i][j].append(([i-1,j-1],DIAG))

            #highest score came from horizontal gap
            elif score == max_horz_score: 
                values = np.where(horz_score==score)[0]
                for v in values:
                    tb_matrix[i][j].append(([i,v], HORZ))

            #highest score came from vertical gap
            elif score == max_vert_score: 
                values = np.where(vert_score==score)[0]
                for v in values:
                    tb_matrix[i][j].append(([v,j], VERT))
    

    # STEP 5: Calculate the traceback and find the best alignment
    bestrow = best_idx[0]
    bestcol = best_idx[1]

    #initialize arrays for the sequences to print
    string1 = ["(",")"]
    string2 = ["(",")"] 
    matching = [" "," "]


    row = bestrow
    col = bestcol

    # add all of the values after the alignment has ended 
    string1.insert(1,seq1[bestcol-1])
    string2.insert(1,seq2[bestrow-1])

    #run the traceback -- until you hit a cell with score 0 keep going
    while score_matrix[row][col] != 0: 


        ## find the traceback value and the label associated with 
        [r,c], VALUE = tb_matrix[row][col][1]

        ## update alignment strings based on the label 

        if VALUE == DIAG: 
            string1.insert(1,seq1[c])
            string2.insert(1,seq2[r])
            if seq2[r] == seq1[c]: 
                matching.insert(1,"|")	
            else: 
                matching.insert(1," ")
        elif VALUE == HORZ:
            for i in range(col-1,c-1,-1):
                string2.insert(1,'-') 
                string1.insert(1,seq1[i])
                matching.insert(1," ")
        elif VALUE == VERT: 
            for j in range(row-1,r-1,-1):
                string1.insert(1,'-')
                string2.insert(1,seq2[j])
                matching.insert(1," ")

        #update row and col to move to the next step
        row = r
        col = c
    
    ## we over-copy twice consistently 
    string1.pop(-2)
    string2.pop(-2)


    # start portion and end portions 
    seq2start = seq2[:row]
    seq1start = seq1[:col]
    seq1end = seq1[bestcol:]	
    seq2end = seq2[bestrow:]

    # find the values to append and do math for both seqs 
    toAppend = len(max([seq2start,seq1start],key=len))

    seq2len = len(seq2start)
    seq1len = len(seq1start)
    seq2gap = toAppend-seq2len
    seq1gap = toAppend-seq1len

    START1 = seq1gap*" " + seq1start
    START2 = seq2gap*" " +seq2start
    matchingStart = toAppend*" "

    sONE = START1+"".join(string1)+seq1end
    sTWO = START2+"".join(string2)+seq2end
    sMATCH = matchingStart+"".join(matching)

    #adding trailing zeros if needed
    maxLen = len(max([sONE,sMATCH,sTWO],key=len))
    if len(sONE) < maxLen:
        diff = maxLen - len(sONE)
        sONE += " "*diff
    if len(sTWO) < maxLen:
        diff = maxLen - len(sTWO)
        sTWO += " "*diff 
    if len(sMATCH) < maxLen:
        diff = maxLen - len(sMATCH)
        sMATCH += " "*diff

    ### STEP 6: Print the results to the output file

    outfile = open('output.txt','w')
    outfile.write("-----------\n")
    outfile.write("|Sequences|\n")
    outfile.write("-----------\n")
    outfile.write("sequence1\n")
    outfile.write(seq1+"\n")
    outfile.write("sequence2\n")
    outfile.write(seq2+"\n")
    outfile.write("--------------\n")
    outfile.write("|Score Matrix|\n")
    outfile.write("--------------\n")

    #write the scoring matrix

    # determine the row and column names
    colNames = [char for char in seq1] 
    rowNames = [char for char in seq2] 	
    colNames.insert(0, "")
    rowNames.insert(0,"")

    outfile.write("\t")
    for char in colNames: 
        outfile.write(char+"\t")
    outfile.write("\n")

    for i in range(len(rowNames)):
        outfile.write(rowNames[i]+"\t")
        for j in range(len(colNames)):
            outfile.write(str(int(score_matrix[i][j]))+"\t")
        outfile.write("\n")		

    outfile.write("----------------------\n")
    outfile.write("|Best Local Alignment|\n")
    outfile.write("----------------------\n")
    outfile.write("Alignment Score:"+str(bestscore)+"\n")
    outfile.write("Alignment Results:\n")

   
    #write the traceback
    outfile.write(sONE+"\n")
    outfile.write(sMATCH+"\n")
    outfile.write(sTWO+"\n")

    outfile.close()	


### Run your Smith-Waterman Algorithm
runSW(args.input, args.score, args.opengap, args.extgap)
