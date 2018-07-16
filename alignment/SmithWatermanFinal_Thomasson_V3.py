'''
Created on Feb 6, 2018

@author: Caley Thomasson
Smith Waterman
(with arg parser)

Finds alignment and max score for two protein sequences
 
data structure: dictionary for reading in blosum matrix txt file
dictionary used to store location, score and direction for scoring matrix
matrix created for scoring matrix
trace back from max score to find alignment

'''
import argparse
import numpy as np


def populate_matrix(matrixFile):
    '''
        reads in blosum matrix txt file and sets key: value pair as tuple(i, j) = score
    '''
    lines = matrixFile.readlines()
    matrixFile.close()
    dictaa = {}
    aminoacidstring = lines[0]
    aminoacidstring = aminoacidstring.split()
    
    i = 1
    while i <= (len(lines)-1):
        row = lines[i]
        row = row.split()

        j = 1
        for character in row[1:25]:
            dictaa[aminoacidstring[i-1],aminoacidstring[j-1]] = character # i,j changes (row, column) to (aa at row, aa at column) for keys
            j+=1
        i+=1

    return(dictaa) 



def scoringMatrix(seq1, seq2, gap_penalty, blosum_dict):
    '''
    dictionary used to store location, score, direction for traceback   
    calculates max score of three directions or 0, records in matrix.
    traceback starts at max score in matrix and traces where score came from until it hits a 0
    gets protein alignment
    calls toString()
    ''' 
    shortestseq = seq1
    longestseq = seq2
    smatrixToDict = {}
    direction = 0 
    runningMaxScore = 0
    maxLocation = 0,0
    gap_penalty = int(gap_penalty)
    
#uses shorter sequence as rows and longer sequence as columns to build matrix
    if(len(seq1) == len(seq2)):
        shortestseq = seq1
        longestseq = seq2
    if(len(seq1) < len(seq2)):
        shortestseq = seq1
        longestseq = seq2
    else:
        shortestseq = seq2
        longestseq = seq1 

    rows = len(shortestseq) 
    columns = len(longestseq)
    smatrix = np.zeros((len(shortestseq)+2, len(longestseq)+2), dtype = np.object)
    
    toprow = smatrix[0]
    longestseq = " "+ " " + longestseq
    shortestseq = " "+ shortestseq
    for i in range(len(longestseq)):
        toprow[i] = longestseq[i] 

    for k in range(len(shortestseq)): 
        smatrix[k+1][0] = shortestseq[k]

    for rowi in range(2, len(shortestseq)+1): 
        
        for colj in range(2, len(longestseq)):  

            blosumScoreKey = (smatrix[rowi][0],smatrix[0][colj])

            smatrix
            scorefromabove = int(smatrix[rowi-1][colj]) +(gap_penalty)

            scorefromleft = int(smatrix[rowi][colj-1]) +(gap_penalty)
            
            scorefromdiag= (smatrix[rowi-1][colj-1]) + int(blosum_dict.get(blosumScoreKey))
            zero = 0
            maxof4scores = max(scorefromabove, scorefromleft, scorefromdiag, zero)

            smatrix[rowi][colj] = maxof4scores

            if(maxof4scores == scorefromabove): # gap goes in seq1 horiz aa
                direction = 1
            if(maxof4scores == scorefromleft): # gap goes in seq2 vert aa 
                direction = 3
            if(maxof4scores == scorefromdiag):
                direction = 2
            if(maxof4scores == 0):
                direction = 0
                
            smatrixToDict[rowi, colj] = direction, maxof4scores 
                    
            if(maxof4scores > runningMaxScore):
                runningMaxScore = max(runningMaxScore, maxof4scores)     
                maxLocation = rowi, colj    
 
    maxScore = smatrix[maxLocation]           
    location = maxLocation       

#TRACEBACK

    directionANDscore = smatrixToDict[location]
    direction, score = (directionANDscore[0], directionANDscore[1])
    shortestSeqString = ""
    longestSeqString = ""

    for r in range(len(shortestseq)):
        for c in range(len(longestseq)):
     
     
            while score >= 1 and direction != 0:      
                if(direction == 1):  #score from above
                    newlocation = ((location[0])-1, location[1])  #location = [row, column] of score
                    longestSeqString += "-"
                    shortestSeqString += smatrix[location[0]][0]  #smatrix[locationtuple(row, col)][0]
                                                         
                if(direction == 3):    #score from left
                    shortestSeqString += "-"   
                    newlocation = (location[0], (location[1])-1)
                    longestSeqString += smatrix[0][location[1]]             

                if(direction == 2): #score from diagonal                      
                    newlocation = ((location[0])-1, (location[1])-1)
                    shortestSeqString += smatrix[location[0]][0]                    
                    longestSeqString += smatrix[0][location[1]]  

                location = newlocation
                directionANDscore = smatrixToDict[location]
                score = (directionANDscore[1]) 
                direction = directionANDscore[0]
   
    print(toString(longestSeqString, shortestSeqString, maxScore, smatrix, smatrixToDict))

   
def toString(longestSeqString, shortestSeqString, maxScore, smatrix, smatrixToDict):    
    longestSeqString = longestSeqString[::-1]
    shortestSeqString = shortestSeqString[::-1]
    
    addbar = ""
    for i in range(len(longestSeqString)):
        if shortestSeqString[i] == longestSeqString[i]:
            addbar += "|"

        else: 
            addbar += " "  
 
    print("Scoring matrix protein sequences:  \n ",smatrix)
    print("Dictionary for Scoring Matrix Traceback \n {(i,j): (direction the score was calculated from, score)} \n",smatrixToDict)
    
    print("Using Gap penalty:", gap_penalty , "MAX SCORE = " , maxScore) 

    print("Alignment: ")
    print("*The longest sequence read in is the top sequence")
    print("The shortest sequence read in is the bottom sequence*")
    print(longestSeqString)   
    print(addbar)
    print(shortestSeqString)    

def getSequence(sequenceFile):
    with open(sequenceFile, 'r') as seqFile:
        seqFileLines = seqFile.readlines()
        
    seqInfoDict = {'ID': None, 'seq': ' '}
    
    for i in range(0, len(seqFileLines)):
        seqFileLines[i] = seqFileLines[i].strip('\n')
        line = seqFileLines[i]
        if line[0] == '>':
            seqInfoDict['ID'] = line.strip('>')
        else:
            seqInfoDict['seq'] += line
 
    return seqInfoDict['seq'] #return the seq

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Count k-mar frequency based on given nucleotide sequence and k')
    parser.add_argument('-s1', dest ='seq_file1', required = True, help='sequence file name(nucleotides only, no headers)')
    parser.add_argument('-s2', dest = 'seq_file2')
    parser.add_argument('-m', dest = 'matrixFile')
    parser.add_argument('-g', dest = 'gap_penalty')
    args = parser.parse_args()
    matrixFile = args.matrixFile

with open(args.seq_file1, 'r') as seqFile:
    seqFileLines = seqFile.readlines()    
    seqInfoDict = {'ID': None, 'seq': ''}
    for i in range(0, len(seqFileLines)):
        seqFileLines[i] = seqFileLines[i].strip('\n')
        line = seqFileLines[i]
        if line[0] == '>':
            seqInfoDict['ID'] = line.strip('>')
        else:
            seqInfoDict['seq'] += line
    seq1 = seqInfoDict['seq'] #return the seq

with open(args.seq_file2, 'r') as seqFile2:
    seqFileLines2 = seqFile2.readlines()    
    seqInfoDict2 = {'ID': None, 'seq': ''}
    for i in range(0, len(seqFileLines2)):
        seqFileLines2[i] = seqFileLines2[i].strip('\n')
        line2 = seqFileLines2[i]
        if line2[0] == '>':
            seqInfoDict2['ID'] = line2.strip('>')
        else:
            seqInfoDict2['seq'] += line2
 
    seq2 = seqInfoDict2['seq'] #return the seq

with open(args.matrixFile, 'r') as matrixFile:
    blosum_dict = populate_matrix(matrixFile)  

gap_penalty = args.gap_penalty
      
scoringMatrix(seq1, seq2, gap_penalty, blosum_dict)

