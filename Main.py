from Bio import Phylo
from io import StringIO
import numpy

#compares sequences and adds to matrix
def sequenceCompare():
    sequenceCJD = 'MANLGCWMLVLFVATWSDLGLCKKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQGGGTHSQWNKPSKPKTNMKHMAGAAAAGAVVGGLGGYVLGSAMSRPIIHFGSDYEDRYYRENMHRYPNQVYYRPMDEYSNQNNFVHNCVNITIKQHTVTTTTKGENFTETDVKMMERVVEQMCITQYERESQAYYQRGSSMVLFSSPPVILLISFLIFL'
    sequenceFFI = 'MANLGCWMLVLFVATWSDLGLCKKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQGGGTHSQWNKPSKPKTNMKHMAGAAAAGAVVGGLGGYMLGSAMSRPIIHFGSDYEDRYYRENMHRYPNQVYYRPMDEYSNQNNFVHNCVNITIKQHTVTTTTKGENFTETDVKMMERVVEQMCITQYERESQAYYQRGSSMVLFSSPPVILLISFLIFL'
    sequenceGSD = 'MANLGCWMLVLFVATWSDLGLCKKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQGGGTHSQWNKLSKPKTNMKHMAGAAAAGAVVGGLGGYMLGSAMSRPIIHFGSDYEDRYYRENMHRYPNQVYYRPMDEYSNQNNFVHDCVNITIKQHTVTTTTKGENFTETDVKMMERVVEQMCITQYERESQAYYQRGSSMVLFSSPPVILLISFLIFL'
    sequencePrion = 'MANLGCWMLVLFVATWSDLGLCKKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQGGGTHSQWNKPSKPKTNMKHMAGAAAAGAVVGGLGGYMLGSAMSRPIIHFGSDYEDRYYRENMHRYPNQVYYRPMDEYSNQNNFVHDCVNITIKQHTVTTTTKGENFTETDVKMMERVVEQMCITQYERESQAYYQRGSSMVLFSSPPVILLISFLIFL'

    #add sequences to array and get length
    sequence_list = [sequenceCJD, sequenceFFI, sequenceGSD, sequencePrion]
    n = len(sequence_list)

    #make a numpy array of size n x n
    upgmaArray = numpy.zeros((n,n))

    #fill array with pairwise distances
    upgmaArray = makeArray(sequence_list, upgmaArray)
    print("UPGMA table for sequences:")
    print(upgmaArray)

    #use the array to get the tree contents
    treeContent = computeTree(upgmaArray)

    #add data to file and draw to tree
    tree = open("tree.dnd","w+")
    tree.write(treeContent)
    drawTree(treeContent)

    #legend
    print("a-CJD Sequence")
    print("b-FFI Sequence")
    print("c-GSD Sequence")
    print("d-Normal Prion Sequence")

    #takes tree file and constructs a tree using biopython
def drawTree(treeFile):
    tree = Phylo.read(StringIO(treeFile),"newick")
    Phylo.draw_ascii(tree)

#makes a upgma array using pairwise distances between sequences
def makeArray(list, array):
    for i, element1 in enumerate(list):
        for j, element2 in enumerate(list):
            if j >= i:
                # since matrix is mirrored, no need to enumerate j past i
                break
            # calculate pairwise distance
            distance = pairwise(element1, element2)
            # add pairwise distances to array
            array[i, j] = distance
            array[j, i] = distance
    return array

#use upgma array to implement upgma algorithm
def computeTree(array):
    #create labels for use in Array and tree
    labels = []
    for i in range(ord("A"), ord("D")+1):
        labels.append(chr(i))

    #get lower triangular matrix
    arrayData = [
        [],
        [array[1][0]],
        [array[2][0],array[2][1]],
        [array[3][0],array[3][1],array[3][2]]
    ]

    #distance array for adding branch lengths
    distances = [0,0,0,0,0]

    #while labels are avaliable to be joined
    while len(labels) > 1:
        print("\nTable being analyzed:")
        for i in range(len(arrayData)):
            print(arrayData[i])
        #find lowest value
        x,y = lowValue(arrayData)

        #merges array entries for x,y by averaging data
        distancesx = (arrayData[x][y]/2)-distances[x]
        distancesy = (arrayData[x][y]/2)-distances[y]
        distances[x] = (arrayData[x][y]/2)
        distances[y] = (arrayData[x][y]/2)
        mergeArray(arrayData,x,y,labels,distances,distancesx,distancesy)

    return labels[0]


#returns pairwise distance for compared sequences
def pairwise(seq1, seq2):
    return sum(x != y for x, y in zip(seq1, seq2))

#finds location of lowest value in array
def lowValue(array):
    #initialize lowest value at infinity
    low = float("inf")
    x,y = -1,-1

    #iterate table for lowest value
    for i in range(len(array)):
        for j in range(len(array[i])):
            if array[i][j] < low:
                low = array[i][j]
                x,y = i,j
    return x,y


# merges array entries for x,y by averaging data
def mergeArray(array,x,y,labels,distance,distancex,distancey):
    #makes sure array is only worked on left side, swap if not
    if y<x:
        x,y = y,x
        tempy = distancey
        distancey = distancex
        distancex = tempy

    distancex = str(distancex)
    distancey = str(distancey)
    #create new label for merged data
    labels[x] = "(" + labels[x] + ":" + distancex + "," + labels[y] + ":" + distancey + ")"

    #recontruct the x row
    row = []
    for i in range(0,x):
        row.append((array[x][i] + array[y][i])/2)
    array[x] = row

    #reconstruct x column
    for i in range(x+1,y):
        array[i][x] = (array[i][x] + array[y][i])/2

    for i in range(y+1, len(array)):
        array[i][x] = (array[i][x] + array[i][y])/2
        #since data merged, delete leftover column
        del array[i][y]

    #since data merged, delete leftover row
    del array[y]

    #delete old label
    del labels[y]

    #delete old distance
    del distance[y]

if __name__ == '__main__':
    print("Comparing sequences of neurodegenerative prion diseases...\n")
    sequenceCompare()
