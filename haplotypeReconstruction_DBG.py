
from helperClasses import Read
import os
from helperClasses import FilteredRead
from stagewiseregression import StagewiseRegression
import numpy as np

"""
This class was used to generate the global haplotypes. 
"""
class HaplotypeReconstruction_DBG:

    def __init__(self):
        self.clusteringTable = []
        self.referenceGenome = ""
        self.variantPositions = []
        self.readList = []
        self.pairedEndReadList = []
        self.successors = {}
        self.badSuccessor = set()
        self.haplotypes = []
        self.specialList = []
        self.rawPositionList = []
        self.positionList = []
        self.haplotypeList = []
        self.nucleotides = []
        self.filteredPositions = []
        self.readTable = []
        self.list = []
        self.edges = {}

    """
    This method generates the clusteringTable by using the corrected reads from 
    the file ClusteringFile.txt
    """
    def createClusteringTable(self):
        file = open('Documents\\ClusteringFile.txt', 'r')
        lines = file.readlines()
        file.close()
        count = 0
        for line in lines:
            parts = line.split(' ')
            if len(parts) > 1:
                parts[2] = parts[2].replace(" ", "")
                parts[2] = parts[2].replace("\n", "")
                read = Read(count, parts[2], int(parts[0]), int(parts[1]))
                self.clusteringTable.append(read)
                count = count + 1

    """
    This method was used to generate the clusteringTable for Experiments. 
    The individual mentioned files were generated using ART.
    @folder Location where the files are
    """
    def createClusteringTableFromSamFile(self, folder):
        for i in range(0, 3):
            # The simulated files should be stored in these files
            if i == 0:
                file = open(folder+'\\paired_dat_errFree1.sam', 'r')
            elif i == 1:
                file = open(folder+'\\paired_dat_errFree1.sam', 'r')
            else:
                file = open(folder+'\\paired_dat_errFree1.sam', 'r')
            lines = file.readlines()
            file.close()
            count = 0
            rawList = []
            for line in lines:
                count += 1
                if count < 4:
                    continue
                spLine = line.split(" ")
                elements = spLine[0].split("\t")
                rawList.append([elements[3], elements[9], elements[0]])
            count = 0
            t = 0
            for i in range(0, len(rawList)-1, 2):
                if int(rawList[i][0]) < int(rawList[i+1][0]):
                    if int(rawList[i][0]) + len(rawList[i][1]) - 1 < int(rawList[i+1][0]):
                        seq = rawList[i][1]
                        for k in range(int(rawList[i][0]) + len(rawList[i][1]), int(rawList[i+1][0])):
                            seq += "."
                        seq += rawList[i+1][1]
                        seq = seq.replace("N", "-")
                        read = Read(count, seq, int(rawList[i][0]))
                    else:
                        pos = int(rawList[i][0]) + len(rawList[i][1]) - int(rawList[i+1][0])
                        seq = rawList[i][1] + rawList[i+1][1][pos:]
                        seq = seq.replace("N", "-")
                        read = Read(count, seq, int(rawList[i][0]))
                else:
                    if int(rawList[i+1][0]) + len(rawList[i+1][1]) - 1 < int(rawList[i][0]):
                        seq = rawList[i+1][1]
                        for k in range(int(rawList[i+1][0]) + len(rawList[i+1][1]), int(rawList[i][0])):
                            seq += "."
                        seq += rawList[i][1]
                        seq = seq.replace("N", "-")
                        read = Read(count, seq, int(rawList[i+1][0]))
                    else:
                        pos = int(rawList[i+1][0]) + len(rawList[i+1][1]) - int(rawList[i][0])
                        seq = rawList[i+1][1] + rawList[i][1][pos:]
                        seq = seq.replace("N", "-")
                        read = Read(count, seq, int(rawList[i+1][0]))
                self.clusteringTable.append(read)
            for i in range(0, len(rawList) - 1):
                read = Read(count, rawList[i][1], int(rawList[i][0]))
                count += 1
                self.clusteringTable.append(read)
            print(count)

    """
    This method generates the readTable containing the raw reads. The reads, which were corrected
    within shannon entropy will be removed. Afterwards, the compact form of the reads are
    computed and stored in filteredReadList.txt. This will be used during the scoring system.
    """
    def createFilteredReadTable(self):
        file = open('Archive\\5V_allPairs.txt', 'r')
        lines = file.readlines()
        file.close()
        nonePositions = set()
        rawTable = []
        count = 0
        print("Reading all reads from csv file")
        for line in lines:
            parts = line.split(' ')
            read = Read(count, parts[7], int(parts[3]))
            rawTable.append(read)
            count = count + 1
        rawTable.sort(key=lambda x: x.startPos, reverse=False)

        # Algorithm to merge duplicate reads
        print("Merging all duplicates")
        currentPos = -1
        sameStartPos = 0
        readTablePos = 0
        self.clusteringTable = []
        for read in rawTable:
            if read.startPos == currentPos:
                duplicate = False
                for i in range(readTablePos - sameStartPos, readTablePos-1):
                    if self.clusteringTable[i].read == read.read:
                        self.clusteringTable[i].count += 1
                        duplicate = True
                        break
                if not duplicate:
                    read.id = readTablePos
                    self.clusteringTable.append(read)
                    sameStartPos = sameStartPos + 1
                    readTablePos = readTablePos + 1
            else:
                currentPos = read.startPos
                sameStartPos = 1
                read.id = readTablePos
                self.clusteringTable.append(read)
                readTablePos = readTablePos + 1

        # if the read was corrected during entropy error correction,
        # store none
        checkCorrection = np.load('Documents\\checkCorrection.npy')
        print(checkCorrection)
        for i in range(0, len(checkCorrection)):
            if checkCorrection[i]:
                self.clusteringTable[i] = None
                nonePositions.add(i)

        if self.filteredPositions == []:
            self.createFilteredList("Documents")

        if len(self.variantPositions) == len(self.filteredPositions):
            print("VariantList has same length as filteredPositions")
            return

        # In this part, we compute the compact version of the raw reads, which were not
        # corrected during entropy correction
        checkPositions = list(set(self.filteredPositions) - set(self.variantPositions))
        checkPositions.sort(key=lambda x: x, reverse=False)
        print("Difference: ", checkPositions)
        del self.filteredPositions
        del checkCorrection

        # store the position of raw reads which were changed during entropy error correction
        for i in range(0, len(self.clusteringTable)):
            if self.clusteringTable[i] is None:
                continue
            for k in range(0, len(checkPositions)):
                if self.clusteringTable[i].startPos > checkPositions[k] + 1:
                    continue
                if self.clusteringTable[i].startPos + len(self.clusteringTable[i].read) - 1 < checkPositions[k] + 1:
                    break
                pos = (checkPositions[k] + 1) - self.clusteringTable[i].startPos
                if self.clusteringTable[i].read[pos] != self.nucleotides[checkPositions[k]] and self.clusteringTable[i].read[pos] != ".":
                    self.clusteringTable[i] = None
                    nonePositions.add(i)
                    break
        id = 0
        rawList = []

        # generate compact form raw reads and store them in rawList
        for i in range(0, len(self.clusteringTable)):
            if self.clusteringTable[i] is None:
                rawList.append(None)
            else:
                filteredStart, filteredEnd = self.getFilteredPositions(self.clusteringTable[i].read,
                                                                       self.clusteringTable[i].startPos)
                if filteredStart == -1:
                    continue
                filteredRead = self.getFilteredRead(self.clusteringTable[i].read, self.clusteringTable[i].startPos,
                                                    filteredStart, filteredEnd)
                read = FilteredRead(id, filteredRead, filteredStart, filteredEnd, self.clusteringTable[i].count)
                rawList.append(read)
            id += 1

        # The compact raw reads are stored in filteredReadTable.txt
        file = open('Documents\\filteredReadTable.txt', 'w')
        for i in range(0, len(rawList)):
            if rawList[i] is None:
                file.write("_\n")
            else:
                file.write(
                    str(rawList[i].startPos) + " " + str(rawList[i].endPos) + " " + str(
                        rawList[i].count) + " " + rawList[i].read + "\n")
        file.close()
        self.clusteringTable = []

        return nonePositions

    """
    This method stores the reference genome from file
    """
    def createReferenceGenome(self):
        # Read reference genome
        file = open('Documents\\5V_ref_seq.fasta', 'r')
        lines = file.readlines()
        self.referenceGenome = lines[1]

    """
    This method stores in variantPositions all positions relevant
    in compact form (after error correction).
    """
    def createVariantPos(self, folder):
        file = open(folder+'\\variantPositions.txt', 'r')
        line = file.readlines()
        file.close()
        numbers = line[0].split(" ")
        for i in range(0, len(numbers)):
            self.variantPositions.append(int(numbers[i]))

    """
    This method stores in filteredPositions all positions relevant
    in compact form (before error correction)
    """
    def createFilteredList(self, folder):
        file = open(folder+'\\filteredPositions.txt', 'r')
        line = file.readlines()
        file.close()
        numbers = line[0].split(" ")
        for i in range(0, len(numbers)):
            self.filteredPositions.append(int(numbers[i]))

    """
    This method computes variantPositions. These are the positions which contains variants
    after error correction. 
    """
    def calculateFilteredPosition(self, folder, mode, val=None):
        self.createReferenceGenome()

        list = []
        if mode == 0:
            for i in range(0, len(self.referenceGenome)):
                list.append([False, False, False, False, False])
        else:
            for i in range(0, val):
                list.append([False, False, False, False, False])

        # go trough each read and compute the statistics
        for i in range(0, len(self.clusteringTable)):
            read = self.clusteringTable[i]
            for k in range(0, len(read.read)):
                if read.read[k] == "A":
                    list[read.startPos + k - 1][0] = True
                elif read.read[k] == "C":
                    list[read.startPos + k - 1][1] = True
                elif read.read[k] == "G":
                    list[read.startPos + k - 1][2] = True
                elif read.read[k] == "T":
                    list[read.startPos + k - 1][3] = True
                elif read.read[k] == "-":
                    list[read.startPos + k - 1][4] = True

        # If there is a position which has more than one detected nucleotide
        # store _ in nucleotides. Otherwise, store detected nucleotide.
        # Store in variantPositions the positions with more than one detected
        # nucleotide
        variantPositions = []
        nucleotideLetters = ["A", "C", "G", "T", "-"]
        for i in range(0, len(list)):
            trueVariants = 0
            index = -1
            for k in range(0, len(list[i])):
                if list[i][k]:
                    trueVariants += 1
                    index = k
            if trueVariants > 1:
                variantPositions.append(i)
                self.nucleotides.append("_")
            else:
                self.nucleotides.append(nucleotideLetters[index])

        # compute for each read the read and its start position in compact form
        for i in range(0, len(self.clusteringTable)):
            filteredRead = ""
            filteredStartPos = -1
            for k in range(0, len(variantPositions)):
                if self.clusteringTable[i].startPos <= variantPositions[k] + 1 <= self.clusteringTable[
                    i].startPos + len(
                    self.clusteringTable[i].read) - 1:
                    if filteredStartPos == -1:
                        filteredStartPos = k
                    filteredRead += self.clusteringTable[i].read[
                        (variantPositions[k] + 1) - self.clusteringTable[i].startPos]
                self.clusteringTable[i].filteredStartPos = filteredStartPos
                self.clusteringTable[i].filteredRead = filteredRead

        # store the relevant positions in variantPositions.txt
        file = open(folder+'\\variantPositions.txt', 'w')
        variantText = ""
        for i in range(0, len(variantPositions)):
            if i == len(variantPositions) - 1:
                variantText += str(variantPositions[i])
            else:
                variantText += str(variantPositions[i]) + " "
        file.write(variantText)
        file.close()

        # Store the shared nucleotides in NucleotideFile.txt
        file = open(folder+'\\NucleotideFile.txt', 'w')
        for i in range(0, len(self.nucleotides)):
            nstr = str(self.nucleotides[i])
            file.write(nstr)
        file.close()

    """
    This method stores the nucleotides from file into nucleotides
    @folder folder where we can find the file
    """
    def getNucleotides(self, folder):
        file = open(folder+'\\NucleotideFile.txt', 'r')
        line = file.readlines()
        file.close()
        nucleotides = line[0].split("\n")
        for i in range(0, len(nucleotides[0])):
            self.nucleotides.append(nucleotides[0][i])
        print(self.nucleotides)

    """
    This method is used to compute for a read its start and end position
    in compact form
    @read nulceotide sequence of the read
    @startPos start position of read in normal form
    """
    def getFilteredPositions(self, read, startPos):
        filteredStart = 0
        filteredEnd = 0
        check1 = False
        check2 = False
        for i in range(0, len(self.variantPositions)):
            if filteredStart == 0 and startPos - 1 <= self.variantPositions[i]:
                filteredStart = i
                check1 = True
            if filteredEnd == 0 and startPos + len(read) - 2 <= self.variantPositions[i]:
                filteredEnd = i - 1
                check2 = True
                break
        if filteredEnd == 0 and not check2:
            filteredEnd = len(self.variantPositions) - 1
        if filteredStart == 0 and not check1:
            filteredStart = -1

        return filteredStart, filteredEnd

    """
    This method computes the compact form of a given read
    @read sequence of the read
    @readStart start position of the read in normal form
    @filteredStart start position of the read in compact form
    @filteredEnd end position of the read in compact form
    """
    def getFilteredRead(self, read, readStart, filteredStart, filteredEnd):
        filteredRead = ""

        for k in range(filteredStart, filteredEnd + 1):
            try:
                filteredRead += read[self.variantPositions[k] + 1 - readStart]
            except Exception as e:
                print(filteredStart, filteredEnd + 1, read, readStart, readStart + len(read), self.variantPositions[k], self.variantPositions[k] + 1 - readStart)
        return filteredRead

    """
    This method is used to compute the readList. The read list contains the splitted parts of the reads
    from the clusteringTable in compact form. It generates a file readList.txt which contains these parts. 
    @folder folder in which the file readList.txt should be stored
    """
    def createFilteredReadList(self, folder):
        id = 0
        rawList = []
        # the parts are stored first in rawList.
        for i in range(0, len(self.clusteringTable)):
            # There were reads in the dataset which contained in the very first as well as in the very last
            # position a dot. To ensure that they don't
            if self.clusteringTable[i].read[0] == ".":
                self.clusteringTable[i].read = self.clusteringTable[i].read[1:]

            if self.clusteringTable[i].read[len(self.clusteringTable[i].read) - 1] == ".":
                self.clusteringTable[i].read = self.clusteringTable[i].read[:len(self.clusteringTable[i].read) - 1]

            nucleotideStrings = self.clusteringTable[i].read.split(".")
            # check whether read contains a gap. If not use entire read. Otherwise, split the read into two parts at
            # the gap
            if len(nucleotideStrings) == 1:
                filteredStart, filteredEnd = self.getFilteredPositions(self.clusteringTable[i].read,
                                                                       self.clusteringTable[i].startPos)
                if filteredStart == -1:
                    continue
                filteredRead = self.getFilteredRead(self.clusteringTable[i].read, self.clusteringTable[i].startPos,
                                                    filteredStart, filteredEnd)
                if filteredRead != "":
                    read = FilteredRead(id, filteredRead, filteredStart, filteredEnd, self.clusteringTable[i].count)
                    rawList.append(read)
                    id += 1
            else:
                count = 1

                filteredStart, filteredEnd = self.getFilteredPositions(
                    self.clusteringTable[i].read[:len(nucleotideStrings[0])],
                    self.clusteringTable[i].startPos)

                if filteredStart == -1:
                    continue

                filteredRead = self.getFilteredRead(self.clusteringTable[i].read[:len(nucleotideStrings[0])],
                                                    self.clusteringTable[i].startPos,
                                                    filteredStart, filteredEnd)

                if len(self.clusteringTable[i].read[:len(nucleotideStrings[0])]) == 0:
                    print(self.clusteringTable[i].read)

                if filteredRead != "":
                    read1 = FilteredRead(id, filteredRead, filteredStart, filteredEnd, self.clusteringTable[i].count)
                    rawList.append(read1)
                    id += 1

                try:
                    while nucleotideStrings[count] == "":
                        count += 1
                except Exception as e:
                    print(self.clusteringTable[i].read)
                    continue

                filteredStart, filteredEnd = self.getFilteredPositions(
                    self.clusteringTable[i].read[count + len(nucleotideStrings[0]):],
                    self.clusteringTable[i].startPos + count + len(nucleotideStrings[0]))

                if filteredStart == -1:
                    continue

                filteredRead = self.getFilteredRead(self.clusteringTable[i].read[count + len(nucleotideStrings[0]):],
                                                    self.clusteringTable[i].startPos + count + len(
                                                        nucleotideStrings[0]),
                                                    filteredStart, filteredEnd)
                if len(self.clusteringTable[i].read[count + len(nucleotideStrings[0]):]) == 0:
                    print(self.clusteringTable[i].read)
                if filteredRead != "":
                    read2 = FilteredRead(id, filteredRead, filteredStart, filteredEnd, self.clusteringTable[i].count)
                    id += 1
                    rawList.append(read2)

        print("Merging all duplicates")
        # merge duplicates
        rawList.sort(key=lambda x: (x.startPos, x.endPos), reverse=False)
        currentPos = -1
        sameStartPos = 0
        readTablePos = 0
        self.readList = []
        for read in rawList:
            if read.startPos == currentPos:
                duplicate = False
                for i in range(readTablePos - sameStartPos, readTablePos):
                    if self.readList[i].read == read.read:
                        self.readList[i].count += read.count
                        duplicate = True
                        break
                if not duplicate:
                    read.id = readTablePos
                    self.readList.append(read)
                    sameStartPos = sameStartPos + 1
                    readTablePos = readTablePos + 1
            else:
                currentPos = read.startPos
                sameStartPos = 1
                read.id = readTablePos
                self.readList.append(read)
                readTablePos = readTablePos + 1

        # the parts in compact form are stored in ReadList.txt
        file = open(folder+'\\ReadList.txt', 'w')
        for i in range(0, len(self.readList)):
            file.write(str(self.readList[i].startPos) + " " + str(self.readList[i].endPos) + " " + str(
                self.readList[i].count) + " " + self.readList[i].read + "\n")
        file.close()

    """
    This method is used to compute the pairedEndReadList. The pairedEndReadList contains the reads
    from the clusteringTable in compact form. It generates a file pairedEndReadList.txt which contains these reads.
    Since compact form of different reads can be merged into one, we keep track where a read in clusterinTable is 
    stored in pairedEndReadList. This information will be essential for the scoring system. Because there, we need
    to find the corresponding raw reads.  
    @nonePositions position of the reads in table which were corrected with shannon entropy
    @folder folder in which the file pairedEndReadList.txt should be stored
    """
    def createFilteredPairedEndReadList(self, nonePositions, folder):
        id = 0
        rawList = []
        # compute for each read its compact form and store it in rawList
        for i in range(0, len(self.clusteringTable)):
            filteredStart, filteredEnd = self.getFilteredPositions(self.clusteringTable[i].read,
                                                                   self.clusteringTable[i].startPos)
            if filteredStart == -1:
                continue

            filteredRead = self.getFilteredRead(self.clusteringTable[i].read, self.clusteringTable[i].startPos,
                                                filteredStart, filteredEnd)
            read = FilteredRead(id, filteredRead, filteredStart, filteredEnd, self.clusteringTable[i].count)
            rawList.append(read)

        print("Merging all duplicates")
        self.pairedEndReadList = []
        # we keep track at which position a read is stored/merged in pairedEndReadList
        # This will be important for the scoring system.
        positionPointer = []
        tracker, rawList = zip(*sorted(enumerate(rawList), key=lambda x: (x[1].startPos, x[1].endPos), reverse=False))
        currentPos = -1
        sameStartPos = 0
        readTablePos = 0
        counter = 0
        for read in rawList:
            if read.startPos == currentPos:
                duplicate = False
                for i in range(readTablePos - sameStartPos, readTablePos):
                    if self.pairedEndReadList[i].read == read.read:
                        self.pairedEndReadList[i].count += read.count
                        duplicate = True
                        positionPointer[i].append(tracker[counter])
                        break
                if not duplicate:
                    read.id = readTablePos
                    self.pairedEndReadList.append(read)
                    positionPointer.append([tracker[counter]])
                    sameStartPos = sameStartPos + 1
                    readTablePos = readTablePos + 1
            else:
                currentPos = read.startPos
                sameStartPos = 1
                read.id = readTablePos
                self.pairedEndReadList.append(read)
                positionPointer.append([tracker[counter]])
                readTablePos = readTablePos + 1
            counter += 1
        del rawList
        # store compact reads in file
        file = open(folder+'\\PairedEndReadList.txt', 'w')
        for i in range(0, len(self.pairedEndReadList)):
            file.write(
                str(self.pairedEndReadList[i].startPos) + " " + str(self.pairedEndReadList[i].endPos) + " " + str(
                    self.pairedEndReadList[i].count) + " " + self.pairedEndReadList[i].read + "\n")
        file.close()
        self.clusteringTable = []
        # store the position where a read from clusteringTable is located in pairedEndReadList
        # in tracker. This is only used for the 5-virus mix.
        if folder == "Documents":
            file = open('Documents\\tracker.txt', 'w')
            for i in range(0, len(positionPointer)):
                pos = ""
                for k in range(0, len(positionPointer[i])):
                    if positionPointer[i][k] in nonePositions:
                        continue
                    if k != len(positionPointer[i]) - 1:
                        pos += str(positionPointer[i][k]) + " "
                    else:
                        pos += str(positionPointer[i][k])
                if pos == "":
                    pos = "_"
                file.write(pos + "\n")

    """
    This method reads from file ReadList.txt the splitted compact reads and
    stores them into readList
    @folder folder where the file is
    """
    def getReadList(self, folder):
        file = open(folder+'\\ReadList.txt', 'r')
        lines = file.readlines()
        file.close()
        count = 0
        for line in lines:
            parts = line.split(' ')
            if len(parts) > 1:
                parts[3] = parts[3].replace(" ", "")
                parts[3] = parts[3].replace("\n", "")
                read = FilteredRead(count, parts[3], int(parts[0]), int(parts[1]), int(parts[2]))
                self.readList.append(read)
                count = count + 1

    """
    This method reads from file PairedEndReadList.txt the compact reads 
    and stores them into pairedEndReadList
    @folder folder where the file is
    """
    def getPairedEndReadList(self, folder):
        file = open(folder+'\\PairedEndReadList.txt', 'r')
        lines = file.readlines()
        file.close()
        count = 0
        for line in lines:
            parts = line.split(' ')
            if len(parts) > 1:
                parts[3] = parts[3].replace(" ", "")
                parts[3] = parts[3].replace("\n", "")
                read = FilteredRead(count, parts[3], int(parts[0]), int(parts[1]), int(parts[2]))
                self.pairedEndReadList.append(read)
                count = count + 1

    """
    This method reads from file FilteredReadTable.txt the compact raw reads and
    stores them into readTable
    @folder folder where the file is
    """
    def getFilteredReadTable(self, folder):
        file = open(folder+'\\FilteredReadTable.txt', 'r')
        lines = file.readlines()
        file.close()
        count = 0
        for line in lines:
            parts = line.split(' ')
            if len(parts) == 1:
                parts[0] = parts[0].replace("\n", "")
                if parts[0] == "_":
                    self.readTable.append(None)
                    count = count + 1
            if len(parts) > 1:
                parts[3] = parts[3].replace(" ", "")
                parts[3] = parts[3].replace("\n", "")
                read = FilteredRead(count, parts[3], int(parts[0]), int(parts[1]), int(parts[2]))
                self.readTable.append(read)
                count = count + 1

    """
    This method was used to build the read graph. The read graph was build using an adapted de Bruijn
    Graph algorithm. 
    """
    def deBruijnAlgorithm(self, nlist, k, border1, border2, folder):
        # store here the parent-child relationship
        self.successors = {}
        # store here the positions where splits in the graph was observed
        self.specialList = []
        # store the nodes which were created
        nodes = set()
        print("Number of reads: ", len(nlist))
        print("k = ", k)
        print("Border: ", border1, border2)
        # In self.edges, we store for each parent-child the number of appearance

        # Store here the roots and how many times they appeared in the dataset
        root = {}
        for i in range(0, len(nlist)):
            for j in range(0, len(self.readList[nlist[i]].read) - k + 1):
                if border1 <= j + self.readList[nlist[i]].startPos <= border2 and border1 <= j + self.readList[
                    nlist[i]].startPos + k <= border2:
                    # Check if predecessor node exists. If not, then it must be a root
                    if (self.readList[nlist[i]].read[j:j + k - 1], self.readList[nlist[i]].startPos + j) not in nodes:
                        # detected node for the first time. Add it as a root
                        root[(self.readList[nlist[i]].read[j:j + k - 1], self.readList[nlist[i]].startPos + j)] = self.readList[nlist[i]].count
                    elif (self.readList[nlist[i]].read[j:j + k - 1], self.readList[nlist[i]].startPos + j) in root:
                        # if node is already a root, update the number of times the root has appeared in the dataset
                        current = root[(self.readList[nlist[i]].read[j:j + k - 1], self.readList[nlist[i]].startPos + j)]
                        root.update({(self.readList[nlist[i]].read[j:j + k - 1], self.readList[nlist[i]].startPos + j): current + self.readList[nlist[i]].count})

                    # Check if predecessor is in dictionary. If not, add it.
                    if (self.readList[nlist[i]].read[j:j + k - 1],
                        self.readList[nlist[i]].startPos + j) not in self.successors:
                        self.successors[
                            (self.readList[nlist[i]].read[j:j + k - 1], self.readList[nlist[i]].startPos + j)] = [
                            (self.readList[nlist[i]].read[j + 1:j + k], self.readList[nlist[i]].startPos + j + 1)]
                        # Store the count for them in edges
                        self.edges[((self.readList[nlist[i]].read[j:j + k - 1], self.readList[nlist[i]].startPos + j),
                                    (self.readList[nlist[i]].read[j + 1:j + k], self.readList[nlist[i]].startPos + j + 1))] \
                            = self.readList[nlist[i]].count
                    # add successor to predecessor in dictionary.
                    else:
                        # check whether successor already exists for predecessor
                        values = self.successors[
                            (self.readList[nlist[i]].read[j:j + k - 1], self.readList[nlist[i]].startPos + j)]
                        if values.count((self.readList[nlist[i]].read[j + 1:j + k],
                                         self.readList[nlist[i]].startPos + j + 1)) == 0:
                            # if not, add it as successor
                            self.successors[(
                                self.readList[nlist[i]].read[j:j + k - 1],
                                self.readList[nlist[i]].startPos + j)].append(
                                (self.readList[nlist[i]].read[j + 1:j + k], self.readList[nlist[i]].startPos + j + 1))
                            # generate new entry to compute number of appearance
                            self.edges[
                                ((self.readList[nlist[i]].read[j:j + k - 1], self.readList[nlist[i]].startPos + j),
                                 (self.readList[nlist[i]].read[j + 1:j + k], self.readList[nlist[i]].startPos + j + 1))] \
                                = self.readList[nlist[i]].count
                        else:
                            # if yes, just update number of appearance
                            current = self.edges[
                                ((self.readList[nlist[i]].read[j:j + k - 1], self.readList[nlist[i]].startPos + j),
                                 (self.readList[nlist[i]].read[j + 1:j + k], self.readList[nlist[i]].startPos + j + 1))]

                            self.edges.update({((self.readList[nlist[i]].read[j:j + k - 1], self.readList[nlist[i]].startPos + j),
                                 (self.readList[nlist[i]].read[j + 1:j + k], self.readList[nlist[i]].startPos + j + 1)): current + self.readList[nlist[i]].count})

                    # add node
                    nodes.add((self.readList[nlist[i]].read[j:j + k - 1], self.readList[nlist[i]].startPos + j))
                    nodes.add(
                        (self.readList[nlist[i]].read[j + 1:j + k], self.readList[nlist[i]].startPos + j + 1))

        rootBorder = 0
        count = 0
        # At this point, we compute the number of times the roots have appeared in the
        # dataset. With the mean of it, we will calculate a min threshold. It the root count
        # is smaller than a min threshold, no depth first search is applied on that root.
        for item in root.items():
            if int(item[1]) > 1:
                rootBorder += int(item[1])
                count += 1
        rootBorder = rootBorder / count * 0.20
        print("Root Border: ", rootBorder)
        self.haplotypes = []
        badEdges = []
        badSuccessor = []
        print("Number of edges before singleton removal: ", len(self.edges))

        # Here, we already remove edges which have an edge count of less than 3.
        # Is used to reduce the tree complexity
        for edge in self.edges.items():
            if edge[1] < 3:
                values = self.successors[edge[0][0]]
                values.remove(edge[0][1])
                if not values:
                    badSuccessor.append(edge[0][0])
                else:
                    self.successors.update({edge[0][0]: values})
                badEdges.append(edge[0])

        for i in range(0, len(badSuccessor)):
            self.successors.pop(badSuccessor[i], None)

        for i in range(0, len(badEdges)):
            self.edges.pop(badEdges[i], None)

        print("Applying DFS")
        # For each root having a count of > min threshold, the dfs is applied.
        if folder == "Documents":
            for item in root.items():
                if item[1] >= rootBorder:
                    print("Root ", item[0], item[1], rootBorder)
                    self.depthFirstSearch(item[0], "", k-1, border2 - border1)
        else:
            for item in root.items():
                if item[1] >= rootBorder:
                    print("Root ", item[0], item[1], rootBorder)
                    self.depthFirstSearch(item[0], "", k-1, border2 - border1-1)

        print("Num of local haplotypes: ", len(self.haplotypes))
        print("SpecialList ", self.specialList)
        self.badSuccessor = set()
        self.successors = {}

        return len(self.haplotypes)


    """
    This method is used to generate the global haplotypes by using depth first search. 
    It traverses the read graph. If a complete haplotype is found, it stores it in the 
    list haplotypes. It uses the tip removal algorithm derived from Velvet. Therefore, 
    if we have more than one successor node, we check for minority count to remove possible
    tips. If the path to a leaf is smaller than the window, it is considered as a tip.
    @node current node
    @haplotype current sequence of haplotype which is build up during dfs
    @currentDepth the current depth of the node in the tree
    @maxDepth the max depth of the tree resp max length of a haplotype/window
    """
    def depthFirstSearch(self, node, haplotype, currentDepth, maxDepth):
        # Check whether leaf node
        if node not in self.successors:
            # If it is a tip, return
            if currentDepth < maxDepth:
                return
            # Otherwise, add haplotype in haplotypes
            haplotype += node[0]
            self.haplotypes.append(haplotype)
            return

        # extend haplotype by sequence of node
        haplotype += node[0][0]
        successors = self.successors[node]

        minThreshold = 0
        if len(successors) > 1:
            # If more than one successor is present, compute largest outgoing arc value.
            # Calculate then the minimum threshold.
            for i in range(0, len(successors)):
                if self.edges[(node, successors[i])] > minThreshold:
                    minThreshold = self.edges[(node, successors[i])]
            minThreshold = minThreshold * 0.01
            i = 0
            while i < len(successors):
                outgoing = self.edges[(node, successors[i])]
                # minority count and singleton removal
                # If number of outgoing arcs to this successor node is larger than
                # minThreshold, apply recursively depth first search on it.
                # Otherwise, we have a tip and remove it from edges and successors.
                if outgoing >= minThreshold:
                    self.depthFirstSearch(successors[i], haplotype, currentDepth + 1, maxDepth)
                    i = i + 1
                else:
                    del self.edges[(node, successors[i])]
                    successors.remove(successors[i])
                    self.successors.update({node: successors})

            # Store the position of the split in the read graph in special list
            if self.specialList.count(node[1] + len(node[0])) == 0:
                self.specialList.append(node[1] + len(node[0]))
        else:
            # only one successor. Directly apply depth first search
            self.depthFirstSearch(successors[0], haplotype, currentDepth + 1, maxDepth)
        return

    """
    This method is used to compute whether given read support a haplotype at a given position.
    @readId position of the read in list
    @haploId position of the haplotype in haplotypes
    @border1 start position of haplotype on genome
    @border2 end position of haplotype on genome
    @specialListPos position which needs to be expanded 
    """
    def similarity(self, readId, haploId, border1, border2, specialListPos):

        start = 0
        if self.list[readId].startPos < border1:
            start = border1
        else:
            start = self.list[readId].startPos

        comparedNucleotides = 0
        leftExtension = 0
        rightExtension = 0
        crucialPosition = 0

        for i in range(start, len(self.variantPositions)):
            if self.list[readId].endPos < i or border2 - 1 < i:
                break
            if self.list[readId].startPos > i or border1 > i:
                continue
            if border1 == 0:
                posOfHaplo = i - border1 - 1
            else:
                posOfHaplo = i - border1
            posOfRead = i - self.list[readId].startPos

            # actual comparison
            if self.list[readId].read[posOfRead] != self.haplotypes[haploId][posOfHaplo]:
                if self.list[readId].read[posOfRead] != ".":
                    return False
            else:
                comparedNucleotides += 1

            # if compared position doesn't contain a dot, increase to which side the read extended
            # the haplotype
            if i < self.positionList[specialListPos] and self.list[readId].read[posOfRead] != ".":
                leftExtension += 1
            elif i > self.positionList[specialListPos] and self.list[readId].read[posOfRead] != ".":
                rightExtension += 1
            elif i == self.positionList[specialListPos] and self.list[readId].read[posOfRead] != ".":
                crucialPosition += 1

        # This would mean that the read didn't extend haplotype to the right resp to the left. Then,
        # don't use it as a supporting read
        if rightExtension == 0 and crucialPosition == 0 or leftExtension == 0 and crucialPosition == 0:
            return False

        # If actually compared nucleotides are at max 4 (out of 13),
        # return false (not good enough). In the Experiments, it was not used

        if comparedNucleotides <= 4:
            return False

        # Was used for score
        #if crucialPosition == 0:
        #    return False
        # Read extends position
        return True


    """
    This method can be used to remove inconsistent/false haplotypes. For the graph, 
    it was only used for local haplotypes. For the consistency check, it was used to 
    remove false global haplotypes. An haplotype is removed if it doesn't have many 
    reads supporting a position compared to others.
    @border1 start position of the haplotype on the genome
    @border2 end position of the haplotype on the genome
    @mode 0 if it is used to remove local haplotypes (read graph), 1 if it is used to 
    remove global (consistency check)
    @nList only used for consistency check. Shows the position of the raw reads which needs
    to be considered in this check
    """
    def removeInconsistentHaplotypes(self, border1, border2, mode, nList=None):

        list = []
        haplotypes = []
        haploStatistics = []
        totalExtension = []
        if mode == 0:
            # In this mode, we only remove local haplotypes. As a consequence,
            # we consider all reads from the paired end read list and we look at
            # all positions where splits were detected
            self.list = self.pairedEndReadList
            self.positionList = self.specialList
        else:
            # In this modes, we remove global haplotypes. We store all raw reads,
            # which needs to be considered. These are the reads, which were not corrected
            # during any correction algorithm. We only use the positions mentioned in the
            # consistency check
            self.list = []
            self.positionList = []
            for i in range(0, len(nList)):
                self.list.append(self.readTable[nList[i]])
            self.positionList = self.rawPositionList

        for i in range(0, len(self.positionList)):
            list.append(0)
            # we store the supporting statistics here
            haploStatistics.append(0)
            # we store the total number of detected supporting reads here
            totalExtension.append(0)
        for i in range(0, len(self.haplotypes)):
            haplotypes.append(list.copy())

        # number of compact positions by which a read needs to expand a position
        # of a haplotype to the left and to the right
        minLeftExtension = 5
        minRightExtension = 5

        # We go through each read. We calculate for each haplotype at the mentioned
        # position the number of extending reads. The information will be stored then
        # in haploStatistics
        for i in range(0, len(self.list)):
            a = 0
            b = 0

            if self.list[i].startPos >= border2:
                break
            if self.list[i].endPos <= border1:
                continue
            if self.list[i].startPos > border1:
                a = self.list[i].startPos - border1
            if self.list[i].startPos + len(self.list[i].read) < border2:
                b = border2 - (self.list[i].startPos + len(self.list[i].read))

            # If read doesn't overlap window by this number of compact positions,
            # skip it. This can be chosen by the user.
            if (border2 - border1 - a - b) < 11:
                continue

            # If needed, we can only focus on the reads containing gaps.
            #if self.list[i].read.count(".") == 0 and mode == 0:
            #    continue

            for t in range(0, len(self.positionList)):

                # Constraints on length of a supporting read. If the read doesn't cover
                # the positions which needs to be covered at least to be considered as supporting,
                # skip it.
                check = False
                if self.list[i].startPos <= self.positionList[t] - minLeftExtension and \
                    self.list[i].endPos >= self.positionList[t] + minRightExtension:
                    for j in range(0, len(haplotypes)):
                        # Check whether read expands haplotype
                        if self.similarity(i, j, border1, border2, t):
                            # If yes, increase number of supporting reads for haplotype j at
                            # position t
                            haplotypes[j][t] += self.list[i].count
                            if not check:
                                # If read didn't extended a haplotype yet, increase the number
                                # of total expanding reads (used to avoid duplicate increase)
                                totalExtension[t] += self.list[i].count
                                check = True

        if mode == 0:
            for i in range(0, len(haplotypes)):
                for j in range(0, len(haplotypes[i])):
                    # compute total number of extending reads at each position
                    haploStatistics[j] += haplotypes[i][j]

            for i in range(0, len(haploStatistics)):
                # Compute min threshold
                haploStatistics[i] = haploStatistics[i] / len(haplotypes) * 0.01
            print(haplotypes[0])
            deleted = False
            i = 0
            while i < len(haplotypes):
                deleted = False
                for k in range(0, len(haplotypes[i])):
                    # If number of reads extending haplotype i at position k is lower than
                    # the computed threshold, remove it
                    if haplotypes[i][k] < haploStatistics[k]:
                        del haplotypes[i]
                        del self.haplotypes[i]
                        deleted = True
                        break
                if not deleted:
                    i = i + 1
        elif mode == 1:
            del haploStatistics
        else:
            return totalExtension, haplotypes

        return haplotypes[0]

    """
    This method is used to apply the construction of the read graphs. It collects
    all the reads which should be used. Then it calls the de Bruijn graph method.
    If necessary, it applies consistency check to remove clear false Haplotypes. 
    It also used to call the regression to reduce the number of local. At the end,
    it returns the number of local/global haplotypes in given region.
    haplotypes
    @border1 start position of local haplotype
    @border2 end position of local haplotype
    @roh current roh value
    @mode 0 if we want to create the read graph. 1 if we want to use regression
    @folder location where the necessary documents are stored
    @haploFolder folder containing the information about haplotypes
    """
    def getLocalHaplotypes(self, border1, border2, roh, mode, folder, haploFolder):
        nlist = []
        pos = 0
        t = []
        for i in range(0, 200):
            t.append(0)
        # Selects reads from readList which overlap with window
        for read in self.readList:
            a = 0
            b = 0
            # include all reads overlapping window
            if read.startPos > border1:
                a = read.startPos - border1
            if read.startPos + len(read.read) < border2:
                b = border2 - (read.startPos + len(read.read))

            if (border2 - border1 - a - b) > 1:
                nlist.append(pos)
            pos += 1

        if mode == 0:
            # Construct read graph. NumOfHaplo stores the number of generated
            # global/local haplotypes
            numOfHaplo = self.deBruijnAlgorithm(nlist, 11, border1, border2, folder)
            # If this condition is met, apply consistency check to remove false haplotypes
            # Don't further increase this value. Otherwise, the execution time will be
            # very large
            if numOfHaplo <= 1500:
                # For this check, only paired end reads are used
                if not self.pairedEndReadList:
                    self.getPairedEndReadList(folder)
                self.specialList.sort(key=lambda x: x, reverse=False)
                self.removeInconsistentHaplotypes(border1, border2, 0)
                numOfHaplo = len(self.haplotypes)
        else:
            # Collect paired end reads overlapping window. Use them for regression
            numOfHaplo = len(self.haplotypes)
            if numOfHaplo <= 1000:
                if not self.pairedEndReadList:
                    self.getPairedEndReadList(folder)
                pairedList = []
                pos = 0
                for read in self.pairedEndReadList:
                    a = 0
                    b = 0
                    if read.startPos > border1:
                        a = read.startPos - border1
                    if read.startPos + len(read.read) < border2:
                        b = border2 - (read.startPos + len(read.read))
                        # TODO change
                    if (border2 - border1 - a - b) >= 1:
                        pairedList.append(pos)
                    pos += 1
                # Call regression
                regression = StagewiseRegression(self, pairedList, roh)
                numOfHaplo = regression.pipe(border1, border2, haploFolder)
        # return number of local/global haplotypes after 1. generating read graph or after 2.
        # reducing the number of local haplotypes
        return numOfHaplo

    """
    Before the above method is called to reduce number of local haplotypes in a window,
    this method is called. It recalls all the stored haplotypes from the given window.
    @border1 start position of local haplotype
    @border2 end position of local haplotype
    @roh current roh value
    @folder location where the necessary documents are stored
    @haploFolder location where the file containing the haplotypes needs to be stored
    """
    def reduceLocalHaplotypes(self, border1, border2, roh, folder, haploFolder):
        self.haplotypes = []
        if border1 < 100:
            string1 = "0" + str(border1)
        else:
            string1 = str(border1)
        if border2 < 100:
            string2 = "0" + str(border2)
        else:
            string2 = str(border2)
        # Open the file containing the haplotypes which were already computed for
        # the given window
        file = open(haploFolder+'\\LSR_' + string1 + '_' + string2 + '.txt', 'r')
        lines = file.readlines()
        file.close()
        for line in lines:
            if len(line) > 1:
                line = line.replace(" ", "")
                line = line.replace("\n", "")
                # Store haplotypes here.
                self.haplotypes.append(line)
        return self.getLocalHaplotypes(border1, border2, roh, 1, folder, haploFolder)

    """
    This method returns for any position the closest position in compact form
    @originalPos original position
    """
    def getVariantPosition(self, originalPos):
        for i in range(0, len(self.variantPositions)):
            if originalPos <= self.variantPositions[i]:
                return i
        return len(self.variantPositions) - 1

    """
    Once all local haplotypes from all local subregions are created, this method is called.
    It generates then the global haplotypes using depth first search
    @depth current subregion. The first subregion has as depth 0. All of its haplotypes are
    stored in haplotypeList at the 0th. position. 
    """
    def recursiveGlobalHaplo(self, depth, haplotype):
        if depth == len(self.haplotypeList) - 1:
            for i in range(0, len(self.haplotypeList[depth])):
                self.haplotypes.append(haplotype + self.haplotypeList[depth][i])
            return
        for i in range(0, len(self.haplotypeList[depth])):
            self.recursiveGlobalHaplo(depth + 1, haplotype + self.haplotypeList[depth][i])

"""
This method contains the entire pipeline to produce the global haplotypes.
"""
def pipe():
    dbg = HaplotypeReconstruction_DBG()
    # Documents folder was used for the 5 virus mix dataset.
    # A folder called Experiment was used for the experiments.
    folder = "Documents"
    # In this variable, we store the location where the files produces during
    # graph construction is stored
    haplotypeFolder = "Haplotypes"
    if folder == "Documents":
        # If we focus on the virus mix, then we create the clustering table
        # and the compact positions
        dbg.createClusteringTable()
        dbg.createReferenceGenome()
        dbg.calculateFilteredPosition(folder, 0)
    else:
        print("Applying Experiment")
        # For the experiments, the clustering tables are created using the sam
        # files produced by ART
        dbg.createClusteringTableFromSamFile(folder)
        dbg.calculateFilteredPosition(folder, 0, 4306)
    # generate datastructrures variantPositions, readList, nucleotides
    dbg.createVariantPos(folder)
    dbg.createFilteredReadList(folder)
    dbg.getNucleotides(folder)
    if folder == "Documents":
        # In the experiments, we compute the compact form reads of the raw data
        nonePositions = dbg.createFilteredReadTable()
        dbg.createClusteringTable()
        # Compute compact form paired end read List. Use nonPositions to
        # keep track of the movements. Further details can be found in the
        # description of the method
        dbg.createFilteredPairedEndReadList(nonePositions, folder)
    else:
        # compute compact form reads. Raw reads are irrelevant for the experiments
        dbg.createFilteredPairedEndReadList([], folder)
    # get read list
    dbg.getReadList(folder)

    intervals = []

    # Store in the intervals the compact start and end position of the entire region.
    # For the virus mix, we only focus on the gag-pol
    if folder == "Documents":
        intervals.append([dbg.getVariantPosition(790), dbg.getVariantPosition(5096)])
    else:
        intervals.append([dbg.getVariantPosition(0), dbg.getVariantPosition(4306)])

    numOfGlobalHaplo = 0
    acceptedIntervals = []

    checkPositions = []
    # Entire process of generating the global haplotypes
    while len(intervals) > 0:
        # pop pending region interval
        border1, border2 = intervals.pop(0)
        # generate read graph
        numOfHaplo = dbg.getLocalHaplotypes(border1, border2, 0, 0, folder, haplotypeFolder)
        # if number of generated haplotypes is too large, split region into subregions
        if numOfHaplo > 1000:
            midBorder = int((border1 + border2) / 2)
            intervals.append([border1, midBorder])
            intervals.append([midBorder, border2])
            dbg.haplotypes = []
        else:
            # Used to compute total number of global haplotypes
            if numOfGlobalHaplo == 0:
                numOfGlobalHaplo = numOfHaplo
            else:
                numOfGlobalHaplo *= numOfHaplo
            # Store interval as accepted.
            acceptedIntervals.append([border1, border2])
            # Store its haplotypes in file. File name contains information
            # about the regions interval
            if border1 < 100:
                string1 = "0" + str(border1)
            else:
                string1 = str(border1)
            if border2 < 100:
                string2 = "0" + str(border2)
            else:
                string2 = str(border2)
            file = open(haplotypeFolder+'\\LSR_' + string1 + '_' + string2 + '.txt', 'w')
            for i in range(0, len(dbg.haplotypes)):
                file.write(dbg.haplotypes[i] + "\n")
            file.close()
            print("special ", dbg.specialList)
            # Store positions where splits happened
            for i in range(0, len(dbg.specialList)):
                if checkPositions.count(dbg.specialList[i]) == 0:
                    checkPositions.append(dbg.specialList[i])

    print("accepted Intervals: ", acceptedIntervals)
    print("number of global Haplos: ", numOfGlobalHaplo)
    checkPositions.sort(key=lambda x: x, reverse=False)
    dbg.readList = []

    # Store all positions where split happens
    file = open(folder+'\\checkPositions.txt', 'w')
    positions = ""
    for i in range(0, len(checkPositions)):
        if i != len(checkPositions) - 1:
            positions += str(checkPositions[i]) + " "
        else:
            positions += str(checkPositions[i])
    file.write(positions)

    # If number of global haplotypes is too large, use regression to reduce it
    roh = 0.001
    while numOfGlobalHaplo >= 1000:
        print("roh ", roh)
        numOfGlobalHaplo = 0
        for i in range(0, len(acceptedIntervals)):
            numOfHaplo = dbg.reduceLocalHaplotypes(acceptedIntervals[i][0], acceptedIntervals[i][1], roh, folder, haplotypeFolder)
            if numOfGlobalHaplo == 0:
                numOfGlobalHaplo = numOfHaplo
            else:
                numOfGlobalHaplo *= numOfHaplo
        roh *= 2

    dbg.haplotypeList = []
    dbg.haplotypes = []
    count = 0
    # Generate global haplotypes with the local
    for file in os.listdir(haplotypeFolder):
        if file.startswith("LSR"):
            # Files are ordered in folder according to ascending border interval
            helperList = []
            f = open(haplotypeFolder+"\\" + file, 'r')
            lines = f.readlines()
            f.close()
            for line in lines:
                parts = line.split(' ')
                if len(parts[0]) > 1 and parts[0] != " ":
                    parts[0] = parts[0].replace("\n", "")
                    helperList.append(parts[0])
            dbg.haplotypeList.append(helperList)
            count += 1
    # If count is 1, then only one region, namely the entire region was
    # considered. Otherwise, we had subregions.
    if count == 1:
        file = open(haplotypeFolder+"\\GlobalHaplotypes.txt", 'w')
        for i in range(0, len(dbg.haplotypeList[0])):
            file.write(dbg.haplotypeList[0][i] + "\n")
        file.close()
    else:
        # Compute global haplotypes using depth first search and store it
        # in file
        dbg.recursiveGlobalHaplo(0, "")
        file = open(haplotypeFolder+"\\GlobalHaplotypes.txt", 'w')
        for i in range(0, len(dbg.haplotypes)):
            file.write(dbg.haplotypes[i] + "\n")
        file.close()

# Call main function to execute pipe line to compute global haplotypes
if __name__ == '__main__':
    pipe()
