import math
import numpy as np
import statistic
from helperClasses import Read

class EntropyErrorCorrection:

    """
    This class serves to correct reads with Shannon Entropy. We calculate for each position the entropy of the
    nucleotides A, C, G, T and -. If a position has a lower entropy than a certain threshold, then it is corrected
    in all overlapping reads using majority rule. Afterwards, the reads are stored in the file EntropyFile.txt
    """
    def __init__(self):
        self.readTable = []
        self.referenceGenome = ""
        self.checkCorrection = np.zeros(0)

    """
    Reads from 5V_allPairs.txt the raw reads with given format
    and merges duplicates.
    """
    def createReadTable(self):
        # Read all reads and store them into rawTable
        file = open('Documents\\5V_allPairs.txt', 'r')
        lines = file.readlines()
        file.close()

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
        self.readTable = []
        for read in rawTable:
            if read.startPos == currentPos:
                duplicate = False
                for i in range(readTablePos - sameStartPos, readTablePos-1):
                    if self.readTable[i].read == read.read:
                        self.readTable[i].count += 1
                        duplicate = True
                        break
                if not duplicate:
                    read.id = readTablePos
                    self.readTable.append(read)
                    sameStartPos = sameStartPos + 1
                    readTablePos = readTablePos + 1
            else:
                currentPos = read.startPos
                sameStartPos = 1
                read.id = readTablePos
                self.readTable.append(read)
                readTablePos = readTablePos + 1
        # In checkCorrection, we store which reads were corrected in this step.
        # The corrected reads will be excluded in the computation of the scores.
        self.checkCorrection = np.full(len(self.readTable), False)

    """
    Reads and stores the reference genome from given file
    """
    def createReferenceGenome(self):
        # Read reference genome
        file = open('Documents\\5V_ref_seq.fasta', 'r')
        lines = file.readlines()
        self.referenceGenome = lines[1]

    """
    This method computed the hamming distance between a read and a haplotype
    without counting a dot as a mistake
    @read read object
    @haplotype haplotype sequence
    """
    def similarity(self, read, haplotype):
        if read.startPos > len(haplotype):
            return None
        distance = 0
        for i in range(0, len(read.read)):
            if read.startPos + i < len(haplotype):
                if read.read[i] != haplotype[read.startPos + i - 1]:
                    if read.read[i] != ".":
                        distance += 1
            else:
                break
        return distance

    """
    This method applies the error correction using Shannon entropy
    """
    def entropyCorrection(self):
        self.createReadTable()
        self.createReferenceGenome()

        print("Calculating statistics")

        list = []
        dot = []
        coverage = []
        for i in range(0, len(self.referenceGenome)):
            list.append([0, 0, 0, 0, 0])
            dot.append(0)
            coverage.append(0)

        # Compute the statistics, how many times a nucleotide in {A, C, G, T, -} appears
        # at each position
        for i in range(0, len(self.readTable)):
            read = self.readTable[i]
            for k in range(0, len(read.read)):
                if read.read[k] == "A":
                    list[read.startPos + k - 1][0] += read.count
                    coverage[read.startPos + k - 1] += read.count
                elif read.read[k] == "C":
                    list[read.startPos + k - 1][1] += read.count
                    coverage[read.startPos + k - 1] += read.count
                elif read.read[k] == "G":
                    list[read.startPos + k - 1][2] += read.count
                    coverage[read.startPos + k - 1] += read.count
                elif read.read[k] == "T":
                    list[read.startPos + k - 1][3] += read.count
                    coverage[read.startPos + k - 1] += read.count
                elif read.read[k] == "-":
                    list[read.startPos + k - 1][4] += read.count
                    coverage[read.startPos + k - 1] += read.count
                else:
                    dot[read.startPos + k - 1] += read.count

        # compute the entropy at each position with the frequencies which
        # were calculated above
        entropy = []
        for i in range(0, len(list)):
            probA = 0
            probC = 0
            probG = 0
            probT = 0
            probS = 0
            if coverage[i] != 0:
                probA = list[i][0] / coverage[i]
                probC = list[i][1] / coverage[i]
                probG = list[i][2] / coverage[i]
                probT = list[i][3] / coverage[i]
                probS = list[i][4] / coverage[i]
            eA = 0
            eC = 0
            eG = 0
            eT = 0
            eS = 0
            if probA != 0:
                eA = probA * math.log2(probA)
            if probC != 0:
                eC = probC * math.log2(probC)
            if probG != 0:
                eG = probG * math.log2(probG)
            if probT != 0:
                eT = probT * math.log2(probT)
            if probS != 0:
                eS = probS * math.log2(probS)
            # normalize entropy
            entropy.append(-(eA + eC + eG + eT + eS) / 2.32192809)

        truePos = []
        nucleotides = ["A", "C", "G", "T", "-"]
        # Decide using a threshold, which position is a non-variant call position.
        # Store for them the nucleotide which appeared the most
        for i in range(0, len(self.referenceGenome)):
            # If normalized with 5, use 0.016 as threshold
            if entropy[i] < 0.03445:
                index = self.highestValue(list[i])
                truePos.append(nucleotides[index])
            else:
                truePos.append("")

        # Correct reads with the majority rule resp the above stored nucleotides
        for i in range(0, len(self.readTable)):
            for k in range(0, len(self.readTable[i].read)):
                if self.readTable[i].read[k] != "." and entropy[k + self.readTable[i].startPos - 1] < 0.03445:
                    if self.readTable[i].read[k] != truePos[k + self.readTable[i].startPos - 1]:
                        self.readTable[i].read = self.readTable[i].read[:k] + truePos[k + self.readTable[i].startPos - 1] + \
                                                 self.readTable[i].read[k + 1:]
                        self.checkCorrection[i] = True

        # Store the position which have variant call
        filteredList = []
        for i in range(0, len(entropy)):
            if entropy[i] >= 0.03445:
                filteredList.append(i)

        # Store reads in EntropyFile.txt
        file = open('Documents\\EntropyFile.txt', 'w')
        for i in range(0, len(self.readTable)):
            filteredRead = ""
            for k in range(0, len(filteredList)):
                if self.readTable[i].startPos <= filteredList[k] + 1 <= self.readTable[i].startPos + len(self.readTable[i].read) - 1:
                    filteredRead += self.readTable[i].read[(filteredList[k] + 1) - self.readTable[i].startPos]
            file.write(str(self.readTable[i].startPos) + " " + str(self.readTable[i].count) + " " + self.readTable[i].read + " " + filteredRead + "\n")
        file.close()

        # Store variant call positions in filteredPositions.txt
        file = open('Documents\\filteredPositions.txt', 'w')
        filteredPositions = ""
        for i in range(0, len(filteredList)):
            if i != len(filteredList) - 1:
                filteredPositions += str(filteredList[i]) + " "
            else:
                filteredPositions += str(filteredList[i])
        file.write(filteredPositions)
        file.close()
        print("Number of filtered Positions = ", len(self.referenceGenome) - len(filteredList))

        # Store which reads were corrected in this step in checkCorrection.npy
        np.save('Documents\\checkCorrection.npy', self.checkCorrection)
        return entropy, list, dot

    """
    This method returns the index at which the highest value is located
    in the list
    @list list containing values
    """
    def highestValue(self, list):
        maxValue = list[0]
        index = 0
        for i in range(1, len(list)):
            if list[i] is None:
                continue
            if list[i] > maxValue:
                maxValue = list[i]
                index = i
        return index

    """
    This method returns the index at which the lowest value is located
    in the list
    @list list containing values
    """
    def leastValue(self, list):
        minValue = list[0]
        index = 0
        for i in range(1, len(list)):
            if list[i] is None:
                continue
            if list[i] < minValue:
                minValue = list[i]
                index = i
        if minValue > 10:
            return index, True
        else:
            return index, False

    """
    This method was used to locate positions at which the sequencing machine made 
    non-uniform errors.
    @entropy a list containing the entropy values of the positions
    @list a list containing how many reads overlapped a position
    @dot a list contatining number of detected dots at each position
    """
    def test(self, entropy, list, dot):
        haplotypes = statistic.getHaplotypes()

        variantList = []
        nucleotideCount = []
        for i in range(0, len(haplotypes[0])):
            nucleotideCount.append([0, 0, 0, 0, 0])

        for i in range(0, len(nucleotideCount[0])):
            for j in range(0, len(haplotypes)):
                if haplotypes[j][i] == "A":
                    nucleotideCount[i][0] += 1
                elif haplotypes[j][i] == "C":
                    nucleotideCount[i][1] += 1
                elif haplotypes[j][i] == "G":
                    nucleotideCount[i][2] += 1
                elif haplotypes[j][i] == "T":
                    nucleotideCount[i][3] += 1
                else:
                    nucleotideCount[i][4] += 1

        for i in range(0, len(haplotypes[0])):
            str = ""
            if nucleotideCount[i][0] >= 1:
                str += "A"
            if nucleotideCount[i][1] >= 1:
                str += "C"
            if nucleotideCount[i][2] >= 1:
                str += "G"
            if nucleotideCount[i][3] >= 1:
                str += "T"
            if nucleotideCount[i][4] >= 1:
                str += "-"
            variantList.append(str)

        list1 = []
        for i in range(0, len(self.referenceGenome)):
            list1.append([False, False, False, False, False])

        for i in range(0, len(self.readTable)):
            read = self.readTable[i]
            for k in range(0, len(read.read)):
                if read.read[k] == "A":
                    list1[read.startPos + k - 1][0] = True
                elif read.read[k] == "C":
                    list1[read.startPos + k - 1][1] = True
                elif read.read[k] == "G":
                    list1[read.startPos + k - 1][2] = True
                elif read.read[k] == "T":
                    list1[read.startPos + k - 1][3] = True
                elif read.read[k] == "-":
                    list1[read.startPos + k - 1][4] = True

        suggestedVariants = []
        nucleotides = ["A", "C", "G", "T", "-"]
        for i in range(0, len(list1)):
            var = ""
            for k in range(0, 5):
                if list1[i][k] > 0:
                    var += nucleotides[k]
            suggestedVariants.append(var)

        for i in range(1, 9300):
            hapCount = [0, 0, 0, 0, 0]
            table = [0, 0, 0, 0, 0]
            dots = [0, 0, 0, 0, 0]
            problemCases = 0
            # print("-------------------------------------------------------------------------------------------------")
            print("i = ", i, ", count = ", list[i], ", entropy = ", entropy[i], ", variants = ", variantList[i],
                  ", dots = ", dot[i])
            for k in range(0, len(haplotypes)):
                if len(haplotypes[k]) > i:
                    print("hap", k, " has here ", haplotypes[k][i])
            for k in range(0, len(self.readTable)):
                if self.readTable[k].startPos > i + 1:
                    break
                distance = [0, 0, 0, 0, 0]
                if self.readTable[k].startPos <= i + 1 <= self.readTable[k].startPos + len(
                        self.readTable[k].read) - 1:
                    for t in range(0, len(haplotypes)):
                        distance[t] = self.similarity(self.readTable[k], haplotypes[t])

                    index, check = self.leastValue(distance)
                    hapCount[index] += self.readTable[k].count
                    if check:
                        table[index] += self.readTable[k].count
                        problemCases += 1
                    if self.readTable[k].read[(i + 1) - self.readTable[k].startPos] == ".":
                        dots[index] += self.readTable[k].count

            print("HapCount = ", hapCount)
            print("DotCount = ", dots)
            print("problematic Cases = ", problemCases)
            print("Looking at classes = ", table)


if __name__ == '__main__':
    correction = EntropyErrorCorrection()
    entropy, list, dot = correction.entropyCorrection()
    #correction.test(entropy, list, dot)
