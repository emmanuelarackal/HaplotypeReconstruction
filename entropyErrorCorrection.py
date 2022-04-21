import math

import statistic
from helperClasses import Read

class EntropyErrorCorrection:

    def __init__(self):
        self.readTable = []
        self.referenceGenome = ""

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
                for i in range(readTablePos - sameStartPos, readTablePos):
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
        print(len(self.readTable))

    def createReferenceGenome(self):
        # Read reference genome
        file = open('Documents\\5V_ref_seq.fasta', 'r')
        lines = file.readlines()
        self.referenceGenome = lines[1]

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
            entropy.append(-(eA + eC + eG + eT + eS) / 5)

        truePos = []
        nucleotides = ["A", "C", "G", "T", "-"]
        for i in range(0, len(self.referenceGenome)):
            if entropy[i] < 0.016:
                index = self.highestValue(list[i])
                truePos.append(nucleotides[index])
            else:
                truePos.append("")

        for i in range(0, len(self.readTable)):
            for k in range(0, len(self.readTable[i].read)):
                if self.readTable[i].read[k] != "." and entropy[k + self.readTable[i].startPos - 1] < 0.016:
                    self.readTable[i].read = self.readTable[i].read[:k] + truePos[k + self.readTable[i].startPos - 1] + \
                                             self.readTable[i].read[k + 1:]

        filteredList = []
        for i in range(0, len(entropy)):
            if entropy[i] >= 0.016:
                filteredList.append(i)


        file = open('Documents/EntropyFile.txt', 'w')

        # TODO: merge
        for i in range(0, len(self.readTable)):
            filteredRead = ""
            for k in range(0, len(filteredList)):
                if self.readTable[i].startPos <= filteredList[k] + 1 <= self.readTable[i].startPos + len(self.readTable[i].read) - 1:
                    filteredRead += self.readTable[i].read[(filteredList[k] + 1) - self.readTable[i].startPos]
            file.write(str(self.readTable[i].startPos) + " " + str(self.readTable[i].count) + " " + self.readTable[i].read + " " + filteredRead + "\n")
        file.close()

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
        return entropy, list, dot


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

        for i in range(7397, 7432):
            #if entropy[i] < 0.016 and len(variantList[i]) > 1:
            #if len(variantList[i]) > len(suggestedVariants[i]):
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
