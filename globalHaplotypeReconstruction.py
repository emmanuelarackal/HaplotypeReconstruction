import numpy as np
from stagewiseregression import StagewiseRegression
from lagrangianRegression import LagrangianRegression
from haplotypeReconstruction_DBG import HaplotypeReconstruction_DBG
from helperClasses import FilteredRead
import os
from threading import Thread, Lock
import time

"""
This class was used to compute the optimal global haplotypes.
"""


class GlobalHaplotypeReconstruction:
    frequencies = np.zeros(0)
    finishedThreads = []

    def __init__(self):
        self.haplotypes = []
        self.frequencies = np.zeros(0)
        self.pairedEndReadList = []
        self.nList = []
        self.roh = [0.1, 1, 2, 3.5, 5, 7, 10]
        self.variantPositions = []

    """
    This method is used to store the global haplotypes in self.haplotypes, 
    to generate the pairedEndReadList, the variantPositions containing the 
    positions in compact form, the nList containing the 
    positions of the read in pairedEndReadList relevant for optimization
    @border1 start position of global region in compact form
    @border2 end position of global region in compact form
    @folder folder where pairedEndReadList file is stored
    @gFolder folder where the global/local haplotype files are stored
    """

    def initializeVariables(self, border1, border2, folder, gFolder):
        # TODO: As user, define which file you use.
        # Either use the global haplotypes or the consistent haplotypes. The consistent are those
        # who passed the consistency check.
        file = open(gFolder + "\\GlobalHaplotypes.txt", 'r')
        # file = open(gFolder + "\\ConsistentGlobalHaplotypes.txt", 'r')
        lines = file.readlines()
        file.close()
        for line in lines:
            if len(line) > 1:
                line = line.replace("\n", "")
                self.haplotypes.append(line)

        self.frequencies = np.zeros((len(self.haplotypes), 700))

        file = open(folder + '\\PairedEndReadList.txt', 'r')
        lines = file.readlines()
        file.close()
        count = 0
        self.pairedEndReadList = []
        for line in lines:
            parts = line.split(' ')
            if len(parts) > 1:
                parts[3] = parts[3].replace(" ", "")
                parts[3] = parts[3].replace("\n", "")
                read = FilteredRead(count, parts[3], int(parts[0]), int(parts[1]), int(parts[2]))
                self.pairedEndReadList.append(read)
                count = count + 1

        pos = 0
        for read in self.pairedEndReadList:
            a = 0
            b = 0

            if read.startPos > border1:
                a = read.startPos - border1
            if read.startPos + len(read.read) < border2:
                b = border2 - (read.startPos + len(read.read))
            # As user, you can specifiy by how many positions a read should
            # overlap the region
            if (border2 - border1 - a - b) >= 6:
                self.nList.append(pos)
            pos += 1

        file = open(folder + '\\variantPositions.txt', 'r')
        line = file.readlines()
        file.close()
        numbers = line[0].split(" ")
        for i in range(0, len(numbers)):
            self.variantPositions.append(int(numbers[i]))

    """
    Given a set, find the best fit. Is used to find the best fit in 
    a bucket K. Returns then the index containing best fit.
    @set containing the costs of the fits
    """

    def findBestFitInSet(self, set):
        bestFit = 100000
        index = -1
        for i in range(0, len(set)):
            if set[i][1] < bestFit:
                bestFit = set[i][1]
                index = i
        return set[index]

    """
    Returns the compact position of a given position
    @originalPos position for which we want to compute the compact form
    """

    def getVariantPosition(self, originalPos):
        for i in range(0, len(self.variantPositions)):
            if originalPos <= self.variantPositions[i]:
                return i
        return len(self.variantPositions) - 1

    """
    This method prints out the scores for the haplotypes in the result. It uses the 
    file Solution.txt which contains the haplotypes and their frequency estimates. Then,
    it uses minimum credibility to compute the score.
    @border1 start position of global region in compact form
    @border2 end position of global region in compact form
    @folder folder where pairedEndReadList file is stored
    @gFolder folder where the global/local haplotype files are stored
    @sFolder folder where Solution.txt is stored
    """

    def defineScore(self, border1, border2, folder, gFolder, sFolder):
        file = open(sFolder + '\\Solution.txt')
        lines = file.readlines()
        solutions = []
        for line in lines:
            elements = line.split(" ")
            solutions.append([elements[1], float(elements[0])])
        dbg = HaplotypeReconstruction_DBG()
        dbg.getFilteredReadTable(gFolder)
        dbg.createVariantPos(gFolder)
        # We need to use the raw reads. The information about them are stored
        # in tracker.txt
        file = open(folder + '\\tracker.txt')
        lines = file.readlines()
        file.close()
        tracker = []
        for line in lines:
            line = line.replace("\n", "")
            List = []
            elements = line.split(" ")
            if elements[0] == "_":
                tracker.append(List)
                continue
            for k in range(0, len(elements)):
                if elements[k] == "":
                    continue
                List.append(int(elements[k]))
            tracker.append(List)
        p = np.load(gFolder + '\\p.npy')

        # Either, we can use the positions where splits where observed in the
        # read graphs
        """file = open(folder + '\\checkPositions.txt', 'r')
        line = file.readlines()
        file.close()
        numbers = line[0].split(" ")

        for i in range(0, len(numbers)):
            dbg.rawPositionList.append(int(numbers[i]))"""

        # Or we can just look at all positions
        for i in range(border1, border2):
            dbg.rawPositionList.append(i)
        dbg.rawPositionList.sort(key=lambda x: x, reverse=False)

        # In nList, we store the reads which had an entry in vector y
        self.nList = []
        file = open(gFolder + '\\nList.txt', 'r')
        line = file.readlines()
        file.close()
        numbers = line[0].split(" ")
        for i in range(0, len(numbers)):
            self.nList.append(int(numbers[i]))

        nList = []
        row, col = np.shape(p)
        # Go trough each y entry. If it had a haplotype which mached with it,
        # store its raw reads in nList. We can check that by looking if the
        # row of p at the location of the read has any 1 entry.
        for i in range(0, row):
            if np.count_nonzero(p[i, :]) != 0:
                nList += tracker[self.nList[i]]
        # Store the haplotypes mentioned in Solution.txt in haplotypes
        for i in range(0, len(solutions)):
            dbg.haplotypes.append(solutions[i][0])

        # Use consistency check to compute how many reads extend the haplotypes at each position
        totalExtension, haploStatistics = dbg.removeInconsistentHaplotypes(border1, border2, 2, nList)

        print("Total Extension: ", totalExtension)
        for i in range(0, len(haploStatistics)):
            print(i, " ", haploStatistics[i])
        score = []
        # Apply scoring system
        for i in range(0, len(solutions)):
            worstCase = 1
            scoreNew = []
            count = 0
            for j in range(0, len(dbg.rawPositionList)):
                # If expected value differ from actual value, compute the difference.
                # Store the least value of 1-difference as score
                if solutions[i][1] * totalExtension[j] > haploStatistics[i][j]:
                    difference = (haploStatistics[i][j] - solutions[i][1] * totalExtension[j]) / (
                                solutions[i][1] * totalExtension[j])
                    scoreNew.append(difference)
                    count += 1
                    print(i, j, solutions[i][1], totalExtension[j], solutions[i][1] * totalExtension[j],
                          haploStatistics[i][j], difference)
                    if worstCase > 1 + difference:
                        worstCase = 1 + difference
            score.append(worstCase)
        print("worst case score: ", score)

    """
    This method applies the consistency check method.
    @border1 start position of global region in compact form
    @border2 end position of global region in compact form
    @folder folder where pairedEndReadList file is stored
    @gFolder folder where the global/local haplotype files are stored
    """

    def consistencyCheckWithRawData(self, folder, gFolder, border1, border2):
        dbg = HaplotypeReconstruction_DBG()
        dbg.getFilteredReadTable(folder)
        dbg.getPairedEndReadList(folder)
        dbg.createVariantPos(folder)

        dbg.positionList = []

        # Either, use only the split positions
        """file = open(folder+'\\checkPositions.txt', 'r')
        line = file.readlines()
        file.close()
        numbers = line[0].split(" ")
        for i in range(0, len(numbers)):
            dbg.specialList.append(int(numbers[i]))

        for fileName in os.listdir("H"):
            if fileName.startswith("LSR"):
                rep = fileName.replace("LSR_", "")
                rep = rep.replace(".txt", "")
                borders = rep.split("_")
                if borders[1][0] == "0":
                    borders[1] = borders[1][1:]
                if int(borders[1]) not in dbg.specialList:
                    dbg.specialList.append(int(borders[1]))"""

        # Or use all positions
        for i in range(border1, border2):
            dbg.specialList.append(i)

        haploStatistics = []

        # Apply consistency check and compute for each haplotype, how many reads
        # extends its positions
        for i in range(0, len(self.haplotypes)):
            dbg.haplotypes = [self.haplotypes[i]]
            dbg.pairedEndReadList = self.pairedEndReadList
            stat = dbg.removeInconsistentHaplotypes(border1, border2, 0, self.nList)
            haploStatistics.append(stat)

        threshold = []

        # Compute for each position the largest found number of extending reads
        for j in range(0, len(dbg.specialList)):
            value = 0
            for i in range(0, len(haploStatistics)):
                if value < haploStatistics[i][j]:
                    value = haploStatistics[i][j]
            threshold.append(value)

        print("threshold: ", threshold)
        # If position in haplotype is supported by at max 1 read, remove it
        i = 0
        while i < len(self.haplotypes):
            check = True
            for k in range(0, len(haploStatistics[i])):
                if haploStatistics[i][k] < 2:
                    del self.haplotypes[i]
                    del haploStatistics[i]
                    check = False
                    break
            if check:
                i = i + 1
        print("Length of haplotypes after removing the ones only supported at max by 1 raw read: ",
              len(self.haplotypes))
        i = 0
        # Apply consistency check
        while i < len(self.haplotypes):
            check = True
            for k in range(1, len(haploStatistics[i])):
                # If there is a decrease in the number of reads extending two neighbouring position [left of current pos],
                # that is larger than 99%, remove haplotype
                if (haploStatistics[i][k] - haploStatistics[i][k - 1]) / haploStatistics[i][k - 1] < -0.99:
                    print("Deleting Big Decrease in Flow ", i, ": ", haploStatistics[i][k], haploStatistics[i][k - 1])
                    del haploStatistics[i]
                    del self.haplotypes[i]
                    check = False
                    break
                # If multiple inferior, remove haplotype
                if check and (haploStatistics[i][k] - threshold[k]) / threshold[k] < -0.99:
                    print("Deleting Big Decrease in Max ", k, ": ", haploStatistics[i][k], threshold[k], k)
                    del haploStatistics[i]
                    del self.haplotypes[i]
                    check = False
                    break
                # If there is a decrease in the number of reads extending two neighbouring position [right of current pos],
                # that is larger than 99%, remove haplotype
                if check and k < len(haploStatistics[i]) - 1 and \
                        (haploStatistics[i][k] - haploStatistics[i][k + 1]) / haploStatistics[i][k + 1] < -0.99:
                    print("Deleting Big Jump ", i, ": ", haploStatistics[i][k], haploStatistics[i][k + 1])
                    del haploStatistics[i]
                    del self.haplotypes[i]
                    check = False
                    break
            if check:
                i = i + 1
        for i in range(0, len(haploStatistics)):
            print(haploStatistics[i])

        # Store in this file the haplotypes which passed consistency check
        file = open(gFolder + '\\ConsistentGlobalHaplotypes.txt', 'w')
        for i in range(0, len(self.haplotypes)):
            file.write(self.haplotypes[i] + "\n")
        file.close()

    """
    After we computed all frequency estimates, this method is applied. It computes
    the best fit.
    @gFolder folder where the global/local haplotype files are stored
    @fFolder folder where frequency files are stored
    @sFolder folder where best fit file needs to be stored
    """
    def calculateBestFit(self, gFolder, fFolder, sFolder):
        set = []
        sk = []
        reg = LagrangianRegression(self, self.nList, 0, 0, 2, 0.0001, 1.0001, 0.005)
        reg.p = np.load(gFolder + '\\p.npy')
        reg.y = np.load(gFolder + '\\y.npy')
        i = 0
        self.frequencies = np.zeros((len(self.haplotypes), 700))
        # Go through each file and store the frequency estimates in frequencies
        for file in os.listdir(fFolder):
            f = open(fFolder + "\\" + file, 'r')
            firstsplit = file.split("_")
            number = firstsplit[1].split(".")
            freqpos = int(number[0])
            lines = f.readlines()
            f.close()
            count = 0
            for line in lines:
                parts = line.replace("\n", "")
                try:
                    self.frequencies[count][freqpos] = float(parts)
                except Exception as e:
                    print("Exception")
                count += 1
            print(count)

        # we can have at max number of haplotypes non zero values
        for i in range(0, len(self.frequencies[0])):
            set.append(sk)

        # Compute the bins according to the number of non-zero entries
        for i in range(0, len(self.frequencies)):
            k = np.count_nonzero(self.frequencies[:, i])
            vector = self.frequencies[:, i]
            if np.count_nonzero(vector) != 0:
                set[k].append([vector, reg.calculateFit(vector)])

        solution = []
        # Find for each bin the best fit
        for i in range(0, len(set)):
            if len(set[i]) > 0:
                solution.append(self.findBestFitInSet(set[i]))

        previousFit = 100000
        index = -1
        # Find best fit
        for i in range(0, len(solution)):
            currentFit = solution[i][1]
            if currentFit > previousFit * 0.9:
                index = i
                break
            previousFit = currentFit

        print("Most suitable Haplotypes: ")
        file = open(sFolder + "\\Solution.txt", 'w')
        # Store best fit containing frequency and haplotype in Solution.txt
        for i in range(0, len(self.haplotypes)):
            if solution[index][0][i] > 0.001:
                print(solution[index][0][i], self.haplotypes[i])
                s = str(solution[index][0][i]) + " " + self.haplotypes[i]
                file.write(s + "\n")


"""
This method controls the computation of the fits.
@border1 start position of global region in compact form
@border2 end position of global region in compact form
@gbr object of class
@folder folder where pairedEndReadList file is stored
@gFolder folder where the global/local haplotype files are stored
@fFolder folder where frequencies needs to be stored
"""
def regressionPipe(border1, border2, gbr, folder, gFolder, fFolder):
    roh = [0.1, 1, 2, 3.5, 5, 7, 10]
    GlobalHaplotypeReconstruction.frequencies = np.zeros((len(gbr.haplotypes), 700))
    # Compute the matrix P and vector y
    reg = StagewiseRegression(gbr, gbr.nList, 0)
    reg.calculateP(border1, border2, gFolder)
    reg.calculateY(gFolder)
    # We can either compute all 700 fits. In my case, I just used some fits for each rho value
    for i in range(0, len(roh)):
        for k in range(0, 100):
            # Use 5 Threads to make computation parallel
            print("Starting ", k)
            o1 = GlobalHaplotypeReconstruction()
            o2 = GlobalHaplotypeReconstruction()
            o3 = GlobalHaplotypeReconstruction()
            o4 = GlobalHaplotypeReconstruction()
            o5 = GlobalHaplotypeReconstruction()

            w1 = FCThread(o1, roh[i], i * 100 + k, border1, border2, folder, gFolder, fFolder)
            w2 = FCThread(o2, roh[i], i * 100 + k, border1, border2, folder, gFolder, fFolder)
            w3 = FCThread(o3, roh[i], i * 100 + k, border1, border2, folder, gFolder, fFolder)
            w4 = FCThread(o4, roh[i], i * 100 + k, border1, border2, folder, gFolder, fFolder)
            w5 = FCThread(o5, roh[i], i * 100 + k, border1, border2, folder, gFolder, fFolder)

            w1.start()
            w2.start()
            w3.start()
            w4.start()
            w5.start()

            while len(GlobalHaplotypeReconstruction.finishedThreads) < 5:
                time.sleep(1000)

            print("All threads of ", i, k, " are finished")
            w1.join()
            w2.join()
            w3.join()
            w4.join()
            w5.join()
            GlobalHaplotypeReconstruction.finishedThreads = []

    GlobalHaplotypeReconstruction.frequencies = []


def pipe():
    # TODO: As user, you need to define these variables
    # Start position of region
    border1 = 32
    # End position of region
    border2 = 358
    # folder where files with reads are
    folder = "Documents"
    # folder where the haplotye files are stored
    gFolder = "Haplotypes"
    # folder where the solution is stored
    sFolder = "Solution"
    # folder where the frequency estimates are stored
    fFolder = "Frequencies"
    # If 0, no consistency check. If 1, consistency check
    mode = 0

    gbr = GlobalHaplotypeReconstruction()
    gbr.initializeVariables(border1, border2, folder, gFolder)
    if mode == 0:
        # without consistency check. Compute the first and cacluclate optimal solution
        regressionPipe(border1, border2, gbr, folder, gFolder, fFolder)
        gbr.calculateBestFit(gFolder, fFolder, sFolder)
    else:
        # With consistency check. Apply it, compute the fits and apply scoring system
        gbr.consistencyCheckWithRawData(folder, gFolder, border1, border2)
        regressionPipe(border1, border2, gbr, folder, gFolder, fFolder)
        gbr.calculateBestFit(gFolder, fFolder, sFolder)
        gbr.defineScore(border1, border2, folder, gFolder, sFolder)


"""
This is a thread class. It is used to compute a fit.
"""
class FCThread(Thread):
    __lock = Lock()
    """
    @gbr object of class
    @roh rho value which should be used in the regression
    @freqPos name of file
    @border1 start position of region
    @border2 end position of region
    @folder location where the reads are stored
    @gFolder folder where the global/local haplotype files are stored
    @fFolder folder where frequencies needs to be stored
    """
    def __init__(self, gbr, roh, freqPos, border1, border2, folder, gFolder, fFolder):
        Thread.__init__(self)
        self.gbr = gbr
        self.roh = roh
        self.freqPos = freqPos
        self.border1 = border1
        self.border2 = border2
        self.folder = folder
        self.gFolder = gFolder
        self.fFolder = fFolder

    """
    Run a computation for a fit
    """
    def run(self):
        try:
            index = []
            self.gbr.initializeVariables(self.border1, self.border2, self.folder, self.gFolder)
            reg = LagrangianRegression(self.gbr, self.gbr.nList, self.roh, 0, 2, 0.0001, 1.0001, 0.005)
            # load p and y
            reg.p = np.load(self.gFolder + "\\p.npy")
            reg.y = np.load(self.gFolder + "\\y.npy")
            # compute fit h
            reg.calculateH(self.gFolder)
            # remove all entries containing the value 0
            for j in range(0, len(reg.h)):
                if reg.h[j] != 0:
                    index.append(j)
            reg.removeRedundancyAfterH()
            reg.roh = 0
            # Compute h#
            reg.calculateH(self.gFolder)
            # Store frequency estimates in file
            file = open(self.fFolder + "\\Frequency_" + str(self.freqPos) + ".txt", 'w')
            count = 0
            for i in range(0, len(self.gbr.haplotypes)):
                if i in index:
                    file.write(str(reg.h[count]) + "\n")
                    count += 1
                else:
                    file.write("0\n")
            file.close()
        except Exception as e:
            print("Problem occured for ", self.freqPos)
        GlobalHaplotypeReconstruction.finishedThreads.append(1)

if __name__ == '__main__':
    pipe()
