import numpy
import numpy as np
import os


class StagewiseRegression:

    """
    This class is used to solve the regression by using the generalized monotone
    incremental forward stagewise regression concept of Hastie et al.
    """

    DEBUG = False

    def __init__(self, dbg, nlist, roh):
        self.p = np.zeros(0)
        self.h = np.zeros(0)
        self.y = np.zeros(0)
        self.dbg = dbg
        # nlist contains the positions of the reads in pairedEndReadList
        # which are allowed in the optimization.
        self.nlist = nlist
        self.roh = roh

    def setVariablesToInitialConfiguration(self):
        self.p = np.zeros(0)
        self.y = np.zeros(0)
        self.h = np.zeros(0)

    """
    This method computed the matrix P
    @border1 start position of haplotype/region
    @border2 end position of haplotype/region
    @haploFolder folder where p.npy needs to be stored.
    """
    def calculateP(self, border1, border2, haploFolder):
        self.p = np.zeros((len(self.nlist), len(self.dbg.haplotypes)))

        for i in range(0, len(self.nlist)):
            for j in range(0, len(self.dbg.haplotypes)):
                # If haplotype j and read i match, store the value 1
                # in the matrix P
                if self.similarity(self.nlist[i], j, border1, border2):
                    self.p[i, j] = 1
        np.save(haploFolder+'\\p.npy', self.p)
        return

    """
    This method computes whether a given read matches a given haplotype in the
    overlapping region
    @readId position of read in pairedEndReadList
    @haploId position of haplotypen in haplotypes list
    @border1 start position of haplotype/region
    @border2 end position of haplotype/region
    """
    def similarity(self, readId, haploId, border1, border2):
        if self.dbg.pairedEndReadList[readId].startPos < border1:
            start = self.dbg.pairedEndReadList[readId].startPos
        else:
            start = border1
        count = 0
        for i in range(start, len(self.dbg.variantPositions)):
            if self.dbg.pairedEndReadList[readId].endPos < i or border2 - 1 < i:
                break
            if self.dbg.pairedEndReadList[readId].startPos > i or border1 > i:
                continue

            if border1 == 0:
                posOfRead = i - self.dbg.pairedEndReadList[readId].startPos
                posOfHaplo = i - border1 - 1
            else:
                posOfRead = i - self.dbg.pairedEndReadList[readId].startPos
                posOfHaplo = i - border1
            if self.dbg.pairedEndReadList[readId].read[posOfRead] != self.dbg.haplotypes[haploId][posOfHaplo]:
                # If mismatch is found, return false indicating that they don't match
                if self.dbg.pairedEndReadList[readId].read[posOfRead] != ".":
                    return False
            count += 1
        # If only 0 or 1 comparisons have happened, return false
        if count < 2:
            return False
        return True

    """
    This method computes the error/cost of the regression
    @h fit vector h
    """
    def calculateFit(self, h):
        res = self.y - self.p.dot(h)
        return res.dot(res)

    """
    This method computes the y vector    
    @haploFolder folder where p.npy needs to be stored.
    @folder folder where documents are stored
    """
    def calculateY(self, haploFolder, folder):
        self.y = np.zeros(len(self.nlist))

        # Create the partition lists
        partitionSets = []
        a = self.dbg.pairedEndReadList[self.nlist[0]].startPos
        b = self.dbg.pairedEndReadList[self.nlist[0]].endPos
        set = []
        count = 0
        # serves for gaining further details on the partition sets
        list = []
        list1 = []
        list2 = []
        # computes the individual sets and stores them in partitionSets
        for i in range(0, len(self.nlist)):
            if self.dbg.pairedEndReadList[self.nlist[i]].startPos == a and \
                    self.dbg.pairedEndReadList[self.nlist[i]].endPos == b:
                if np.count_nonzero(self.p[i, :]) != 0:
                    count += 1
                    set.append(self.dbg.pairedEndReadList[self.nlist[i]].count)
            if self.dbg.pairedEndReadList[self.nlist[i]].startPos != a or \
                    self.dbg.pairedEndReadList[self.nlist[i]].endPos != b:
                if np.count_nonzero(self.p[i, :]) != 0:
                    count += 1
                    partitionSets.append(set)
                    list1.append(len(set))
                    if len(set) == 3:
                        list.append(len(set))
                    if len(set) > 3:
                        list2.append(len(set))
                    set = []
                    a = self.dbg.pairedEndReadList[self.nlist[i]].startPos
                    b = self.dbg.pairedEndReadList[self.nlist[i]].endPos
                    set.append(self.dbg.pairedEndReadList[self.nlist[i]].count)
        if set:
            partitionSets.append(set)
            list1.append(len(set))
            if len(set) == 3:
                list.append(len(set))
            if len(set) > 3:
                list2.append(len(set))
        if StagewiseRegression.DEBUG:
            print("Number of partition sets ", len(list1))
            print("Having length 3: ", len(list))
            print("Having length >3: ", len(list2))
            print("Average length in > 3: ", sum(list2) / len(list2))
            print("Having length <3: ", len(list1) - len(list) - len(list2))
            print("Average length ", sum(list1) / len(list1))

        pos = 0
        self.y = np.zeros(count)
        # normalizes a partition set and stores the values in vector y
        for i in range(0, len(partitionSets)):
            for j in range(0, len(partitionSets[i])):
                self.y[pos] = partitionSets[i][j] / sum(partitionSets[i])
                pos += 1

        i = 0
        # if there are situations in which a read in y doesn't match to any
        # haplotype remove them. The corresponding algorithm for y is already
        # included in the calculation of y.
        while i < len(self.p):
            if np.count_nonzero(self.p[i]) == 0:
                self.p = np.delete(self.p, i, 0)
                del self.nlist[i]
            else:
                i = i + 1
        np.save(haploFolder+'\\y.npy', self.y)
        np.save(haploFolder+'\\p.npy', self.p)

        # stores the positions of the reads in pairedEndReadList
        # which were actually used in optimization in nList.txt.
        # This is relevant for the extension.
        file = open(folder+'\\nList.txt', 'w')
        nText = ""
        for i in range(0, len(self.nlist)):
            if i == len(self.nlist) - 1:
                nText += str(self.nlist[i])
            else:
                nText += str(self.nlist[i]) + " "
        file.write(nText)
        file.close()

    """
    After calculating the first fit h, this method removes
    all entries in h having the value 0. Accordingly, all
    columns in P corresponding to the are also removed.
    """
    def removeRedundancyAfterH(self):
        i = 0
        while i < len(self.h):
            if self.h[i] == 0:
                self.p = np.delete(self.p, i, 1)
                self.h = np.delete(self.h, i, 0)
            else:
                i = i + 1

    """
    generalized monotone incremental forward stagewise regression concept
    @haploFolder
    """
    def calculateH(self, haploFolder):
        eps = 0.00001
        if np.count_nonzero(self.h) == 0:
            self.h = np.zeros(len(self.dbg.haplotypes))
        else:
            self.h = np.zeros(len(self.h))
        constraint = 0

        if np.all((self.p == 0)) and np.all((self.y == 0)):
            self.p = np.load(haploFolder+'\\p.npy')
            self.y = np.load(haploFolder+'\\y.npy')

        # helper matrix m to compute the derivative of the penalty
        # in faster way
        m = np.ones((len(self.h), len(self.h)))
        np.fill_diagonal(m, 0)

        while constraint < 1:
            # gradient parts
            firstDerivative = 2 * (np.transpose(self.p)).dot(self.p.dot(self.h) - self.y)
            secondDerivative = 2 * self.roh * (m.dot(self.h))
            maxValue = -1000000
            index = -1
            # select always highest gradient and increase at the end corresponding
            # entry in h
            for i in range(0, len(self.h)):
                value = -(firstDerivative[i] + secondDerivative[i])
                if value > maxValue:
                    maxValue = value
                    index = i
            self.h[index] += eps
            constraint = sum(self.h)
        if StagewiseRegression.DEBUG:
            print("h: ", self.h)
            print("cost: ", self.calculateFit(self.h))
        return 0

    """
    The pipe is used during read graph generation to reduce local
    haplotypes
    @border1 start position of region
    @border2 end position of region
    @haploFolder folder where the files about haplotypes needs to be stored resp are stored
    @folder folder where the documents are stored
    
    """
    def pipe(self, border1, border2, haploFolder, folder):
        self.calculateP(border1, border2, haploFolder)
        print("Calculated P")
        self.calculateY(haploFolder, folder)
        print("Calculated y")
        if border1 < 100:
            string1 = "0" + str(border1)
        else:
            string1 = str(border1)
        if border2 < 100:
            string2 = "0" + str(border2)
        else:
            string2 = str(border2)

        if os.path.exists(haploFolder+"\\LSR_"+string1+"_"+string2+".txt"):
            os.remove(haploFolder+"\\LSR_"+string1+"_"+string2+".txt")

        self.calculateH(haploFolder)
        file = open(haploFolder+'\\LSR_' + string1 + '_' + string2 + '.txt', 'w')
        count = 0
        for i in range(0, len(self.dbg.haplotypes)):
            if self.h[i] >= 0.001:
                file.write(self.dbg.haplotypes[i] + "\n")
                count += 1
        file.close()
        return count
