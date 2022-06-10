import math
import random
import distance
import scipy.special
import os
from collections import Counter

from helperClasses import Read
from helperClasses import Window
from helperClasses import HaplotypeClass
import numpy
import time
from threading import Thread, Lock

"""
This class is used to apply the clustering error correction. In the methods
pipe() and ECThreads.run() are explained, how it was used to correct the reads
"""
class ClusteringErrorCorrection:
    finishedThreads = []

    def __init__(self):
        self.readTable = []
        self.referenceGenome = ""
        self.correctTable = []
        self.filteredPosition = []

    """
    This method reads all reads from EntropyFile and stores their compact form in readTable
    """
    def createReadTable(self):
        # Read all reads and store them into rawTable
        file = open('Documents/EntropyFile.txt', 'r')
        lines = file.readlines()
        file.close()

        count = 0
        for line in lines:
            parts = line.split(' ')
            if len(parts) > 2:
                read = Read(count, parts[2], int(parts[0]), int(parts[1]))
                self.readTable.append(read)
                count = count + 1

        id = 0
        # Transform each read into compact form and stores them in readTable
        for i in range(0, len(self.readTable)):
            filteredStart, filteredEnd = self.getFilteredPositions(self.readTable[i].read,
                                                                   self.readTable[i].startPos)
            filteredRead = self.getFilteredRead(self.readTable[i].read, self.readTable[i].startPos,
                                                filteredStart, filteredEnd)
            self.readTable[i].filteredRead = filteredRead
            self.readTable[i].filteredStart = filteredStart
            # If -1 is returned, then these values won't be considered during similarity computation
            a = filteredRead.find(".")
            if a != -1:
                self.readTable[i].dotStart = filteredStart + a
                self.readTable[i].dotEnd = filteredStart + filteredRead.rfind(".")
            else:
                self.readTable[i].dotStart = -1
                self.readTable[i].dotEnd = -1
            id += 1

        for read in self.readTable:
            self.correctTable.append([read, [], ""])

    """
    This method stores the relevant positions in compact form from file
    """
    def createFilteredPositionList(self):
        file = open('Documents\\filteredPositions.txt', 'r')
        line = file.readlines()
        file.close()
        numbers = line[0].split(" ")
        for i in range(0, len(numbers)):
            self.filteredPosition.append(int(numbers[i]))

    """
    This method generates the reference genome from file
    """
    def createReferenceGenome(self):
        # Read reference genome
        file = open('Documents\\5V_ref_seq.fasta', 'r')
        lines = file.readlines()
        self.referenceGenome = lines[1]

    """
    This method returns the start and end position of a given read in compact form
    @read read sequence
    @startPos start position of the read
    """
    def getFilteredPositions(self, read, startPos):
        filteredStart = 0
        filteredEnd = 0
        for i in range(0, len(self.filteredPosition)):
            if filteredStart == 0 and startPos - 1 <= self.filteredPosition[i]:
                filteredStart = i
            if filteredEnd == 0 and startPos + len(read) - 2 <= self.filteredPosition[i]:
                filteredEnd = i - 1
                break
        if filteredEnd == 0:
            filteredEnd = len(self.filteredPosition) - 1
        return filteredStart, filteredEnd

    """
    This method return the compact form of a given read
    @read read object
    @readStart start position of the read
    @filteredStart start position of the read in compact form
    @filteredEnd end position of the read in compact form
    """
    def getFilteredRead(self, read, readStart, filteredStart, filteredEnd):
        filteredRead = ""
        for k in range(filteredStart, filteredEnd + 1):
            try:
                filteredRead += read[self.filteredPosition[k] + 1 - readStart]
            except Exception as e:
                print(read, readStart, readStart + len(read), self.filteredPosition[k], self.filteredPosition[k] + 1 - readStart)

        return filteredRead

    """
    This method creates a window
    @id id of window
    @startPos start position of window
    @endPos end position of window
    @startFilteredList start position in compact form
    @endFilteredList end position in compact form
    """
    def createWindow(self, id, startPos, endPos, startFilteredList, endFilteredList):
        print("Creating window in [", startPos, ", ", endPos, "]")
        window = Window(id, startPos, endPos, startFilteredList, endFilteredList)
        for read in self.readTable:
            if read.startPos >= window.endPos:
                break
            readEndPos = read.startPos + len(read.read)

            # A read is considered in a window, if 75% of read is inside window
            fraction = len(read.read)
            if read.startPos < window.startPos:
                fraction = fraction - (window.startPos - read.startPos)
            if readEndPos > window.endPos:
                fraction = fraction - (readEndPos - window.endPos)
            if fraction / len(read.read) > 0.75:
                window.list.append([read, -1, ""])
                continue
        return window

    """
    This method samples a class
    @classListPos position of class in windows classList
    @w1 window object
    @theta current theta
    """
    def sampleHaplotype(self, classListPos, w1, theta):
        hapStart = w1.startFilteredList
        hapEnd = w1.endFilteredList
        haplotype = w1.classList[classListPos]
        bigB = 5
        nucleotides = ["A", "C", "T", "G", "-"]

        haplotype.haplotype = ""
        # For each position, compute the frequency of the nucleotides given the read set of the class
        # Calculate the probabilities and sample nucleotide.
        for pos in range(hapStart, hapEnd + 1):
            nucleotideCount = [0, 0, 0, 0, 0]
            nucleotideProb = [0, 0, 0, 0, 0]
            overlappingReads = 0
            # Compute frequency of each nucleotide at position pos
            for i in range(0, len(haplotype.list)):
                if haplotype.list[i].startPos > self.filteredPosition[pos] + 1 or (
                        haplotype.list[i].startPos + len(haplotype.list[i].read) - 1) < self.filteredPosition[pos] + 1:
                    continue
                else:
                    nucleotide = haplotype.list[i].read[self.filteredPosition[pos] + 1 - haplotype.list[i].startPos]
                    if nucleotide == "A":
                        nucleotideCount[0] += haplotype.list[i].count
                        overlappingReads += haplotype.list[i].count
                    elif nucleotide == "C":
                        nucleotideCount[1] += haplotype.list[i].count
                        overlappingReads += haplotype.list[i].count
                    elif nucleotide == "T":
                        nucleotideCount[2] += haplotype.list[i].count
                        overlappingReads += haplotype.list[i].count
                    elif nucleotide == "G":
                        nucleotideCount[3] += haplotype.list[i].count
                        overlappingReads += haplotype.list[i].count
                    elif nucleotide == "-":
                        nucleotideCount[4] += haplotype.list[i].count
                        overlappingReads += haplotype.list[i].count
            if theta != 1:
                if overlappingReads == 0:
                    haplotype.haplotype += "*"
                else:
                    # Compute the log probabilities
                    for i in range(0, 5):
                        nucleotideProb[i] = nucleotideCount[i] * math.log(theta) + \
                                            (overlappingReads - nucleotideCount[i]) * math.log((1 - theta) / (bigB - 1))
                    # If you want, you can just use nucleotide with highest probability
                    #i = self.highestValue(nucleotideProb)

                    # We used logsumexp to normalize the log probabilities
                    # After that, we created a uniform distribution with b = sum of
                    # normalized vector. With this, we avoid that the highest log probability
                    # is always chosen. Instead, log probabilities which are a bit smaller
                    # then the highest can be chosen as well. Therefore, some randomness is added
                    # in the sampling

                    array = numpy.zeros(len(nucleotideProb))
                    for i in range(0, len(nucleotideProb)):
                        array[i] = nucleotideProb[i]
                    logSumExp = scipy.special.logsumexp(array)
                    pValues = numpy.exp(array - logSumExp)
                    indices = pValues.argsort()[::-1][:len(pValues)]
                    sum = math.fsum(pValues)
                    x = numpy.random.uniform(0, sum)
                    for t in range(0, len(indices)):
                        x -= pValues[indices[t]]
                        if x < 0:
                            break
                    haplotype.haplotype += nucleotides[t]

            if theta == 1:
                # If theta is 1 (only at the beginning) use majority rule to sample
                i = self.highestValue(nucleotideCount)
                # The case at the very beginning
                haplotype.haplotype += nucleotides[i]
        return

    """
    This class assigns a given read to a class
    @currentPos position of the read in list of window
    @w1 window object
    @theta current theta
    @gamma current gamma
    @alpha alpha
    """
    def assignToClass(self, currentPos, w1, theta, gamma, alpha):
        read = w1.list[currentPos][0]
        currentClassId = w1.list[currentPos][1]
        hapStart = w1.startFilteredList
        hapEnd = w1.endFilteredList
        bigB = 5
        currentIndex = -1
        logClass = []
        check = False
        # Goes through each class and computes log probability. They are stored in
        # logClass
        for i in range(0, len(w1.classList)):
            logP = 0
            # compute m and m'
            (similarity, hamDistance) = self.similarity(read, w1.classList[i].haplotype, hapStart, hapEnd)
            # If current class is viewed
            if w1.classList[i].id == currentClassId:
                currentIndex = i
                # If current class only contains of examined read, remove it and replace it with
                # the new class. Set check = true to remark that a class was already newly instantiated
                # Ensures, that class id doesn't get lost.
                if (w1.classList[i].count - read.count) == 0:
                    check = True
                    w1.classList[i].haplotype = self.getFilteredReferenceGenome(hapStart, hapEnd)
                    (similarity, hamDistance) = self.similarity(read, w1.classList[i].haplotype, hapStart, hapEnd)
                    if theta == 1 and gamma == 1:
                        print("theta and gamme are 1")
                    else:
                        logP = math.log(alpha) + math.log(read.count) + \
                               similarity * \
                               math.log(theta * gamma + (1 - gamma) * (1 - theta) / (bigB - 1)) + \
                               hamDistance * \
                               math.log((theta + gamma + bigB * (1 - gamma * theta) - 2) / math.pow((bigB - 1), 2))
                        logClass.append(logP)
                else:
                    # Current class contains further reads
                    try:
                        logP = math.log(w1.classList[i].count - read.count) + \
                               similarity * math.log(theta) + \
                               hamDistance * math.log((1 - theta) / (bigB - 1))
                    except Exception as e:
                        print("LOG ", w1.classList[i].count, ", ", read.count)
                    logClass.append(logP)
            else:
                # For other classes
                if theta != 1:
                    logP = math.log(w1.classList[i].count) + \
                           similarity * math.log(theta) + \
                           hamDistance * math.log((1 - theta) / (bigB - 1))
                    logClass.append(logP)
                else:
                    print("theta is 1")

        # If a class was not instantiated yet, do it now
        if not check:
            newClass = HaplotypeClass(w1.haploCount)
            newClass.haplotype = self.getFilteredReferenceGenome(hapStart, hapEnd)
            w1.classList.append(newClass)
            w1.haploCount += 1
            if theta == 1 and gamma == 1:
                print("theta and gamme are 1")
            else:
                (similarity, hamDistance) = self.similarity(read, newClass.haplotype, hapStart, hapEnd)
                logP = math.log(alpha) + math.log(read.count) + \
                       similarity * \
                       math.log(theta * gamma + (1 - gamma) * (1 - theta) / (bigB - 1)) + \
                       hamDistance * \
                       math.log((theta + gamma + bigB * (1 - gamma * theta) - 2) / math.pow((bigB - 1), 2))
                logClass.append(logP)

        # Use logSumExp and inversive method to decide to which class
        # the read is assigned to
        array = numpy.zeros(len(logClass))
        for i in range(0, len(logClass)):
            array[i] = logClass[i]
        logSumExp = scipy.special.logsumexp(array)
        pValues = numpy.exp(array - logSumExp)
        indices = pValues.argsort()[::-1][:len(pValues)]
        sum = math.fsum(pValues)
        x = numpy.random.uniform(0, sum)
        for t in range(0, len(indices)):
            x -= pValues[indices[t]]
            if x < 0:
                break
        index = indices[t]

        if currentIndex != -1:
            newID = w1.classList[currentIndex].id
        else:
            newID = -1
        # If class has changed, remove read from old class and add it to new class
        if index != currentIndex:
            if currentIndex != -1:
                w1.classList[currentIndex].list.remove(read)
                w1.classList[currentIndex].count -= read.count
            w1.classList[index].list.append(read)
            w1.classList[index].count += read.count
            newID = w1.classList[index].id
            # If old class only contained that read, delete it
            if check:
                if currentIndex != -1:
                    del w1.classList[currentIndex]
        # Delete newly instantiated class if it was not selected
        if not check and index != len(w1.classList) - 1:
            del w1.classList[len(w1.classList) - 1]
        # return id of new class
        return newID

    """
    This method computes and returns the compact form of the reference genome
    @hapStart start position of window/haplotype
    @hapEnd end position of window/haplotype
    """
    def getFilteredReferenceGenome(self, hapStart, hapEnd):
        filteredReference = ""
        for i in range(hapStart, hapEnd + 1):
            filteredReference += self.referenceGenome[self.filteredPosition[i]]
        return filteredReference

    """
    This method computes and returns the index at which the highest value is located 
    in the list
    @list list
    """
    def highestValue(self, list):
        maxValue = list[0]
        index = 0
        for i in range(1, len(list)):
            if list[i] > maxValue:
                maxValue = list[i]
                index = i
        return index

    """
    This method computes the number of matches and mismatches between a read and a haplotype.
    We use here no libraries. Therefore, not fast.
    @read read object
    @haplo sequence of haplotype
    @hapStart start position of window/haplotype
    @hapEnd end position of window/haplotype
    """
    def similarityNaive(self, read, haplo, hapStart, hapEnd):
        match = 0
        mismatch = 0
        count = 0
        for pos in range(hapStart, hapEnd + 1):
            if read.startPos - 1 > self.filteredPosition[pos]:
                count += 1
                continue
            elif (read.startPos + len(read.read) - 1) < self.filteredPosition[pos] + 1:
                break
            else:
                if read.read[self.filteredPosition[pos] + 1 - read.startPos] == ".":
                    count += 1
                    continue
                if read.read[self.filteredPosition[pos] + 1 - read.startPos] == haplo[count]:
                    match = match + 1
                else:
                    mismatch = mismatch + 1
            count += 1
        return match, mismatch

    """
    This method computes the number of matches and mismatches between a read and a haplotype.
    We use the library distance which makes the computation faster
    @read read object
    @haplo sequence of haplotype
    @hapStart start position of window/haplotype
    @hapEnd end position of window/haplotype
    """
    def similarity(self, read, haplo, hapStart, hapEnd):
        start = hapStart
        if read.filteredStart > hapStart:
            start = read.filteredStart
        end = hapEnd
        if read.filteredStart + len(read.filteredRead) - 1 < hapEnd:
            end = read.filteredStart + len(read.filteredRead) - 1

        # TODO: remove except part since not necessary
        try:
            mismatch = distance.hamming(read.filteredRead[start - read.filteredStart:end - read.filteredStart + 1],
                                    haplo[start - hapStart:end-hapStart + 1])
        except Exception as e:
            print("Exception")
            mismatch = distance.hamming(read.filteredRead[start - read.filteredStart:end - read.filteredStart + 1],
                                        haplo[start - hapStart:end - hapStart])
        a = 0
        # In distance.hamming, dots are seen as mistakes. Therefore, we need to make
        # sure that this doesn't influence number of matches and mismatches. This piece
        # of code ensures that dots are not considered in the computation.
        if read.dotStart != -1:
            if read.dotStart < start:
                a = start - read.dotStart
            b = 0
            if read.dotEnd > hapEnd:
                b = read.dotEnd - hapEnd
            match = end + 1 - start - mismatch
            mismatch = mismatch - (read.dotEnd + 1 - read.dotStart - a - b)
        else:
            match = end + 1 - start - mismatch
        return match, mismatch

    """
    This method updates theta
    @w1 window object
    """
    def updateTheta(self, w1):
        hapStart = w1.startFilteredList
        hapEnd = w1.endFilteredList
        match = 0
        numOfPos = 0
        readCount = 0
        for i in range(0, len(w1.classList)):
            for j in range(0, len(w1.classList[i].list)):
                count = 0
                readCount += 1
                for pos in range(hapStart, hapEnd + 1):
                    if w1.classList[i].list[j].startPos > self.filteredPosition[pos] + 1:
                        count += 1
                        continue
                    elif (w1.classList[i].list[j].startPos + len(
                            w1.classList[i].list[j].read) - 1) < self.filteredPosition[pos] + 1:
                        break
                    else:
                        numOfPos += w1.classList[i].list[j].count
                        if w1.classList[i].list[j].read[self.filteredPosition[pos] + 1 - w1.classList[i].list[j].startPos] == ".":
                            count += 1
                            continue
                        if w1.classList[i].list[j].read[self.filteredPosition[pos] + 1 - w1.classList[i].list[j].startPos] == w1.classList[i].haplotype[count]:
                            match += 1
                        count += 1
        # Since we compute in compact form, all other positions are automatically considered as correct
        return (match + (w1.endPos - w1.startPos - (hapEnd - hapStart)) * readCount) / ((w1.endPos - w1.startPos) * readCount)


    """
    This method updates gamma
    @w1 window object
    """
    def updateGamma(self, w1):
        hapStart = w1.startFilteredList
        hapEnd = w1.endFilteredList
        match = 0
        total = 0
        for i in range(0, len(w1.classList)):
            count = 0
            for pos in range(hapStart, hapEnd + 1):
                if w1.classList[i].haplotype[count] == self.referenceGenome[self.filteredPosition[pos]]:
                    match = match + 1
                count += 1
                total += 1
        # Since we compute in compact form, all other positions are automatically considered as correct
        return (match + (w1.endPos - w1.startPos - (hapEnd - hapStart)) * len(w1.classList)) / (
                    (w1.endPos - w1.startPos) * len(w1.classList))

    """
    This method serves to merge classes which have the same haplotype sequence
    @w1 window object
    """
    def mergeHaplotypes(self, w1):
        i = 0
        nameList = []
        for k in range(0, len(w1.classList)):
            nameList.append(w1.classList[k].haplotype)
        while True:
            if i >= len(w1.classList):
                break
            if nameList.count(w1.classList[i].haplotype) > 1:
                k = i + 1
                while True:
                    if k >= len(w1.classList):
                        break
                    # Detected two classes with same haplotype. Add the reads
                    # of the second class to the first class and remove second
                    if w1.classList[i].haplotype == w1.classList[k].haplotype:
                        for t in range(0, len(w1.classList[k].list)):
                            w1.classList[i].list.append(w1.classList[k].list[t])
                            w1.classList[i].count += w1.classList[k].list[t].count
                        for t in range(0, len(w1.list)):
                            if w1.list[t][1] == w1.classList[k].id:
                                w1.list[t][1] = w1.classList[i].id
                        del w1.classList[k]
                        nameList.remove(w1.classList[i].haplotype)
                    else:
                        k += 1
            else:
                i = i + 1

"""
The pipe controls the execution of the threads. In our case, 
we executed together 64 threads parallely (2 Laptops, each 4 main threads, each
main thread executed 8 ECThreads). The ranges of k and the initial value of i
were set accordingly. In the first step, the error correction was applied with window length 500 
on the reads from the file EntropyFile (EntropyErrorCorrection). In the second step, the output
file of the first iteration (ClusteringFile) was renamed to EntropyFile.txt. Then, the error 
correction was applied again with windows of length 600. In the third iteration, we used 
windows of length 700. The output file of it, ClusteringFile.txt, was used for the read graph
construction
"""
def pipe():

    # Here, we can set the starting position of the first window and the window length.
    i = 400
    windowLength = 500
    for k in range(0, 6):

        o1 = ClusteringErrorCorrection()
        o2 = ClusteringErrorCorrection()
        o3 = ClusteringErrorCorrection()
        o4 = ClusteringErrorCorrection()
        o5 = ClusteringErrorCorrection()
        o6 = ClusteringErrorCorrection()
        o7 = ClusteringErrorCorrection()
        o8 = ClusteringErrorCorrection()
        o9 = ClusteringErrorCorrection()

        o1.createFilteredPositionList()
        o1.createReadTable()
        o1.createReferenceGenome()

        o2.createFilteredPositionList()
        o2.createReadTable()
        o2.createReferenceGenome()

        o3.createReadTable()
        o3.createReferenceGenome()
        o3.createFilteredPositionList()

        o4.createReadTable()
        o4.createReferenceGenome()
        o4.createFilteredPositionList()

        o5.createReadTable()
        o5.createReferenceGenome()
        o5.createFilteredPositionList()

        o6.createReadTable()
        o6.createReferenceGenome()
        o6.createFilteredPositionList()

        o7.createReadTable()
        o7.createReferenceGenome()
        o7.createFilteredPositionList()

        o8.createReadTable()
        o8.createReferenceGenome()
        o8.createFilteredPositionList()

        o9.createReadTable()
        o9.createReferenceGenome()
        o9.createFilteredPositionList()

        w1 = ECThread(i, i, i+windowLength, 150, o1)
        i += 100

        w2 = ECThread(i, i, i+windowLength, 150, o2)
        i += 100

        w3 = ECThread(i, i, i+windowLength, 150, o3)
        i += 100

        w4 = ECThread(i, i, i+windowLength, 150, o4)
        i += 100

        w5 = ECThread(i, i, i+windowLength, 150, o5)
        i += 100

        w6 = ECThread(i, i, i+windowLength, 150, o6)
        i += 100

        w7 = ECThread(i, i, i+windowLength, 150, o7)
        i += 100

        w8 = ECThread(i, i, i+windowLength, 150, o8)
        i += 100

        w9 = ECThread(i, i, i+windowLength, 150, o9)
        i += 100

        w1.start()
        w2.start()
        w3.start()
        w4.start()
        w5.start()
        w6.start()
        w7.start()
        w8.start()
        w9.start()
        # Adapt here according to the number of threads
        while len(ClusteringErrorCorrection.finishedThreads) < 9:
            time.sleep(1000)
        print("All threads of ", k, " are finished")
        w1.join()
        w2.join()
        w3.join()
        w4.join()
        w5.join()
        w6.join()
        w7.join()
        w8.join()
        w9.join()

        ClusteringErrorCorrection.finishedThreads = []

    o = ClusteringErrorCorrection()
    o.createReadTable()

    # This part serves to apply error correction based on majority rule
    # and to create the file ClusteringFile.txt
    for file in os.listdir("Data"):
        if file.startswith("ClusteringFile"):
            f = open("Data\\"+file, 'r')
            lines = f.readlines()
            f.close()
            count = 0
            for line in lines:
                parts = line.split(' ')
                if len(parts[2]) > 1:
                    o.correctTable[count][1].append(parts[2])
                count += 1

    for elem in o.correctTable:
        for i in range(0, len(elem[0].read)):
            if elem[0].read[i] != ".":
                nucleotideCount = [0, 0, 0, 0, 0]
                nucleotides = ["A", "C", "T", "G", "."]
                # Compute statistics for a position of a read
                for h in range(0, len(elem[1])):
                    if len(elem[1][h]) - 1 > i:
                        if elem[1][h][i] == "A":
                            nucleotideCount[0] += 1
                        elif elem[1][h][i] == "C":
                            nucleotideCount[1] += 1
                        elif elem[1][h][i] == "T":
                            nucleotideCount[2] += 1
                        elif elem[1][h][i] == "G":
                            nucleotideCount[3] += 1
                        elif elem[1][h][i] == "-":
                            nucleotideCount[4] += 1
                        else:
                            continue
                maxValue = -1
                index = -1
                # majority rule
                for t in range(0, len(nucleotideCount)):
                    if nucleotideCount[t] > maxValue:
                        maxValue = nucleotideCount[t]
                        index = t
                if index != -1 and maxValue > 0:
                    elem[2] += nucleotides[index]
                else:
                    elem[2] += elem[0].read[i]
            else:
                elem[2] += "."

    # Store corrected reads in this file
    file = open('Documents\\ClusteringFile.txt', 'w')

    for i in range(0, len(o.correctTable)):
        file.write(str(o.correctTable[i][0].startPos) + " " + str(o.correctTable[i][0].count) + " " + o.correctTable[i][2] + "\n")
    file.close()
    return

"""
This thread serves to apply the error correction
"""
class ECThread(Thread):
    __lock = Lock()

    def __init__(self, id, windowStart, windowEnd, iteration, rc):
        Thread.__init__(self)
        self.id = id
        self.windowStart = windowStart
        self.windowEnd = windowEnd
        self.iteration = iteration
        self.rc = rc

    """
    The windowing and the iterative approach are executed in the run.
    """
    def run(self):
        try:
            # If necessary, change alpha.
            alpha = 0.04
            theta = 1

            filteredStart = 0
            filteredEnd = 0
            for i in range(0, len(self.rc.filteredPosition)):
                if filteredStart == 0 and self.windowStart - 1 <= self.rc.filteredPosition[i]:
                    filteredStart = i
                if filteredEnd == 0 and self.windowEnd - 1 <= self.rc.filteredPosition[i]:
                    filteredEnd = i - 1
                    break
            if filteredEnd == 0:
                filteredEnd = self.rc.filteredPosition[len(self.rc.filteredPosition) - 1]
            # Create window
            w1 = self.rc.createWindow(1, self.windowStart, self.windowEnd, filteredStart, filteredEnd)
            k = 10

            print("Total Number of reads: ", len(self.rc.readTable))
            print("Number of reads: ", len(w1.list))
            # create k classes/haplotypes
            for w1.haploCount in range(0, k):
                w1.classList.append(HaplotypeClass(w1.haploCount))
            w1.haploCount += 1

            # assign the reads randomly to classes
            print("Assigning reads randomly to classes")
            for i in range(0, len(w1.list)):
                randomClassId = random.randint(0, k - 1)
                w1.classList[randomClassId].list.append(w1.list[i][0])
                w1.classList[randomClassId].count += w1.list[i][0].count
                w1.list[i][1] = randomClassId

            # remove all classes which didn't get any read
            breaker = False
            i = 0
            while not breaker:
                if i >= len(w1.classList):
                    break
                if len(w1.classList[i].list) == 0:
                    del w1.classList[i]
                else:
                    i = i + 1

            # sample haplotypes
            for i in range(0, len(w1.classList)):
                self.rc.sampleHaplotype1(i, w1, theta)

            # Update theta and gamma
            theta = self.rc.updateTheta(w1)
            gamma = self.rc.updateGamma(w1)

            for k in range(1, self.iteration):
                # Iteratively, assign reads to classes, sample haplotypes, update theta and gamma
                for i in range(0, len(w1.list)):
                    newClassId = self.rc.assignToClass(i, w1, theta, gamma, alpha)
                    w1.list[i][1] = newClassId
                    # After a burn-in phase, store for each read id of the class to
                    # which it was assigned to.
                    if k >= 100:
                        w1.list[i][2] += " " + str(newClassId)
                for i in range(0, len(w1.classList)):
                    self.rc.sampleHaplotype1(i, w1, theta)
                self.rc.mergeHaplotypes(w1)
                theta = self.rc.updateTheta(w1)
                gamma = self.rc.updateGamma(w1)

            for i in range(0, len(w1.list)):
                correctedRead = ""
                classId = -1
                list = w1.list[i][2].split(" ")
                # Find position of current classId
                for k in range(0, len(w1.classList)):
                    if w1.list[i][1] == w1.classList[k].id:
                        classId = k
                        break
                # Check how many times the current classId was mentioned in the records
                currentCount = list.count(str(w1.classList[classId].id))
                # If it is lower than a half, check whether another class was choosen more often
                if currentCount < self.iteration * 0.5:
                    c = Counter(list)
                    # find most common chosen class
                    mostCommon = c.most_common(1)
                    possibleClassId = int(mostCommon[0][0])
                    if mostCommon[0][1] > currentCount:
                        # Check whether that class still exists
                        for k in range(0, len(w1.classList)):
                            if possibleClassId == w1.classList[k].id:
                                # If so, choose that class for error correction
                                w1.classList[k].count += w1.list[i][0].count
                                w1.classList[classId].count -= w1.list[i][0].count
                                classId = k
                                break
                # For each read in the window, we store the haplotype in the following way:
                # If there is a nucleotide in the read, which positionwise doesn't overlap
                # with the haplotype, then _ is stored in the file.
                # If there is a nucleotide in the read, which positionwise does overlap
                # with the haplotype, then the nucleotide of the haplotype is stored in the file.
                # If however the read has a dot at this position, a _ is stored.
                if w1.list[i][0].startPos < w1.startPos:
                    for base in range(0, w1.startPos - w1.list[i][0].startPos):
                        correctedRead += "_"
                count = 0
                for base in range(0, w1.endPos - w1.startPos):
                    if w1.list[i][0].startPos > base + w1.startPos:
                        if base + w1.startPos - 1 in self.rc.filteredPosition:
                            count += 1
                        continue
                    elif (w1.list[i][0].startPos + len(w1.list[i][0].read) - 1) < base + w1.startPos:
                        break
                    else:
                        if base + w1.startPos - 1 in self.rc.filteredPosition:
                            if w1.classList[classId].haplotype[count] != "*" and w1.classList[classId].haplotype[count] != ".":
                                correctedRead += w1.classList[classId].haplotype[count]
                            else:
                                correctedRead += "_"
                            count += 1
                        else:
                            correctedRead += "_"
                self.rc.correctTable[w1.list[i][0].id][2] = correctedRead

            print("Steps are over for ", self.windowStart)
            for i in range(0, len(w1.classList)):
                print("Class has ", w1.classList[i].count, " reads")
            for i in range(0, len(w1.classList)):
                print("Class has ", w1.classList[i].haplotype)

            # Store the above generated correction string into file
            fileName = "Data\\ClusteringFile400_" + str(self.windowStart) + ".txt"
            file = open(fileName, 'w')
            for i in range(0, len(self.rc.correctTable)):
                file.write(str(self.rc.correctTable[i][0].startPos) + " " + str(self.rc.correctTable[i][0].count) + " " +
                           self.rc.correctTable[i][2] + "\n")
            file.close()
        except Exception as e:
            print("Problem occurred in thread ", self.id)
        ClusteringErrorCorrection.finishedThreads.append(0)


if __name__ == '__main__':
    pipe()
