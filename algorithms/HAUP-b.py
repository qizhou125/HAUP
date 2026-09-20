import math
import sys
import time
import json
from collections import defaultdict, OrderedDict
import Pdata
from memory_profiler import memory_usage
def Min_FP(SeqNum, S, sort_item, minsup):
    P = {}
    LP = {}
    FP = []
    FP1=[]
    CanNum = 0
    count = 0
    for i in sort_item:
        CanNum += 1
        for j in range(SeqNum):
            count += len(S[i][j])
        if count >= minsup:
            FP1.append(i)
            FP.append(i)
            compute_utility(i, count)
            P[i] = [[] for k in range(SeqNum)]
            for k in range(SeqNum):
                for pos in S[i][k]:
                    P[i][k].append((pos,))
            LP[i] = [[] for k in range(SeqNum)]
            for k in range(SeqNum):
                for t in range(len(P[i][k])):
                    LP[i][k].append(P[i][k][t][-1] - 1)
        count = 0
    flag = 0
    for fp in FP:
        for item in FP1:
            pattern = fp + " -1 " + item
            CanNum += 1
            LP[pattern] = [[] for k in range(SeqNum)]
            LP[pattern] = [[pos[-1] for pos in P[fp][k]] for k in range(SeqNum)]
            P[pattern] = [[] for i in range(SeqNum)]
            for i in range(SeqNum):
                for j in range(len(P[fp][i])):
                    last_pos = P[fp][i][j][-1]
                    if flag == len(S[item][i]):
                        break
                    for k in range(flag, len(S[item][i])):
                        if S[item][i][k] > last_pos:
                            P[pattern][i].append(P[fp][i][j] + (S[item][i][k],))
                            count += 1
                            flag = k + 1
                            break
                        if k == len(S[item][i]) - 1:
                            flag = len(S[item][i])
                flag = 0
            if count >= minsup:
                FP.append(pattern)
                compute_utility(pattern, count)
            else:
                del P[pattern]
                del LP[pattern]
            count = 0
            if fp.split()[len(fp.split())-1] < item:
                CanNum += 1
                pattern = fp + ' ' + item
                LP[pattern] = [[] for i in range(SeqNum)]
                LP[pattern] = LP[fp]
                P[pattern] = [[] for i in range(SeqNum)]
                unit = S[item][:]
                pattern_array = fp.strip().split(' ')
                f = -1
                for i in range(len(pattern_array)):
                    if pattern_array[i] == '-1':
                        f = i
                for i in range(f+1, len(pattern_array)):
                    for j in range(len(unit)):
                        unit[j] = sorted(list(set(unit[j]) & set(S[pattern_array[i]][j])))
                for i in range(SeqNum):
                    for j in range(len(LP[fp][i])):
                        if flag == len(unit[i]):
                            break
                        for k in range(flag, len(unit[i])):
                            if unit[i][k] > LP[fp][i][j]:
                                P[pattern][i].append(P[fp][i][j][:-1] + (unit[i][k],))
                                count += 1
                                flag = k + 1
                                break
                            if k == len(unit[i]) - 1:
                                flag = len(unit[i])
                    flag = 0
                if count >= minsup:
                    FP.append(pattern)
                    compute_utility(pattern, count)
                else:
                    del P[pattern]
                    del LP[pattern]
            count = 0
        del P[fp]
        del LP[fp]
    print ("Number of frequent patterns: " + str(len(FP)))
    print ("Number of candidate patterns: " + str(CanNum))
    print("Number of AUNP:" + str(len(AUNP)))
Utility = {}
def getUtility(utilityFileName):
    global Utility
    global minsup
    pdata1 = Pdata.processingData()
    maxau = 0
    lines = pdata1.read_file(utilityFileName)
    for line in lines:
        str = line.split(' ')
        Utility[str[0]] = int(str[1])
        if int(str[1]) > maxau:
            maxau = int(str[1])
    minsup = math.ceil(int(minau)/maxau)
def compute_utility(pattern,count):
    global AUNP
    au = 0
    len = 0
    p = pattern.split(" -1 ")
    for i in p:
        j = i.split(" ")
        for K in j:
            len += 1
            au += Utility[str(K)]
    au = au*count/len
    if au>=int(minau):
        AUNP.append(pattern)
if __name__ == '__main__':
    try:
        readFileName = sys.argv[1]
    except Exception as e:
        print(e)
        exit(0)
    pdata = Pdata.processingData()
    minsup = 0
    AUNP = []
    S = {}
    last_index = readFileName.rindex(".")
    dataFileName =readFileName[0:last_index]
    last_index = dataFileName.rindex("/")
    utilityFileName = dataFileName[0:last_index] + "/utility" + dataFileName[last_index:]+"_utility"+".txt"
    SeqNum, S, sort_item = pdata.datap(readFileName, S)
    del pdata
    for minau in sys.argv[2:]:
        getUtility(utilityFileName)
        print('AUNP-B:', readFileName, 'minau=', minau, ':')
        starttime = time.time()
        Min_FP(SeqNum, S, sort_item, int(minsup))
        endtime = time.time()
        print ("Running time: " + str(int(round(endtime * 1000)) - int(round(starttime * 1000))) + "ms")
