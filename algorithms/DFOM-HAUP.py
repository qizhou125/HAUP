import math
import sys
import time
import json
from collections import defaultdict, OrderedDict
import Pdata
pdata = Pdata.processingData()
from memory_profiler import memory_usage
CanNum = 0
def _last(x):
    return x[-1] if isinstance(x, tuple) else x
def Matching_I(list1, list2):
    list3 = [[] for k in range (SeqNum)]
    count = 0
    for i in range (SeqNum):
        set1 = set(_last(t) for t in list1[i])
        set2 = set(_last(t) for t in list2[i])
        common = sorted(set1 & set2)
        list3[i] = [(pos,) for pos in common]
        count += len(list3[i])
    return count, list3
def DFOM(pattern, ItemS):
    count = 0
    Nettree = [[[] for i in range(SeqNum)] for k in range(len(pattern))]
    unit = [[] for k in range(len(pattern))]
    for i in range(len(pattern)):
        unit[i] = ItemS[str(pattern[i])]
    list3 = [[] for i in range(SeqNum)]
    for i in range(SeqNum):
        for m in range(len(unit[0][i])):
            bbb = 0
            Nettree[0][i].append(unit[0][i][m])
            for j in range(1, len(unit)):
                n = 1
                if unit[j][i] == []:
                    bbb = 1
                    break
                else:
                    for k in range(len(unit[j][i])):
                        if Nettree[j][i] == []:
                            aaa = -1
                        else:
                            aaa = Nettree[j][i][-1]
                        if _last(unit[j][i][k]) > _last(Nettree[j - 1][i][-1]) and _last(unit[j][i][k]) > _last(aaa):
                            Nettree[j][i].append(unit[j][i][k])
                            if j == len(unit) - 1:
                                count += 1
                                full_path = tuple(_last(Nettree[t][i][-1]) for t in range(len(unit)))
                                list3[i].append(full_path)
                            n = 0
                            break
                    if n:
                        bbb = 1
                        break
            if bbb:
                break
    return count, list3
def Join_I(FP, ExpSet, ItemS):
    global CanNum
    while ExpSet != []:
        Temp = []
        for p in ExpSet:
            sufP = p[1:]
            for q in ExpSet:
                preQ = q[:len(q)-1]
                if sufP == preQ:
                    pattern = p[:]
                    pattern.append(q[-1])
                    CanNum += 1
                    qq = []
                    qq.append(q[-1])
                    count, ItemS[str(pattern)] = Matching_I(ItemS[str(p)], ItemS[str(qq)])
                    if count >= int(minsup):
                        compute_utility([pattern], count)
                        FP.append(pattern)
                        Temp.append(pattern)
                    else:
                        del ItemS[str(pattern)]
        ExpSet = Temp[:]
def Gen_I(FP, ItemS):
    global CanNum
    Item = FP[:]
    ExpSet = []
    for i in range (len(Item)):
        for j in range (i+1, len(Item)):
            pre = Item[i]
            suf = Item[j]
            p = []
            p.append(pre[0])
            p.append(suf[0])
            CanNum += 1
            count, ItemS[str(p)] = Matching_I(ItemS[str(pre)], ItemS[str(suf)])
            if count >= int(minsup):
                compute_utility([p], count)
                FP.append(p)
                ExpSet.append(p)
            else:
                del ItemS[str(p)]
    Join_I(FP, ExpSet, ItemS)
def Mine_ItemS(FP, ItemS):
    global CanNum
    for i in sort_item:
        CanNum += 1
        count = 0
        for j in range(SeqNum):
            count += len(S[i][j])
        if count >= int(minsup):
            p = []
            p.append(i)
            compute_utility([p], count)
            FP.append(p)
            ItemS[str(p)] = [[] for k in range(SeqNum)]
            for j in range(SeqNum):
                ItemS[str(p)][j] = [(pos,) for pos in S[i][j]]
    Gen_I(FP, ItemS)
def Join_S(FP, ExpSet, ItemS):
    global CanNum
    while ExpSet != []:
        Temp = []
        for p in ExpSet:
            for q in ExpSet:
                sufP = p[1:]
                preQ = q[:len(q)-1]
                if sufP == preQ:
                    pattern = p[:]
                    pattern.append(q[-1])
                    CanNum += 1
                    count, itemS_val = DFOM(pattern, ItemS)
                    if count >= int(minsup):
                        ItemS[str(pattern)] = itemS_val
                        compute_utility(pattern, count)
                        FP.append(pattern)
                        Temp.append(pattern)
        ExpSet = Temp[:]
def Mine_Pattern(FP, ItemS):
    global CanNum
    ExpSet = []
    temp = FP[:]
    for pre in temp:
        for suf in temp:
            pattern = [pre, suf]
            CanNum += 1
            count, itemS_val = DFOM(pattern, ItemS)
            if count >= int(minsup):
                ItemS[str(pattern)] = itemS_val
                FP.append(pattern)
                compute_utility(pattern, count)
                ExpSet.append(pattern)
    Join_S(FP, ExpSet, ItemS)
def Miner():
    FP = []
    ItemS = {}
    Mine_ItemS(FP, ItemS)
    Mine_Pattern(FP, ItemS)
    print("Number of frequent patterns:" + str(len(FP)))
    print("Number of candidate patterns:" + str(CanNum))
Utility = {}
def getUtility(utilityFileName):
    global Utility
    global minsup
    maxau = 0
    pdata = Pdata.processingData()
    lines = pdata.read_file(utilityFileName)
    for line in lines:
        str = line.split(' ')
        Utility[str[0]] = int(str[1])
        if int(str[1]) > maxau:
            maxau = int(str[1])
    minsup = math.ceil(int(minau) / maxau)
def compute_utility(pattern, count):
    global AUNP
    au = 0
    len = 0
    for i in pattern:
        for j in i:
            len += 1
            au += Utility[str(j)]
    au = au * count / len
    if au >= int(minau):
        AUNP.append(pattern)
if __name__ == '__main__':
    try:
        readFileName = sys.argv[1]
    except Exception as e:
        print(e)
        exit(0)
    minsup = 0
    AUNP = []
    S = {}
    last_index = readFileName.rindex(".")
    dataFileName = readFileName[0:last_index]
    last_index = dataFileName.rindex("/")
    utilityFileName = dataFileName[0:last_index] + "/utility" + dataFileName[last_index:] + "_utility" + ".txt"
    SeqNum, S, sort_item = pdata.datap(readFileName, S)
    del pdata
    for minau in sys.argv[2:]:
        getUtility(utilityFileName)
        print('DFOM-AUNP:', readFileName, 'minau=', minau, ':')
        starttime = time.time()
        Miner()
        endtime = time.time()
        print("Running time: " + str(int(round(endtime * 1000)) - int(round(starttime * 1000))) + "ms")
