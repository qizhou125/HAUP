import copy
import math
import sys
import time
import json
from collections import defaultdict, OrderedDict
import pandas as pd
import Pdata
pdata = Pdata.processingData()
from memory_profiler import memory_usage
CanNum = 0
SItem = {}
SItems = {}
def campute_PreSuf(pattern):
    global SItem, SItems
    if len(pattern[0]) > 1:
        if str(pattern[0][0]) not in SItem:
            SItem[str(pattern[0][0])] = {}
        SItem[str(pattern[0][0])][str(pattern[0][1])] = ''
    else:
        if str(pattern[0][0]) not in SItems:
            SItems[str(pattern[0][0])] = {}
        SItems[str(pattern[0][0])][str(pattern[1][0])] = ''
def Matching_I(list1, list2):
    list3 = [[] for k in range(SeqNum)]
    count = 0
    for i in range(SeqNum):
        list3[i] = sorted(list(set(list1[i]) & set(list2[i])))
        count += len(list3[i])
    return count, list3
def Matching_S(list1, list2):
    list3 = [[] for k in range(SeqNum)]
    flag = 0
    count = 0
    for i in range(SeqNum):
        for j in range(len(list1[i])):
            if flag == len(list2[i]):
                break
            for k in range(flag, len(list2[i])):
                if list2[i][k] > list1[i][j]:
                    list3[i].append(list2[i][k])
                    count += 1
                    flag = k + 1
                    break
                if k == len(list2[i]) - 1:
                    flag = len(list2[i])
        flag = 0
    return count, list3
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
            FP.append([p])
            compute_utility([p], count)
            ItemS[str([p])] = [[] for k in range(SeqNum)]
            ItemS[str([p])] = S[i]
def Miner():
    global CanNum
    FP = []
    ItemS = {}
    Mine_ItemS(FP, ItemS)
    twoLenPattern = two_len(FP, ItemS)
    more_len(FP, ItemS, twoLenPattern)
    print("Number of AUNP:" + str(len(AUNP)))
    print("Number of candidate patterns:" + str(CanNum))
    CanNum = 0
def two_len(FP, ItemS):
    global CanNum
    twoLenPattern = []
    Item = copy.deepcopy(FP)
    for pre in range(len(Item)):
        for suf in range(pre + 1, len(Item)):
            t = copy.deepcopy(Item[pre][0])
            t.append(Item[suf][0][0])
            p = [t]
            CanNum += 1
            count, ItemS[str(p)] = Matching_I(ItemS[str(Item[pre])], ItemS[str(Item[suf])])
            if count >= int(minsup):
                FP.append(p)
                campute_PreSuf(p)
                twoLenPattern.append(p)
                compute_utility(p, count)
            else:
                del ItemS[str(p)]
    for m in Item:
        pre = copy.deepcopy(m)
        for n in Item:
            suf = copy.deepcopy(n)
            p = [pre[0], suf[0]]
            CanNum += 1
            count, ItemS[str(p)] = Matching_S(ItemS[str(pre)], ItemS[str(suf)])
            if count >= int(minsup):
                FP.append(p)
                campute_PreSuf(p)
                twoLenPattern.append(p)
                compute_utility(p, count)
            else:
                del ItemS[str(p)]
    return twoLenPattern
def more_len(FP, ItemS, ExpSet):
    global CanNum
    while ExpSet != []:
        temp = ExpSet[:]
        ExpSet = []
        for m in temp:
            suf = copy.deepcopy(m)
            suf[0].pop(0)
            if suf[0] == []:
                sufstr = str(suf[1:])
            else:
                sufstr = str(suf)
            for n in temp:
                pre = copy.deepcopy(n)
                pre[-1].pop(-1)
                if pre[-1] == []:
                    prestr = str(pre[:-1])
                else:
                    prestr = str(pre)
                if sufstr == prestr:
                    if len(n[-1]) == 1:
                        if str(m[0][0]) in SItems.keys():
                            if str(n[-1][0]) in SItems[str(m[0][0])].keys():
                                pattern = copy.deepcopy(m)
                                pattern.append(n[-1])
                                CanNum += 1
                                count, ItemS[str(pattern)] = Matching_S(ItemS[str(pattern[:-1])],
                                                                        ItemS[str([pattern[-1]])])
                                if count >= int(minsup):
                                    FP.append(pattern)
                                    compute_utility(pattern, count)
                                    ExpSet.append(pattern)
                                else:
                                    del ItemS[str(pattern)]
                    else:
                        pattern = copy.deepcopy(m)
                        pattern[-1].append(n[-1][-1])
                        if len(pattern) > 1:
                            if str(pattern[-1][-1]) in SItems[str(pattern[0][0])].keys():
                                CanNum += 1
                                count, ItemS[str(pattern)] = Matching_S(ItemS[str(pattern[:-1])],
                                                                        ItemS[str([pattern[-1]])])
                                if count >= int(minsup):
                                    FP.append(pattern)
                                    ExpSet.append(pattern)
                                    compute_utility(pattern, count)
                                else:
                                    del ItemS[str(pattern)]
                        else:
                            if str(pattern[-1][-1]) in SItem[str(pattern[0][0])].keys():
                                CanNum += 1
                                count, ItemS[str(pattern)] = Matching_I(ItemS[str(m)], S[str(n[-1][-1])])
                                if count >= int(minsup):
                                    FP.append(pattern)
                                    ExpSet.append(pattern)
                                    compute_utility(pattern, count)
                                else:
                                    del ItemS[str(pattern)]
def count(FP, ItemS):
    SeqSup = {}
    for i in FP:
        SeqSup[str(i)] = []
        for j in ItemS[str(i)]:
            SeqSup[str(i)].append(len(j))
    df = pd.DataFrame(SeqSup)
    strlist = readFileName.split('/')
    while strlist != []:
        str1 = strlist.pop(0)
    str1 = str1.split('.')[0]
    df.to_excel(str1 + 'SeqSup.xlsx', sheet_name='sheet1')
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
        print('AUNP-Miner:', readFileName, 'minau=', minau, ':')
        starttime = time.time()
        Miner()
        endtime = time.time()
        print("Running time: " + str(int(round(endtime * 1000)) - int(round(starttime * 1000))) + "ms")
