import copy
import math
import sys
import time
import json
import bisect
from collections import defaultdict, Counter
import Pdata
from memory_profiler import memory_usage
class PatternMiner:
    def __init__(self):
        self.CanNum = 0
        self.AuubPruned = 0
        self.BfPsubPruned = 0
        self.SItem = {}
        self.SItems = {}
        self.S = {}
        self.sort_item = []
        self.SeqNum = 0
        self.Utility = {}
        self.maxUtil = 0
        self.AUNP = []
        self.minsup = 0
        self.minau = 0
        self.pdata = Pdata.processingData()
        self.ItemS_map = {}
        self.support = {}
        self.support_all = {}
        self.util_map = {}
        self.mine_result_file = None
        self.mining_time_ms = None
    @staticmethod
    def _norm_key(pat):
        if len(pat) == 1 and len(pat[0]) == 1:
            return str([pat[0][0]])
        return str(pat)
    @staticmethod
    def _flatten_pattern(p):
        return [list(itemset) for itemset in p]
    def campute_PreSuf(self, pattern):
        if len(pattern[0]) > 1:
            if str(pattern[0][0]) not in self.SItem:
                self.SItem[str(pattern[0][0])] = {}
            self.SItem[str(pattern[0][0])][str(pattern[0][1])] = ''
        else:
            if str(pattern[0][0]) not in self.SItems:
                self.SItems[str(pattern[0][0])] = {}
            self.SItems[str(pattern[0][0])][str(pattern[1][0])] = ''
    def Matching_I(self, list1, list2):
        list3 = [[] for _ in range(self.SeqNum)]
        count = 0
        for i in range(self.SeqNum):
            b2 = list2[i]
            if not b2:
                continue
            b_starts = {
                (pos_tuple2[0] if isinstance(pos_tuple2, tuple) else pos_tuple2)
                for pos_tuple2 in b2
            }
            seq_out = list3[i]
            for pos_tuple1 in list1[i]:
                if pos_tuple1[-1] in b_starts:
                    seq_out.append(pos_tuple1)
            count += len(seq_out)
        return count, list3
    def Matching_S(self, list1, list2):
        list3 = [[] for _ in range(self.SeqNum)]
        count = 0
        for i in range(self.SeqNum):
            a = list1[i]
            b = list2[i]
            len_a = len(a)
            len_b = len(b)
            seq_out = list3[i]
            if len_a == 0 or len_b == 0:
                continue
            b_is_singleton = (len(b[0]) == 1)
            k = 0
            for j in range(len_a):
                if k == len_b:
                    break
                pos_tuple1 = a[j]
                last_pos1 = pos_tuple1[-1]
                while k < len_b:
                    pos_tuple2 = b[k]
                    pos2 = pos_tuple2[0] if b_is_singleton else pos_tuple2[-1]
                    if pos2 > last_pos1:
                        seq_out.append(pos_tuple1 + (pos2,))
                        count += 1
                        k += 1
                        break
                    k += 1
        return count, list3
    def _support_upper_bound(self, list1, list2):
        return sum(min(len(list1[i]), len(list2[i])) for i in range(self.SeqNum))
    def _pos_upper_bound_S(self, list1, list2):
        total = 0
        for i in range(self.SeqNum):
            a = list1[i]
            b = list2[i]
            if not a or not b:
                continue
            min_a_last = a[0][-1]
            b_vals = [x[0] for x in b]
            idx = bisect.bisect_right(b_vals, min_a_last)
            valid_b_count = len(b) - idx
            total += min(len(a), valid_b_count)
        return total
    def _pass_auub_S(self, list1, list2):
        return self._pos_upper_bound_S(list1, list2) >= int(self.minsup)
    def Mine_ItemS(self, FP, ItemS):
        for i in self.sort_item:
            self.CanNum += 1
            count = 0
            ItemS[str([i])] = [[] for k in range(self.SeqNum)]
            for j in range(self.SeqNum):
                for pos in self.S[i][j]:
                    ItemS[str([i])][j].append((pos,))
                count += len(ItemS[str([i])][j])
            if count >= int(self.minsup):
                p = [i]
                FP.append([p])
                self.support_all[str([i])] = count
                self.compute_utility([p], count)
            else:
                del ItemS[str([i])]
    def two_len(self, FP, ItemS):
        twoLenPattern = []
        Item = copy.deepcopy(FP)
        for pre in range(len(Item)):
            for suf in range(pre + 1, len(Item)):
                t = Item[pre][0] + [Item[suf][0][0]]
                p = [t]
                pre_key = str([Item[pre][0][0]])
                suf_key = str([Item[suf][0][0]])
                self.CanNum += 1
                count, ItemS[str(p)] = self.Matching_I(ItemS[pre_key], ItemS[suf_key])
                if count >= int(self.minsup):
                    FP.append(p)
                    self.campute_PreSuf(p)
                    twoLenPattern.append(p)
                    self.support_all[str(p)] = count
                    self.compute_utility(p, count)
                else:
                    del ItemS[str(p)]
        for m in Item:
            for n in Item:
                p = [m[0], n[0]]
                pre_key = str([m[0][0]])
                suf_key = str([n[0][0]])
                if not self._pass_auub_S(ItemS[pre_key], ItemS[suf_key]):
                    self.AuubPruned += 1
                    self.BfPsubPruned += 1
                    continue
                self.CanNum += 1
                count, ItemS[str(p)] = self.Matching_S(ItemS[pre_key], ItemS[suf_key])
                if count >= int(self.minsup):
                    FP.append(p)
                    self.campute_PreSuf(p)
                    twoLenPattern.append(p)
                    self.support_all[str(p)] = count
                    self.compute_utility(p, count)
                else:
                    del ItemS[str(p)]
        return twoLenPattern
    def more_len(self, FP, ItemS, ExpSet):
        while ExpSet != []:
            temp = ExpSet[:]
            ExpSet = []
            pre_index = defaultdict(list)
            for n in temp:
                last_iset = n[-1]
                if len(last_iset) == 1:
                    prestr = str(n[:-1])
                else:
                    prestr = str(n[:-1] + [last_iset[:-1]])
                pre_index[prestr].append(n)
            for m in temp:
                first_iset = m[0]
                if len(first_iset) == 1:
                    sufstr = str(m[1:])
                else:
                    sufstr = str([first_iset[1:]] + m[1:])
                for n in pre_index.get(sufstr, []):
                    if len(n[-1]) == 1:
                        if str(m[0][0]) in self.SItems.keys():
                            if str(n[-1][0]) in self.SItems[str(m[0][0])].keys():
                                pattern = m + [n[-1]]
                                last_item_key = str([n[-1][0]])
                                prefix_key = self._norm_key(pattern[:-1])
                                if not self._pass_auub_S(ItemS[prefix_key], ItemS[last_item_key]):
                                    self.AuubPruned += 1
                                    self.BfPsubPruned += 1
                                    continue
                                self.CanNum += 1
                                count, ItemS[str(pattern)] = self.Matching_S(
                                    ItemS[prefix_key], ItemS[last_item_key])
                                if count >= int(self.minsup):
                                    FP.append(pattern)
                                    self.support_all[str(pattern)] = count
                                    self.compute_utility(pattern, count)
                                    ExpSet.append(pattern)
                                else:
                                    del ItemS[str(pattern)]
                    else:
                        pattern = m[:-1] + [m[-1] + [n[-1][-1]]]
                        if len(pattern) > 1:
                            if str(pattern[-1][-1]) in self.SItems[str(pattern[0][0])].keys():
                                last_item_key = self._norm_key([pattern[-1]])
                                prefix_key = self._norm_key(pattern[:-1])
                                if not self._pass_auub_S(ItemS[prefix_key], ItemS[last_item_key]):
                                    self.AuubPruned += 1
                                    self.BfPsubPruned += 1
                                    continue
                                self.CanNum += 1
                                count, ItemS[str(pattern)] = self.Matching_S(
                                    ItemS[prefix_key], ItemS[last_item_key])
                                if count >= int(self.minsup):
                                    FP.append(pattern)
                                    ExpSet.append(pattern)
                                    self.support_all[str(pattern)] = count
                                    self.compute_utility(pattern, count)
                                else:
                                    del ItemS[str(pattern)]
                        else:
                            if str(pattern[-1][-1]) in self.SItem[str(pattern[0][0])].keys():
                                last_item_key = str([n[-1][-1]])
                                m_key = self._norm_key(m)
                                self.CanNum += 1
                                count, ItemS[str(pattern)] = self.Matching_I(
                                    ItemS[m_key], ItemS[last_item_key])
                                if count >= int(self.minsup):
                                    FP.append(pattern)
                                    ExpSet.append(pattern)
                                    self.support_all[str(pattern)] = count
                                    self.compute_utility(pattern, count)
                                else:
                                    del ItemS[str(pattern)]
    def getUtility(self, utilityFileName):
        maxau = 0
        lines = self.pdata.read_file(utilityFileName)
        for line in lines:
            parts = line.split(' ')
            if len(parts) >= 2:
                self.Utility[parts[0]] = int(parts[1])
                if int(parts[1]) > maxau:
                    maxau = int(parts[1])
        if maxau > 0:
            self.minsup = math.ceil(int(self.minau) / maxau)
        self.maxUtil = maxau
    def compute_utility(self, pattern, count):
        pattern_str = str(pattern)
        au = 0
        pattern_len = 0
        for itemset in pattern:
            for item in itemset:
                pattern_len += 1
                au += self.Utility.get(str(item), 0)
        if pattern_len > 0:
            au = au * count / pattern_len
            if au >= int(self.minau):
                self.AUNP.append(pattern)
                self.util_map[pattern_str] = au
                self.support[pattern_str] = int(count)
    def Miner(self):
        start = time.time()
        FP = []
        ItemS = {}
        self.Mine_ItemS(FP, ItemS)
        twoLenPattern = self.two_len(FP, ItemS)
        print("Frequent patterns :" + str(FP))
        self.more_len(FP, ItemS, twoLenPattern)
        print("Frequent patterns with size=1:" + str(FP))
        print("Number of frequent patterns with size=1:" + str(len(FP)))
        self.mining_time_ms = int((time.time() - start) * 1000)
        self.ItemS_map = {}
        for key, value in ItemS.items():
            try:
                p_obj = eval(key)
                if isinstance(p_obj, list) and len(p_obj) == 1 and not isinstance(p_obj[0], list):
                    new_key = key
                else:
                    new_key = str(self._flatten_pattern(p_obj))
                self.ItemS_map[new_key] = value
            except Exception as e:
                print(f"Warning: Could not process key {key}: {e}")
                self.ItemS_map[key] = value
    def get_stats(self):
        return {
            "candidate_count": self.CanNum,
            "auub_pruned_count": self.AuubPruned,
            "bfpsub_pruned_count": self.BfPsubPruned,
            "aunp_count": len(self.AUNP),
            "mining_time_ms": self.mining_time_ms,
        }
    def mine_patterns(self, dataset_path, min_utility_threshold, utility_file_path=None, save_mine_path=None):
        flatten_pattern = self._flatten_pattern
        self.__init__()
        self.minau = min_utility_threshold
        self.SeqNum, self.S, self.sort_item = self.pdata.datap(dataset_path, self.S)
        if utility_file_path is None:
            last_index = dataset_path.rindex(".")
            dataFileName = dataset_path[0:last_index]
            last_index = dataFileName.rindex("/")
            utility_file_path = dataFileName[0:last_index] + "/utility" + dataFileName[last_index:] + "_utility.txt"
        self.getUtility(utility_file_path)
        self.Miner()
        if save_mine_path is None:
            base = dataset_path
            try:
                last_index = base.rindex(".")
                base_no_ext = base[:last_index]
            except ValueError:
                base_no_ext = base
            save_mine_path = base_no_ext + "_mine_results.json"
            patterns_info = []
            pattern_ItemS = {}
            for p in self.AUNP:
                flat_p = flatten_pattern(p)
                p_str = str(flat_p)
                info = {
                    "pattern": p_str,
                    "support": int(self.support.get(str(p), 0)),
                    "utility": float(self.util_map.get(str(p), 0.0))
                }
                patterns_info.append(info)
                occ_info = []
                pattern_positions = self.ItemS_map.get(p_str, [])
                for seq_id in range(self.SeqNum):
                    if seq_id < len(pattern_positions) and pattern_positions[seq_id]:
                        occ_info.append({
                            "seq": seq_id,
                            "positions": pattern_positions[seq_id]
                        })
                pattern_ItemS[p_str] = occ_info
            flattened_util_map = {}
            flattened_support_map = {}
            for p_str, val in self.util_map.items():
                try:
                    p_obj = eval(p_str)
                    flat_p = flatten_pattern(p_obj)
                    flattened_util_map[str(flat_p)] = val
                except Exception:
                    flattened_util_map[p_str] = val
            for p_str, val in self.support.items():
                try:
                    p_obj = eval(p_str)
                    flat_p = flatten_pattern(p_obj)
                    flattened_support_map[str(flat_p)] = val
                except Exception:
                    flattened_support_map[p_str] = val
            mine_results = {
                "seq_num": self.SeqNum,
                "S": self.S,
                "sort_item": self.sort_item,
                "Utility": self.Utility,
                "minau": self.minau,
                "minsup": self.minsup,
                "patterns": patterns_info,
                "pattern_ItemS": pattern_ItemS,
                "util_map": flattened_util_map,
                "support_map": flattened_support_map,
                "candidate_count": self.CanNum,
                "auub_pruned_count": self.AuubPruned,
                "bfpsub_pruned_count": self.BfPsubPruned,
                "mining_time_ms": self.mining_time_ms,
            }
            with open(save_mine_path, "w", encoding="utf-8") as f:
                json.dump(mine_results, f, ensure_ascii=False, indent=2)
            self.mine_result_file = save_mine_path
        self.AUNP = [flatten_pattern(p) for p in self.AUNP]
        return {
            'aunp_count': len(self.AUNP),
            'aunp_patterns': self.AUNP,
            'mine_result_file': self.mine_result_file,
            'candidate_count': self.CanNum,
            'auub_pruned_count': self.AuubPruned,
            'bfpsub_pruned_count': self.BfPsubPruned,
            'mining_time_ms': self.mining_time_ms,
        }
    def mine_patterns_from_S(self, S, seq_num, sort_item, utility_dict, min_utility_threshold):
        self.__init__()
        self.S = S
        self.SeqNum = seq_num
        self.sort_item = sort_item
        self.Utility = utility_dict
        self.minau = min_utility_threshold
        maxau = max(utility_dict.values()) if utility_dict else 1
        self.minsup = math.ceil(int(min_utility_threshold) / maxau) if maxau > 0 else 0
        self.maxUtil = maxau
        self.Miner()
        self.after_AUNP = [self._flatten_pattern(p) for p in self.AUNP]
        return {
            'aunp_count': len(self.after_AUNP),
            'aunp_patterns': self.after_AUNP,
            'candidate_count': self.CanNum,
            'auub_pruned_count': self.AuubPruned,
            'bfpsub_pruned_count': self.BfPsubPruned,
            'mining_time_ms': self.mining_time_ms,
        }
def run_pattern_mining(dataset_path, min_utility_threshold, utility_file_path=None):
    miner = PatternMiner()
    return miner.mine_patterns(dataset_path, min_utility_threshold, utility_file_path)
if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python HAUP_v2.py <dataset_path> <min_utility_threshold> [utility_file_path]")
        sys.exit(1)
    dataset_path = sys.argv[1]
    last_index = dataset_path.rindex(".")
    dataFileName = dataset_path[0:last_index]
    last_index = dataFileName.rindex("/")
    utility_file = dataFileName[0:last_index] + "/utility" + dataFileName[last_index:] + "_utility" + ".txt"
    for min_utility in sys.argv[2:]:
        print('HAUP-Miner:', dataset_path, 'minau=', min_utility, ':')
        results = run_pattern_mining(dataset_path, min_utility, utility_file)
        print(f"Number of AUNP: {results['aunp_count']}")
        print(f"AUNP patterns: {results['aunp_patterns']}")
        print(f"Candidates evaluated: {results['candidate_count']}, AUUB pruned (total): {results['auub_pruned_count']} "
          f"Time: {results['mining_time_ms']}ms")
