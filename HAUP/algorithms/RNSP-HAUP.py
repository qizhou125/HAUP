import math
import sys
import time
from bisect import bisect_right
from typing import Dict, List, Set, Tuple, Iterator
import os
import gc
import psutil
import threading

# ============================================================
# Data Preprocessing (inlined from Pdata.py, minimal & aligned)
#   Output structure (RNP-Miner style):
#     - SeqNum: int (number of sequences)
#     - S: Dict[item, List[List[int]]]   item -> positions per sequence
#     - sort_item: List[item]            global item vocabulary, sorted
#     - itemcount: int                   total number of itemsets in DB
# ============================================================

class DataProcessor:
    """Pdata-style preprocessing using *position lists* (no bitmaps).

    Parsing rules (as requested):
      - Sequence boundary: newline OR token '-2'
      - Itemset boundary: token '-1'
    Outputs:
      - SeqNum: number of sequences
      - S: item -> positions per sequence (OccType = List[List[int]])
      - sort_item: global items, sorted
      - itemcount: total number of itemsets in DB
    """

    @staticmethod
    def _parse_sequences(readfilename: str) -> List[List[List[str]]]:
        """Return sequences as list of itemsets, each itemset is list of items (strings)."""
        sequences: List[List[List[str]]] = []
        cur_seq: List[List[str]] = []
        cur_iset: List[str] = []

        def close_itemset():
            nonlocal cur_iset, cur_seq
            # In standard SPM format, every itemset ends with -1.
            # We still close on newline/end for robustness.
            cur_seq.append(cur_iset)
            cur_iset = []

        def close_sequence():
            nonlocal cur_seq, cur_iset, sequences
            if cur_iset or (cur_seq and cur_seq[-1] is not cur_iset):
                # If there's an unclosed itemset (no trailing -1), close it.
                if cur_iset:
                    close_itemset()
            if cur_seq:
                sequences.append(cur_seq)
            cur_seq = []
            cur_iset = []

        with open(readfilename, "r", encoding="utf-8") as f:
            for raw_line in f:
                line = raw_line.strip()
                if not line:
                    # treat empty line as a sequence boundary
                    close_sequence()
                    continue

                tokens = [t for t in line.replace("  ", " ").split(" ") if t != ""]
                ended_by_minus2 = False

                for tok in tokens:
                    if tok == "-1":
                        close_itemset()
                    elif tok == "-2":
                        ended_by_minus2 = True
                        close_sequence()
                    else:
                        cur_iset.append(tok)

                # newline is also a sequence boundary (unless we already closed by -2)
                if not ended_by_minus2:
                    close_sequence()

        # file end
        close_sequence()
        return sequences

    @staticmethod
    def _item_sorted(items: List[str]) -> List[str]:
        """Dedupe by first appearance, then global sort."""
        seen = set()
        no_repeat: List[str] = []
        for it in items:
            if it not in seen:
                seen.add(it)
                no_repeat.append(it)
        return sorted(no_repeat)

    def datap(self, readFileName: str) -> Tuple[int, Dict[str, List[List[int]]], List[str], int]:
        seq_itemsets = self._parse_sequences(readFileName)
        seqnum = len(seq_itemsets)

        # collect items in first-appearance order (Pdata semantics)
        items: List[str] = []
        seen_global: Set[str] = set()
        for seq in seq_itemsets:
            for iset in seq:
                for it in iset:
                    if it not in seen_global:
                        seen_global.add(it)
                        items.append(it)

        sort_item = self._item_sorted(items)

        # init S[item][sid] = []
        S: Dict[str, List[List[int]]] = {it: [[] for _ in range(seqnum)] for it in sort_item}

        # fill positions
        itemcount = 0
        for sid, seq in enumerate(seq_itemsets):
            itemcount += len(seq)
            for k, iset in enumerate(seq):
                # match Pdata's `if item in itemset`: duplicates inside one itemset count once
                for it in set(iset):
                    if it in S:
                        S[it][sid].append(k)

        return seqnum, S, sort_item, itemcount


# ============================================================
# Pattern + occurrence types (position lists, NOT bitmaps)
# ============================================================

PatternType = Tuple[Tuple[str, ...], ...]               # e.g., ((a,), (b,c)) represents (a)(bc)
OccType = List[List[int]]                              # positions per sequence


def support_of_occ(occ: OccType) -> int:
    return sum(len(lst) for lst in occ)


def _intersect_two_sorted(a: List[int], b: List[int]) -> List[int]:
    """Intersection of two increasing integer lists (unique positions)."""
    i = j = 0
    out: List[int] = []
    la, lb = len(a), len(b)
    while i < la and j < lb:
        va, vb = a[i], b[j]
        if va == vb:
            out.append(va)
            i += 1
            j += 1
        elif va < vb:
            i += 1
        else:
            j += 1
    return out


def matching_I(list1: OccType, list2: OccType) -> OccType:
    """RNP Matching_I: per-sequence intersection."""
    return [_intersect_two_sorted(list1[sid], list2[sid]) for sid in range(len(list1))]


def matching_S(list1: OccType, list2: OccType) -> OccType:
    """RNP Matching_S: greedy non-overlapping sequential match (positions strictly increasing)."""
    seqnum = len(list1)
    out: OccType = [[] for _ in range(seqnum)]
    for sid in range(seqnum):
        A = list1[sid]
        B = list2[sid]
        if not A or not B:
            continue
        flag = 0
        for a_pos in A:
            # advance in B until find first b_pos > a_pos
            while flag < len(B) and B[flag] <= a_pos:
                flag += 1
            if flag >= len(B):
                break
            out[sid].append(B[flag])
            flag += 1
    return out


# ============================================================
# Utility model (added for AUNP comparison against HUP-Miner)
#   - Utility: item -> utility value, read from a utility file
#     with lines "item utility_value" (same file format/role as
#     HUP-Miner's getUtility)
#   - getUtility(): fills the global Utility dict and returns
#     maxau, the maximum single-item utility value found; this is
#     used to convert minau -> minsup via
#     minsup = ceil(minau / maxau), exactly as in HUP-Miner
#   - compute_utility(): average-utility computation, ported from
#     HUP-Miner's compute_utility and adapted only to RNSP's
#     PatternType (tuple-of-tuples instead of list-of-lists). It is
#     a pure post-hoc filter applied to patterns RNSP already found
#     frequent — it introduces no new candidate-pruning strategy.
# ============================================================

Utility: Dict[str, int] = {}


def getUtility(utility_file_name: str) -> int:
    """Read a utility file (format: 'item utility_value' per line).

    Fills the module-level Utility dict and returns maxau, the
    maximum utility value found across all items.
    """
    global Utility
    Utility = {}
    maxau = 0
    with open(utility_file_name, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(' ')
            item, value = parts[0], int(parts[1])
            Utility[item] = value
            if value > maxau:
                maxau = value
    return maxau


# ============================================================
# CSC-NIPFs (static) — now fully position-structure based
#   - removed all bitmap logic
#   - support counting follows RNP-Miner <-> Pdata.py interface
# ============================================================

class NIPFStaticMiner:

    def __init__(self, min_support_rate: float = 0.4, minau=None, min_sup_count: int = None):
        self.lambda_rate = min_support_rate
        self.generated_candidates_count = 0

        # --- utility-mining additions ---
        # minau: the average-utility threshold (mirrors HUP-Miner's minau)
        # min_sup_count: an absolute support count for B_min, computed
        #   upstream as minsup = ceil(minau / maxau); when given it
        #   overrides the lambda_rate-based B_min below, so the two
        #   algorithms are compared under the same minsup.
        self.minau = minau
        self.min_sup_count = min_sup_count
        self.AUNP: List[Tuple[PatternType, int]] = []

        self.SeqNum: int = 0
        self.S: Dict[str, OccType] = {}
        self.B_min: int = 0

        self.freq_patterns: Dict[PatternType, int] = {}

        # caches
        self.pat_occ: Dict[PatternType, OccType] = {}
        self.itemset_occ: Dict[Tuple[str, ...], OccType] = {}

        # pruning structures
        self.F_t: Set[PatternType] = set()
        self.I_adj: Dict[str, Set[str]] = {}
        self.S_adj: Dict[str, Set[str]] = {}

    # ---------- utility (added; no new pruning strategy) ----------
    def compute_utility(self, pattern: "PatternType", count: int) -> None:
        """Average-utility filter, ported from HUP-Miner's compute_utility.

        au = (sum of each item's utility across the pattern) * count / (number of items in pattern)

        A pattern already found frequent by RNSP (sup >= B_min) is
        additionally kept in self.AUNP when au >= minau. This does
        not change which patterns RNSP explores or how it prunes
        candidates — it only labels/filters the frequent patterns
        RNSP already produces, same as HUP-Miner does for its own
        frequent patterns.
        """
        if self.minau is None:
            return
        total_utility = 0
        length = 0
        for itemset in pattern:
            for item in itemset:
                length += 1
                total_utility += Utility[str(item)]
        au = total_utility * count / length
        if au >= int(self.minau):
            self.AUNP.append((pattern, count))

    # ---------- public ----------
    def load_and_mine(self, filename: str) -> Dict[PatternType, int]:
        t0 = time.time()
        dp = DataProcessor()
        self.SeqNum, self.S, items, total_itemsets = dp.datap(filename)

        # if an absolute minsup (derived from minau) was supplied, use it
        # directly instead of the lambda_rate-based density threshold,
        # so RNSP and HUP-Miner run under the same minsup.
        if self.min_sup_count is not None:
            self.B_min = self.min_sup_count
        else:
            self.B_min = total_itemsets * self.lambda_rate

        self.AUNP = []
        self._mine(items)

        t1 = time.time()
        print(f"--- Loading: {filename} ---")
        print(f"minsup: {self.B_min}")
        print(f"Frequent Patterns: {len(self.freq_patterns)}")
        print(f"Candidates Generated: {self.generated_candidates_count}")
        if self.minau is not None:
            print(f"minau: {self.minau}")
            print(f"Number of AUNP: {len(self.AUNP)}")
        print(f"Total Time: {t1 - t0:.2f}s")
        return self.freq_patterns

    # ---------- mining ----------
    def _mine(self, items: List[str]) -> None:
        self.freq_patterns.clear()
        self.pat_occ.clear()
        self.itemset_occ.clear()
        self.I_adj.clear()
        self.S_adj.clear()
        self.F_t.clear()
        self.generated_candidates_count = 0

        # item support
        item_sup = {it: support_of_occ(self.S[it]) for it in items}

        # keep only frequent items
        items = [it for it in items if item_sup[it] >= self.B_min]
        for it in list(self.S.keys()):
            if it not in item_sup or item_sup[it] < self.B_min:
                del self.S[it]

        # --- level 1 ---
        level: List[PatternType] = []
        level_sup: Dict[PatternType, int] = {}

        for it in items:
            pat = ((it,),)
            sup = item_sup[it]
            level.append(pat)
            level_sup[pat] = sup
            self.freq_patterns[pat] = sup
            self.pat_occ[pat] = self.S[it]
            self.compute_utility(pat, sup)

        # --- level k (k >= 2) ---
        k = 2
        while level:
            self.F_t = set(level)
            if k == 3:
                self._build_f2_adjacency()

            next_level: List[PatternType] = []
            next_sup: Dict[PatternType, int] = {}

            for cand in self._candidates(level, k, level_sup, item_sup, items):
                self.generated_candidates_count += 1
                sup, occ = self._support(cand)

                if sup >= self.B_min:
                    next_level.append(cand)
                    next_sup[cand] = sup
                    self.freq_patterns[cand] = sup
                    self.pat_occ[cand] = occ
                    self.compute_utility(cand, sup)

            self._shrink_caches(next_level)
            level, level_sup = next_level, next_sup
            k += 1

    # ---------- support ----------
    def _support(self, pat: PatternType) -> Tuple[int, OccType]:
        # (1) one itemset pattern: (a) or (ab)
        if len(pat) == 1:
            occ = self._itemset_occ(pat[0])
            return support_of_occ(occ), occ

        # (2) sequential pattern: prefix + last itemset
        prefix = pat[:-1]
        last = pat[-1]

        prefix_occ = self.pat_occ.get(prefix)
        if prefix_occ is None:
            # only allow fallback for length-1 prefix
            if len(prefix) == 1:
                prefix_occ = self._itemset_occ(prefix[0])
            else:
                return 0, [[] for _ in range(self.SeqNum)]

        tail_occ = self._itemset_occ(last)
        occ = matching_S(prefix_occ, tail_occ)
        return support_of_occ(occ), occ

    def _itemset_occ(self, itset: Tuple[str, ...]) -> OccType:
        """Get occurrence lists for an itemset; cache intersections for multi-item itemsets."""
        if len(itset) == 1:
            return self.S[itset[0]]

        cached = self.itemset_occ.get(itset)
        if cached is not None:
            return cached

        occ = self.S[itset[0]]
        for it in itset[1:]:
            occ = matching_I(occ, self.S[it])
        self.itemset_occ[itset] = occ
        return occ

    # ---------- candidate generation ----------
    def _candidates(
        self,
        prev_level: List[PatternType],
        k: int,
        prev_sup: Dict[PatternType, int],
        item_sup: Dict[str, int],
        items_sorted: List[str],
    ) -> Iterator[PatternType]:
        for p in prev_level:
            last_itemset = p[-1]
            last_item = last_itemset[-1]
            p_sup = prev_sup[p]
            anchor = p[0][0]

            # I-extension: extend last itemset with x > last_item
            start = bisect_right(items_sorted, last_item)
            allow_I = self._allow_I(anchor, last_itemset, p) if (k > 2 and self.I_adj) else None

            for x in items_sorted[start:]:
                if allow_I is not None and x not in allow_I:
                    continue
                # SHUP is disabled in this ablation file (kept as-is)
                cand = p[:-1] + (last_itemset + (x,),)
                if k <= 2 or self._suffix_in_Ft(cand):    # FS
                    yield cand

            # S-extension: append new singleton itemset (y)
            allow_S = self._allow_S(anchor, last_itemset) if (k > 2 and self.S_adj) else None
            y_iter = allow_S if allow_S is not None else items_sorted

            for y in y_iter:
                cand = p + ((y,),)
                if k <= 2 or self._suffix_in_Ft(cand):     # FS
                    yield cand

    def _suffix_in_Ft(self, pat: PatternType) -> bool:
        first = pat[0]
        suffix = ((first[1:],) + pat[1:]) if len(first) > 1 else pat[1:]
        return suffix in self.F_t

    def _build_f2_adjacency(self) -> None:
        self.I_adj.clear()
        self.S_adj.clear()

        for p in self.F_t:
            if len(p) == 1 and len(p[0]) == 2:  # (ab)
                a, b = p[0]
                self.I_adj.setdefault(a, set()).add(b)
            elif len(p) == 2 and len(p[0]) == 1 and len(p[1]) == 1:  # (a)(b)
                a, b = p[0][0], p[1][0]
                self.S_adj.setdefault(a, set()).add(b)

    def _allow_I(self, anchor: str, last_itemset: Tuple[str, ...], pat: PatternType) -> Set[str]:
        anchor_map = self.I_adj if len(pat) == 1 else self.S_adj
        pbm = anchor_map.get(anchor)
        if not pbm:
            return set()

        tail = self._intersect_adjs(self.I_adj, last_itemset)
        return pbm & tail if tail else set()

    def _allow_S(self, anchor: str, last_itemset: Tuple[str, ...]) -> Set[str]:
        pbm = self.S_adj.get(anchor)
        if not pbm:
            return set()

        tail = self._intersect_adjs(self.S_adj, last_itemset)
        return pbm & tail if tail else set()

    @staticmethod
    def _intersect_adjs(adj: Dict[str, Set[str]], keys: Tuple[str, ...]) -> Set[str]:
        if not keys:
            return set()
        it = iter(keys)
        first = next(it, None)
        if first is None:
            return set()
        s = set(adj.get(first, ()))
        if not s:
            return set()
        for k in it:
            s &= adj.get(k, set())
            if not s:
                return set()
        return s

    # ---------- cache shrink ----------
    def _shrink_caches(self, next_level: List[PatternType]) -> None:
        if not next_level:
            self.pat_occ.clear()
            self.itemset_occ.clear()
            return

        needed = set(next_level)
        needed.update(p[:-1] for p in next_level if len(p) > 1)
        self.pat_occ = {p: self.pat_occ[p] for p in needed if p in self.pat_occ}

        needed_itemsets = {p[-1] for p in next_level if len(p[-1]) > 1}
        self.itemset_occ = {s: self.itemset_occ[s] for s in needed_itemsets if s in self.itemset_occ}


# =========================
# Peak RSS helper (unchanged)
# =========================

_MB = 1024 * 1024

def _rss_bytes() -> int:
    return psutil.Process(os.getpid()).memory_info().rss

def run_with_peak_rss(func, *, interval: float = 0.02):
    baseline = _rss_bytes()
    peak = baseline

    stop = threading.Event()

    def sampler():
        nonlocal peak
        while not stop.is_set():
            rss = _rss_bytes()
            if rss > peak:
                peak = rss
            stop.wait(interval)

    t = threading.Thread(target=sampler, daemon=True)
    t.start()
    try:
        result = func()
    finally:
        stop.set()
        t.join()

    end = _rss_bytes()
    return result, baseline, peak, end


if __name__ == "__main__":
    print("=== RNSP-Miner (AUNP comparison mode) ===")
    try:
        readFileName = sys.argv[1]
        # minsup = int(sys.argv[2])
    except Exception as e:
        print(e)
        exit(0)

    #f = "../data/1.txt"
    last_index = readFileName.rindex(".")
    dataFileName = readFileName[0:last_index]
    last_index = dataFileName.rindex("/")
    utilityFileName = dataFileName[0:last_index] + "/utility" + dataFileName[last_index:] + "_utility" + ".txt"

    # utility file: one line per item, "item utility_value"
    # (same format/role as HUP-Miner's utility file, e.g. SDB9_utility.txt)
    #utility_file = "../data/1_utility.txt"

    # min_support_rate is kept only as a fallback for plain RNSP runs
    # (used when no minau is supplied); it is NOT used in AUNP mode.
    min_support_rate = 0.25



    repeats = 1
    sample_interval = 0.05

    def run_once(minau):
        # convert minau -> minsup, identical to HUP-Miner:
        #   minsup = ceil(minau / maxau)
        maxau = getUtility(utilityFileName)
        minsup = math.ceil(int(minau) / maxau)

        gc.collect()
        miner = NIPFStaticMiner(min_support_rate=min_support_rate,
                                 minau=minau, min_sup_count=minsup)
        t0 = time.perf_counter()
        _, base, peak, _ = run_with_peak_rss(lambda: miner.load_and_mine(readFileName), interval=sample_interval)
        t1 = time.perf_counter()
        return (t1 - t0), base, peak, miner

    for minau in sys.argv[2:]:
        print(f"\nAUNP-Miner (RNSP-based): {readFileName}, minau={minau}:")
        results = [run_once(minau) for _ in range(repeats)]
        avg_time = sum(t for t, _, _, _ in results) / repeats
        avg_peak_mb = sum(peak for _, _, peak, _ in results) / repeats / _MB

        print(f"Runs={repeats} | Time(avg)={avg_time:.4f}s | PeakRSS(avg)={avg_peak_mb:.2f}MB")

