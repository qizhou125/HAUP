import ast
import struct
import time
from array import array
from collections import OrderedDict
from itertools import accumulate
from operator import sub

_CODE = {1: "B", 2: "H", 4: "I"}
assert array("I").itemsize == 4, "本实验假设 array('I') 为 4 字节"
_HDR = struct.Struct("<BBBI")


def _width(maxv):
    if maxv <= 0xFF:
        return 1
    if maxv <= 0xFFFF:
        return 2
    if maxv <= 0xFFFFFFFF:
        return 4
    raise ValueError("数值超出 32 位，当前编码不支持")


def encode_entry(value, m):

    sid_deltas = []
    counts = []
    vals = []
    prev = 0
    for sid, seq in enumerate(value):
        if not seq:
            continue
        if len(seq[0]) != m:
            raise ValueError("路径长度 %d 与项集数 %d 不一致" % (len(seq[0]), m))
        sid_deltas.append(sid - prev)
        prev = sid
        k = len(seq)
        counts.append(k)
        for col in zip(*seq):
            vals.append(col[0])
            if k > 1:
                diffs = list(map(sub, col[1:], col))
                if min(diffs) <= 0:
                    raise ValueError(
                        "同一列在路径之间不是严格递增，违反引理1，无法做列差分: %r" % (col,))
                vals.extend(diffs)
    n = len(sid_deltas)
    w_sid = _width(max(sid_deltas, default=0))
    w_cnt = _width(max(counts, default=0))
    w_val = _width(max(vals, default=0))
    return b"".join((
        _HDR.pack(w_sid, w_cnt, w_val, n),
        array(_CODE[w_sid], sid_deltas).tobytes(),
        array(_CODE[w_cnt], counts).tobytes(),
        array(_CODE[w_val], vals).tobytes(),
    ))


def decode_entry(blob, m, seq_num):
    w_sid, w_cnt, w_val, n = _HDR.unpack_from(blob, 0)
    mv = memoryview(blob)
    off = _HDR.size
    sids = array(_CODE[w_sid])
    sids.frombytes(mv[off:off + n * w_sid])
    off += n * w_sid
    counts = array(_CODE[w_cnt])
    counts.frombytes(mv[off:off + n * w_cnt])
    off += n * w_cnt
    vals = array(_CODE[w_val])
    vals.frombytes(mv[off:])
    out = [[] for _ in range(seq_num)]
    pos = 0
    for sid, k in zip(accumulate(sids), counts):
        cols = []
        for _ in range(m):
            cols.append(list(accumulate(vals[pos:pos + k])))
            pos += k
        out[sid] = list(zip(*cols))
    return out


def pattern_of_key(key):

    obj = ast.literal_eval(key)
    if obj and not isinstance(obj[0], list):
        return [obj]
    return obj


class CheckpointStore:


    def __init__(self, miner, d, encoding="packed", cache_values=200000):
        if d < 1:
            raise ValueError("d 必须 >= 1（用很大的整数表示 inf）")
        if encoding not in ("raw", "packed"):
            raise ValueError("encoding 只能是 raw 或 packed")
        self.miner = miner
        self.d = d
        self.encoding = encoding
        self.cache_budget = cache_values
        self.phase = "mining"

        self._pending = None
        self.meta = {}
        self.packed = {}
        self.cache = OrderedDict()
        self.cache_size = 0


        self.total_values = 0
        self.ckpt_values = 0
        self.payload_bytes = 0
        self.rebuilt_keys = set()
        self.cnt = {"hits": 0, "misses": 0, "decode_mining": 0, "decode_export": 0,
                    "rebuild_mining": 0, "rebuild_export": 0, "rebuild_steps": 0,
                    "encode": 0}
        self.tm = {"encode": 0.0, "decode_mining": 0.0, "decode_export": 0.0,
                   "rebuild_mining": 0.0, "rebuild_export": 0.0}
        self.rebuilt_by_phase = {"mining": set(), "export": set()}


    def is_checkpoint_m(self, m):
        return m == 1 or m % self.d == 0


    def __setitem__(self, key, value):
        self._flush()
        self._pending = (key, value)

    def __delitem__(self, key):
        if self._pending is not None and self._pending[0] == key:
            self._pending = None
            return
        raise KeyError(key)

    def __contains__(self, key):
        self._flush()
        return key in self.meta

    def __len__(self):
        self._flush()
        return len(self.meta)

    def __getitem__(self, key):
        self._flush()
        if self.encoding == "raw":
            v = self.packed.get(key)
            if v is not None:
                return v
        ent = self.cache.get(key)
        if ent is not None:
            self.cache.move_to_end(key)
            self.cnt["hits"] += 1
            return ent[0]
        self.cnt["misses"] += 1
        info = self.meta.get(key)
        if info is None:
            raise KeyError(key)
        m, values = info
        blob = self.packed.get(key)
        if blob is not None:
            value = self._decode(blob, m)
        else:
            value = self._rebuild(key)
        self._cache_put(key, value, values)
        return value

    def flush(self):
        self._flush()


    def _flush(self):
        if self._pending is None:
            return
        key, value = self._pending
        self._pending = None
        if key in self.meta:
            self._remove(key)
        m = len(pattern_of_key(key))
        values = m * sum(map(len, value))
        self.meta[key] = (m, values)
        self.total_values += values
        if self.is_checkpoint_m(m):
            if self.encoding == "packed":
                t0 = time.perf_counter()
                blob = encode_entry(value, m)
                self.tm["encode"] += time.perf_counter() - t0
                self.cnt["encode"] += 1
                self.packed[key] = blob
                self.payload_bytes += len(blob)
            else:
                self.packed[key] = value
            self.ckpt_values += values
            if self.encoding == "packed":
                self._cache_put(key, value, values)
        else:
            self._cache_put(key, value, values)

    def _remove(self, key):
        m, values = self.meta.pop(key)
        self.total_values -= values
        old = self.packed.pop(key, None)
        if old is not None:
            self.ckpt_values -= values
            if self.encoding == "packed":
                self.payload_bytes -= len(old)
        ent = self.cache.pop(key, None)
        if ent is not None:
            self.cache_size -= ent[1]

    def _cache_put(self, key, value, size):
        if self.cache_budget <= 0:
            return
        old = self.cache.pop(key, None)
        if old is not None:
            self.cache_size -= old[1]
        self.cache[key] = (value, size)
        self.cache_size += size
        while self.cache_size > self.cache_budget and self.cache:
            _, (_, s) = self.cache.popitem(last=False)
            self.cache_size -= s

    def _decode(self, blob, m):
        t0 = time.perf_counter()
        v = decode_entry(blob, m, self.miner.SeqNum)
        self.tm["decode_" + self.phase] += time.perf_counter() - t0
        self.cnt["decode_" + self.phase] += 1
        return v

    def _rebuild(self, key):

        t0 = time.perf_counter()
        pattern = pattern_of_key(key)
        m = len(pattern)
        c = (m // self.d) * self.d
        base_j = c if c >= 1 else 1
        norm = self.miner._norm_key
        cur = self[norm(pattern[:base_j])]
        for j in range(base_j + 1, m + 1):
            t_paths = self[norm([pattern[j - 1]])]
            _, cur = self.miner.Matching_S(cur, t_paths)
            self.cnt["rebuild_steps"] += 1
        self.cnt["rebuild_" + self.phase] += 1
        self.rebuilt_keys.add(key)
        self.rebuilt_by_phase[self.phase].add(key)
        self.tm["rebuild_" + self.phase] += time.perf_counter() - t0
        return cur


    def begin_export(self):

        self._flush()
        self.cache.clear()
        self.cache_size = 0
        self.phase = "export"

    def fetch_export(self, key):

        if self.encoding == "raw":
            v = self.packed.get(key)
            if v is not None:
                return v
        info = self.meta[key]
        blob = self.packed.get(key)
        if blob is not None:
            return self._decode(blob, info[0])
        return self._rebuild(key)


    def stats(self):
        F = len(self.meta)
        C = len(self.packed)
        R = len(self.rebuilt_keys)
        s = {
            "F": F, "C": C,
            "C_over_F": (C / F) if F else 0.0,
            "values_total": self.total_values,
            "values_ckpt": self.ckpt_values,
            "payload_bytes": self.payload_bytes if self.encoding == "packed" else None,
            "bytes_per_value_ckpt": (self.payload_bytes / self.ckpt_values
                                     if self.encoding == "packed" and self.ckpt_values else None),
            "bytes_per_value_F": (self.payload_bytes / self.total_values
                                  if self.encoding == "packed" and self.total_values else None),
            "R": R,
            "q": (R / F) if F else 0.0,
            "q_mining": (len(self.rebuilt_by_phase["mining"]) / F) if F else 0.0,
            "q_export": (len(self.rebuilt_by_phase["export"]) / F) if F else 0.0,
            "cache_hits": self.cnt["hits"], "cache_misses": self.cnt["misses"],
            "rebuild_calls_mining": self.cnt["rebuild_mining"],
            "rebuild_calls_export": self.cnt["rebuild_export"],
            "rebuild_steps": self.cnt["rebuild_steps"],
            "decode_calls_mining": self.cnt["decode_mining"],
            "decode_calls_export": self.cnt["decode_export"],
        }
        for k, v in self.tm.items():
            s[k + "_ms"] = v * 1000.0
        return s
