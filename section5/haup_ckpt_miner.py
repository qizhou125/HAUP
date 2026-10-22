import importlib.util
import os
import sys
import time

from ckpt_store import CheckpointStore


def load_miner_module(path):

    path = os.path.abspath(path)
    d = os.path.dirname(path)
    if d not in sys.path:
        sys.path.insert(0, d)
    name = os.path.splitext(os.path.basename(path))[0].replace("-", "_")
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    if not hasattr(mod, "PatternMiner"):
        raise ImportError("%s 中没有找到 PatternMiner" % path)
    return mod


def build_miners(miner_mod):
    Base = miner_mod.PatternMiner

    class BaselineMiner(Base):

        phase_callback = None
        export_paths = True

        def Miner(self):
            t0 = time.perf_counter()
            FP = []
            ItemS = {}
            self.Mine_ItemS(FP, ItemS)
            two = self.two_len(FP, ItemS)
            self.more_len(FP, ItemS, two)
            self.mining_time_s = time.perf_counter() - t0
            self.mining_time_ms = int(self.mining_time_s * 1000)
            self.path_store = ItemS
            self.store_stats = {"F": len(ItemS)}
            if self.phase_callback:
                self.phase_callback("mining_done")
            if not self.export_paths:
                self.export_time_s = 0.0
                return
            t1 = time.perf_counter()
            self.ItemS_map = {}
            for p in self.AUNP:
                key = self._norm_key(p)
                self.ItemS_map[key] = ItemS[key]
            self.export_time_s = time.perf_counter() - t1

    class CheckpointMiner(Base):

        phase_callback = None
        export_paths = True
        ckpt_d = 1
        encoding = "packed"
        cache_values = 200000

        def configure(self, d, encoding="packed", cache_values=200000):
            self.ckpt_d = d
            self.encoding = encoding
            self.cache_values = cache_values
            return self

        def Miner(self):
            store = CheckpointStore(self, self.ckpt_d, self.encoding, self.cache_values)
            self.store = store
            t0 = time.perf_counter()
            FP = []
            self.Mine_ItemS(FP, store)
            two = self.two_len(FP, store)
            self.more_len(FP, store, two)
            store.flush()
            self.mining_time_s = time.perf_counter() - t0
            self.mining_time_ms = int(self.mining_time_s * 1000)
            if self.phase_callback:
                self.phase_callback("mining_done")
            if not self.export_paths:
                self.export_time_s = 0.0
                self.store_stats = store.stats()
                return

            t1 = time.perf_counter()
            store.begin_export()
            self.ItemS_map = {}
            for p in self.AUNP:
                key = self._norm_key(p)
                self.ItemS_map[key] = store.fetch_export(key)
            self.export_time_s = time.perf_counter() - t1
            self.store_stats = store.stats()

    return BaselineMiner, CheckpointMiner
