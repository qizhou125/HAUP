import argparse
import csv
import gc
import hashlib
import json
import os
import platform
from pathlib import Path
if os.name == "nt":
    import psutil
else:
    import resource
import statistics
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
INF = 10 ** 9
MARK = "@@RESULT@@"


def parse_d(tok):
    tok = tok.strip().lower()
    return INF if tok in ("inf", "infty", "max") else int(tok)


def d_label(d):
    return "inf" if d >= INF else str(d)


def default_utility_path(dataset):

    path = Path(dataset)
    candidates = [path.parent / "utility" / (path.stem + "_utility.txt"),
                  path.with_name(path.stem + "_utility.txt")]
    for candidate in candidates:
        if candidate.is_file():
            return str(candidate)
    raise FileNotFoundError("请使用 --utility 指定效用文件")


def _fingerprint(m):
    h = hashlib.sha256()
    for p in sorted(m.AUNP, key=str):
        key = m._norm_key(p)
        h.update(repr((str(p), int(m.support[str(p)]), float(m.util_map[str(p)]),
                       (m.ItemS_map[key] if m.export_paths else
                        m.path_store[key] if hasattr(m, "path_store") else m.store[key]))).encode())
    return h.hexdigest()


def _maxrss_mb():
    if os.name == "nt":
        return psutil.Process().memory_info().peak_wset / 1048576.0
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / (1024.0 * 1024.0) if sys.platform == "darwin" else r / 1024.0


def worker(mode, cfg):
    from haup_ckpt_miner import load_miner_module, build_miners
    mod = load_miner_module(cfg["miner_path"])
    Baseline, Ckpt = build_miners(mod)
    if cfg["kind"] == "baseline":
        m = Baseline()
    else:
        m = Ckpt().configure(cfg["d"], cfg["encoding"], cfg["cache_values"])
    m.export_paths = cfg.get("workload", "mining") == "export"
    m.minau = cfg["minau"]
    m.SeqNum, m.S, m.sort_item = m.pdata.datap(cfg["dataset"], m.S)
    m.getUtility(cfg["utility"])
    gc.collect()

    mem = {}
    memory_backend = cfg.get("memory_backend", "memory_profiler")
    if mode == "mem" and memory_backend == "tracemalloc":
        import tracemalloc
        tracemalloc.start()

        def cb(_):
            mem["peak_mining_mb"] = tracemalloc.get_traced_memory()[1] / 1048576.0
            tracemalloc.reset_peak()
        m.phase_callback = cb

    if mode == "mem" and memory_backend == "memory_profiler":
        from memory_profiler import memory_usage
        samples = memory_usage((m.Miner, (), {}), interval=cfg.get("memory_interval", 0.1),
                               max_iterations=1, backend="psutil")
    else:
        m.Miner()

    res = {"kind": cfg["kind"], "minsup": m.minsup,
           "aunp": len(m.AUNP), "candidates": m.CanNum, "auub_pruned": m.AuubPruned,
           "mining_ms": m.mining_time_s * 1000.0, "export_ms": m.export_time_s * 1000.0}
    res.update(m.store_stats)
    if mode == "mem" and memory_backend == "memory_profiler":
        res.update(memory_samples_mb=samples, memory_sample_count=len(samples),
                   sampled_peak_rss_mb=max(samples), sampled_min_rss_mb=min(samples),
                   memory_mb=max(samples)-min(samples), memory_backend=memory_backend,
                   memory_interval=cfg.get("memory_interval", 0.1))
    elif mode == "mem":
        cur, peak = tracemalloc.get_traced_memory()
        res["peak_export_mb"] = peak / 1048576.0
        res["peak_mining_mb"] = mem["peak_mining_mb"]
        res["peak_total_mb"] = max(mem["peak_mining_mb"], res["peak_export_mb"])
        res["end_traced_mb"] = cur / 1048576.0
        res["memory_mb"] = res["peak_total_mb"]
        res["memory_backend"] = memory_backend
        tracemalloc.stop()
    else:
        res["peak_rss_mb"] = _maxrss_mb()
    res["total_rebuild_ms"] = (res.get("rebuild_mining_ms", 0) + res.get("rebuild_export_ms", 0))
    res["end_to_end_ms"] = res["mining_ms"] + res["export_ms"]

    res["fingerprint"] = _fingerprint(m)
    print(MARK + json.dumps(res))


def run_worker(mode, cfg, pythonpath, workdir, timeout):
    env = dict(os.environ)
    if pythonpath:
        env["PYTHONPATH"] = os.pathsep.join([pythonpath, env.get("PYTHONPATH", "")])
    p = subprocess.run(
        [sys.executable, os.path.abspath(__file__), "--worker", mode, "--cfg", json.dumps(cfg)],
        cwd=workdir, env=env, capture_output=True, text=True, timeout=timeout)
    if p.returncode:
        raise RuntimeError(p.stderr[-4000:] or p.stdout[-4000:])
    for line in p.stdout.splitlines():
        if line.startswith(MARK):
            return json.loads(line[len(MARK):])
    raise RuntimeError("子进程失败 (mode=%s cfg=%s)\n--- stdout ---\n%s\n--- stderr ---\n%s"
                       % (mode, cfg, p.stdout[-2000:], p.stderr[-4000:]))


TIMING_KEYS = ("mining_ms", "export_ms", "encode_ms", "decode_mining_ms", "decode_export_ms",
               "rebuild_mining_ms", "rebuild_export_ms", "peak_rss_mb",
               "total_rebuild_ms", "end_to_end_ms")


def agg_time(results):
    out = dict(results[0])
    for k in TIMING_KEYS:
        if k in out and out[k] is not None:
            out[k] = statistics.median(r[k] for r in results)
    out["fp_consistent"] = len({r["fingerprint"] for r in results}) == 1
    return out


def make_row(args, minau, name, cfg, t, mm):
    r = {"dataset": os.path.basename(args.dataset), "minau": minau, "config": name,
         "kind": cfg["kind"], "d": d_label(cfg["d"]) if cfg["kind"] != "baseline" else "-",
         "encoding": cfg.get("encoding", "-"), "cache_values": cfg.get("cache_values", "-"),
         "minsup": t["minsup"], "aunp": t["aunp"], "candidates": t["candidates"],
         "auub_pruned": t["auub_pruned"], "F": t["F"], "fingerprint": t["fingerprint"][:16]}
    keys = ["C", "C_over_F", "values_total", "values_ckpt", "payload_bytes",
            "bytes_per_value_ckpt", "bytes_per_value_F", "R", "q", "q_mining", "q_export",
            "cache_hits", "cache_misses", "rebuild_calls_mining", "rebuild_calls_export",
            "rebuild_steps", "decode_calls_mining", "decode_calls_export",
            "encode_ms", "decode_mining_ms", "decode_export_ms",
            "rebuild_mining_ms", "rebuild_export_ms"]
    for k in keys:
        r[k] = t.get(k)
    r["mining_ms"] = t["mining_ms"]
    r["export_ms"] = t["export_ms"]
    r["total_rebuild_ms"] = t["total_rebuild_ms"]
    r["end_to_end_ms"] = t["end_to_end_ms"]
    r["workload"] = args.workload
    r["peak_rss_mb"] = t.get("peak_rss_mb")
    r["peak_mining_mb"] = mm.get("peak_mining_mb")
    r["peak_total_mb"] = mm.get("peak_total_mb")
    for key in ("memory_mb", "memory_backend", "memory_interval", "sampled_peak_rss_mb", "sampled_min_rss_mb", "memory_sample_count"):
        r[key] = mm.get(key)
    return r


def fmt(v, nd=2):
    if v is None:
        return "-"
    if isinstance(v, float):
        return ("%." + str(nd) + "f") % v
    return str(v)


def print_table(rows):
    cols = [("config", "config"), ("F", "|F|"), ("C", "|C|"), ("C_over_F", "|C|/|F|"), ("q", "q"),
            ("bytes_per_value_ckpt", "B/pos(C)"), ("bytes_per_value_F", "B/pos(F)"),
            ("memory_mb", "memoryMiB"), ("peak_rss_mb", "rssMB"), ("mining_ms", "mine ms"),
            ("total_rebuild_ms", "rebuild ms"), ("end_to_end_ms", "e2e ms"), ("equiv", "same")]
    widths = [max(len(h), *(len(fmt(r.get(k), 3 if k in ("C_over_F", "q") else 2)) for r in rows))
              for k, h in cols]
    line = "  ".join(h.rjust(w) for (k, h), w in zip(cols, widths))
    print(line)
    print("-" * len(line))
    for r in rows:
        print("  ".join(fmt(r.get(k), 3 if k in ("C_over_F", "q") else 2).rjust(w)
                        for (k, _), w in zip(cols, widths)))


def plot(rows, outdir, dataset, minau, encoding):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print("（未安装 matplotlib，跳过画图：%s）" % e)
        return
    base = next(r for r in rows if r["kind"] == "baseline")
    ck = [r for r in rows if r["kind"] == "ckpt" and r["encoding"] == encoding]
    if not ck:
        return
    xs = [r["d"] for r in ck]
    fig, ax = plt.subplots(2, 2, figsize=(11, 7.5))
    a = ax[0][0]
    a.plot(xs, [r["memory_mb"] for r in ck], "o-", label="HAUP-ckpt")
    a.axhline(base["memory_mb"], color="gray", ls="--", label="baseline (plain dict)")
    a.set_title("RSS max-min (MiB)" if base["memory_backend"] == "memory_profiler" else "tracemalloc peak (MiB)")
    a.set_xlabel("d")
    a.legend()
    a = ax[0][1]
    a.plot(xs, [r["mining_ms"] for r in ck], "o-", label="mining")
    a.plot(xs, [r["total_rebuild_ms"] for r in ck], "s-", label="rebuild (total)")
    a.plot(xs, [r["end_to_end_ms"] for r in ck], "^-", label="mining+export")
    a.axhline(base["mining_ms"], color="gray", ls="--", label="baseline mining")
    a.set_title("time (ms)")
    a.set_xlabel("d")
    a.legend()
    a = ax[1][0]
    a.plot(xs, [r["C_over_F"] for r in ck], "o-", label="|C|/|F|")
    a.plot(xs, [r["q"] for r in ck], "s-", label="q = |R|/|F|")
    a.set_ylim(0, 1.05)
    a.set_title("checkpoint ratio and rebuild ratio")
    a.set_xlabel("d")
    a.legend()
    a = ax[1][1]
    if encoding == "packed":
        a.plot(xs, [r["bytes_per_value_ckpt"] for r in ck], "o-", label="bytes/position over C")
        a.plot(xs, [r["bytes_per_value_F"] for r in ck], "s-", label="bytes/position over F")
        a.legend()
    a.set_title("compressed bytes per position")
    a.set_xlabel("d")
    fig.suptitle("%s  minau=%s  encoding=%s" % (dataset, minau, encoding))
    fig.tight_layout()
    path = os.path.join(outdir, "tradeoff_%s_minau%s_%s.png" % (
        os.path.splitext(dataset)[0], minau, encoding))
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print("曲线图：", path)


def main():
    if "--worker" in sys.argv:
        wp = argparse.ArgumentParser()
        wp.add_argument("--worker", choices=["time", "mem"])
        wp.add_argument("--cfg")
        wa = wp.parse_args()
        worker(wa.worker, json.loads(wa.cfg))
        return
    ap = argparse.ArgumentParser(description="HAUP 检查点/按需重建 时间-空间权衡实验")
    ap.add_argument("--miner-path", required=True, help="你的 HAUP_Miner3.py 路径")
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--utility", default=None, help="效用文件；默认按 HAUP 原规则推断")
    ap.add_argument("--minau", nargs="+", required=True, help="一个或多个 minau 阈值")
    ap.add_argument("--d-list", default="1,2,3,4,inf", help="检查点间隔，逗号分隔，inf 表示只存单项集模式")
    ap.add_argument("--encodings", default="packed", help="raw 和/或 packed，逗号分隔")
    ap.add_argument("--cache-values", type=int, default=0,
                    help="LRU 缓存上限（以位置值个数计），0=不缓存")
    ap.add_argument("--workload", choices=["mining", "export"], default="mining",
                    help="mining: 仅挖掘；export: 含 AUNP 路径导出")
    ap.add_argument("--memory-backend", choices=["memory_profiler", "tracemalloc"], default="memory_profiler")
    ap.add_argument("--memory-interval", type=float, default=0.1)
    ap.add_argument("--repeat", type=int, default=3, help="计时运行的重复次数（取中位数）")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--pythonpath", default=None,
                    help="附加到子进程 PYTHONPATH 的目录（放 Pdata/memory_profiler 的位置）")
    ap.add_argument("--timeout", type=int, default=3600)
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    if args.repeat < 1 or args.cache_values < 0 or args.memory_interval <= 0:
        ap.error("repeat 必须 >=1，cache-values 必须 >=0")
    args.dataset = os.path.abspath(args.dataset)
    args.miner_path = os.path.abspath(args.miner_path)
    utility = os.path.abspath(args.utility) if args.utility else default_utility_path(args.dataset)
    outdir = os.path.abspath(args.outdir)
    work = os.path.join(outdir, "_work")
    os.makedirs(work, exist_ok=True)
    ds = [parse_d(x) for x in args.d_list.split(",") if x.strip()]
    encs = [e.strip() for e in args.encodings.split(",") if e.strip()]

    if not ds or any(d < 1 for d in ds) or not encs or any(e not in ("raw", "packed") for e in encs):
        ap.error("d 必须 >=1；encodings 只能为 raw,packed")
    if os.path.exists(os.path.join(outdir, "manifest.json")):
        ap.error("输出目录已有实验，请指定新的 --outdir 以保留原始结果")
    sources = [args.dataset, utility, args.miner_path,
               os.path.join(os.path.dirname(args.miner_path), "Pdata.py")]
    sources += [os.path.join(HERE, f) for f in ("run_experiments.py", "haup_ckpt_miner.py", "ckpt_store.py")]
    manifest = {"args": vars(args), "python": sys.version, "executable": sys.executable,
                "platform": platform.platform(), "utility": utility,
                "sha256": {p: hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in sources}}
    Path(outdir, "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    all_rows = []
    for minau in args.minau:
        base_cfg = {"miner_path": args.miner_path, "dataset": args.dataset, "utility": utility,
                    "minau": minau, "workload": args.workload,
                    "memory_backend": args.memory_backend, "memory_interval": args.memory_interval}
        configs = [("baseline", dict(base_cfg, kind="baseline"))]
        for enc in encs:
            for d in ds:
                configs.append(("ckpt-%s-d%s" % (enc, d_label(d)),
                                dict(base_cfg, kind="ckpt", d=d, encoding=enc,
                                     cache_values=args.cache_values)))
        rows = []
        base_row = None
        for name, cfg in configs:
            print("[minau=%s] %s ..." % (minau, name), flush=True)
            runs = [run_worker("time", cfg, args.pythonpath, work, args.timeout)
                    for _ in range(args.repeat)]
            t = agg_time(runs)
            mm = run_worker("mem", cfg, args.pythonpath, work, args.timeout)
            raw_path = Path(outdir, "raw_minau%s_%s.json" % (minau, name))
            raw_path.write_text(json.dumps({"config": cfg, "timing_runs": runs, "memory_run": mm},
                                           ensure_ascii=False, indent=2), encoding="utf-8")
            row = make_row(args, minau, name, cfg, t, mm)
            if base_row is None:
                base_row = row
                base_fingerprint = t["fingerprint"]
            row["equiv"] = bool(t["fp_consistent"] and mm["fingerprint"] == t["fingerprint"]
                                and t["fingerprint"] == base_fingerprint
                                and row["F"] == base_row["F"]
                                and row["candidates"] == base_row["candidates"]
                                and row["aunp"] == base_row["aunp"])
            row["mem_ratio_vs_base"] = (row["memory_mb"] / base_row["memory_mb"]
                                        if base_row["memory_mb"] else None)
            row["mining_time_ratio_vs_base"] = (row["mining_ms"] / base_row["mining_ms"]
                                                if base_row["mining_ms"] else None)
            row["e2e_time_ratio_vs_base"] = (row["end_to_end_ms"] / base_row["end_to_end_ms"]
                                             if base_row["end_to_end_ms"] else None)
            rows.append(row)
        print()
        print("== %s  minau=%s  minsup=%s  |AUNP|=%s ==" % (
            os.path.basename(args.dataset), minau, rows[0]["minsup"], rows[0]["aunp"]))
        print_table(rows)
        print()
        bad = [r["config"] for r in rows if not r["equiv"]]
        if bad:
            print("!! 以下配置与基线结果不一致（请检查）：", bad)
        all_rows.extend(rows)
        if not args.no_plot:
            for enc in encs:
                plot(rows, outdir, os.path.basename(args.dataset), minau, enc)

    csv_path = os.path.join(outdir, "results.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
        w.writeheader()
        w.writerows(all_rows)
    with open(os.path.join(outdir, "results.json"), "w", encoding="utf-8") as f:
        json.dump(all_rows, f, ensure_ascii=False, indent=2)
    print("已写入：", csv_path)
    if any(not r["equiv"] for r in all_rows):
        raise SystemExit("结果不一致，请勿用于论文")


if __name__ == "__main__":
    main()
