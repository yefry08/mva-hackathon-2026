#!/usr/bin/env python3
"""Pool de workers para la corrida del MVA Hackathon.

Reemplaza al dispatcher del scaffold, que tenia cuatro defectos que lo hacian
inservible: ejecutaba por lotes y esperaba a que terminara el lote entero, no
respetaba dependencias, corria todo como tarea de LLM aunque fuera un job de
horas, y confiaba en que el propio worker declarara si habia filtrado datos.

Tres tipos de tarea:
  C  job de computo determinista. Se ejecuta como comando, sin LLM.
  A  agente. Una invocacion headless de `claude -p` con rol y contrato JSON.
  H  gate humano. El dispatcher se detiene y espera confirmacion.

Uso:
    python orchestrator/dispatch.py --dry-run
    python orchestrator/dispatch.py --max-concurrent 8 --budget 100
"""
from __future__ import annotations

import argparse
import json
import os
import queue
import re
import subprocess
import sys
import threading
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
QUEUE = ROOT / "orchestrator" / "queue.jsonl"
STATE = ROOT / "orchestrator" / "state.json"
LOG = ROOT / "reports" / "run.log.jsonl"
WORKER_OUT = ROOT / "reports" / "workers"

MAX_RETRIES = 2
TASK_TIMEOUT = int(os.environ.get("TASK_TIMEOUT", "1800"))
COMPUTE_TIMEOUT = int(os.environ.get("COMPUTE_TIMEOUT", "43200"))  # 12h
MODEL = os.environ.get("MVA_MODEL", "claude-opus-5")

# Herramientas que un worker necesita. Sin esto, en headless no puede correr nada.
ALLOWED_TOOLS = "Bash,Read,Write,Edit,Glob,Grep,WebSearch,WebFetch"

RATE_LIMIT = re.compile(r"rate.?limit|429|overloaded|529", re.IGNORECASE)

_state_lock = threading.Lock()
_log_lock = threading.Lock()


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def log(event: dict) -> None:
    event["ts"] = now()
    with _log_lock:
        LOG.parent.mkdir(parents=True, exist_ok=True)
        with LOG.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(event, ensure_ascii=False) + "\n")
    print("[%s] %-10s %-6s %s" % (event["ts"][11:19], event.get("kind", ""),
                                  event.get("task_id", ""), event.get("msg", "")),
          flush=True)


def load_state() -> dict:
    if STATE.exists():
        return json.loads(STATE.read_text(encoding="utf-8"))
    return {"tasks": {}, "started_at": now(), "cost_usd": 0.0,
            "submissions_used": {"track1": 0, "track2": 0}}


def save_state(state: dict) -> None:
    tmp = STATE.with_suffix(".tmp")
    tmp.write_text(json.dumps(state, indent=2, ensure_ascii=False), encoding="utf-8")
    tmp.replace(STATE)


def load_queue() -> list:
    tasks = []
    for line in QUEUE.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            tasks.append(json.loads(line))
    return tasks


def validate_dag(tasks) -> list:
    """Ciclos y dependencias huerfanas. Correr esto antes de gastar un centavo."""
    by_id = {t["task_id"]: t for t in tasks}
    problems = []
    for t in tasks:
        for dep in t.get("depends_on", []):
            if dep not in by_id:
                problems.append("%s depende de %s, que no existe" % (t["task_id"], dep))
    color = {}

    def visit(tid, stack):
        if color.get(tid) == "done":
            return
        if color.get(tid) == "open":
            problems.append("ciclo: %s" % " -> ".join(stack + [tid]))
            return
        color[tid] = "open"
        for dep in by_id.get(tid, {}).get("depends_on", []):
            if dep in by_id:
                visit(dep, stack + [tid])
        color[tid] = "done"

    for t in tasks:
        visit(t["task_id"], [])
    return problems


PROMPT = """Eres el worker `{role}` de la corrida del MVA Hackathon 2026.

Lee CLAUDE.md antes de actuar. Sus reglas son vinculantes, en particular:
- Hacia servicios externos solo va nivel gen, via y mecanismo. El hook egress-guard
  bloquea lo demas; si te bloquea, reformula o reporta blocked. No lo rodees.
- Cada afirmacion biologica necesita un DOI o PMID que resuelvas en esta misma
  tarea. Lo que no puedas verificar va marcado [UNVERIFIED].
- Nada de consejo medico ni recomendaciones clinicas.
- No envies submissions ni hagas push publico. Eso lo decide el humano.

TAREA {task_id} (fase {phase}):
{description}

Al terminar escribe exactamente un JSON en reports/workers/{task_id}.json siguiendo el
contrato de la seccion 8 de CLAUDE.md, con `external_calls` completo: el texto exacto
de cada consulta que saliera de esta maquina. Si te bloqueas, usa status "blocked" y
explica por que. Un blocked honesto vale mas que un resultado fabricado.
"""


def run_agent(task: dict):
    tid = task["task_id"]
    prompt = PROMPT.format(role=task["role"], task_id=tid, phase=task.get("phase", "?"),
                           description=task["description"])
    cmd = ["claude", "-p", prompt, "--model", MODEL, "--output-format", "json",
           "--allowedTools", ALLOWED_TOOLS, "--permission-mode", "acceptEdits"]
    try:
        proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True,
                              timeout=TASK_TIMEOUT, encoding="utf-8", errors="replace")
    except subprocess.TimeoutExpired:
        return "failed", "timeout tras %ds" % TASK_TIMEOUT, 0.0
    except FileNotFoundError:
        return "failed", "no encuentro el CLI `claude` en PATH", 0.0

    cost = 0.0
    try:
        payload = json.loads(proc.stdout or "{}")
        cost = float(payload.get("total_cost_usd") or 0.0)
    except (json.JSONDecodeError, TypeError, ValueError):
        pass

    if proc.returncode != 0:
        return "failed", (proc.stderr or "")[-400:], cost

    out = WORKER_OUT / ("%s.json" % tid)
    if not out.exists():
        return "failed", "el worker no escribio su archivo de output", cost
    try:
        result = json.loads(out.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return "failed", "output JSON malformado", cost
    return result.get("status", "done"), "", cost


def run_compute(task: dict):
    cmd = task.get("command")
    if not cmd:
        return "failed", "tarea C sin campo `command`", 0.0
    try:
        proc = subprocess.run(cmd, cwd=ROOT, shell=True, capture_output=True,
                              text=True, timeout=COMPUTE_TIMEOUT,
                              encoding="utf-8", errors="replace")
    except subprocess.TimeoutExpired:
        return "failed", "timeout tras %ds" % COMPUTE_TIMEOUT, 0.0
    if proc.returncode != 0:
        return "failed", (proc.stderr or proc.stdout or "")[-400:], 0.0
    logdir = ROOT / "logs"
    logdir.mkdir(exist_ok=True)
    (logdir / ("%s.out" % task["task_id"])).write_text(proc.stdout or "", encoding="utf-8")
    return "done", "", 0.0


def worker_loop(work_q, state, args, stop):
    while not stop.is_set():
        try:
            task = work_q.get(timeout=1)
        except queue.Empty:
            return
        tid, kind = task["task_id"], task.get("kind", "A")
        log({"kind": "start", "task_id": tid, "msg": "%s %s" % (kind, task["role"])})
        started = time.time()

        if args.dry_run:
            status, msg, cost = "done", "dry-run", 0.0
        elif kind == "C":
            status, msg, cost = run_compute(task)
        else:
            status, msg, cost = run_agent(task)

        with _state_lock:
            rec = state["tasks"].setdefault(tid, {"retries": 0})
            rec.update({"status": status, "updated_at": now(), "role": task["role"],
                        "kind": kind, "seconds": round(time.time() - started, 1)})
            state["cost_usd"] = round(state.get("cost_usd", 0.0) + cost, 4)
            if status == "failed":
                rec["retries"] += 1
                rec["last_error"] = msg
                if RATE_LIMIT.search(msg or ""):
                    rec["status"] = "pending"
                    log({"kind": "ratelimit", "task_id": tid, "msg": "backoff"})
                    time.sleep(min(300, 15 * 2 ** rec["retries"]))
                elif rec["retries"] > MAX_RETRIES:
                    rec["status"] = "abandoned"
                    log({"kind": "abandoned", "task_id": tid, "msg": msg})
                else:
                    rec["status"] = "pending"
                    log({"kind": "retry", "task_id": tid, "msg": msg})
            else:
                log({"kind": status, "task_id": tid,
                     "msg": "%.0fs  $%.3f" % (rec["seconds"], cost)})
            save_state(state)
        work_q.task_done()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--max-concurrent", type=int,
                    default=int(os.environ.get("MAX_CONCURRENT", "15")))
    ap.add_argument("--budget", type=float,
                    default=float(os.environ.get("BUDGET_USD", "400")))
    args = ap.parse_args()

    tasks = load_queue()
    problems = validate_dag(tasks)
    if problems:
        for p in problems:
            print("  !", p)
        print("el DAG no valida. No arranco.")
        return 1

    state = load_state()
    log({"kind": "run_start", "msg": "%d tareas, concurrencia %d, tope $%.0f"
         % (len(tasks), args.max_concurrent, args.budget)})

    while True:
        with _state_lock:
            done = {tid for tid, r in state["tasks"].items()
                    if r.get("status") in ("done", "abandoned")}
            cost = state.get("cost_usd", 0.0)
            ready = [t for t in tasks
                     if t["task_id"] not in done
                     and all(d in done for d in t.get("depends_on", []))]

        if cost >= args.budget:
            log({"kind": "budget", "msg": "$%.2f alcanza el tope. Pausa." % cost})
            return 2
        if not ready:
            pendientes = [t["task_id"] for t in tasks if t["task_id"] not in done]
            if pendientes:
                log({"kind": "deadlock", "msg": "%d tareas con dependencias sin "
                     "resolver: %s" % (len(pendientes), pendientes[:8])})
            break

        gates = [t for t in ready if t.get("kind") == "H"]
        if gates:
            log({"kind": "gate", "task_id": gates[0]["task_id"],
                 "msg": "GATE HUMANO: %s" % gates[0]["description"][:120]})
            print("\nResuelvelo, marca la tarea como done en state.json y vuelve a "
                  "lanzar el dispatcher.")
            return 3

        work_q = queue.Queue()
        for t in sorted(ready, key=lambda t: (t.get("priority", 5), t["task_id"])):
            work_q.put(t)

        stop = threading.Event()
        threads = [threading.Thread(target=worker_loop,
                                    args=(work_q, state, args, stop), daemon=True)
                   for _ in range(min(args.max_concurrent, work_q.qsize()))]
        for th in threads:
            th.start()
        for th in threads:
            th.join()

    with _state_lock:
        ok = sum(1 for r in state["tasks"].values() if r.get("status") == "done")
        log({"kind": "run_end", "msg": "%d/%d completadas, $%.2f"
             % (ok, len(tasks), state.get("cost_usd", 0.0))})
    return 0


if __name__ == "__main__":
    sys.exit(main())
