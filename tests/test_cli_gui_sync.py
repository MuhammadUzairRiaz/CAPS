#!/usr/bin/env python3
"""The command line and CAPS Studio in step.

cli/gui_sync.json places every caps command and option in the Studio (or says why it has none). This test fails when:
  * caps usage() or the command dispatch gains a command the file does not place;
  * a command's usage lines name an option the file does not place for that command;
  * known_options() holds an option no command entry or hidden_options entry places;
  * an anchor FILE|TEXT no longer holds (the Studio control was renamed or removed);
  * a C API entry point named for a command is not declared in caps_c.h, or not bound by the Studio (Native.cs) when the
    command is in the Studio;
  * the Studio's "Copy command" builders emit an option the command does not take;
  * a React run through the C API (as the Studio runs it) and the same run through the command line disagree.

usage: test_cli_gui_sync.py CAPS_EXECUTABLE SOURCE_ROOT
"""
import json
import os
import re
import subprocess
import sys
import tempfile

caps, root = sys.argv[1], sys.argv[2]
studio = os.path.join(root, "studio", "CapsStudio")
fails = []


def fail(msg):
    fails.append(msg)
    print("FAIL", msg)


src = open(os.path.join(root, "cli", "main.cpp"), encoding="utf-8").read()
k0 = src.index("const std::set<std::string>& known_options() {")
known = set(re.findall(r'"(--[a-z0-9-]+)"', src[k0:src.index("return k;", k0)]))
dispatch = set(re.findall(r'cmd == "([a-z0-9-]+)"', src))

# the usage text as caps prints it
usage = subprocess.run([caps], capture_output=True, text=True).stderr
per = {}
cur = None
for line in usage.splitlines():
    m = re.match(r"\s+caps ([a-z][a-z0-9-]*)(?: ([a-z][a-z0-9-]*))?", line)
    if m:
        cur = m.group(1)
        if cur in ("ff", "ml", "job") and m.group(2):
            cur += " " + m.group(2)
        per.setdefault(cur, set())
    if cur:
        per[cur].update(re.findall(r"(--[a-z0-9][a-z0-9-]*)", line))
# commands named in passing on a usage line (caps vasp-conv setup|collect CASE · vasp-scan make|fit · …)
for line in usage.splitlines():
    for c in re.findall(r"·\s+(vasp-[a-z]+)", line):
        per.setdefault(c, set())

sync = json.load(open(os.path.join(root, "cli", "gui_sync.json"), encoding="utf-8"))
cmds = sync["commands"]
header = open(os.path.join(root, "capi", "include", "caps_c.h"), encoding="utf-8").read()
native = open(os.path.join(studio, "Interop", "Native.cs"), encoding="utf-8").read()


def anchor_ok(where, what):
    if "|" not in where:
        fail(f"{what}: anchor '{where}' is not FILE|TEXT")
        return
    path, text = where.split("|", 1)
    full = os.path.join(studio, path)
    if not os.path.exists(full):
        fail(f"{what}: {path} does not exist")
    elif text not in open(full, encoding="utf-8").read():
        fail(f"{what}: '{text}' not found in {path}")


# every command placed
for c in sorted(per):
    if c not in cmds:
        fail(f"caps {c}: in usage() but not in cli/gui_sync.json")
for c in sorted(dispatch):
    if not any(k == c or k.startswith(c + " ") for k in cmds):
        fail(f"caps {c}: dispatched in main() but not in cli/gui_sync.json")

placed = set()
for c, e in cmds.items():
    if e["status"] in ("gui", "partial"):
        anchor_ok(e["where"], f"caps {c}")
    for fn in e.get("capi", []):
        if not re.search(r"\b" + re.escape(fn) + r"\s*\(", header):
            fail(f"caps {c}: C API entry point {fn} not declared in caps_c.h")
        elif e["status"] in ("gui", "partial") and f'"{fn}"' not in native:
            fail(f"caps {c}: the Studio does not bind {fn} (Interop/Native.cs)")
    opts = e["options"]
    if opts == "shared":
        continue
    for o in sorted(per.get(c, set())):
        if o not in opts:
            fail(f"caps {c} {o}: in usage() but not placed in cli/gui_sync.json")
    for o, v in opts.items():
        placed.add(o)
        if o not in known and o != "--help":
            fail(f"caps {c} {o}: placed in cli/gui_sync.json but caps does not know it")
        if v == "page" or v.startswith("partial:") or v.startswith("cli:"):
            if v.startswith(("partial:", "cli:")) and len(v.split(":", 1)[1].strip()) < 10:
                fail(f"caps {c} {o}: say why ({v})")
            continue
        anchor_ok(v, f"caps {c} {o}")
for o in sorted(known - placed):
    if o not in sync["hidden_options"]:
        fail(f"{o}: a known option no command entry or hidden_options places")

# the Studio's Copy command builders emit only options their command takes
copy = open(os.path.join(studio, "ViewModels", "MainViewModel.CopyCommand.cs"), encoding="utf-8").read()
field = open(os.path.join(studio, "ViewModels", "FieldViewModel.Sync.cs"), encoding="utf-8").read()
cli_args = field[field.index("public string? CliArgs()"):field.index("public string ApplyCommand")]
field_opts = set(re.findall(r"(--[a-z][a-z-]+)", cli_args))
for m in re.finditer(r"public string (\w+Command)\(\)\s*\{(.*?)\n    \}", copy, re.S):
    body = m.group(2)
    c = re.search(r'caps ([a-z]+) \{', body)
    if not c:
        fail(f"{m.group(1)}: names no caps command")
        continue
    name = c.group(1)
    emitted = set(re.findall(r"(--[a-z][a-z-]+)", body))
    if "CliField(" in body:
        emitted |= field_opts
    for o in sorted(emitted - per.get(name, set())):
        fail(f"Studio {m.group(1)}: emits {o}, which caps {name} does not list")
    print(f"ok   Studio {m.group(1)} → caps {name}: {len(emitted)} options, all taken")
apply = field[field.index("public string ApplyCommand"):]
for o in sorted(set(re.findall(r"(--[a-z][a-z-]+)", apply)) - per.get("ff apply", set())):
    fail(f"Studio ApplyCommand: emits {o}, which caps ff apply does not list")

counts = {}
for e in cmds.values():
    counts[e["status"]] = counts.get(e["status"], 0) + 1
print(f"ok   {len(cmds)} commands placed ({', '.join(f'{k} {v}' for k, v in sorted(counts.items()))}); {len(known)} options known")

# the same React run through the C API (the Studio's path) and through the command line
lib = os.environ.get("CAPS_LIB")
if lib and os.path.exists(lib):
    sys.path.insert(0, os.path.join(root, "data", "python"))
    import caps as capsy   # noqa: E402
    with tempfile.TemporaryDirectory() as td:
        cell = os.path.join(td, "pe.data")
        r = subprocess.run([caps, "grow", "-o", cell, "--chains", "4", "--dp", "10", "--density", "0.5", "--units", "*CC*", "--seed", "3"],
                           capture_output=True, text=True)
        if r.returncode:
            fail("caps grow: " + r.stderr)
        opts = ["--template", "cc_crosslink", "--cycles", "6", "--per-cycle", "2", "--seed", "4", "--no-relax", "--between-chains",
                "--crosslinks", "per-chain:1", "--auto-capture", "--capture-max", "9"]
        summary = os.path.join(td, "net.json")
        r = subprocess.run([caps, "react", cell, "-o", os.path.join(td, "net.data"), *opts, "--summary", summary, "--quiet"], capture_output=True, text=True)
        if r.returncode:
            fail("caps react: " + r.stderr)
        cli = json.load(open(summary))
        doc = capsy.open(cell)
        doc.react("cc_crosslink", cycles=6, per_cycle=2, seed=4, relax=False, between_chains=True, crosslinks=("per_chain", 1.0),
                  auto_capture=True, capture_max=9.0, default_field=True, temperature=300.0, md_ps=0.0)
        api = doc.react_summary()
        keys = ["reactions", "crosslinks", "target", "chains", "byproducts", "initial_sites"]
        same = all(cli[k] == api[k] for k in keys) and [(l["chain_a"], l["unit_a"], l["chain_b"], l["unit_b"]) for l in cli["links"]] == \
            [(l["chain_a"], l["unit_a"], l["chain_b"], l["unit_b"]) for l in api["links"]]
        if not same:
            fail("React through the C API and the command line differ: " + json.dumps({k: (cli[k], api[k]) for k in keys}))
        else:
            print(f"ok   React: caps react and the C API give the same network ({cli['reactions']} reactions, {cli['crosslinks']} links)")
else:
    print("skip React parity: CAPS_LIB not set")

print(f"{len(fails)} failures" if fails else "all sync checks passed")
sys.exit(1 if fails else 0)
