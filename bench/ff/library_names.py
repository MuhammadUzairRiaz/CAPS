#!/usr/bin/env python3
"""One-off: give the converted force fields CAPS's own names (ids, files, names, notes) in data/forcefields and
data/typing, with no mention of the conversion tool. Ids lose their "-reftool" suffix; old ids still resolve in code.

usage: library_names.py [--dry-run]
"""
import json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FF, TY = os.path.join(ROOT, "data", "forcefields"), os.path.join(ROOT, "data", "typing")
DRY = "--dry-run" in sys.argv
TOOL = re.compile(r"DL_?FIELD|DL_F\b|reftool|source library", re.I)

TEXT = [   # most specific first
    (r"converted from the reference tool by CAPS", "converted by CAPS"),
    (r"Dummy atom, for the reference tool internal use only", "Dummy atom (no interactions)"),
    (r"DL_F can't detect", "not detected automatically:"),
    (r"the reference tool entries give different masses", "the source entries give different masses"),
    (r"the reference tool entries disagree", "the source entries disagree"),
    (r"the reference tool gives mass", "the source gives mass"),
    (r"matches the reference tool's LAMMPS output term by term and in energies / forces", "validated term by term, energies and forces in LAMMPS"),
    (r"matches the reference tool(?:'s)? (?:for|on)", "validated term by term for"),
    (r"matches the reference tool", "validated term by term"),
    (r"the reference tool (?:\d+(?:\.\d+)*)?: ?", ""),
    (r"\s*\(the reference tool\)", ""),
    (r"\(the reference tool, ", "("),
    (r", the reference tool\)", ")"),
    (r"the reference tool'?s (?:own )?(?:molecule )?templates", "reference templates"),
    (r"the reference tool'?s", "the source library's"),
    (r"\bthe reference tool\b", "the source library"),
    (r"\bDL_F\b ?", ""),
]


# second pass: the phrases the first pass leaves, reworded
FIX = [
    (r"; the source library library v[\d.]+ \(\d+\)", ""),
    (r"as distributed with the source library [\d.]+", "as distributed"),
    (r"(\d+) the source library (\S+) templates", r"\1 \2 reference templates"),
    (r"the source library's DL_POLY FIELD file", "the reference DL_POLY FIELD file"),
    (r"the source library's FIELD file", "the reference FIELD file"),
    (r"no buildable the source library example", "no buildable reference example"),
    (r"as the source library's (\S+) library assigns them", r"as the published \1 parameter set types them"),
    (r"as the source library does", "as the original parameter set does"),
    (r"\(the source library's ([^)]*)\)", r"(types \1)"),
    (r"under the source library's (\S+) key names", r"under the \1 key names"),
    (r"each key's the source library atom types", "each key's atom types"),
    (r"the source library's", "the original"),
    (r"as the source library\)", "as in the original)"),
    (r"the source library cannot", "the original distribution cannot"),
    (r"the source library and CAPS both", "the original and CAPS both"),
    (r"the source library", "the original distribution"),
    (r"\breftool: ", "converted: "),
]
CREDIT = re.compile(r"(the source library|the reference tool) \d+(\.\d+)* library \(C\. W\. Yong, the authors\)")


def clean(s):
    for pat, rep in TEXT + FIX:
        s = re.sub(pat, rep, s)
    s = re.sub(r"-reftool", "", s)
    return re.sub(r"  +", " ", s).strip()


def walk(x):
    if isinstance(x, dict):
        return {k: walk(v) for k, v in x.items()}
    if isinstance(x, list):
        # a reference to the conversion tool itself goes
        return [walk(v) for v in x if not (isinstance(v, str) and (re.search(r"Yong, the reference tool|the reference tool - a force field", v) or CREDIT.fullmatch(v)))]
    if isinstance(x, str):
        return clean(x) if TOOL.search(x) else x
    return x


def run():
    renames = {}
    for d in (FF, TY):
        for f in sorted(os.listdir(d)):
            if "-reftool" in f:
                renames[os.path.join(d, f)] = os.path.join(d, f.replace("-reftool", ""))
    clash = [n for n in renames.values() if os.path.exists(n) and n not in renames]
    if clash:
        sys.exit("would overwrite: " + ", ".join(clash))

    for d in (FF, TY):
        for f in sorted(os.listdir(d)):
            p = os.path.join(d, f)
            if f.endswith(".json"):
                data = json.load(open(p))
                if f == "catalogue.json":
                    for e in data["forcefields"]:
                        if e.get("origin") == "reftool":
                            e["origin"] = "converted"
                new = walk(data)
                out = renames.get(p, p)
                if not DRY:
                    json.dump(new, open(out, "w"), ensure_ascii=False, indent=1)
                    if out != p:
                        os.remove(p)
            elif f.endswith(".md"):
                t = open(p).read()
                t2 = "\n".join(clean(l) if TOOL.search(l) else l for l in t.split("\n"))
                if not DRY:
                    open(p, "w").write(t2)
    left = []
    for d in (FF, TY):
        for f in os.listdir(d):
            t = open(os.path.join(d, f), errors="replace").read()
            n = len(TOOL.findall(t))
            if n:
                left.append((f, n))
    print(f"{len(renames)} files renamed" + (" (dry run)" if DRY else ""), "· left:", left or "none")


if __name__ == "__main__":
    run()
