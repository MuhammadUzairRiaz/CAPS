#!/usr/bin/env python3
"""Write release.json (format caps-release/1): what an installed CAPS Studio reads to tell whether a newer CAPS is out
and what changed.

    python3 tools/release/make_release_json.py -o dist/release.json --sums dist/SHA256SUMS.txt
    python3 tools/release/make_release_json.py -o <app>/data/release.json            (the copy bundled in an app)

The commit and its date come from git (HEAD, or --commit); the version from CMakeLists.txt (or --version); the title
and highlights from tools/release/highlights.json (or --highlights FILE, or --from-whatsnew for the newest section of
tools/docs/src/whatsnew.md). With --sums, the assets are the installers SHA256SUMS.txt lists, each with its SHA-256;
without, the file names CAPS's packaging scripts write for this version, without checksums (the copy inside an app is
written before its installer exists). Standard library only.
"""
import argparse
import datetime as dt
import json
import os
import re
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
REPO = "https://github.com/MuhammadUzairRiaz/CAPS"


def version_from_cmake():
    text = open(os.path.join(ROOT, "CMakeLists.txt"), encoding="utf-8").read()
    m = re.search(r"project\(CAPS VERSION ([0-9.]+)", text)
    if not m:
        sys.exit("make_release_json: no 'project(CAPS VERSION …)' in CMakeLists.txt")
    return m.group(1)


def git(*args):
    return subprocess.run(["git", "-C", ROOT, *args], capture_output=True, text=True, check=True).stdout.strip()


def utc_iso(text):
    """'2026-10-22T11:14:03+02:00' → '2026-10-22T09:14:03Z'."""
    t = dt.datetime.fromisoformat(text.strip().replace("Z", "+00:00"))
    if t.tzinfo is None:
        t = t.replace(tzinfo=dt.timezone.utc)
    return t.astimezone(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def asset_key(name, version):
    """The platform key the Studio asks for, from an installer's file name (None for other files)."""
    v = re.escape(version)
    rules = [
        (rf"^CAPS-{v}-macos-(arm64|x86_64)\.dmg$", lambda m: f"macos-{m.group(1)}"),
        (rf"^CAPS-{v}-windows-x64-setup\.exe$", lambda m: "windows-x64"),
        (rf"^CAPS-{v}-windows-x64-portable\.zip$", lambda m: "windows-x64-portable"),
        (rf"^CAPS-{v}-linux-(x86_64|aarch64)\.AppImage$", lambda m: f"linux-{m.group(1)}"),
        (rf"^CAPS-{v}-linux-(x86_64|aarch64)\.tar\.gz$", lambda m: f"linux-{m.group(1)}-tar"),
        (rf"^caps-studio_{v}_(amd64|arm64)\.deb$", lambda m: "linux-" + {"amd64": "x86_64", "arm64": "aarch64"}[m.group(1)] + "-deb"),
    ]
    for pattern, key in rules:
        m = re.match(pattern, name)
        if m:
            return key(m)
    return None


def default_assets(version):
    names = [f"CAPS-{version}-macos-arm64.dmg", f"CAPS-{version}-windows-x64-setup.exe", f"CAPS-{version}-windows-x64-portable.zip",
             f"CAPS-{version}-linux-x86_64.AppImage", f"CAPS-{version}-linux-x86_64.tar.gz", f"caps-studio_{version}_amd64.deb"]
    return {asset_key(n, version): {"file": n} for n in names}


def assets_from_sums(path, version):
    assets = {}
    for line in open(path, encoding="utf-8"):
        parts = line.split()
        if len(parts) < 2:
            continue
        sha, name = parts[0].lower(), os.path.basename(parts[-1].lstrip("*"))
        if not re.fullmatch(r"[0-9a-f]{64}", sha):
            sys.exit(f"make_release_json: {path}: '{line.strip()}' is not '<sha256>  <file>'")
        key = asset_key(name, version)
        if key:
            assets[key] = {"file": name, "sha256": sha}
    if not assets:
        sys.exit(f"make_release_json: {path} lists no CAPS {version} installer")
    return assets


def plain(md):
    """Markdown inline text as plain text."""
    md = re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", md)
    return re.sub(r"[*`]", "", md).strip()


def from_whatsnew(path):
    """The newest '## ' section: its heading is the title, each '### ' heading a highlight with its first line as detail."""
    title, items, cur = None, [], None
    for line in open(path, encoding="utf-8"):
        line = line.rstrip()
        if line.startswith("## "):
            if title is not None:
                break
            title = plain(line[3:])
        elif title is not None and line.startswith("### "):
            cur = {"title": plain(line[4:])}
            items.append(cur)
        elif cur is not None and "detail" not in cur and line.strip():
            text = plain(line.lstrip("- "))
            cur["detail"] = text if len(text) <= 240 else text[:text.rfind(" ", 0, 237)] + " …"
    if title is None:
        sys.exit(f"make_release_json: no '## ' section in {path}")
    return {"title": title, "highlights": items}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--output", required=True, help="where to write release.json")
    ap.add_argument("--sums", help="SHA256SUMS.txt of the release's files (assets with checksums)")
    ap.add_argument("--highlights", default=os.path.join(ROOT, "tools", "release", "highlights.json"))
    ap.add_argument("--from-whatsnew", action="store_true", help="highlights from the newest section of tools/docs/src/whatsnew.md")
    ap.add_argument("--commit", help="the build commit (default: git HEAD)")
    ap.add_argument("--date", help="its date, ISO 8601 (default: the commit's date from git)")
    ap.add_argument("--version", help="default: CMakeLists.txt")
    a = ap.parse_args()

    version = a.version or version_from_cmake()
    try:
        commit = git("rev-parse", a.commit or "HEAD")
        date = a.date or git("show", "-s", "--format=%cI", commit)
    except (subprocess.CalledProcessError, FileNotFoundError):
        if not (a.commit and a.date):
            sys.exit("make_release_json: not a git checkout; give --commit and --date")
        commit, date = a.commit, a.date
    if not re.fullmatch(r"[0-9a-fA-F]{7,40}", commit):
        sys.exit(f"make_release_json: '{commit}' is not a commit SHA")
    src = from_whatsnew(os.path.join(ROOT, "tools", "docs", "src", "whatsnew.md")) if a.from_whatsnew else \
        json.load(open(a.highlights, encoding="utf-8"))
    highlights = []
    for h in src.get("highlights", []):
        item = {"title": h["title"]}
        for k in ("detail", "page"):
            if h.get(k):
                item[k] = h[k]
        highlights.append(item)

    manifest = {
        "format": "caps-release/1",
        "version": version,
        "commit": commit.lower(),
        "date": utc_iso(date),
        "title": src.get("title", f"CAPS {version}"),
        "highlights": highlights,
        "notes_url": f"{REPO}/releases/latest",
        "assets": assets_from_sums(a.sums, version) if a.sums else default_assets(version),
    }
    os.makedirs(os.path.dirname(os.path.abspath(a.output)), exist_ok=True)
    with open(a.output, "w", encoding="utf-8", newline="\n") as f:
        json.dump(manifest, f, indent=1, ensure_ascii=False)
        f.write("\n")
    print(f"wrote {a.output}: CAPS {version} · {commit[:10]} · {manifest['date']} · {len(highlights)} highlights · {len(manifest['assets'])} files")


if __name__ == "__main__":
    main()
