#!/usr/bin/env bash
# Installs the freshly built CAPS Studio.app (build_dmg.sh) into /Applications, replacing the one there.
# Refuses while CAPS Studio is running, so an open session is never pulled from under the user.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
ARCH="${1:-$(uname -m)}"
APP="$ROOT/build-pkg/macos-$ARCH/CAPS Studio.app"
DEST="/Applications/CAPS Studio.app"
[ -d "$APP" ] || { echo "no build at $APP: run packaging/macos/build_dmg.sh first"; exit 1; }
if pgrep -f "$DEST/Contents/MacOS/CapsStudio" > /dev/null; then echo "CAPS Studio is running: quit it, then run this again"; exit 3; fi
rm -rf "$DEST.new"
ditto "$APP" "$DEST.new"
rm -rf "$DEST"
mv "$DEST.new" "$DEST"
codesign --verify --deep --strict "$DEST"
echo "installed $DEST ($(defaults read "$DEST/Contents/Info" CFBundleShortVersionString 2>/dev/null || echo "?"), built $(stat -f %Sm "$APP"))"
