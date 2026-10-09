#!/usr/bin/env bash
# CAPS Studio for Linux: native core, self-contained .NET publish (linux-x64 or linux-arm64), then
#   dist/CAPS-<version>-linux-<arch>.AppImage   (any distribution)
#   dist/caps-studio_<version>_<debarch>.deb    (Debian, Ubuntu)
#   dist/CAPS-<version>-linux-<arch>.tar.gz     (portable)
# The app lives in /opt/caps (deb) or the AppImage's usr/lib/caps, with data/ and samples/ next to it.
#   packaging/linux/build.sh [x86_64|aarch64]
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
ARCH="${1:-$(uname -m)}"
case "$ARCH" in x86_64) RID=linux-x64; DEBARCH=amd64 ;; aarch64|arm64) ARCH=aarch64; RID=linux-arm64; DEBARCH=arm64 ;; *) echo "arch?"; exit 2 ;; esac
VERSION="$(sed -n 's/^project(CAPS VERSION \([0-9.]*\).*/\1/p' "$ROOT/CMakeLists.txt")"
ID=io.github.muhammaduzairriaz.CAPS
WORK="$ROOT/build-pkg/linux-$ARCH"
OUT="$ROOT/dist"
APP="$WORK/app"
export DOTNET_CLI_TELEMETRY_OPTOUT=1
rm -rf "$WORK" && mkdir -p "$WORK" "$OUT" "$APP"

echo "== native core"
cmake -S "$ROOT" -B "$WORK/native" -DCMAKE_BUILD_TYPE=Release -DCAPS_BUILD_TESTS=OFF > "$WORK/cmake.log"
cmake --build "$WORK/native" -j "$(nproc)" > "$WORK/build.log"

echo "== Studio (.NET, self-contained, $RID)"
dotnet publish "$ROOT/studio/CapsStudio/CapsStudio.csproj" -c Release -r "$RID" --self-contained true \
       -p:CapsNativeDir="$WORK/native/capi" -p:Version="$VERSION" -o "$APP" > "$WORK/publish.log"
cp "$WORK/native/cli/caps" "$APP/"
mkdir -p "$APP/data"
cp -R "$ROOT/data/." "$APP/data/"   # every data folder
find "$APP/data" -name "__pycache__" -type d -prune -exec rm -rf {} +
cp -R "$ROOT/licenses" "$APP/licenses"
cp "$ROOT/LICENSE" "$APP/LICENSE"
cp -R "$ROOT/samples" "$APP/"
find "$APP" -maxdepth 1 -name "*.pdb" -delete
strip "$APP/caps" "$APP/libcaps.so" 2>/dev/null || true

echo "== self-test of the staged app"
"$APP/CapsStudio" --selftest "$APP/samples" "$WORK" | tee "$WORK/selftest.log" | tail -1
grep -q "all checks passed" "$WORK/selftest.log" || { grep "^FAIL" "$WORK/selftest.log"; exit 1; }

icons() {   # hicolor icons into $1/share/icons
  for s in 16 32 48 64 128 256 512; do
    mkdir -p "$1/share/icons/hicolor/${s}x${s}/apps"
    cp "$ROOT/packaging/icon/caps_${s}.png" "$1/share/icons/hicolor/${s}x${s}/apps/$ID.png"
  done
}

echo "== AppImage"
AD="$WORK/AppDir"
mkdir -p "$AD/usr/lib" "$AD/usr/bin" "$AD/usr/share/applications" "$AD/usr/share/metainfo"
cp -R "$APP" "$AD/usr/lib/caps"
ln -s ../lib/caps/CapsStudio "$AD/usr/bin/caps-studio"
ln -s ../lib/caps/caps "$AD/usr/bin/caps"
cp "$ROOT/packaging/linux/$ID.desktop" "$AD/usr/share/applications/"
cp "$ROOT/packaging/linux/$ID.metainfo.xml" "$AD/usr/share/metainfo/$ID.appdata.xml"
install -Dm644 "$ROOT/packaging/linux/$ID.mime.xml" "$AD/usr/share/mime/packages/$ID.xml"
icons "$AD/usr"
cp "$ROOT/packaging/linux/$ID.desktop" "$AD/"
cp "$ROOT/packaging/icon/caps_256.png" "$AD/$ID.png"
cat > "$AD/AppRun" <<'RUN'
#!/bin/sh
HERE="$(dirname "$(readlink -f "$0")")"
# "caps …" runs the command line; anything else opens the Studio
if [ "$(basename "$ARGV0")" = "caps" ] || [ "$1" = "--cli" ]; then [ "$1" = "--cli" ] && shift; exec "$HERE/usr/lib/caps/caps" "$@"; fi
exec "$HERE/usr/lib/caps/CapsStudio" "$@"
RUN
chmod +x "$AD/AppRun"
TOOL="${APPIMAGETOOL:-$WORK/appimagetool}"
if [ ! -x "$TOOL" ]; then
  curl -sSfL -o "$TOOL" "https://github.com/AppImage/appimagetool/releases/download/continuous/appimagetool-$ARCH.AppImage"
  chmod +x "$TOOL"
fi
# the AppStream metadata checked here without the network (the homepage cannot be reached while the repository is
# private), then appimagetool told not to check it again online
if command -v appstreamcli > /dev/null; then
  appstreamcli validate --no-net "$ROOT/packaging/linux/$ID.metainfo.xml" > "$WORK/appstream.log" 2>&1 || { cat "$WORK/appstream.log"; exit 1; }
fi
ARCH="$ARCH" "$TOOL" --appimage-extract-and-run -n "$AD" "$OUT/CAPS-$VERSION-linux-$ARCH.AppImage" > "$WORK/appimage.log" 2>&1 || { tail -20 "$WORK/appimage.log"; exit 1; }

echo "== deb"
DEB="$WORK/deb"
mkdir -p "$DEB/DEBIAN" "$DEB/opt" "$DEB/usr/bin" "$DEB/usr/share/applications" "$DEB/usr/share/metainfo"
cp -R "$APP" "$DEB/opt/caps"
ln -s /opt/caps/CapsStudio "$DEB/usr/bin/caps-studio"
ln -s /opt/caps/caps "$DEB/usr/bin/caps"
cp "$ROOT/packaging/linux/$ID.desktop" "$DEB/usr/share/applications/"
cp "$ROOT/packaging/linux/$ID.metainfo.xml" "$DEB/usr/share/metainfo/"
install -Dm644 "$ROOT/packaging/linux/$ID.mime.xml" "$DEB/usr/share/mime/packages/$ID.xml"   # .capsproj: a CAPS project
icons "$DEB/usr"
cat > "$DEB/DEBIAN/control" <<CTL
Package: caps-studio
Version: $VERSION
Section: science
Priority: optional
Architecture: $DEBARCH
Depends: libc6 (>= 2.35), libstdc++6, zlib1g, libfontconfig1, libice6, libsm6, libx11-6, libicu70 | libicu72 | libicu74 | libicu-dev
Maintainer: Muhammad Uzair Riaz <274318895+MuhammadUzairRiaz@users.noreply.github.com>
Homepage: https://github.com/MuhammadUzairRiaz/CAPS
Description: CAPS Studio - build, simulate and analyse amorphous polymer cells
 Chain Assembly and Packing Suite: grow and pack polymer cells, assign force fields,
 relax, equilibrate, run dynamics and crosslinking, and analyse the results.
CTL
dpkg-deb --root-owner-group --build "$DEB" "$OUT/caps-studio_${VERSION}_${DEBARCH}.deb" > /dev/null

echo "== tarball"
tar -C "$WORK" -czf "$OUT/CAPS-$VERSION-linux-$ARCH.tar.gz" --transform "s,^app,CAPS-$VERSION," app
ls -la "$OUT"
