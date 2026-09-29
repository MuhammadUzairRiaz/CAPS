#!/usr/bin/env bash
# CAPS Studio for macOS: native core (macOS 14+), self-contained .NET publish, CAPS Studio.app, CAPS-<version>-<arch>.dmg.
#   packaging/macos/build_dmg.sh [arm64|x86_64]
# Signing: ad hoc by default (runs on this Mac; Gatekeeper warns elsewhere). With CAPS_CODESIGN_ID="Developer ID
# Application: …" the app is signed with the hardened runtime; with CAPS_NOTARY_PROFILE (xcrun notarytool
# store-credentials) the .dmg is also notarised and stapled.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
ARCH="${1:-$(uname -m)}"
case "$ARCH" in arm64) RID=osx-arm64 ;; x86_64) RID=osx-x64 ;; *) echo "arch must be arm64 or x86_64"; exit 2 ;; esac
VERSION="$(sed -n 's/^project(CAPS VERSION \([0-9.]*\).*/\1/p' "$ROOT/CMakeLists.txt")"
OUT="$ROOT/dist"
WORK="$ROOT/build-pkg/macos-$ARCH"
APP="$WORK/CAPS Studio.app"
export DOTNET_ROOT="${DOTNET_ROOT:-/opt/homebrew/opt/dotnet/libexec}"
export DOTNET_CLI_TELEMETRY_OPTOUT=1
rm -rf "$WORK" && mkdir -p "$WORK" "$OUT"

echo "== native core ($ARCH, macOS 14+)"
cmake -S "$ROOT" -B "$WORK/native" -DCMAKE_BUILD_TYPE=Release -DCAPS_BUILD_TESTS=OFF \
      -DCMAKE_OSX_ARCHITECTURES="$ARCH" -DCMAKE_OSX_DEPLOYMENT_TARGET=14.0 > "$WORK/cmake.log"
cmake --build "$WORK/native" -j "$(sysctl -n hw.ncpu)" > "$WORK/build.log"

echo "== Studio (.NET, self-contained, $RID)"
dotnet publish "$ROOT/studio/CapsStudio/CapsStudio.csproj" -c Release -r "$RID" --self-contained true \
       -p:CapsNativeDir="$WORK/native/capi" -p:Version="$VERSION" -o "$WORK/publish" > "$WORK/publish.log"

echo "== app bundle"
mkdir -p "$APP/Contents/MacOS" "$APP/Contents/Resources"
cp -R "$WORK/publish/." "$APP/Contents/MacOS/"
cp "$WORK/native/cli/caps" "$APP/Contents/MacOS/caps"
cp "$ROOT/packaging/icon/caps.icns" "$APP/Contents/Resources/"
mkdir -p "$APP/Contents/Resources/data"
# every data folder (force fields, typing, polymers, crystals, reactions, manual, python …): a new one is never left out
cp -R "$ROOT/data/." "$APP/Contents/Resources/data/"
find "$APP/Contents/Resources/data" -name "__pycache__" -type d -prune -exec rm -rf {} +
cp -R "$ROOT/licenses" "$APP/Contents/Resources/licenses"
cp "$ROOT/LICENSE" "$APP/Contents/Resources/LICENSE"
cp -R "$ROOT/samples" "$APP/Contents/Resources/"
sed "s/@VERSION@/$VERSION/g" "$ROOT/packaging/macos/Info.plist.in" > "$APP/Contents/Info.plist"
# debug symbols are not shipped
find "$APP/Contents/MacOS" -name "*.pdb" -delete

echo "== signing"
if [ -n "${CAPS_CODESIGN_ID:-}" ]; then
  # inside out: every Mach-O file, then the bundle
  find "$APP/Contents/MacOS" -type f \( -perm -u+x -o -name "*.dylib" \) -print0 |
    xargs -0 codesign --force --timestamp --options runtime --entitlements "$ROOT/packaging/macos/entitlements.plist" --sign "$CAPS_CODESIGN_ID"
  codesign --force --timestamp --options runtime --entitlements "$ROOT/packaging/macos/entitlements.plist" --sign "$CAPS_CODESIGN_ID" "$APP"
else
  codesign --force --deep --sign - "$APP"
fi
codesign --verify --deep --strict "$APP"

echo "== disk image"
STAGE="$WORK/dmg"
mkdir -p "$STAGE"
cp -R "$APP" "$STAGE/"
ln -s /Applications "$STAGE/Applications"
DMG="$OUT/CAPS-$VERSION-macos-$ARCH.dmg"
rm -f "$DMG"
hdiutil create -quiet -volname "CAPS $VERSION" -srcfolder "$STAGE" -ov -format UDZO "$DMG" > /dev/null
if [ -n "${CAPS_CODESIGN_ID:-}" ]; then codesign --force --sign "$CAPS_CODESIGN_ID" "$DMG"; fi
if [ -n "${CAPS_NOTARY_PROFILE:-}" ]; then
  xcrun notarytool submit "$DMG" --keychain-profile "$CAPS_NOTARY_PROFILE" --wait
  xcrun stapler staple "$DMG"
fi
echo "wrote $DMG ($(du -h "$DMG" | cut -f1))"
