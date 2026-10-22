# Building the CAPS installers

| Product | Built on | Command | Output (dist/) |
|---|---|---|---|
| macOS (.dmg) | a Mac | `bash packaging/macos/build_dmg.sh` | `CAPS-<v>-macos-arm64.dmg` |
| Windows (setup .exe, portable .zip) | Windows 10/11 | see below | `CAPS-<v>-windows-x64-setup.exe`, `…-portable.zip` |
| Linux (AppImage, .deb, .tar.gz) | Linux or WSL 2 | `bash packaging/linux/build.sh x86_64` | `CAPS-<v>-linux-x86_64.AppImage`, `caps-studio_<v>_amd64.deb`, `…tar.gz` |

Every build runs the Studio self-test on the staged app and stops unless it ends with "all checks passed".

## Windows (native)

One-time setup:

1. MSYS2 (https://www.msys2.org). In the **MSYS2 UCRT64** shell:
   `pacman -S --needed git mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-gcc-fortran mingw-w64-ucrt-x86_64-cmake mingw-w64-ucrt-x86_64-ninja`
2. .NET 10 SDK (https://dotnet.microsoft.com/download), PowerShell 7 (`winget install Microsoft.PowerShell`).
3. Inno Setup 6 (`winget install JRSoftware.InnoSetup`).
4. Python 3 (`winget install Python.Python.3.12`): the build writes `data\release.json` with it.

Build (the native core in the UCRT64 shell, the rest in PowerShell 7 from the repository root):

    cmake -S . -B build-pkg/native -G Ninja -DCMAKE_BUILD_TYPE=Release -DCAPS_BUILD_TESTS=OFF
    cmake --build build-pkg/native
    pwsh packaging/windows/build.ps1 -Native build-pkg/native

## Windows core cross-compiled on a Mac

With Homebrew's `mingw-w64` (gcc, g++, gfortran for x86_64-w64-mingw32):

    cmake -S . -B build-pkg/win-native -DCMAKE_BUILD_TYPE=Release -DCAPS_BUILD_TESTS=OFF \
          -DCMAKE_TOOLCHAIN_FILE=$PWD/packaging/windows/mingw-cross.cmake
    cmake --build build-pkg/win-native --parallel 8

`caps.dll` and `caps.exe` link the MinGW runtimes statically and need only Windows' own DLLs. The Studio publishes for
win-x64 on any OS (`dotnet publish … -r win-x64 --self-contained true -p:CapsNativeDir=<win-native>/capi`); the
self-test and the Inno Setup installer need Windows, so copy the staged app (or the portable zip) there and run
`CapsStudio.exe --selftest samples <a scratch folder>`.

## Linux (or WSL 2 on Windows)

Ubuntu 22.04 or later: `sudo apt-get install -y build-essential cmake gfortran zlib1g-dev libfontconfig1 libice6 libsm6 file`
and the .NET 10 SDK (`sudo apt-get install -y dotnet-sdk-10.0`, or Microsoft's install script). Then
`bash packaging/linux/build.sh x86_64`. Arch: `makepkg` with `packaging/linux/PKGBUILD`.

## Releasing: release.json

An installed CAPS Studio offers a newer CAPS by reading one file, `release.json`, from the latest GitHub release
(`https://github.com/MuhammadUzairRiaz/CAPS/releases/latest/download/release.json`). It holds the version, the build
commit and its date, a title, a few highlights (each may name the Studio page to open), the release page and each
installer's file name and SHA-256 (format `caps-release/1`, written by `tools/release/make_release_json.py`). The
version decides first; while it stays the same, a release whose commit differs from the installed build's and whose
date is later counts as newer. Every packaging script also puts this build's own manifest in the app as
`data/release.json` (no checksums: the installer does not exist yet), so the Studio shows What's new after an update
without the network.

To make a release:

1. Write the highlights for users in `tools/release/highlights.json` (title and a handful of lines; `page` is a Studio
   page such as `semicrystalline` or a palette command id such as `module.open thermal`), and the full changelog in
   `tools/docs/src/whatsnew.md`. Commit both.
2. Build every installer from that same commit (the sections above). Each one bundles `data/release.json` for it.
3. Checksums and the manifest, in the repository root:

       (cd dist && shasum -a 256 CAPS-* caps-studio_* > SHA256SUMS.txt)
       python3 tools/release/make_release_json.py -o dist/release.json --sums dist/SHA256SUMS.txt

   (`sha256sum` instead of `shasum -a 256` on Linux.) The script prints the version, commit, date and the number of
   highlights and files; check that the commit is the one the installers were built from.
4. Upload `release.json` as a release asset beside the installers and `SHA256SUMS.txt`, on the release marked
   Latest: `gh release upload <tag> dist/release.json dist/SHA256SUMS.txt dist/CAPS-* dist/caps-studio_* --clobber`
   (or drag the files into the release on GitHub). Without `release.json` on the latest release, installed Studios
   see no update and stay quiet.
5. Check from any installed or source build: `caps --check-update` prints the latest release's commit and date and
   exits 0 (this build is the latest or newer), 1 (a newer release is out) or 2 (cannot tell: offline, turned off,
   no manifest).

The CI release job (a `v*` tag) does steps 3 and 4 by itself. `CAPS_NO_UPDATE_CHECK=1` turns the check off on a
computer (offline or managed installs); the self-test and the `--screenshot` and `--figures` tools never use the
network.
