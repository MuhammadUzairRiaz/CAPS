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
