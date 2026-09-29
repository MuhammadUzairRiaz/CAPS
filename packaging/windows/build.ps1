# CAPS Studio for Windows: native core (MinGW-w64: gcc, g++, gfortran; built beforehand, see -Native), self-contained
# .NET publish (win-x64), staged app folder, Inno Setup installer dist\CAPS-<version>-windows-x64-setup.exe and a
# portable .zip.
#   (in an MSYS2 UCRT64 shell) cmake -S . -B build-pkg/native -G Ninja -DCMAKE_BUILD_TYPE=Release -DCAPS_BUILD_TESTS=OFF
#                              cmake --build build-pkg/native
#   pwsh packaging/windows/build.ps1 -Native build-pkg/native
param([string]$Native = "build-pkg/native")
$ErrorActionPreference = "Stop"
$Root = (Resolve-Path "$PSScriptRoot\..\..").Path
$Version = (Select-String -Path "$Root\CMakeLists.txt" -Pattern 'project\(CAPS VERSION ([0-9.]+)').Matches[0].Groups[1].Value
$Work = "$Root\build-pkg\windows"
$App = "$Work\app"
$Dist = "$Root\dist"
$Native = (Resolve-Path (Join-Path $Root $Native)).Path
Remove-Item -Recurse -Force $App -ErrorAction SilentlyContinue
New-Item -ItemType Directory -Force $Work, $App, $Dist | Out-Null
$env:DOTNET_CLI_TELEMETRY_OPTOUT = "1"

if (-not (Test-Path "$Native\capi\caps.dll")) { throw "no native core at $Native (build it with MinGW-w64 first)" }

Write-Host "== Studio (.NET, self-contained, win-x64)"
dotnet publish "$Root\studio\CapsStudio\CapsStudio.csproj" -c Release -r win-x64 --self-contained true `
  -p:CapsNativeDir="$Native\capi" -p:Version=$Version -o $App | Out-Null
Copy-Item "$Native\cli\caps.exe" $App
New-Item -ItemType Directory -Force "$App\data" | Out-Null
# every data folder (a new one is never left out), without Python caches
Get-ChildItem "$Root\data" -Directory | ForEach-Object { Copy-Item -Recurse $_.FullName "$App\data" }
Get-ChildItem "$App\data" -Recurse -Directory -Filter "__pycache__" | Remove-Item -Recurse -Force
Copy-Item -Recurse "$Root\licenses" "$App\licenses"
Copy-Item "$Root\LICENSE" "$App\LICENSE"
Copy-Item -Recurse "$Root\samples" $App
Get-ChildItem $App -Filter *.pdb -File | Remove-Item   # debug symbols (the samples folder keeps its .pdb structures)

Write-Host "== self-test of the staged app"
& "$App\CapsStudio.exe" --selftest "$App\samples" $env:TEMP | Tee-Object -Variable log | Out-Null
if (-not ($log -match "all checks passed")) { $log; throw "self-test failed" }

Write-Host "== installer"
$iscc = (Get-Command iscc -ErrorAction SilentlyContinue).Source
if (-not $iscc) { $iscc = "${env:ProgramFiles(x86)}\Inno Setup 6\ISCC.exe" }
& $iscc /Q "/DVersion=$Version" "/DSource=$App" "/DOutDir=$Dist" "$Root\packaging\windows\caps.iss"
Compress-Archive -Path "$App\*" -DestinationPath "$Dist\CAPS-$Version-windows-x64-portable.zip" -Force
Get-ChildItem $Dist
