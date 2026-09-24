# CAPS Studio for Windows: native core (MSVC), self-contained .NET publish (win-x64), staged app folder, Inno Setup
# installer dist\CAPS-<version>-windows-x64-setup.exe and a portable .zip.
#   pwsh packaging/windows/build.ps1
$ErrorActionPreference = "Stop"
$Root = (Resolve-Path "$PSScriptRoot\..\..").Path
$Version = (Select-String -Path "$Root\CMakeLists.txt" -Pattern 'project\(CAPS VERSION ([0-9.]+)').Matches[0].Groups[1].Value
$Work = "$Root\build-pkg\windows"
$App = "$Work\app"
$Dist = "$Root\dist"
Remove-Item -Recurse -Force $Work -ErrorAction SilentlyContinue
New-Item -ItemType Directory -Force $Work, $App, $Dist | Out-Null
$env:DOTNET_CLI_TELEMETRY_OPTOUT = "1"

Write-Host "== native core (MSVC)"
cmake -S $Root -B "$Work\native" -DCAPS_BUILD_TESTS=OFF | Out-Null
cmake --build "$Work\native" --config Release --parallel | Out-Null

Write-Host "== Studio (.NET, self-contained, win-x64)"
dotnet publish "$Root\studio\CapsStudio\CapsStudio.csproj" -c Release -r win-x64 --self-contained true `
  -p:CapsNativeDir="$Work\native\capi\Release" -p:Version=$Version -o $App | Out-Null
Copy-Item "$Work\native\cli\Release\caps.exe" $App
New-Item -ItemType Directory -Force "$App\data" | Out-Null
Copy-Item -Recurse "$Root\data\forcefields", "$Root\data\typing", "$Root\data\reference" "$App\data"
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
