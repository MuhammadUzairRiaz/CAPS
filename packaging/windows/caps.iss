; CAPS Studio installer (Inno Setup 6). Built by packaging/windows/build.ps1:
;   iscc /DVersion=0.1.0 /DSource=<staged app folder> /DOutDir=<dist> caps.iss
#ifndef Version
  #define Version "0.1.0"
#endif
#ifndef Source
  #define Source "..\..\build-pkg\windows\app"
#endif
#ifndef OutDir
  #define OutDir "..\..\dist"
#endif

[Setup]
AppId={{6E4B7B2E-0B7C-4C8A-9E0A-CA95C0A5F001}
AppName=CAPS Studio
ChangesAssociations=yes
AppVersion={#Version}
AppVerName=CAPS Studio {#Version}
AppPublisher=CAPS
AppPublisherURL=https://github.com/MuhammadUzairRiaz/CAPS
AppCopyright=Copyright (c) 2026 Muhammad Uzair Riaz
LicenseFile=..\..\LICENSE
DefaultDirName={autopf}\CAPS
DefaultGroupName=CAPS
DisableProgramGroupPage=yes
OutputDir={#OutDir}
OutputBaseFilename=CAPS-{#Version}-windows-x64-setup
SetupIconFile=..\icon\caps.ico
UninstallDisplayIcon={app}\CapsStudio.exe
Compression=lzma2/max
SolidCompression=yes
WizardStyle=modern
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
PrivilegesRequiredOverridesAllowed=dialog
ChangesEnvironment=yes

[Tasks]
Name: "desktopicon"; Description: "Create a desktop shortcut"; GroupDescription: "Shortcuts:"
Name: "addtopath"; Description: "Add the caps command line to PATH"; GroupDescription: "Command line:"; Flags: unchecked

[Files]
Source: "{#Source}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\CAPS Studio"; Filename: "{app}\CapsStudio.exe"
Name: "{group}\Uninstall CAPS"; Filename: "{uninstallexe}"
Name: "{autodesktop}\CAPS Studio"; Filename: "{app}\CapsStudio.exe"; Tasks: desktopicon

[Registry]
; .capsproj opens in CAPS Studio (a CAPS project: its structures and sessions)
Root: HKA; Subkey: "Software\Classes\.capsproj"; ValueType: string; ValueName: ""; ValueData: "CAPS.Project"; Flags: uninsdeletevalue
Root: HKA; Subkey: "Software\Classes\CAPS.Project"; ValueType: string; ValueName: ""; ValueData: "CAPS project"; Flags: uninsdeletekey
Root: HKA; Subkey: "Software\Classes\CAPS.Project\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}\CapsStudio.exe,0"
Root: HKA; Subkey: "Software\Classes\CAPS.Project\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\CapsStudio.exe"" ""%1"""
; PATH: the user's (HKCU\Environment) for a per-user install, the machine's for an install for all users — HKLM has no
; "Environment" key of its own (writing there fails with code 87)
Root: HKCU; Subkey: "Environment"; ValueType: expandsz; ValueName: "Path"; ValueData: "{olddata};{app}"; Tasks: addtopath; Check: (not IsAdminInstallMode) and NeedsAddPath(ExpandConstant('{app}'))
Root: HKLM; Subkey: "SYSTEM\CurrentControlSet\Control\Session Manager\Environment"; ValueType: expandsz; ValueName: "Path"; ValueData: "{olddata};{app}"; Tasks: addtopath; Check: IsAdminInstallMode and NeedsAddPath(ExpandConstant('{app}'))

[Run]
Filename: "{app}\CapsStudio.exe"; Description: "Start CAPS Studio"; Flags: nowait postinstall skipifsilent

[Code]
function NeedsAddPath(Dir: string): boolean;
var Paths: string;
begin
  if IsAdminInstallMode then begin
    if not RegQueryStringValue(HKEY_LOCAL_MACHINE, 'SYSTEM\CurrentControlSet\Control\Session Manager\Environment', 'Path', Paths) then begin Result := True; exit; end;
  end else
    if not RegQueryStringValue(HKEY_CURRENT_USER, 'Environment', 'Path', Paths) then begin Result := True; exit; end;
  Result := Pos(';' + Uppercase(Dir) + ';', ';' + Uppercase(Paths) + ';') = 0;
end;
