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
AppVersion={#Version}
AppVerName=CAPS Studio {#Version}
AppPublisher=CAPS
AppPublisherURL=https://github.com/MuhammadUzairRiaz/CAPS
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
Root: HKA; Subkey: "Environment"; ValueType: expandsz; ValueName: "Path"; ValueData: "{olddata};{app}"; Tasks: addtopath; Check: NeedsAddPath(ExpandConstant('{app}'))

[Run]
Filename: "{app}\CapsStudio.exe"; Description: "Start CAPS Studio"; Flags: nowait postinstall skipifsilent

[Code]
function NeedsAddPath(Dir: string): boolean;
var Paths: string;
begin
  if not RegQueryStringValue(HKEY_CURRENT_USER, 'Environment', 'Path', Paths) then begin Result := True; exit; end;
  Result := Pos(';' + Uppercase(Dir) + ';', ';' + Uppercase(Paths) + ';') = 0;
end;
