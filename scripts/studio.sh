#!/usr/bin/env bash
# Launch CAPS Studio. Optional: a file to open, and a LAMMPS data file as its topology.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export DOTNET_ROOT="${DOTNET_ROOT:-/opt/homebrew/opt/dotnet/libexec}"
exec "$ROOT/studio/CapsStudio/bin/Release/net10.0/CapsStudio" "$@"
