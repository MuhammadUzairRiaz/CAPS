#!/usr/bin/env bash
# Build the C++ core, CLI, C ABI and tests, then the Avalonia Studio. Run from anywhere.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export DOTNET_ROOT="${DOTNET_ROOT:-/opt/homebrew/opt/dotnet/libexec}"
export DOTNET_CLI_TELEMETRY_OPTOUT=1
cmake -S "$ROOT" -B "$ROOT/build" -DCMAKE_BUILD_TYPE=Release
cmake --build "$ROOT/build" -j "$(sysctl -n hw.ncpu 2>/dev/null || nproc)"
ctest --test-dir "$ROOT/build" --output-on-failure
dotnet build -c Release "$ROOT/studio/CapsStudio/CapsStudio.csproj"
"$ROOT/studio/CapsStudio/bin/Release/net10.0/CapsStudio" --selftest "$ROOT/samples" "${TMPDIR:-/tmp}"
