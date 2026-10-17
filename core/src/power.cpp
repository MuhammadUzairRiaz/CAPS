// Windows: run at full speed when Windows takes the process for background work (see parallel.hpp)
#include "parallel.hpp"

#ifdef _WIN32
#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0A00
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

namespace caps {

void prefer_fast_cores() {
#ifdef _WIN32
  static const bool done = [] {
    // Opt out of execution-speed power throttling (EcoQoS). Without it Windows 11 puts a process with no foreground window
    // (caps md from a terminal over SSH, a Studio job) on the efficiency cores of a hybrid processor at a low clock.
    PROCESS_POWER_THROTTLING_STATE s{};
    s.Version = PROCESS_POWER_THROTTLING_CURRENT_VERSION;
    s.ControlMask = PROCESS_POWER_THROTTLING_EXECUTION_SPEED;
    s.StateMask = 0;
    return SetProcessInformation(GetCurrentProcess(), ProcessPowerThrottling, &s, sizeof s) != 0;
  }();
  (void)done;
#endif
}

}  // namespace caps
