#pragma once
// The coarse-graining workflow as commands (one implementation behind `caps <command>`, the C API, Python and the Studio's
// Coarse-grain page, as the DFT workbench's): each takes its options as a Json object (the CLI's --flags with '-' → '_',
// positional arguments as "inputs"), returns a Json report with "ok", "command" (the equivalent command line) and what it
// wrote, and records provenance beside what it writes.
#include <string>
#include <vector>

#include "caps/json.hpp"

namespace caps {

// cgmap (more as the workflow grows)
const std::vector<std::string>& cg_commands();
bool is_cg_command(const std::string& cmd);
bool cg_is_switch(const std::string& cmd, std::string name);
// whether the command takes this option (--name, -name; - and _ alike)
bool cg_is_option(const std::string& cmd, std::string name);
std::string cg_help(const std::string& cmd);
// Throws std::invalid_argument / std::runtime_error with a message for the user on failure.
Json cg_run(const std::string& cmd, const Json& args, const std::string& data_dir);

}  // namespace caps
