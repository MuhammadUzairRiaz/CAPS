#pragma once
// The DFT surface & adsorption workbench as commands (one implementation behind `caps <command>`, the C API, Python and
// the Studio pages, so defaults and results are identical everywhere): each takes its options as a Json object (the
// CLI's --flags with '-' → '_', positional arguments as "inputs"), returns a Json report with "ok", "command" (the
// equivalent command line) and what it wrote, and records provenance beside the structures it writes.
#include <string>
#include <vector>

#include "caps/json.hpp"

namespace caps {

// sheet, terminate, validate, adsorb-dft, vasp-set, vasp-conv, vasp-scan, vasp-derived, vasp-jobs, vasp-check,
// vasp-progress, vasp-health, vasp-bind, vasp-analyze
const std::vector<std::string>& dft_commands();
bool is_dft_command(const std::string& cmd);
// an option that takes no value (--janus, --force, --json …)
bool dft_is_switch(const std::string& cmd, std::string name);
// --help text: what it does, every option with its default and the reason, examples
std::string dft_help(const std::string& cmd);
// Runs a command; throws std::invalid_argument / std::runtime_error with a message for the user on failure.
Json dft_run(const std::string& cmd, const Json& args, const std::string& data_dir);
// The CLI's flags and positionals as the args object
Json dft_args_from_cli(const std::vector<std::string>& positional, const std::vector<std::pair<std::string, std::string>>& flags);
// The command line for args (copyable, shown in the Studio)
std::string dft_command_line(const std::string& cmd, const Json& args);

}  // namespace caps
