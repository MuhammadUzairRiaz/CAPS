// A YAML subset for CAPS recipes and pipelines: block mappings and sequences by indentation, "- " items (which may
// open a mapping), flow mappings { a: 1, b: x } and sequences [a, b], quoted and plain scalars, true/false/null and
// numbers, '#' comments. No anchors, tags or multi-line strings. Errors name the line.
#pragma once
#include <string>

#include "caps/json.hpp"

namespace caps {

Json yaml_parse(const std::string& text);

}  // namespace caps
