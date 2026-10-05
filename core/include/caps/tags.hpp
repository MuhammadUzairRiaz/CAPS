#pragma once
// Tags (design/boards/Tags): named, coloured sets of atoms. An atom carries any of a structure's tags (up to 32: bit k
// of Atom::tags is System::tags[k]), so edits that keep, copy or delete atoms keep their tags with them. Tags are saved
// beside a structure file (NAME.tags.json) and become groups everywhere: LAMMPS group lines, GROMACS index groups,
// Analyze of a tag, the Python package.
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

constexpr size_t kMaxTags = 32;

// The tag's index, or -1.
int tag_index(const System& s, const std::string& name);
// Atoms (0-based) carrying tag k.
std::vector<size_t> tag_atoms(const System& s, int k);
// op: "set" (exactly these atoms), "add", "remove"; the tag is made when it does not exist ("remove" then does nothing).
// Returns the tag's index. Throws std::runtime_error past kMaxTags or for an empty name.
int tag_edit(System& s, const std::string& name, const std::string& colour, const std::vector<size_t>& atoms, const std::string& op);
// Deletes tag k: the atoms lose it and the later tags move down a bit.
void tag_delete(System& s, int k);
// A name LAMMPS and GROMACS take as a group name: letters, digits and _ (others become _), not starting with a digit.
std::string tag_group_name(const std::string& name);
// "1:20 25 30:40" (1-based, LAMMPS group id ranges) for these 0-based atoms.
std::string id_ranges(const std::vector<size_t>& atoms);
// {"format":"caps-tags","atoms":N,"tags":[{"name","colour","atoms":"1:20 25"}]} and back (a different atom count: none read).
std::string tags_json(const System& s);
bool tags_from_json(System& s, const std::string& json);
std::string tags_sidecar(const std::string& structure_path);   // PATH.tags.json
void write_tags(const System& s, const std::string& structure_path);   // nothing to write: an old sidecar removed
bool read_tags(System& s, const std::string& structure_path);
// "group NAME id …" lines for every tag (long lists continued with &).
std::string lammps_tag_groups(const System& s);

}  // namespace caps
