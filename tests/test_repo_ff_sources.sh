#!/bin/bash
# CAPS carries nothing from the force-field conversion program whose licence forbids redistributing derived data: no
# tracked file may name it. The name is built from pieces so this file does not match itself.
# usage: test_repo_ff_sources.sh SOURCE_DIR   (exit 77, skipped, outside a git checkout)
cd "$1" || exit 1
git rev-parse --is-inside-work-tree > /dev/null 2>&1 || { echo "not a git checkout: skipped"; exit 77; }
a="dl"
b="field"
hits=$(git grep -l -i -E "${a}_${b}|${a}${b}|${a}-${b}|${a} ${b}" -- . ':!docs/paper')
if [ -n "$hits" ]; then
  echo "files naming the program:"
  echo "$hits"
  exit 1
fi
echo "no tracked file names the program"
