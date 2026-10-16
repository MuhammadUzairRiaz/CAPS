#!/bin/bash
# CAPS carries nothing from the force-field conversion program whose licence forbids redistributing derived data: no
# tracked file may name it, and no tracked file may name the OPLS library converted from it as a library id (its
# ignore patterns in .gitignore aside, which keep local copies out of the repository). Names are built from pieces so
# this file does not match itself.
# usage: test_repo_ff_sources.sh SOURCE_DIR   (exit 77, skipped, outside a git checkout)
cd "$1" || exit 1
git rev-parse --is-inside-work-tree > /dev/null 2>&1 || { echo "not a git checkout: skipped"; exit 77; }
status=0
a="dl"
b="field"
hits=$(git grep -l -i -E "${a}_${b}|${a}${b}|${a}-${b}|${a} ${b}" -- . ':!docs/paper')
if [ -n "$hits" ]; then
  echo "files naming the program:"
  echo "$hits"
  status=1
fi
# Still pointing at the library until their pending rewrites land (other branches): the moltemplate OPLS copies'
# charge_increments_from, and the Studio Field page's Fill / companion-charges sources. Drop these exclusions then.
id="opls""2005"
hits=$(git grep -l -i -E "${id}" -- . ':!docs/paper' ':!.gitignore' ':!data/forcefields/oplsaa2008-moltemplate.json' \
  ':!data/forcefields/oplsaa2024-moltemplate.json' ':!studio/CapsStudio/ViewModels/FieldViewModel.cs')
if [ -n "$hits" ]; then
  echo "files naming the converted OPLS library as a library id:"
  echo "$hits"
  status=1
fi
[ $status = 0 ] && echo "no tracked file names the program or its converted OPLS library"
exit $status
