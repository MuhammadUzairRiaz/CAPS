---
title: Draw and build molecules
lede: Draw styrene by hand, fuse rings into naphthalene and anthracene, type a curative as SMILES, edit a molecule in 3D with the Studio's tools, and send the result to Packing and the Polymer builder.
time: 15 minutes
level: first steps
uses: Studio · Molecule
---

## What you will make

Four small molecules, each built a different way: styrene drawn atom by atom, naphthalene and anthracene from fused benzene rings, the rubber accelerator CBS typed as SMILES, and chlorobenzene made by editing benzene in 3D. At the end, styrene goes to the Polymer builder as a repeat unit and CBS goes to Packing.

Everything happens in the **Molecule** module of CAPS Studio. Open it from the Start page (the **Molecule** card) or with **Build › Molecule builder…**.

![The Molecule module: the 2D sketch on the left, the Studio view with every shelf in the middle, the molecule's SMILES, conformers and buttons on the right.](03_molecule.jpg)

## 1. Find your way around

The module has three parts.

- **The 2D sketch, on the left.** Its tools sit above the drawing: Draw, Bond, Ring, Charge, Move, Cut and Erase, then the bond orders (single, double, triple), the ring list, Undo, Redo, Tidy and Clear, and the element keys with **…** for the periodic table.
- **The Studio view, in the middle.** It shows the molecule in 3D as soon as there is one, with every Studio shelf above it: select and lasso, translate and rotate, place atom, draw bond, cut, **+H**, **Fragments**, rings, and the **Modify** row.
- **The molecule's panel, on the right.** The SMILES, the conformers with their energies, the properties (formula, mass, stereocentres, valence check), the 3D embedding options and the buttons that send the molecule on.

Hover over any tool to read what it does. Over the drawing, the label next to an atom says what the next click will do.

> CAPS adds hydrogens for you from each atom's valence. You draw the heavy atoms; the hover label shows the hydrogens each one carries (CH₃, CH₂, OH).

## 2. Draw styrene by hand

Styrene is a benzene ring with a vinyl group, C₆H₅–CH=CH₂.

1. Click **Ring** (the hexagon) and check that the ring list says **Benzene**. Click in the middle of the sketch. A benzene ring appears and the panel reads C₆H₆, valence *ok*.
2. Click **Draw** (the pencil) and make sure **C** is the chosen element. Click one carbon of the ring. A new carbon appears, bonded to it: you now have toluene, C₇H₈.
3. Click the new carbon. Another carbon is added in a zig-zag: ethylbenzene, C₈H₁₀.
4. Click **Bond** (second tool), then the **double** key (two lines). Click the bond between the two chain carbons. It becomes a double bond: styrene, C₈H₈.

![The Bond tool with double bonds chosen, on styrene.](03c_molecule_bonds.jpg)

About half a second after each change, the view in the middle builds the molecule in 3D and the panel lists its conformers.

If you click the wrong place, press **Undo** (or ⌘Z on a Mac, Ctrl+Z on Windows and Linux, while the sketch has the focus).

> CAPS keeps every atom within its valence. Try a double bond from a ring carbon that already has one: the bond is not drawn, and the status line says *Not drawn: C takes 4 bonds*.

## 3. Fuse rings: naphthalene, anthracene, phenanthrene

1. Press **Clear** (the ×). With **Ring** and **Benzene**, click empty space: benzene.
2. Click one of the ring's bonds. A second benzene ring is fused onto that bond: naphthalene, C₁₀H₈.
3. Click a bond on the far edge of the new ring for anthracene (three rings in a line), or a bond on its angled edge for phenanthrene. Either way the panel reads C₁₄H₁₀, valence *ok*.

CAPS places the double bonds over the whole ring system each time, shifting the older ring's double bonds where needed. Every CH keeps exactly one hydrogen, so the result is fully conjugated.

![Phenanthrene from three fused benzene rings: C₁₄H₁₀, one H on every CH.](03b_molecule_rings.jpg)

The ring list also has saturated rings of 3 to 8 atoms. Click an atom instead of a bond to put a ring through that atom (a spiro ring).

## 4. Type a SMILES: the accelerator CBS

For a larger molecule, typing is quicker than drawing. CBS (N-cyclohexyl-2-benzothiazolesulfenamide) is a common sulfur-cure accelerator.

1. Click the SMILES box under **Identifiers** on the right and type `C1CCC(CC1)NSc1nc2ccccc2s1`.
2. The sketch draws the molecule and the view builds it: C₁₃H₁₆N₂S₂, 264.4 g/mol, 33 atoms.
3. The **Conformers** list shows five conformers minimised with GAFF2, lowest energy first. Click a row to see that conformer in the view.

To change how the 3D structure is made, open **3D embedding** in the panel: the method (distance bounds, or with a rotor search over the rotatable bonds), the number of conformers (1 to 50), the force field for the clean-up (GAFF2, or UFF for any element) and how hydrogens are added.

A SMILES and the sketch are the same molecule. Edit the drawing and the SMILES changes; edit the SMILES and the drawing is redrawn.

## 5. Any element: the periodic table

The element row has the common elements. For anything else:

1. Click **…** at the end of the element row. The periodic table opens over the sketch.
2. Click an element, for example **Se**. It gets its own key next to the others and becomes the element **Draw** places.
3. Click empty space: a selenium atom with its two hydrogens, written `[SeH2]`.

![The periodic table, opened from the … key.](03d_molecule_table.jpg)

## 6. Edit in 3D with the Studio's tools

The view in the middle is a real Studio view, so its tools work on the molecule. Here you turn benzene into chlorobenzene in 3D.

1. Clear the sketch and type `c1ccccc1` in the SMILES box. Benzene appears in the view.
2. In the shelf above the view, pick the **Select** tool (the arrow) and click one hydrogen in the view.
3. In the **Modify** row (the second row of the shelf), click **Cl**.

The hydrogen becomes chlorine in 3D. The SMILES box changes to `c1(ccccc1)Cl` and the sketch redraws as chlorobenzene. The molecule is not built again, so the 3D geometry you edited stays.

![Benzene with one hydrogen made Cl in the view: the SMILES and the sketch follow.](03e_molecule_3dedit.jpg)

Every other tool works the same way: **+H** adds missing hydrogens, **Fragments** attaches a group to a picked atom, and **draw bond** and **cut** work on the 3D structure. After an edit in 3D, the conformer list shows one entry, *edited in the Studio view*. That edited structure is what the buttons on the right use.

## 7. Send the molecules on

- **A repeat unit for a polymer.** Type `*CC(*)c1ccccc1` (styrene's repeat unit: `*` marks the head and the tail) and click **Use as repeat unit**. The Polymer builder opens with it as unit A.
- **Molecules to pack.** With CBS built, click **Use in Pack**. CAPS asks where to save it, then opens Packing with CBS on its list of molecules.
- **Into a structure you had open.** **Insert into document** puts the molecule into the structure that was open in the Studio before you came to the Molecule module. Click in the view where it should go; the mouse wheel turns it.
- **On its own.** **Open in Studio** keeps the molecule as a structure of the project and goes to the Studio. The save button writes mol2, PDB, XYZ or LAMMPS data.

## The same from the command line

The Studio uses the same core as the `caps` command, so these steps can be scripted:

```bash
caps build 'C=Cc1ccccc1' -o styrene.mol2 --conformers 5 --ff uff
caps build 'C1CCC(CC1)NSc1nc2ccccc2s1' -o cbs.mol2 --conformers 5
caps build 'c1ccccc1' -o benzene.mol2
caps edit benzene.mol2 --ops 'ring 1 2 6 aromatic' -o naphthalene.mol2
```

The last line fuses an aromatic six-ring onto the bond between atoms 1 and 2 of benzene, which gives naphthalene.

## Next

- [Molecule module guide](../modules/molecule.html): every tool and option on one page.
- [Tutorial 01](first-cell.html): grow a polymer cell from a repeat unit like the one above.
