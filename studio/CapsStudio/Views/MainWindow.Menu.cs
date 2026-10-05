using CapsStudio.ViewModels;
using System.Linq;
using System;
using Avalonia.Controls;
using Avalonia.Input;

namespace CapsStudio.Views;

/// <summary>The menu bar (the system menu bar on macOS, in the window elsewhere), organised as desktop modelling suites
/// are (Materials Studio's File · Edit · View · Build · Modules · Tools · Help): every item runs a palette command or opens
/// a page, so the menu, the palette, the rail and the shortcuts never disagree.</summary>
public partial class MainWindow
{
    private NativeMenu BuildMenu()
    {
        var mac = OperatingSystem.IsMacOS();
        NativeMenuItem Item(string header, Action run, string? key = null)
        {
            var it = new NativeMenuItem(header);
            // on macOS the menu owns its key equivalents; elsewhere the window's own key handling does
            if (mac && key != null) it.Gesture = KeyGesture.Parse(key);   // only keys the window does not handle itself
            it.Click += (_, _) => run();
            return it;
        }
        NativeMenuItem Cmd(string header, string id, string? key = null) => Item(header, () => { if (!_vm.RunCommand(id)) _vm.Status = header.TrimEnd('…') + " is not available now"; }, key);
        NativeMenuItem Page(string header, int m, string? key = null) => Item(header, () => _vm.GoModule(m), key);
        NativeMenuItem Sub(string header, params NativeMenuItemBase[] items)
        {
            var sub = new NativeMenu();
            foreach (var i in items) sub.Add(i);
            return new NativeMenuItem(header) { Menu = sub };
        }
        NativeMenuItemSeparator Sep() => new();
        NativeMenuItem Top(string header, params NativeMenuItemBase[] items) => Sub(header, items);

        var menu = new NativeMenu();
        menu.Add(Top("File",
            Item("New project…", () => _vm.OpenNewProject()),
            Item("Open project…", () => _ = OpenProjectFileDialog()),   // its folder
            Item("Save project", () => _vm.Status = _vm.HasCapsProject ? _vm.SaveCapsProject() : "No project is open: New project… makes one"),
            Item("Close project", () => _vm.Status = _vm.HasCapsProject ? _vm.CloseCapsProject() : "No project is open"),
            Sep(),
            Cmd("Open…", "document.open"),
            Cmd("Open with preview…", "document.open preview"),
            Item("Open most recent", OpenMostRecent),
            Item("Read a log or table…", () => _ = ReadTextDialog()),
            Page("Project home…", 42),
            Item("Open in Jupyter notebook…", () => _vm.OpenInNotebook()),
            Sep(),
            Sub("Save as",
                Cmd("LAMMPS data…", "document.save data"),
                Cmd("PDB…", "document.save pdb"),
                Cmd("Tripos mol2…", "document.save mol2"),
                Cmd("Materials Studio .car / .mdf…", "document.save car"),
                Cmd("CIF…", "document.save cif"),
                Cmd("VASP POSCAR…", "document.save poscar"),
                Cmd("Trajectory (LAMMPS dump)…", "trajectory.save")),
            Sub("Export",
                Page("Simulation files for LAMMPS / GROMACS…", 68),
                Page("Structure files…", 21),
                Page("Image…", 90),
                Page("Movie…", 91),
                Page("Figure…", 18),
                Page("Figure bundle…", 26)),
            Sep(),
            Cmd("Close structure", "document.close"),
            Cmd("Start page", "start.open")));
        menu.Add(Top("Edit",
            Item("Undo", () => _vm.UndoEdit(false)),
            Item("Redo", () => _vm.UndoEdit(true)),
            Sep(),
            Cmd("Select by query…", "select.query"),
            Item("Selection & stereo…", () => { _vm.SetModule(8); _vm.SelectionOpen = true; }),
            Item("Brush to select…", () => _vm.OpenBrush()),
            Item("Bond rules…", () => _vm.OpenRules()),
            Item("H autopilot (on / off)", () => _vm.ToggleHAutopilot()),
            Item("History & snapshots…", () => { _vm.SetModule(8); _vm.HistoryOpen = true; }),
            Item("Compare states…", () => { _vm.SetModule(8); _vm.StatesOpen = true; }),
            Sep(),
            Page("Unit cell & supercell…", 58),
            Page("Periodic box (wrap / unwrap)…", 51),
            Page("Add hydrogens…", 66)));
        menu.Add(Top("View",
            Cmd("Reset view", "view.reset"),
            Cmd("Frame selection", "view.frame"),
            Cmd("Perspective / orthographic", "view.projection"),
            Sep(),
            Item("Project panel", () => Panes.Toggle(Body, 0)),
            Item("Analysis dock", () => Panes.Toggle(Centre, 2)),
            Item("Inspector", () => Panes.Toggle(Body, 2)),
            Sub("Shelves",
                Item("Sketch & edit", () => ToggleShelf("tools")),
                Item("View & panels", () => ToggleShelf("view")),
                Item("Modify", () => ToggleShelf("modify")),
                Sep(),
                Item("Make a shelf…", () => _ = MakeShelfDialog()),
                Item("Lock shelves in place", () => { _vm.Settings.ShelvesLocked = !_vm.Settings.ShelvesLocked; _vm.Settings.Save(); _vm.Status = _vm.Settings.ShelvesLocked ? "Shelves locked" : "Shelves unlocked"; }),
                Item("Reset this workspace", ResetWorkspace)),
            Sub("Workspace",
                Item("Sketch", () => UseWorkspace("Sketch")), Item("Assemble", () => UseWorkspace("Assemble")), Item("Analyse", () => UseWorkspace("Analyse")),
                Item("Present", () => UseWorkspace("Present")), Item("Mine", () => UseWorkspace("Mine"))),
            Sep(),
            Page("Display styles…", 65),
            Page("Colour by…", 24),
            Item("Appearance…", () => { _vm.SetModule(8); _vm.AppearanceOpen = true; }),
            Page("Four views", 25),
            Page("Split view", 34),
            Sep(),
            Page("Trajectory player", 32),
            Page("Visualize pipeline…", 20),
            Page("Render…", 19),
            Sep(),
            Cmd("Theme: Graphite (dark)", "theme dark"),
            Cmd("Theme: Paper (light)", "theme light")));
        menu.Add(Top("Build",
            Page("Polymer cell (grow chains)…", 0),
            Page("Polymer builder…", 13),
            Cmd("Repeat unit from the picked head and tail", "polymer.unit from picked"),
            Page("Copolymer builder…", 60),
            Page("Blend builder…", 16),
            Page("Coarse-grained melt…", 44),
            Sep(),
            Page("Packing (molecules into a box or cell)…", 5),
            Page("Solvation…", 31),
            Sep(),
            Page("Molecule builder…", 9),
            Page("Fragment library…", 35),
            Page("Crystal…", 29),
            Sub("Crystal tools (the open structure)",
                Item("Find primitive cell", () => _vm.FindPrimitiveCell()),
                Item("Niggli-reduced cell", () => _vm.NiggliCell()),
                Item("Conventional cell", () => _vm.ConventionalCell()),
                Sep(),
                Page("Redefine lattice…", 29),
                Page("Vacuum slab…", 29),
                Page("Nanowire…", 29),
                Page("Lattice parameters & supercell…", 58)),
            Page("Surface / slab…", 14),
            Page("Nanostructure…", 15),
            Page("Biomolecule…", 30),
            Sep(),
            Page("Density calculator…", 56),
            Page("Model resolution (united atom, CG)…", 67)));
        menu.Add(Top("Modify",
            Sub("Element (picked atoms)", _vmElements(Item)),
            Sub("Bond (two picked atoms)",
                Item("Single", () => _vm.BondOrderPicked(1)),
                Item("Double", () => _vm.BondOrderPicked(2)),
                Item("Triple", () => _vm.BondOrderPicked(3)),
                Item("Break", () => _vm.BreakBondPicked())),
            Sub("Geometry (picked atom)", MainViewModel.Geometries.Select(g => (NativeMenuItemBase)Item(g.Name, () => _vm.GeometryPicked(g.Id))).ToArray()),
            Sep(),
            Page("Add hydrogens…", 66),
            Page("Partial charges…", 50),
            Cmd("Clean up (UFF)", "edit.clean"),
            Sep(),
            Page("Unit cell & supercell…", 58),
            Page("Periodic box (wrap / unwrap)…", 51)));
        menu.Add(Top("Modules",
            Page("Force field…", 7),
            Page("Partial charges…", 50),
            Page("Electrostatics…", 64),
            Sep(),
            Page("Minimise…", 2),
            Page("Equilibrate…", 4),
            Page("Dynamics…", 3),
            Page("React (crosslink, cure)…", 6),
            Page("Reaction templates…", 45),
            Page("Torsion scan…", 33),
            Page("Mesoscale (DPD)…", 72),
            Sep(),
            Sub("Analysis",
                Page("Properties…", 1),
                Page("Mechanics…", 38),
                Page("Glass transition…", 47),
                Page("Diffusion…", 49),
                Page("Scattering…", 39),
                Page("Free volume…", 40),
                Page("Chain statistics…", 55),
                Page("Tacticity…", 62),
                Page("Polydispersity…", 59),
                Page("Orientation…", 52),
                Page("Surface area…", 57),
                Page("Interface…", 48),
                Page("Adsorption locator…", 70),
                Page("Sorption (gas uptake)…", 71),
                Page("Solvent screen…", 61),
                Page("Blend phase diagram…", 63),
                Page("File checks…", 17))));
        menu.Add(Top("Tools",
            Item("Command palette…", TogglePalette),
            Page("Jobs", 11),
            Page("Provenance…", 37),
            Sep(),
            Page("Recipes…", 53),
            Page("Batch runs…", 22),
            Page("Parameter sweep…", 43),
            Page("Macro recorder…", 36),
            Page("Compare cells…", 23),
            Page("Figure composer…", 54),
            Sep(),
            Page("Validation bench…", 12),
            Page("Colour-vision check…", 46),
            Page("Settings…", 10)));
        menu.Add(Top("Help",
            Page("Theory manual…", 41),
            Page("Studio tour", 92),
            Cmd("Check for updates…", "app.update")));
        return menu;
    }

    private NativeMenuItemBase[] _vmElements(Func<string, Action, string?, NativeMenuItem> item) =>
        MainViewModel.QuickElements.Select(e => (NativeMenuItemBase)item(e, () => _vm.ModifyElementPicked(e), null)).ToArray();
}
