import sys, time
import screen_app, screen_builders, screen_main, screen_row11, screen_row12, screen_row14, screen_row15b, screen_row16, screen_row16b, screen_row17b
import screen_row18, screen_row18b, screen_row19, screen_row19b, screen_row20, screen_row21, screen_row21b, screen_row7, screen_row8, screen_studio2
J = {"Grow": screen_app.grow, "GrowAllAtom": lambda: screen_app.grow("All atoms").replace("<title>CAPS Grow — amorphous cell</title>", "<title>CAPS Grow — all-atom view</title>"),
     "ForceField": screen_app.forcefield, "PolymerBuilder": screen_builders.polymer, "SolvationBuilder": screen_builders.solvation, "Main": screen_main.build,
     "CompactLayout": screen_row11.compact, "ProjectHome": screen_row11.project_home, "MillionAtoms": screen_row12.large, "JupyterNotebook": screen_row12.notebook,
     "Motion": screen_row14.motion, "AccessibilityMap": screen_row14.a11y, "DownloadPage": screen_row14.download_page, "History": screen_row15b.history, "Compare": screen_row15b.compare,
     "SmartSelect": screen_row16.smart_select, "FigureComposer": screen_row16.figure, "Charges": screen_row16b.charges, "SurfaceArea": screen_row17b.sasa_board,
     "Polydispersity": screen_row18.polydispersity, "Tacticity": screen_row18b.tacticity, "DisplayStyles": screen_row19.display_styles, "AddHydrogens": screen_row19b.add_hydrogens,
     "LensView": screen_row19b.lens_view, "VisPipeline": screen_row20.pipeline, "ColourBy": screen_row20.colour_by, "OpenLammps": screen_row20.import_formats, "OpenGromacs": screen_row20.import_gromacs,
     "ClusterAnalysis": screen_row21.clusters, "RdfCoordination": screen_row21.rdf_board, "SliceBinning": screen_row21.slice_bin, "ExpressionSelect": screen_row21.expression_select,
     "ComputeProperty": screen_row21b.compute_property, "CreateBonds": screen_row21b.create_bonds, "Replicate": screen_row21b.replicate,
     "SelectionStereo": screen_row7.selection, "LightSplit": screen_row7.light_split, "ExportDialog": screen_row8.export_dialog, "Start": screen_studio2.welcome, "Appearance": screen_studio2.visual}
for n, f in J.items():
    t = time.time()
    h = f()
    if isinstance(h, tuple): h = h[0]
    open(f"stage21/project/{n}.dc.html", "w").write(h)
    print(n, len(h), round(time.time() - t, 1), flush=True)
