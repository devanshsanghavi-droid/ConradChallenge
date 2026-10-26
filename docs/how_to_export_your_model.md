# How to export your model for ResidualMap

ResidualMap reads an EPANET model (`.inp`). Almost every water modelling package can write one, and that export is the route for every package below: ResidualMap never parses a vendor's own database, whose layout is not published. Send the export as it comes out of the package. The app identifies a file by its content, not its name, fixes the common problems in a copy (a byte-order mark, a Windows code page, reaction orders, per-pipe bulk lines), converts Darcy-Weisbach headloss to Hazen-Williams when that moves the water age little (and says so when it cannot), and says in its readiness report what it fixed, what it assumed and what is missing.

Two checks after any export, whatever the package:

1. Open the `.inp` in a text editor and look at `[CONTROLS]` and `[RULES]`. Pump and valve rules are the part exports lose most often, and a pump with no control runs all week.
2. Look at the `[TIMES]` section's `START CLOCKTIME`. ResidualMap converts your grab samples' times of day with it, so it must be the time of day your model's hour 0 stands for.

The app's Get started page (`pages/1_Get_started.py`) shows the section for your package. Its third step, "Check a file", reads any file you have and says what it is, without building the map. Its fourth step builds a skeleton model from a GIS layer of the water mains, for a system that has no hydraulic model: the pipes and where they meet, with no source, tank, pump, valve, control or demand until they are added. Its fifth step does the same from a CAD drawing of the mains saved as DXF: topology and sizes only.

Each section below starts with a line the Get started page shows as the choice.

## WaterGEMS or WaterCAD (Bentley)

1. Open the model and make the scenario you want active: the export writes the active scenario only, with its own demands, controls and initial settings.
2. Convert hydropneumatic tanks first. EPANET has no hydropneumatic tank, so the export cannot carry one as it is: model it as an ordinary tank, or as a reservoir with a head pattern, and check the result.
3. File > Export > EPANET. Choose IDs, not labels, so every junction keeps an ID EPANET accepts (31 characters at most, each one different).
4. Check that the active scenario's controls came through (`[CONTROLS]` and `[RULES]` in the `.inp`).

The native files, `.wtg`, `.wtg.sqlite` and the older `.wtg.mdb`, are not read: their layout is not published, and the supported way to read them needs the licensed program on Windows. The app recognises them and asks for the export.

## InfoWater Pro (Autodesk)

1. Model Explorer > Command Center > Exchange > Export EPANET. This menu path comes from a third-party guide and has not yet been checked against Autodesk's own help; if your menus differ, ask your modeller for the EPANET export.
2. Pump and valve rules are often lost in this export: check `[CONTROLS]` and `[RULES]` in the `.inp`, and add what is missing.

The native project database (`.IWDB`) is not read.

## InfoWorks WS Pro

1. Attach a control object to the network first: the export takes its pump and valve settings from it.
2. Right-click the network > Export > To EPANET v2 File.
3. The export drops rules, the water quality sections and tank mixing. ResidualMap sets the water quality itself, so the quality sections do not matter; check the rules and the tank mixing (`[MIXING]`) by hand.

## KYPipe

1. KYPipe writes an EPANET file, but its menu path is not public: ask your modeller for the EPANET export.
2. Check `[CONTROLS]` and the `[TIMES]` section after the export.

## EPANET (.net or .inp)

1. An `.inp` is ready: upload it as it is.
2. A `.net` file is EPANET's own binary project: open it in EPANET 2 and use File > Export > Network to write the `.inp`.

## AutoCAD or Civil 3D

A drawing holds lines, not a hydraulic model. If a hydraulic model of the system exists, its EPANET export is always better.

1. Save As > DXF, R2000 or a later version, which keeps the drawing's unit (an R12 DXF carries none, so the unit is then asked); a binary DXF is read too.
2. For Civil 3D pressure networks, run AECTOACAD (Export Civil 3D Drawing to AutoCAD) first, so the pipes become plain lines a DXF reader can see.
3. Bind external references first (XREF > Bind): an xref that points to another drawing is not read.
4. A `drawing.bak` is an AutoCAD backup copy: rename it to `.dwg`, open it, and save it as DXF.
5. Without AutoCAD, the free ODA File Converter turns a DWG into a DXF: https://www.opendesign.com/guestfiles/oda_file_converter (a link only: ResidualMap does not install or run it).

Step 5 of the Get started page builds a skeleton model from the DXF: topology and sizes only. It reads the lines, polylines, arcs and splines of model space (mains inside blocks too, where AutoCAD draws them) and each pipe's size from its layer name (WATER-8IN, W_8_PVC) or from the nearest label ('8" DIP', '200 mm'); a number alone in a layer name (WATER-ZONE-2, PVC-ZONE-3) is asked, since it may be a zone, not a size, and so is one the labels on several of its lines dispute. Leave out the layer styles that are not mains in service (abandoned mains drawn dashed on the mains layer, streets, lot lines), and check the drawing's unit. The model has no source, tank, pump, valve, control, demand or elevation, and a drawing cannot say which crossing lines connect, so it needs the operator form (or EPANET, on the downloaded file, until the form is offered on the page; it is built, but not offered yet) before it can run.

## MicroStation

1. File > Export > DWG/DXF and choose DXF, then use step 5 of the Get started page (previous section), or send the water mains as a GIS layer (next sections).
2. A DGN V8 file needs MicroStation itself to convert. In a WaterGEMS for MicroStation project, the drawing is not the model: send the WaterGEMS EPANET export instead (first section).

## ArcGIS Pro or ArcGIS Online

1. ArcGIS Pro: Export Features, to a GeoPackage or a shapefile (zip its `.shp`, `.shx`, `.dbf`, `.prj` and `.cpg` together: the `.prj` says the coordinates' unit and the `.cpg` the text encoding).
2. ArcGIS Online: open the layer's item page and download it as a GeoPackage, a shapefile or GeoJSON.
3. A File Geodatabase (`.gdb`) is not read: export the water mains layer as above.

Step 4 of the Get started page builds a skeleton model from the layer: pipes, their sizes and lengths, and nodes where they meet. It has no source, tank, pump, valve, control or demand, so it is not ready to run until those are added (in EPANET, on the downloaded file, until a form for them is built). If a hydraulic model exists, its EPANET export is always better.

## QGIS

1. Right-click the water mains layer > Export > Save Features As > GeoPackage.

Step 4 of the Get started page builds a skeleton model from the layer: pipes, their sizes and lengths, and nodes where they meet. It has no source, tank, pump, valve, control or demand, so it is not ready to run until those are added (in EPANET, on the downloaded file, until a form for them is built). If a hydraulic model exists, its EPANET export is always better.

## Google Earth

1. Right-click the folder of lines > Save Place As > KMZ.
2. Lines drawn by hand in Google Earth are joined where their ends come within 10 m (about 33 ft) of each other, since hand-drawn ends often miss by several metres; check the table of what a smaller or larger distance does. Write each line's size in its name or description ("8 in", 8" or "200 mm", with no word after "in"), or the size of every pipe is asked.

Step 4 of the Get started page builds a skeleton model from the layer: pipes, their sizes and lengths, and nodes where they meet. It has no source, tank, pump, valve, control or demand, so it is not ready to run until those are added (in EPANET, on the downloaded file, until a form for them is built). If a hydraulic model exists, its EPANET export is always better.

## Paper or PDF maps

A PDF map book cannot be read as a network. Two ways forward:

1. Trace the mains over the scanned map in the QGIS Georeferencer, then save them as a GeoPackage (previous sections).
2. Ask your state's rural water association: many help small systems build a first hydraulic model.

## Not sure

Go to step 3 of the Get started page, "Check a file", and give it any file you have: it says what the file is, by its content, and what to send instead. If the file is a model, it shows the readiness report.

## Your grab samples

On the main page, turn the demo off and upload two files under "2. Your samples": the grab log, with your chlorine readings, and the tap map, which says which junction of the model each sampling tap is at. The templates show the columns; your own export does not have to match them. The Get started page offers them for download, the grab log also as an Excel workbook (in the repository they are in `docs/templates/`):

- `grab_log_template.csv`: date, time, tap_id, free_chlorine_mgL (optional: method, notes).
- `grab_log_total_template.csv`: the same with total_chlorine_mgL, for a chloraminated system. Free and total chlorine are never mixed: a log with both columns asks which one is the reading.
- `tap_map_by_id_template.csv`: tap_id, junction_id (optional: flush_min, notes). An optional role column marks a tap at a zone inlet: write inlet as its role and the inlet's node ID as its junction, and its readings set that inlet's residual instead of being fitted.
- `plant_log_template.csv`: month, temp_C (optional: toc_mgL, dose_mgL). Optional, and read only by the pilot script, not by the app.

Your own export does not need to match the templates. The app guesses the columns from the header (for example "Sample Point", "Sample Date/Time", "Total Chlorine (mg/L)") and shows its guesses for you to check. It reads commas, semicolons or tabs, decimal commas, Windows encodings, a title or "Lab:" lines above the header, Excel's "sep=;" line, month-name dates such as 10-Mar-26, and Excel workbooks with typed dates and times (a cover sheet before the data is passed over). A LIMS export with one row per analyte (an "Analyte" or "Parameter" column and a "Result" column) keeps only the rows of the chlorine you chose and lists the others; the species may sit in brackets ("Chlorine (Free)") or be written "Free Cl", FCL2 or TCL2. When the analytes do not say which rows are the chlorine you chose (a code such as CHLF, or two chlorine analytes that name no species), the app stops with the question: choose the "Analyte value" whose rows are the reading. A reading in µg/L, as its header or units column says, is divided by 1000, with a warning; a reading above 10 mg/L is listed, not used. A reading written "<0.02" or "ND" is kept, listed and marked on the map as below detection, and one written ">2.2" as above the measuring range; neither is fitted. A tap with 12 or more readings a day at a median step of 30 minutes or less on 2 or more days, or with 12 or more readings a day on 7 or more days at any step, is read as an online analyzer: it keeps one reading a day, the one nearest 10:00 when one lies within 60 minutes of it, and a day with none keeps no reading and is listed. Its readings are kept whole instead when you say the taps are a sampling study; an hourly grab-sampling study of a day or two is kept whole anyway. A date that reads two ways (03/04/2026) stops with the question until you choose the date format. Every row that is not used is listed with its reason, the rows dated before the last 3 months included.
