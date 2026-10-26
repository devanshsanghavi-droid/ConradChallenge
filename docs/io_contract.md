# ResidualMap's inputs and outputs

This page lists every file ResidualMap takes and every file and number it gives back, in plain words, with the columns, units and settings. Its machine-readable twin is `residualmap/contract.py`, and the saved check `io_contract_matches_code` holds the two to the code in both directions: every command-line flag, every column the readers take, every template header and every labelled input, title and header of both pages is in the contract and named here, and every name the contract holds is in the code. Every table on this page is the one the contract renders, in its place (`python -m residualmap.contract --write-tables` rewrites them), and so is every sentence of the prose that quotes a number (a limit, a default or a setting), each rendered from the code's own values; the check fails on any other table or row, any other number in the prose, and any phrase in quotes that names no label, section, option or example of the contract. Not listed: captions, expander titles and the pages' other text, and the extra output columns of the pilot's optional modes; the prose's other words are checked only for the names they quote. The Get started page's 'What files it takes' table and its template downloads are read from the same module.

Everything ResidualMap has been tested on is a public model (EPA's examples and the KYPIPE research models bundled with WNTR) with a simulated truth. Its importers are tested as converter self-consistency on those public models, written out by our own writers; no utility's own files and no real grab samples have been used.

## Units

| quantity | unit |
|---|---|
| chlorine | mg/L as Cl2: free chlorine, or total chlorine in the chloramine mode (never both in one run) |
| hour | 0 to 23, a model hour counted from the model's start time (START CLOCKTIME); shown as a time of day when the model does not start at midnight |
| water age | hours |
| bulk decay | per day |
| wall decay | metres per day |
| pipe length | feet or metres on input, as the file says; metres inside the tool |
| pipe size | inches or millimetres on input |
| flow | the model file's own EPANET units |

## What files it takes

The Get started page shows this table in step 3, where files can be given: on the laptop, and in the hosted demo only when its uploads are turned on (they are off by default, and then the hosted demo takes no file).

| What | Formats | Where it goes | Needed? |
|---|---|---|---|
| Your network model | An EPANET model (.inp), as EPANET writes it or as your modelling software exports it; a zip holding one model is fine | Main page, '1. Your network': 'Upload my .inp' | Yes, or build a skeleton in step 4 or step 5 |
| A GIS layer of your water mains | GeoPackage (.gpkg), zipped shapefile (.zip), GeoJSON, or a Google Earth KML or KMZ | This page, step 4 | Only with no hydraulic model: it gives pipes and nodes, and the source, tanks, pumps and demand are added in EPANET before it runs |
| A CAD drawing of your water mains | DXF, text or binary (save a DWG as DXF first) | This page, step 5 | Only with no hydraulic model and no GIS layer: topology and sizes only |
| Grab log | CSV or Excel (.xlsx): the date, the time, the tap ID and one chlorine reading column, free or total chlorine, never both | Main page, '2. Your samples', with the demo off | Yes, or type the samples in the main page's table |
| Tap map | CSV: each tap's ID and the junction ID it is at in the model | Main page, beside the grab log | Yes, with a grab log |
| Inlet residuals | CSV: each inlet's node ID in the model and its measured residual in mg/L | Main page, 'My zone's inlets carry different residuals' | Only for a zone whose inlets carry different residuals; free chlorine only |
| Plant log | CSV: the month and the water temperature (TOC and dose are optional) | The pilot script only; the app does not read it | Optional |
| Operator form (sources, tanks, pumps, valves, demand) | A table, one row per item | Not taken by the app yet: it is built and tested, but not offered here | Until it is offered, add those to a skeleton model in EPANET |

## 1. The network model

**An EPANET model** (`.inp`). It can come from EPANET itself or from any modelling package's EPANET export; `docs/how_to_export_your_model.md` gives the steps for each package, and the Get started page shows them. The file is identified by what is inside it, not by its name, so a renamed file, a backup copy or a zip holding one model is read too, in any common text encoding.

- Read as given: pipes (length, diameter, roughness), junctions (elevation, base demand, pattern), tanks, reservoirs, pumps and their curves, valves, controls and rules, the flow units and headloss formula, the initial quality and the coordinates (used for the maps only). `START CLOCKTIME` is recorded and used to turn clock times into model hours.
- Set by ResidualMap: a 7-day run with hourly hydraulic steps and a 5-minute quality step, of which only the last day is scored; the bulk and wall decay rates, which are calibrated; first-order reactions (a file declaring another order is fixed in a copy, with a note); and the chlorine source: water enters at every reservoir, or else at junctions with an inflow, or else at the tanks, all at one dose (or each inlet at its own residual, below).
- The readiness report says what came from the file, what was derived, what is assumed, what is still needed and why, and what was fixed in a copy (per-pipe reaction lines cleared, a Darcy-Weisbach file converted to Hazen-Williams when the conversion passes its water-age gate, a model with zero elevations next to pressure valves flagged). A model of more than 3,000 junctions needs a confirm before its grid is built: the main page's 'Build the grid anyway: N runs, about M min and S MB (a large model)' box, or the pilot's `--allow-large`.

**A GIS layer of the water mains** (Get started page, step 4): GeoPackage, zipped shapefile, GeoJSON, or a Google Earth KML or KMZ. Pipes are its lines; the size, material, install year, declared length and pipe ID come from its fields when it has them (guessed, and shown for checking). The coordinate system comes from the file; longitude and latitude are handled without extra software; a layer that does not say its unit is asked. Line ends closer than the joining distance (6 ft by default, 10 m for lines drawn by hand in Google Earth) become one node. The result is a skeleton model: pipes and nodes only, with no source, tank, pump, valve, control or demand, so it does not run until those are added. Each pipe's Hazen-Williams C comes from its material alone (new-pipe values), so it cannot see a pipe's age.

**A CAD drawing of the water mains** (Get started page, step 5): DXF, text or binary, or a `.bak` holding one. Pipes are the lines, polylines and arcs of model space on the layer styles kept; sizes come from layer names and nearby labels; the unit comes from the drawing's header or is asked. The result is a skeleton model too: topology and sizes only.

**Recognised but not read.** ResidualMap says what such a file is and what to send instead: a WaterGEMS or WaterCAD project (`.wtg`, which holds its settings, with its `.wtg.sqlite` or `.wtg.mdb` model database), another Access or SQLite database that is not a GeoPackage, or a zip of them: send the EPANET export; a DWG: save it as DXF or export the mains as a GIS layer; a MicroStation DGN: export DXF or a GIS layer; an IFC (BIM) file or a PDF map: send the EPANET export, or else the mains as a GIS layer; an Esri file geodatabase: export a shapefile or GeoPackage; an EPANET results file (`.bin`): send the `.inp` that made it. A file it does not recognise is called unknown, with the same advice to send the EPANET export (EPANET's own binary `.net` project is expected to be one, but no `.net` file was tested).

## 2. Grab samples

**Typed on the main page.** The 'Grab samples' table has three columns: `junction` (a junction ID of the model), `hour` (the model hour, 0 to 23) and `mg/L` (the reading, 0 to 5 mg/L).

**A grab log file** (CSV or Excel `.xlsx`; main page, '2. Your samples', with the demo off). The templates `grab_log_template.csv` (free chlorine) and `grab_log_total_template.csv` (total chlorine, for a chloraminated system) have the columns `date`, `time`, `tap_id`, `free_chlorine_mgL`, `method`, `notes` and `date`, `time`, `tap_id`, `total_chlorine_mgL`, `method`, `notes`. The reader fills these roles from the log's columns:

| role | the template's column | what it is |
|---|---|---|
| `tap` | `tap_id` | the tap or sample point's ID |
| `reading` | `free_chlorine_mgL`, or `total_chlorine_mgL` for a chloraminated system | the chlorine reading in mg/L; exactly one species, never both |
| `date` | `date` | the sample's date (YYYY-MM-DD, or any common form; a date that reads two ways is asked) |
| `time` | `time` | the clock time (HH:MM), or none when the time is in the date column |
| `analyte` | (none) | a long-format export's analyte column, one row per analyte |
| `unit` | (none) | a long-format export's units column (µg/L, ug/L and ppb are divided by 1000) |

The columns are guessed from the header, so a lab or field export's own names work too ('Sample Point', 'Free Chlorine (mg/L)', 'Sample Date/Time', a LIMS code), and the main page shows the guesses for checking. The pilot's `--columns` takes the same roles as keys, `tap`, `reading`, `date`, `time`, `analyte`, `unit`, `analyte_value`, where `analyte_value` names the one analyte whose rows are the reading. The template's `method`, `notes` columns are for the operator and are not read. Any separator, decimal commas and Windows encodings are read. Readings below detection ('<0.02', 'ND') and above the measuring range are shown and listed, not fitted. A tap with 12 or more readings a day at a median step of 30 minutes or less on 2 or more days, or with 12 or more readings a day on 7 or more days at any step, is read as an online analyzer: it keeps one reading a day, the one nearest 10:00 when one lies within 60 minutes of it, and a day with none keeps no reading and is listed. The main page's 'Keep every reading of these taps' box and the pilot's `--keep-every-reading` keep them all. A reading above 10 mg/L is listed, not used. Each clock time becomes a model hour: the time minus the model's start time, with the minutes dropped. The main page fits the last 3 calendar months that hold a reading.

**The tap map** (CSV; `tap_map_by_id_template.csv`, columns `tap_id`, `junction_id`, `flush_min`, `notes`). The reader fills the roles `tap` (the template's `tap_id`) and `junction` (`junction_id`, the model's junction ID); an optional `role` column marks a tap as a zone `inlet`, whose readings then set that inlet's residual and are not fitted. The `flush_min` and `notes` columns are for the operator and are not read. A tap map by coordinates (latitude and longitude, or x and y) is not taken yet: its reader is built, but it missed its pre-registered bar in simulation.

**Inlet residuals** (CSV; free chlorine only; `inlet_residuals_template.csv`, columns `inlet_node_id`, `residual_mgL`, `notes`): the roles `inlet` (the template's `inlet_node_id`, the inlet's node ID in the model) and `residual` (`residual_mgL`, its measured residual in mg/L). Each inlet is then dosed in proportion to its residual, and the dose input means the highest residual entering the zone. The `notes` column is not read. This template is not offered on the Get started page: the pilot's `--inlets` help names it, and the main page's uploader says its two columns.

## 3. The plant log (the pilot script only)

`plant_log_template.csv`, one row per month: `month` (YYYY-MM) and `temp_C` (the plant water's temperature) are required; `dose_mgL` (that month's dose) and `toc_mgL` (organic carbon, carried and not used) are optional. It is read as a grab log is (any separator, decimal commas, any common encoding, or an Excel sheet): a row that leaves out the optional cells at its end is read with them blank, and one with a value beyond the header is refused; a temperature that is blank, not a number or outside 0 to 35 C (the range of every rate law used here; a log whose temperatures are all above it is said to look like Fahrenheit), or a month not written YYYY-MM, is refused with the file's name and its rows. Free chlorine only. It makes the pilot's model temperature-aware, a model that was tested and not adopted as the default.

## 4. The operator form (built and tested, not yet in the app)

A table, CSV or Excel, one row per item, that completes a skeleton model from a GIS layer or a drawing (or fills gaps in an EPANET model). It is built and tested, but it is not yet in the app: the app does not take it, and its own bar on ky4 was missed in simulation. Its command lines are `python -m residualmap.gis_read LAYER --form FORM --out completed.inp`, the same with `residualmap.cad_read`, and `python -m residualmap.operator_form MODEL FORM`.

The `kind` column says what a row is: `source`, `tank`, `pump`, `pump_rule`, `closed_valve`, `check_valve`, `prv`, `psv`, `elevation`, `production`, `large_user`, `pattern`, `pipe_c`. Its columns: `kind`, `id`, `location`, `pipe`, `upstream`, `head_ft`, `dose_mgL`, `residual_mgL`, `base_elev_ft`, `diameter_ft`, `min_level_ft`, `max_level_ft`, `init_level_ft`, `power_hp`, `flow_gpm`, `pump`, `tank`, `on_below_ft`, `off_above_ft`, `on_hours`, `setting_psi`, `elevation_ft`, `diameter_in`, `gallons_per_day`, `multipliers`, `material`, `hw_c`, `note`. A location is a node ID, a point x,y in the drawing's or layer's own coordinates, or latitude,longitude for a layer in degrees. Other column names it takes, and the column each means: `type` (kind), `item` (kind), `name` (id), `node` (location), `at` (location), `node_id` (location), `pipe_id` (pipe), `suction_pipe` (pipe), `from` (upstream), `upstream_node` (upstream), `residual_entering_mgl` (residual_mgL), `residual_mgl` (residual_mgL), `dose_mgl` (dose_mgL), `gpd` (gallons_per_day), `gallons_per_day_gpd` (gallons_per_day), `c` (hw_c), `c_value` (hw_c), `hazen_williams_c` (hw_c), `notes` (note), `comment` (note), `base_elevation_ft` (base_elev_ft), `minimum_level_ft` (min_level_ft), `maximum_level_ft` (max_level_ft), `initial_level_ft` (init_level_ft), `setting_psig` (setting_psi), `flow_gpm_points` (flow_gpm), `head_ft_points` (head_ft). Other words it takes in the `kind` column: `reservoir` (source), `clearwell` (source), `well` (source), `plant` (source), `closed valve` (closed_valve), `boundary_valve` (closed_valve), `check valve` (check_valve), `cv` (check_valve), `pump rule` (pump_rule), `large user` (large_user), `demand` (large_user), `pipe c` (pipe_c), `roughness` (pipe_c). A row it cannot use is refused with its reason, and no model is built until the form is fixed; a model with no water source is never run.

## 5. Templates

| template | header | read as | on the Get started page |
|---|---|---|---|
| `grab_log_template.csv` | date,time,tap_id,free_chlorine_mgL,method,notes | grab log | yes: 'Grab log, free chlorine (CSV)' |
| `grab_log_total_template.csv` | date,time,tap_id,total_chlorine_mgL,method,notes | grab log | yes: 'Grab log, total chlorine for a chloraminated system (CSV)' |
| `tap_map_by_id_template.csv` | tap_id,junction_id,flush_min,notes | tap map | yes: 'Tap map by ID (CSV)' |
| `plant_log_template.csv` | month,temp_C,toc_mgL,dose_mgL | plant log | yes: 'Plant log, optional, read by the pilot script only (CSV)' |
| `inlet_residuals_template.csv` | inlet_node_id,residual_mgL,notes | inlet residuals | no: named in the pilot's --inlets help; the main page's inlet uploader says its two columns |

The Get started page also offers the free-chlorine grab log as an Excel workbook, 'Grab log, free chlorine (Excel)', written on the fly.

## 6. The main page's inputs

In the sidebar:

| label | kind | range and default | what it is |
|---|---|---|---|
| '1. Your network' | header |  | the model |
| 'Model' | radio |  | an example network, your own model file, or a model built on the Get started page (offered while one is kept in this session) |
| 'EPANET .inp file' | file uploader |  | your model: any file, identified by its content |
| 'Example' | selectbox |  | Net3, Net2 or ky4 (Net3 and Net2 only in the hosted demo) |
| 'Disinfectant' | radio |  | free chlorine, or chloramine (total chlorine); never mixed |
| 'Highest residual entering the zone (mg/L)' | number input | 0.05 to 4; default the highest residual given | with inlet residuals: the highest residual entering the zone, treated as uncertain by 10% |
| 'Total chlorine leaving the plant (mg/L as Cl2)' | number input | 0.5 to 4; default 2 | chloramine: total chlorine at the source (an effective dose of 0.65 to 1.10 times this) |
| 'Chlorine dose leaving the plant (mg/L)' | number input | 0.2 to 4; default 1.2 | free chlorine at the source, treated as uncertain by 10% |
| 'Minimum residual, total chlorine (mg/L)' | number input | 0.05 to 2; default 0.5 | chloramine: the threshold a junction is flagged below (a common utility target, not a California rule) |
| 'Minimum residual (mg/L)' | number input | 0.05 to 1; default 0.2 | free chlorine: the threshold a junction is flagged below |
| '2. Your samples' | header |  | the grab samples |
| 'Demo: simulate a hidden truth and draw samples from it' | toggle |  | on: a simulated network with samples drawn from it; off: your grab log, or samples typed in the table |
| 'Grab log (CSV or Excel)' | file uploader |  | your grab-sample log, as your lab or field system exports it |
| 'Tap map (CSV)' | file uploader |  | which junction each tap is at |
| 'My zone's inlets carry different residuals' | checkbox |  | a zone with several inlets, each with its own residual |
| 'Inlet residuals (CSV)' | file uploader |  | each inlet's node ID and its measured residual |
| 'Tap column' | selectbox |  | the grab log's column of tap IDs (guessed; check it) |
| 'Reading column' | selectbox |  | the grab log's chlorine reading column |
| 'Date column' | selectbox |  | the grab log's date column |
| 'Time column' | selectbox |  | the grab log's time column, or the time is in the date column |
| 'Analyte column (a log with one row per analyte)' | selectbox |  | a long-format export's analyte column |
| 'Analyte value (the analyte whose rows are the reading)' | selectbox |  | which analyte's rows are the chlorine reading |
| 'Units column' | selectbox |  | a long-format export's units column |
| 'Date format' | selectbox |  | how the dates are written: guessed, month first, day first, year first, or typed |
| 'Date format, as strptime writes it (%d day, %b month name, %m month number, %y or %Y year, %H:%M time)' | text input |  | a typed date format |
| 'Keep every reading of these taps (a sampling study, not an online analyzer)' | checkbox |  | keep a tap that looks like an online analyzer whole |
| 'Online analyzer: keep the reading nearest this hour (0 to 23)' | number input | 0 to 23; default 10, the reader's hour | an online analyzer keeps one reading a day, the one nearest this clock hour when one lies within 60 minutes of it |
| 'Demo samples (daytime, random junctions)' | slider | 0 to 15; default 8 | the demo's number of samples |
| 'Demo scenario' | number input | 0 to 99; default 0 | the demo's simulated scenario |
| 'Plant water pH (logged)' | number input | 6.5 to 9.5; default 8.0, or in the demo its hidden truth's | chloramine: sets the prior centre of the bulk decay rate |
| 'Chlorine to ammonia-N ratio, by mass (logged)' | number input | 3 to 5; default 4.5, or in the demo its hidden truth's | chloramine: the plant's Cl2:N |
| 'Water temperature this month (C)' | number input | 0 to 35; default 20 | chloramine: used only by the nitrification watch (15 C or warmer) |
| '3. Route' | header |  | next month's sampling route |
| 'Sites on next month's route' | slider | 3 to 12; default 6 | how many sites the route picks |

On the page itself:

| label | kind | range and default | what it is |
|---|---|---|---|
| 'Bring your data: what to export from your software, and 'Check a file'' | page link |  | to the Get started page, under a file that is not a model |
| 'Build the grid anyway: N runs, about N min and N MB (a large model)' | checkbox |  | a model of more than 3,000 junctions: confirm before its runs start |
| 'Show' | radio |  | the map's hour: the daily minimum (the compliance number) or one hour of the day |
| 'Show the true daily-minimum map and score the flags' | toggle |  | the demo only: reveal the hidden truth |
| 'Download the one-page PDF report' | download button |  | the PDF report |
| 'Bring your data: what to export from your software, and templates for the grab log and tap map' | page link |  | to the Get started page (the hosted demo, uploads off) |
| 'Bring your data: what to export from your software, templates for the grab log and tap map, and 'Check a file'' | page link |  | to the Get started page |

## 7. The Get started page's inputs

| label | kind | range and default | what it is |
|---|---|---|---|
| 'Bring your data' | title |  | the Get started page |
| 'Step 1. Which software holds your network?' | header |  | export steps per software |
| 'Which software holds your network?' | radio |  | each choice shows its export steps |
| 'Step 2. Your grab samples: templates' | header |  | the templates |
| 'Grab log, free chlorine (CSV)' | download button |  | a template |
| 'Grab log, total chlorine for a chloraminated system (CSV)' | download button |  | a template |
| 'Tap map by ID (CSV)' | download button |  | a template |
| 'Plant log, optional, read by the pilot script only (CSV)' | download button |  | a template |
| 'Grab log, free chlorine (Excel)' | download button |  | a template, written on the fly |
| 'Step 3. Check a file' | header |  | any file: what it is, by its content, and what to send instead |
| 'What files it takes' | subheader |  | the table of every file ResidualMap takes (above) |
| 'Disinfectant' | radio |  | the species a grab log or a model is checked for |
| 'Check a file' | file uploader |  | any file; nothing is kept |
| 'What this file is' | subheader |  | the check's answer |
| 'Step 4. Build a model from a GIS layer' | header |  | a skeleton model from a mains layer |
| 'GIS layer of your water mains' | file uploader |  | a GeoPackage, zipped shapefile, GeoJSON, KML or KMZ |
| 'Unit of the coordinates' | selectbox |  | asked when the layer does not say, or may be wrong |
| 'Layer' | selectbox |  | which line layer, when the file has several |
| 'Order of the coordinates' | selectbox |  | asked when latitude may be written first |
| 'Pipe size field' | selectbox |  | the layer's field with each pipe's size (guessed) |
| 'Material field' | selectbox |  | the layer's field with each pipe's material (guessed) |
| 'Install year field' | selectbox |  | the layer's field with each pipe's install year (read, not used) |
| 'Declared length field' | selectbox |  | the layer's field with each pipe's length (compared with the drawn length) |
| 'Pipe ID field' | selectbox |  | the layer's field with each pipe's ID |
| 'Join line ends closer than (ft)' | number input | 0.1 to 200; default 6 ft for a GIS layer, 10 m for a hand-drawn KML | the joining distance |
| 'Unit of the pipe sizes' | selectbox |  | asked when sizes are bare numbers |
| 'Size of every pipe (in)' | number input | 1 to 120; default empty until given | asked when no size can be read; empty until given |
| 'Is the layer drawn to scale?' | radio |  | asked for a layer whose unit was chosen |
| 'Pipes whose true length you know (optional): one per line, 'pipe ID, length in ft'' | text area |  | known lengths, to scale a drawing that is not to scale |
| 'Download the model (.inp)' | download button |  | the skeleton model |
| 'Open it on the main page' | button |  | keeps the skeleton in this session and opens the main page |
| 'Step 5. Build a model from a CAD drawing' | header |  | a skeleton model from a DXF drawing |
| 'CAD drawing of your water mains (DXF)' | file uploader |  | a DXF, text or binary |
| 'Layer styles to leave out' | multiselect |  | lines and texts that are not mains in service |
| 'Read the lines inside blocks inserted on frozen layers' | checkbox |  | blocks AutoCAD does not draw |
| 'Read pipe sizes from layer names' | checkbox |  | a size in a layer's name sizes every line on it |
| 'Unit of the drawing' | selectbox |  | the drawing's unit, as it says or chosen |
| 'A label names the nearest line within (text heights)' | number input | 0.5 to 20; default 3 | a label's reach |
| 'Unit of the sizes in layer names' | selectbox |  | asked when layer names carry bare sizes |
| 'Unit of the sizes marked with the diameter symbol' | selectbox |  | asked when labels carry bare sizes |
| 'Numbers in layer names' | selectbox |  | asked: are bare numbers in layer names sizes, and in what unit |
| 'Join line ends closer than (ft)' | number input | 0.1 to 200; default 6 | the joining distance |
| 'Size of every pipe (in)' | number input | 1 to 120; default empty until given | asked when no size can be read; empty until given |
| 'Is the drawing to scale?' | radio |  | asked for a drawing whose unit was chosen |

## 8. Command lines

These are the commands that take a utility's files. The experiments, fixtures, saved checks and numbers scripts are the project's own tools; the README lists them.

`python -m residualmap.pilot`, the real-data path (results for real data go to a git-ignored private folder):

| argument | what it is |
|---|---|
| `--inp` | the utility's model file, any form the importer takes (a bundled network's name runs it by name) |
| `--samples` | the grab log, CSV or Excel |
| `--taps` | the tap map, CSV |
| `--dose` | free chlorine leaving the plant, mg/L (default 1.2); total chlorine with chloramine (default 2.0); with inlet residuals, the highest residual entering the zone |
| `--threshold` | the minimum residual, mg/L: default 0.2 free chlorine, 0.5 total chlorine |
| `--disinfectant` | free_chlorine (the default) or chloramine |
| `--ph` | chloramine only, with --cl2n: the plant water's logged pH |
| `--cl2n` | chloramine only, with --ph: the chlorine to ammonia-N ratio by mass |
| `--window` | months of history to fit on (default 3) |
| `--synthetic` | a bundled network's name: run on a synthetic 6-month log instead of real data |
| `--plant` | the plant log, CSV (free chlorine only): makes the model temperature-aware |
| `--out` | where results go: default the git-ignored private folder for real data, the committed pilot folder for --synthetic |
| `--allow-large` | build the grid of a model with more than 3,000 junctions |
| `--date-format` | how the log's dates are written when they read two ways: MDY, DMY, YMD or a strptime format |
| `--analyzer-hour` | the clock hour an online analyzer's reading is kept nearest (default 10, for 10:00). A tap with 12 or more readings a day at a median step of 30 minutes or less on 2 or more days, or with 12 or more readings a day on 7 or more days at any step, is read as an online analyzer: it keeps one reading a day, the one nearest 10:00 when one lies within 60 minutes of it, and a day with none keeps no reading and is listed. |
| `--keep-every-reading` | keep every reading of a tap that looks like an online analyzer |
| `--columns` | the grab log's columns when the guess is wrong: 'tap=...;reading=...;date=...;time=...', also analyte=, unit= and analyte_value= |
| `--inlets` | the inlet residuals CSV, for a zone whose inlets carry different residuals (free chlorine only) |
| `--holdout-sites` | also leave this many sites out of each held-out month's fit, and predict them |

`python -m residualmap.importers`, what each file is, by its content:

| argument | what it is |
|---|---|
| `files` | one or more files of any kind: each is identified by its content |
| `--prepare` | prepare each EPANET model and print its readiness report |
| `--out` | where a fixed copy goes (default the git-ignored private inputs folder) |
| `--disinfectant` | free_chlorine (the default) or chloramine, for the readiness report |
| `--json` | print JSON instead of plain text |

`python -m residualmap.gis_read`, a skeleton model from a GIS layer:

| argument | what it is |
|---|---|
| `layer` | a GIS layer of the water mains (GeoPackage, zipped shapefile, GeoJSON, KML or KMZ) |
| `--unit` | the coordinates' unit when the layer does not say |
| `--tolerance-ft` | join line ends closer than this, in feet (default 6 ft; 10 m for a KML) |
| `--out` | where the skeleton model (.inp) is written |
| `--form` | an operator form: complete the skeleton with it and write the completed model to --out |

`python -m residualmap.cad_read`, a skeleton model from a DXF drawing:

| argument | what it is |
|---|---|
| `drawing` | a DXF drawing of the water mains (text or binary, or a .bak holding one) |
| `--unit` | the drawing's unit when it does not say |
| `--tolerance-ft` | join line ends closer than this, in feet (default 6 ft) |
| `--out` | where the skeleton model (.inp) is written |
| `--form` | an operator form: complete the skeleton with it and write the completed model to --out |
| `--layer-numbers` | numbers in layer names with no unit beside them: none (not sizes), in, mm or m |
| `--default-in` | the size of pipes no layer name or label sizes, in inches |

`python -m residualmap.operator_form`, the operator form on an EPANET model (built and tested, not yet in the app):

| argument | what it is |
|---|---|
| `model` | an EPANET model whose gaps the form fills |
| `form` | the operator form, CSV or Excel |
| `--out` | where the completed model is written |

## Outputs

**The readiness report**, on the main page under the uploaded model, on the Get started page's 'Check a file', and from the pilot as `import_report_{tag}.json` beside its results. It lists what came from the file (counts of junctions, pipes, tanks, pumps, controls and patterns, the units and the start time), what was derived, what is assumed, what is still needed and why, the fixes made in a copy, a test run of the first hour, and the run cost (the grid's runs, minutes and disk).

**The prepared model**, a copy named by its content, in a git-ignored folder; the original file is never changed.

**The main page**, titled 'ResidualMap', shows these headline numbers: 'Samples used', 'Junctions likely below N mg/L (N)', 'Worst hour of the day', 'Calibrated decay (bulk / wall)', 'Calibrated decay' (the bulk rate per day, the wall rate in m/day, the old-pipe factor, the demand and dose multipliers). Its maps are the median in mg/L, the 90% band's width in mg/L, the probability of being below the threshold and the route's sites. Its sections are 'Grab samples', 'Sample here next', 'Every junction: the daily minimum', 'Nitrification watch (chloramine)', 'Water age', 'Where chlorine is lost', 'When is the network worst?', 'Demo only: reveal the hidden truth', 'What this model accounts for'. 'Sample here next' has the columns 'junction', 'hour', 'P(daily min < N)', 'water age (h)', 'reason'; 'Every junction: the daily minimum' has 'junction', 'median (mg/L)', 'lo90 (mg/L)', 'hi90 (mg/L)', 'P(daily min < N)', 'hour of the daily minimum', with a CSV download from the table's toolbar, and the junctions most likely below the threshold again as a plain table. The one-page PDF report is named after the model file, `residualmap_Net3.inp.pdf` for Net3. In the chloramine mode a nitrification watch lists junctions, and the water-age panels give each junction's age in hours and its range over the grid's settings.

**The pilot's files** (`predictions_{tag}.csv`, `validation_{tag}.csv`, `validation_{tag}.json`, `import_report_{tag}.json`, `synthetic_log_{network}.csv`): `predictions_{tag}.csv` has one row per held-out sample with the columns `junction`, `hour`, `y`, `month`, `date`, `tap_id`, `pred`, `z_mu`, `z_sd`, `p_below`, `train_months`, `in50`, `in80`, `in90`, `in95`, `persistence`, `network_mean`, `seen_tap` (`pred` in mg/L, `z_mu` and `z_sd` in ln mg/L); `validation_{tag}.csv` one row per held-out month and tap subset with `held_out_month`, `taps`, `trained_on`, `n_train`, `n_test`, `rmse`, `mae`, `rmse_persistence`, `rmse_network_mean`, `coverage50`, `coverage80`, `coverage90`, `coverage95`, `band90_width_mgL`, `n_true_below`, `recall_below`, `false_alarms`, `map_kb`, `map_kw`, `map_gamma`, `map_dose` (RMSE and MAE in mg/L); `validation_{tag}.json` holds the same and the means; `synthetic_log_{network}.csv` (with `--synthetic`) has `junction`, `hour`, `y`, `month`, `date`, `tap_id`. The chloramine mode, a plant log, inlet residuals and `--holdout-sites` add columns of their own and their own suffix to the tag.

**Files it writes**, all git-ignored:

| where | what |
|---|---|
| `outputs/app_uploads/` | the main page's prepared copies of uploaded models (git-ignored) |
| `outputs/private/` | real data: the pilot's prepared models and results, the page's checks (git-ignored) |
| `outputs/cache/` | the simulator grids, one per model and dose (git-ignored) |
| temp.inp, temp.rpt, temp.bin | EPANET's scratch files in the working directory, locally (git-ignored) |

## The hosted demo (cloud mode)

The hosted demo runs the same two pages in cloud mode: the example networks only (Net3 and Net2), no file uploader on either page, and a privacy notice at the top of both. Its uploads are off by default. If they are turned on, each file is held to the hosted limits below and read in memory and in temporary files deleted when each run of the page ends; what is made from it is kept in the visitor's session only, and nothing is written under the outputs folder. Real data runs on the laptop (`docs/privacy.md`).

## Limits

| limit | value | what it limits |
|---|---|---|
| large model | 3,000 | junctions above which a grid build needs a confirm (the app's checkbox, the pilot's flag) |
| date window | 3 | the main page fits the last this many calendar months that hold a reading |
| analyzer readings a day | 12 | a day of a tap with this many readings or more can be an online analyzer's |
| analyzer step | 30 | minutes: such a day whose readings are at most this far apart (the median step) is an analyzer day |
| analyzer days | 2 | analyzer days that make a tap an online analyzer |
| analyzer dense days | 7 | days with that many readings, at any step, that make a tap an online analyzer too |
| analyzer hour | 10 | an online analyzer keeps the reading nearest this clock hour |
| analyzer reach | 60 | minutes: a day keeps its reading nearest the analyzer hour only when it lies within this; a day with none keeps none |
| highest reading | 10 | a reading above this many mg/L is listed, not used |
| hosted upload | 5 | the hosted demo's largest file, in MB, if its uploads are turned on |
| hosted model | 300 | the hosted demo's largest model, in junctions |

## What is not built, or not tested

- The importers are tested only as converter self-consistency on public models written out by our own writers, never on a utility's own files; no external layer was used.
- Native vendor databases are not read: send the EPANET export.
- Censored readings ('<0.02', 'ND') are shown, not fitted.
- A tap map by coordinates and the operator form are built and tested, but not offered in the app: each missed its pre-registered bar in simulation.
- The hosted demo shows the example networks only; real data runs on the laptop.
- Windows is untested.
