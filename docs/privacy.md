# Privacy: what runs where

ResidualMap can run in two places: on the laptop (the local version, `streamlit run app.py` and the pilot scripts) or as the hosted demo on Streamlit Community Cloud (`hosted/streamlit_app.py`, see `docs/deploy_streamlit_cloud.md`). A utility's real data runs only on the laptop, unless the utility agrees in writing to something else.

## How a pilot works

Devansh runs ResidualMap on the laptop, on the files a utility emails or shares (the model's EPANET export, the grab log, the tap map), and sends back the monthly PDF and tables. Utility staff use the hosted demo and its Get started page to see the tool and to prepare their files, without sending anything to it.

A utility may want a data-use agreement before it shares anything. Ask before the first file is sent, and keep the answer with the pilot's notes, outside the repository.

## The local version

- **Everything runs on the laptop.** The app's server listens on this computer only (`server.address = "localhost"` in `.streamlit/config.toml`), so no one else on the same Wi-Fi can open it, and Streamlit's usage statistics are off (`browser.gatherUsageStats = false`).
- **Where files go.** An uploaded model is copied, by its content's SHA-1, to `outputs/app_uploads/`, and its simulation grid to `outputs/cache/`; the pilot's copies and results go to `outputs/private/`. All three are git-ignored, so they are never committed or pushed. 'Check a file' on the Get started page works in a temporary folder under `outputs/private/` that is deleted at once.
- **The disk is encrypted.** FileVault is turned on for the laptop that runs the pilot (checked with `fdesetup status` on October 19, 2026), so a lost laptop does not expose the files.
- **Nothing is sent anywhere.** ResidualMap's own code makes no network calls. The repository is public on GitHub, so nothing private may be committed: the names check (`RESIDUALMAP_NAMES_CSV`, see `docs/deploy_streamlit_cloud.md`) runs before every push.
- **Windows** is untested; the local version has been run on macOS only.

## The hosted demo

- **Who runs it.** Streamlit Community Cloud is operated by Snowflake, under Snowflake's privacy policy and Streamlit's terms. Its servers, and the logs they keep, are theirs, not ResidualMap's.
- **What it shows.** By default, the example networks only (EPANET's Net3 and Net2): no file uploader on either page, and a notice at the top of both pages that says so and asks visitors not to use it for a utility's real data unless the utility agrees.
- **What a visitor sends.** The settings they choose and any samples they type into the table go to the server, as with any web app. The typed samples are kept in the visitor's session only. Results computed from the example networks at the settings a visitor picks (the dose, the demo scenario, and the decay rates calibrated on their samples) can stay in the server's shared caches until it restarts, and a dose nobody picked before gets its own grid file on the server's disk; none of them holds anything that says who the visitor is. The example networks' grids are kept on the server's disk until it restarts.
- **Usage statistics.** Community Cloud's documentation lists `browser.gatherUsageStats` among the settings it sets for every app it runs, so the hosted demo's browser sends Streamlit's usage statistics whatever ResidualMap's settings files say.

## If uploads are turned on in the hosted demo

Uploads are off unless `RESIDUALMAP_ALLOW_UPLOADS = "true"` is added to the app's secrets (open decision 8). With it:

- **The file leaves the visitor's computer.** The notice says so.
- **Limits.** Files up to 5 MB (the model, the grab log and the tap map). Every uploader shows the 5 MB limit, so the visitor's browser refuses a larger file; if one still arrives, it is refused once it reaches the server, before ResidualMap reads it. A zip whose model is larger than 5 MB once unpacked is refused by the size its directory declares, before it is unpacked, and an Excel grab log or tap map whose parts add up to more than 20 MB once unpacked, or any of whose parts is packed more than 100 to 1, is refused before it is opened.
- **What a model may cost.** A small file can still be an expensive model, so a model is held to more than its size. Before it is parsed, by counting lines in its text: at most 300 junctions, 600 links (pipes, pumps and valves), 300 tanks and reservoirs, and 20,000 lines of model data (its title and drawing aside). Before EPANET runs it: demand patterns that change no more often than every 15 minutes, and rules checked no more often than every minute. After one timed 7-day simulation: at most 3 s for that run, and at most 10 minutes for the grid it gives on the server's 2 workers. The whole check runs in a process of its own that is stopped after 60 s. A model over any of these is refused, with what it was over and a pointer to the local version. These limits bound the cost of one upload; they do not limit how many uploads arrive.
- **Where it goes.** The file is processed in the server's memory and in a temporary folder that is deleted when each run of the page ends. Every EPANET run on it runs in a process of its own whose working folder, and temporary folder, is that folder: EPANET writes a scratch file of its own to its working folder, and deletes it only when a run ends normally, so it is in that folder too, the grid's runs included. Everything made from the file (the prepared model, its simulation grid, the water age, the demo's hidden truth, the grab log's reading) is kept in that visitor's browser session only, never in the caches the server shares between visitors, and nothing is written under `outputs/`. When the session ends, it is gone. If the server itself is stopped in the middle of a run, the folder stays in the server's temporary directory, not in the app's own folder; whether the platform clears it is not documented. Saved checks (`cloud_mode_uploads_leave_no_files`, `mode_detection_fails_safe`) upload Net3 and check that no file under `outputs/` is new or changed, that nothing is left in the temporary directory, and that EPANET's scratch file appears only in the run's own folder while an upload is processed; and that the oversized and expensive files above are refused.
- **What it does not cover.** The platform's own logs and its handling of the network traffic are Snowflake's, under its policy. That is why a utility's real data stays on the laptop unless the utility agrees in writing.
