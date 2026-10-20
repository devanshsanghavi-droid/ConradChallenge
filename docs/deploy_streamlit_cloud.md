# Deploying the hosted demo on Streamlit Community Cloud

This puts ResidualMap at a web address anyone can open, with no Python and no install, so utility staff can see the tool, read the Get started page and download the templates. It is a demo: it shows the example networks only, and real pilot data still runs on the laptop (see `docs/privacy.md`).

Signing in, accepting Streamlit's terms, giving Streamlit access to the GitHub repository and clicking Deploy are account and consent steps, so Devansh does them himself. Nothing in the repository signs in anywhere.

## What gets deployed

- The main file is `hosted/streamlit_app.py`, not `app.py`. It runs the same two pages as `streamlit run app.py` (the main page, named ResidualMap in the sidebar, and Get started), always in cloud mode: example networks only (Net3 and Net2), no file uploader on either page, and a privacy notice at the top of both pages.
- It has its own settings file, `hosted/.streamlit/config.toml`: `server.address = "0.0.0.0"` (every network interface of the container), a 5 MB upload limit, and usage statistics off. Why: the repository's own setting (`localhost`, task 15) keeps the laptop's app off the Wi-Fi, and Community Cloud's documentation does not say whether visitors' traffic reaches an app that listens on localhost only. A forum report from a Community Cloud app shows the platform checking the app's health at `http://localhost:8501` and failing while the project's config.toml moved the port to 8502, so the platform does read a project's server settings.
- Whether that file applies on Community Cloud is **not verified**. Streamlit 1.64.0's own loader, on the laptop, reads the repository's `.streamlit/config.toml` first and the hosted file after it when `hosted/streamlit_app.py` is the main file (the saved check `hosted_start_settings`). But Community Cloud's documentation on file organization says that when the main file is in a subdirectory, the configuration file must stay at the repository's root. If the platform reads only the root file, the hosted app runs with `server.address = "localhost"` and Streamlit's default upload limit. The platform's health check on `localhost:8501` suggests the app may still be reachable then, and the uploaders carry the 5 MB limit themselves, but only the smoke test settles it (step 4 of "Checking it").
- Packages: `requirements.txt`, every line pinned to the version on the Mac, so the cloud build installs what was tested. WNTR has a ready-built Linux wheel, so nothing needs a compiler.

## Before deploying

1. Decide whether the GitHub repository stays public. It is public today, so every committed file, `outputs/` included, can already be read by anyone. Run the names check before every push, not only before deploying:

   ```bash
   RESIDUALMAP_NAMES_CSV=/path/to/the/private/names.csv python -m residualmap.checks
   ```

   The result is in `outputs/chem/checks_run.json` (git-ignored); the run fails if any listed name is found.

2. Try the hosted demo's behaviour on the laptop first, if you like:

   ```bash
   streamlit run hosted/streamlit_app.py --server.address localhost
   ```

   (The flag keeps this local test on localhost; on the laptop, the hosted settings file would otherwise open it to the network.)

## Deploying

1. Open https://share.streamlit.io and sign in with the GitHub account that owns the repository. Accept Streamlit's terms, and allow Streamlit to read the repository when GitHub asks (an OAuth grant you can revoke later in GitHub's settings, under Applications).
2. Click **Create app**, then choose to deploy from a GitHub repository.
3. Fill in:
   - Repository: `devanshsanghavi-droid/ConradChallenge`
   - Branch: `main`
   - Main file path: `hosted/streamlit_app.py`
4. Open **Advanced settings**:
   - Python version: 3.13 (the version the Mac runs).
   - Secrets: paste the line below. The hosted start sets cloud mode by itself; the secret is a second guard, and the app also treats a checkout under `/mount/src/` as cloud mode when the variable is missing.

     ```toml
     RESIDUALMAP_MODE = "cloud"
     ```

5. Choose the app's address (the part before `.streamlit.app`) and click **Deploy**. The first build installs the packages and takes a few minutes; its log is in the panel at the bottom right of the app.
6. Choose who can open it: in the app's **Share** settings, either public, or only people you invite by email.

## Checking it

1. Open the app's address. The page should show the privacy notice under the title ("This is a public demo on Streamlit Community Cloud, run by Snowflake under its privacy policy. It shows example networks only. ..."), the Model choice with only "Example network", and Net3 and Net2 in the example list.
2. The first visit builds Net3's 675-run grid on the server, so the spinner can stay for minutes. Later visits read it from the server's disk until the app restarts.
3. Click **Get started** in the sidebar: the export guides and the template buttons should be there, and step 3 should say that checking a file is turned off in the public demo.
4. Open the app's log (**Manage app**, bottom right) and find the lines Streamlit prints when it starts: "You can now view your Streamlit app in your browser" and the URL or URLs after it. In Streamlit 1.64.0 a single `URL: http://localhost:8501` means the server listens on localhost only (the hosted settings file did not apply); a `Local URL` followed by a `Network URL` means it listens on every interface (it did). A single URL with any other address means the platform sets the browser's address itself, and the log does not tell. Copy those lines for Claude either way.
5. If the app never loads and the log says the health check failed or the app is unreachable, do not change `.streamlit/config.toml` (that would open the laptop's app to the network). Copy the log's last lines to Claude.
6. Send Claude the address. Claude will make one read-only visit and record it in `outputs/intake/cloud_smoke.json`: the address, the Net3 cold-start and warm-load times, that the Get started page renders, that the notice is visible, and the address the server reported in its log (from the lines you copied in step 4). Claude does not sign in, change settings or upload anything.

## Later decisions (Devansh's)

- **Cold start too slow** (open decision 10): commit the Net3 and Net2 example grids (about 26 MB) under a cloud-only folder with a SHA-256 manifest, so the server reads them instead of building them. Not done now; decide after the smoke test's timing.
- **Uploads** (open decision 8, off by default): only after reading `docs/privacy.md`, add a second secret line:

  ```toml
  RESIDUALMAP_ALLOW_UPLOADS = "true"
  ```

  Write it with the quotes: Streamlit passes text values on to the app, but not a bare `true`. With it, the Model choice offers "Upload my .inp" again, with a 5 MB file limit and a 300-junction model limit, and the notice says that an uploaded file leaves the visitor's computer. A model is also held to what it costs to run (links, tanks, lines of model data, its time steps, and a grid of at most 10 minutes by its own timed run on the server, checked in a process stopped after 60 s), and an Excel grab log to 20 MB once unpacked; `docs/privacy.md` lists them. Remove the line to turn uploads off again.

## What does not work in the cloud, and what the app does instead

- **EPA's chloramine chemistry (EPANET-MSX)**: without a C compiler it runs uncompiled and takes minutes per run (compiled, the app's chloramine demo truth takes about 15 s on the Mac), and whether the default image has a compiler and the OpenMP runtime EPANET-MSX needs is untested. So the chloramine demo goes straight to its reduced truth (EPANET first-order decay at the rate EPA's chemistry gives for the water), with the warning that this flatters the model. The chloramine map itself runs.
- **Speed**: the server is a small shared machine, so grids take longer than on the Mac. The app uses 2 worker processes there, so the container keeps its memory. ky4 is not offered.
- **Cache**: the server's disk is not kept across restarts, so the example grids are built again after each restart (decision 10 above).
- **Sleeping**: an app with no visitors for a while is put to sleep; the next visitor wakes it, which takes longer than a normal load.
- **Usage statistics**: Community Cloud's documentation lists `browser.gatherUsageStats` among the settings it sets for every app, so the hosted demo's browser statistics are on, whatever the settings files say. The local version keeps them off.
- **Platform**: Linux, not macOS. The macOS-only OpenMP preload is not used, and WNTR's Linux EPANET libraries run.
- **Updates**: Community Cloud updates the app when the branch it runs changes, so every push to `main` reaches the public demo. Run the saved checks and the names check before pushing.
