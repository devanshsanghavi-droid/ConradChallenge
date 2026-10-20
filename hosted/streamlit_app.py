"""
hosted/streamlit_app.py: the hosted demo's start on Streamlit Community Cloud (iteration 5, task 20).

    streamlit run hosted/streamlit_app.py     # from the repository's root, as Community Cloud runs it

It is always cloud mode (residualmap/cloud.py), whatever the app's secrets say: example networks only, no uploader, the
privacy notice.  It runs the same two pages as 'streamlit run app.py' (app.py and pages/1_Get_started.py, named
'ResidualMap' and 'Get started').  It lives in its own folder for its settings: Streamlit 1.64.0's loader reads
hosted/.streamlit/config.toml after the repository's .streamlit/config.toml when this file is the main script, so on
this Mac the hosted server listens on every network interface (server.address 0.0.0.0) and takes uploads of up to
5 MB, while 'streamlit run app.py' still listens on localhost only.  Community Cloud's documentation says that an app
whose entrypoint is in a subdirectory takes its settings from the repository's root, so whether the hosted file
applies there is not verified; if it does not, the hosted app runs with the root's localhost setting, which the
platform's own health check on localhost:8501 suggests may still be reachable.  Only the smoke test after the deploy
settles it (docs/deploy_streamlit_cloud.md).
"""
import os
import sys

os.environ["RESIDUALMAP_MODE"] = "cloud"          # this start is only for the hosted demo: never local mode
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:      # Streamlit puts this file's folder on the path, not the repository's root
    sys.path.insert(0, ROOT)

from residualmap.cloud import run_pages  # noqa: E402

run_pages(os.path.join(ROOT, "app.py"), os.path.join(ROOT, "pages", "1_Get_started.py"))
