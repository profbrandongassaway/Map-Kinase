import threading
import time

from MapKinase_WebApp.m5_main_ui import (
    GUI_POPUP,
    HOST,
    PORT,
    _launch_desktop_window,
    app as mapkinase_app,
    uvicorn,
    webview,
)


app = mapkinase_app


def _run_uvicorn_app():
    if uvicorn is None:  # pragma: no cover - runtime guard
        raise RuntimeError("uvicorn is required to run the Shiny app server")
    uvicorn.run(app, host=HOST, port=PORT, log_level="info")


if __name__ == "__main__":
    if GUI_POPUP:
        server_thread = threading.Thread(target=_run_uvicorn_app, daemon=True)
        server_thread.start()
        time.sleep(1.2)
        app_url = f"http://{HOST}:{PORT}"
        try:
            _launch_desktop_window(app_url)
        except KeyboardInterrupt:
            pass
        finally:
            if webview is None:
                server_thread.join()
    else:
        _run_uvicorn_app()
