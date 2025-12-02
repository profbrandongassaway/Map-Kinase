"""start_m5.py

Python runner to start the m5 webapp and display it in a lightweight desktop window.

Usage examples:
    # Run with default data in a desktop window
    python start_m5.py

    # Run in background (headless) with optional logfile
    python start_m5.py --background --logfile m5.log

    # Provide custom files/port and open in desktop window
    python start_m5.py --protein-file "C:\path\to\prot.txt" --ptm-file "C:\path\to\ptm.txt" --port 8004
"""
from __future__ import annotations

import argparse
import os
import socket
import subprocess
import sys
import time
from typing import Iterable, Sequence

DEFAULT_PROT = r"C:\Users\clayt\OneDrive - Brigham Young University\Desktop\Graduate_Documents\data_analysis\BCp2_data_03272025\Phosmap\BCp2_ProtMaphsaanno.txt"
DEFAULT_PTM = r"C:\Users\clayt\OneDrive - Brigham Young University\Desktop\Graduate_Documents\data_analysis\BCp2_data_03272025\Phosmap\BCp2_PhosMaphsanno.txt"
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def _wait_for_server(host: str, port: int, timeout: float = 20.0) -> bool:
    """Poll until the Shiny/uvicorn server accepts connections."""
    deadline = time.time() + timeout
    while time.time() < deadline:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.settimeout(1.0)
            try:
                sock.connect((host, port))
            except OSError:
                time.sleep(0.25)
                continue
            else:
                return True
    return False


def _build_uvicorn_cmd(host: str, port: int, extra_args: str | None) -> Sequence[str]:
    cmd: list[str] = [
        sys.executable,
        '-m',
        'uvicorn',
        'm5_webapp:app.starlette_app',
        '--host',
        host,
        '--port',
        str(port),
        '--log-level',
        'info',
    ]
    if extra_args:
        cmd.extend(extra_args.split())
    return cmd


def _print_launch_summary(cmd: Iterable[str], env: dict[str, str]) -> None:
    print('Starting m5 with:')
    print('  PROTEIN_FILE=', env.get('M5_PROTEIN_FILE', ''))
    print('  PTM_FILE=', env.get('M5_PTM_FILE', ''))
    print('  CMD=', ' '.join(cmd))


def main() -> None:
    parser = argparse.ArgumentParser(
        prog='start_m5.py',
        description='Start the MapKinase m5 webapp (uvicorn) and optionally wrap it in a desktop window.'
    )
    parser.add_argument('--protein-file', '-p', default=DEFAULT_PROT, help='Path to protein file (M5_PROTEIN_FILE)')
    parser.add_argument('--ptm-file', '-t', default=DEFAULT_PTM, help='Path to PTM file (M5_PTM_FILE)')
    parser.add_argument('--port', '-P', type=int, default=8003, help='Port to bind uvicorn')
    parser.add_argument('--background', '-b', action='store_true', help='Start in background (no desktop window)')
    parser.add_argument('--logfile', '-l', help='When backgrounding, redirect stdout/stderr to this logfile')
    parser.add_argument('--host', default='127.0.0.1', help='Host to bind (default 127.0.0.1)')
    parser.add_argument('--uvicorn-args', help='Extra args to pass to uvicorn (quoted string)')
    args = parser.parse_args()

    try:
        import uvicorn  # noqa: F401  # Ensure uvicorn is available before proceeding
    except Exception:
        print('uvicorn is not installed in this Python environment. Install with: pip install uvicorn')
        sys.exit(1)

    env = os.environ.copy()
    env['M5_PROTEIN_FILE'] = args.protein_file
    env['M5_PTM_FILE'] = args.ptm_file

    cmd = list(_build_uvicorn_cmd(args.host, args.port, args.uvicorn_args))
    _print_launch_summary(cmd, env)

    if args.background:
        _launch_background(cmd, env, args.logfile)
        return

    try:
        import webview  # type: ignore
    except Exception:
        print('pywebview is not installed. Install with: pip install pywebview')
        print('Falling back to running the server in this console...')
        try:
            subprocess.run(cmd, env=env, check=True, cwd=SCRIPT_DIR)
        except KeyboardInterrupt:
            print('Interrupted by user')
        except subprocess.CalledProcessError as exc:
            print('uvicorn exited with non-zero status', exc.returncode)
        return

    server_process: subprocess.Popen[str] | None = None
    try:
        server_process = subprocess.Popen(cmd, env=env, cwd=SCRIPT_DIR)
        if not _wait_for_server(args.host, args.port):
            raise RuntimeError('Timed out waiting for uvicorn to become ready')
        window_url = f'http://{args.host}:{args.port}'
        webview.create_window('MapKinase m5', window_url, width=1280, height=800)
        webview.start()
    except KeyboardInterrupt:
        print('Interrupted by user')
    except Exception as exc:
        print('Failed to start desktop wrapper:', exc)
        if server_process and server_process.poll() is None:
            server_process.terminate()
            server_process.wait(timeout=5)
        raise SystemExit(1)
    finally:
        if server_process and server_process.poll() is None:
            server_process.terminate()
            try:
                server_process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                server_process.kill()


def _launch_background(cmd: Sequence[str], env: dict[str, str], logfile: str | None) -> None:
    try:
        creationflags = 0
        if os.name == 'nt':
            creationflags = subprocess.CREATE_NEW_CONSOLE
        if logfile:
            logf = open(logfile, 'a', buffering=1)
            proc = subprocess.Popen(cmd, env=env, stdout=logf, stderr=subprocess.STDOUT, creationflags=creationflags, cwd=SCRIPT_DIR)
            print(f'Started background process PID={proc.pid} (output -> {logfile})')
        else:
            proc = subprocess.Popen(cmd, env=env, creationflags=creationflags, cwd=SCRIPT_DIR)
            print(f'Started background process PID={proc.pid}')
    except Exception as exc:
        print('Failed to start uvicorn in background:', exc)


if __name__ == '__main__':
    main()
