import hashlib
import hmac
import html
import os
import secrets
import time
from typing import Dict, List, Optional, Tuple
from urllib.parse import parse_qs, quote

from MapKinase_WebApp.mk_runtime_env import env_truthy


LOGIN_PATH = "/__login__"
LOGOUT_PATH = "/__logout__"
COOKIE_NAME = "mapkinase_auth"
COOKIE_MAX_AGE_SECONDS = int(os.environ.get("MAPKINASE_AUTH_MAX_AGE_SECONDS", "43200"))
COOKIE_SECURE = str(os.environ.get("MAPKINASE_AUTH_COOKIE_SECURE", "")).strip().lower() in {
    "1",
    "true",
    "yes",
    "on",
}


def _header_lookup(scope) -> Dict[str, str]:
    headers: Dict[str, str] = {}
    for key, value in scope.get("headers", []):
        headers[key.decode("latin-1").lower()] = value.decode("latin-1")
    return headers


def _parse_cookie_header(raw_cookie: str) -> Dict[str, str]:
    cookies: Dict[str, str] = {}
    for part in raw_cookie.split(";"):
        if "=" not in part:
            continue
        name, value = part.split("=", 1)
        cookies[name.strip()] = value.strip()
    return cookies


def _build_html_page(title: str, body: str) -> bytes:
    page = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(title)}</title>
  <style>
    :root {{
      color-scheme: light;
      font-family: "Segoe UI", Tahoma, sans-serif;
      background: #f4f6fb;
      color: #1b2533;
    }}
    body {{
      margin: 0;
      min-height: 100vh;
      display: grid;
      place-items: center;
      background:
        radial-gradient(circle at top, rgba(49, 110, 188, 0.12), transparent 36%),
        linear-gradient(180deg, #f7f9fd 0%, #eef2f9 100%);
    }}
    .card {{
      width: min(420px, calc(100vw - 32px));
      background: rgba(255, 255, 255, 0.96);
      border: 1px solid #d7deea;
      border-radius: 16px;
      box-shadow: 0 22px 60px rgba(15, 23, 42, 0.14);
      padding: 28px;
    }}
    h1 {{
      margin: 0 0 8px;
      font-size: 1.6rem;
    }}
    p {{
      margin: 0 0 18px;
      line-height: 1.45;
      color: #4a5565;
    }}
    label {{
      display: block;
      font-size: 0.92rem;
      font-weight: 600;
      margin-bottom: 6px;
    }}
    input {{
      width: 100%;
      box-sizing: border-box;
      padding: 12px 13px;
      margin: 0 0 14px;
      border: 1px solid #c8d2e2;
      border-radius: 10px;
      font-size: 1rem;
      background: #fbfcfe;
    }}
    input:focus {{
      outline: 2px solid rgba(30, 95, 171, 0.2);
      border-color: #1e5fab;
      background: #fff;
    }}
    button {{
      width: 100%;
      border: 0;
      border-radius: 10px;
      background: linear-gradient(180deg, #1e5fab 0%, #174d8d 100%);
      color: #fff;
      padding: 12px 14px;
      font-size: 1rem;
      font-weight: 700;
      cursor: pointer;
    }}
    .error {{
      margin: 0 0 14px;
      padding: 10px 12px;
      border-radius: 10px;
      background: #fff1f2;
      border: 1px solid #fecdd3;
      color: #9f1239;
      font-size: 0.95rem;
    }}
    .note {{
      margin-top: 14px;
      font-size: 0.86rem;
      color: #64748b;
    }}
  </style>
</head>
<body>
  {body}
</body>
</html>
"""
    return page.encode("utf-8")


class LoginGateASGI:
    def __init__(self, inner_app):
        self.inner_app = inner_app
        self.username = str(os.environ.get("MAPKINASE_LOGIN_USERNAME", "")).strip()
        self.password = str(os.environ.get("MAPKINASE_LOGIN_PASSWORD", ""))
        if not self.username or not self.password:
            raise RuntimeError(
                "MAPKINASE_LOGIN_USERNAME and MAPKINASE_LOGIN_PASSWORD must both be set before starting the protected app."
            )
        self.cookie_secret = os.environ.get("MAPKINASE_AUTH_SECRET") or secrets.token_urlsafe(32)

    async def __call__(self, scope, receive, send):
        scope_type = scope.get("type")
        if scope_type == "lifespan":
            await self.inner_app(scope, receive, send)
            return
        if scope_type == "websocket":
            if self._is_authenticated(scope):
                await self.inner_app(scope, receive, send)
            else:
                await send({"type": "websocket.close", "code": 1008})
            return
        if scope_type != "http":
            await self.inner_app(scope, receive, send)
            return

        path = scope.get("path", "") or "/"
        if path == LOGIN_PATH:
            await self._handle_login(scope, receive, send)
            return
        if path == LOGOUT_PATH:
            await self._handle_logout(send)
            return
        if self._is_authenticated(scope):
            await self.inner_app(scope, receive, send)
            return
        await self._redirect_to_login(scope, send)

    def _is_authenticated(self, scope) -> bool:
        headers = _header_lookup(scope)
        cookie_header = headers.get("cookie", "")
        token = _parse_cookie_header(cookie_header).get(COOKIE_NAME, "")
        return self._validate_token(token)

    def _create_token(self) -> str:
        issued_at = str(int(time.time()))
        payload = f"{self.username}|{issued_at}"
        signature = hmac.new(
            self.cookie_secret.encode("utf-8"),
            payload.encode("utf-8"),
            hashlib.sha256,
        ).hexdigest()
        return f"{payload}|{signature}"

    def _validate_token(self, token: str) -> bool:
        try:
            username, issued_at, signature = token.split("|", 2)
        except ValueError:
            return False
        if not secrets.compare_digest(username, self.username):
            return False
        expected_payload = f"{username}|{issued_at}"
        expected_signature = hmac.new(
            self.cookie_secret.encode("utf-8"),
            expected_payload.encode("utf-8"),
            hashlib.sha256,
        ).hexdigest()
        if not secrets.compare_digest(signature, expected_signature):
            return False
        try:
            issued_at_int = int(issued_at)
        except ValueError:
            return False
        return (int(time.time()) - issued_at_int) <= COOKIE_MAX_AGE_SECONDS

    def _sanitize_next(self, next_path: Optional[str]) -> str:
        candidate = str(next_path or "/").strip()
        if not candidate.startswith("/") or candidate.startswith("//"):
            return "/"
        if candidate in {LOGIN_PATH, LOGOUT_PATH}:
            return "/"
        return candidate

    def _build_set_cookie(self, token: str, clear: bool = False) -> str:
        parts = [
            f"{COOKIE_NAME}={token}",
            "HttpOnly",
            "Path=/",
            "SameSite=Lax",
        ]
        if clear:
            parts.append("Max-Age=0")
            parts.append("Expires=Thu, 01 Jan 1970 00:00:00 GMT")
        else:
            parts.append(f"Max-Age={COOKIE_MAX_AGE_SECONDS}")
        if COOKIE_SECURE:
            parts.append("Secure")
        return "; ".join(parts)

    async def _read_body(self, receive) -> bytes:
        body = bytearray()
        while True:
            message = await receive()
            if message["type"] != "http.request":
                continue
            body.extend(message.get("body", b""))
            if not message.get("more_body", False):
                break
        return bytes(body)

    async def _handle_login(self, scope, receive, send):
        if scope.get("method", "GET").upper() == "POST":
            form_data = parse_qs((await self._read_body(receive)).decode("utf-8", errors="ignore"))
            username = str(form_data.get("username", [""])[0])
            password = str(form_data.get("password", [""])[0])
            next_path = self._sanitize_next(form_data.get("next", ["/"])[0])
            if secrets.compare_digest(username, self.username) and secrets.compare_digest(password, self.password):
                headers = [
                    (b"location", next_path.encode("utf-8")),
                    (b"cache-control", b"no-store"),
                    (b"set-cookie", self._build_set_cookie(self._create_token()).encode("utf-8")),
                ]
                await send({"type": "http.response.start", "status": 303, "headers": headers})
                await send({"type": "http.response.body", "body": b""})
                return
            await self._send_login_page(send, next_path, error_message="Invalid username or password.", status=401)
            return

        if self._is_authenticated(scope):
            headers = [(b"location", b"/"), (b"cache-control", b"no-store")]
            await send({"type": "http.response.start", "status": 303, "headers": headers})
            await send({"type": "http.response.body", "body": b""})
            return

        query_string = scope.get("query_string", b"").decode("utf-8", errors="ignore")
        next_path = self._sanitize_next(parse_qs(query_string).get("next", ["/"])[0])
        await self._send_login_page(send, next_path, status=200)

    async def _handle_logout(self, send):
        headers = [
            (b"location", LOGIN_PATH.encode("utf-8")),
            (b"cache-control", b"no-store"),
            (b"set-cookie", self._build_set_cookie("", clear=True).encode("utf-8")),
        ]
        await send({"type": "http.response.start", "status": 303, "headers": headers})
        await send({"type": "http.response.body", "body": b""})

    async def _redirect_to_login(self, scope, send):
        path = scope.get("path", "") or "/"
        query_string = scope.get("query_string", b"")
        if query_string:
            path = f"{path}?{query_string.decode('utf-8', errors='ignore')}"
        login_target = f"{LOGIN_PATH}?next={quote(path, safe='/?:=&')}"
        headers = [(b"location", login_target.encode("utf-8")), (b"cache-control", b"no-store")]
        await send({"type": "http.response.start", "status": 303, "headers": headers})
        await send({"type": "http.response.body", "body": b""})

    async def _send_login_page(self, send, next_path: str, error_message: str = "", status: int = 200):
        error_block = ""
        if error_message:
            error_block = f'<div class="error">{html.escape(error_message)}</div>'
        body = f"""
<main class="card">
  <h1>MapKinase Login</h1>
  <p>Enter the shared collaborator username and password to open the site.</p>
  {error_block}
  <form method="post" action="{LOGIN_PATH}">
    <input type="hidden" name="next" value="{html.escape(next_path, quote=True)}">
    <label for="username">Username</label>
    <input id="username" name="username" type="text" autocomplete="username" required>
    <label for="password">Password</label>
    <input id="password" name="password" type="password" autocomplete="current-password" required>
    <button type="submit">Log In</button>
  </form>
  <p class="note">This temporary gate protects both page loads and the app's websocket connection.</p>
</main>
"""
        payload = _build_html_page("MapKinase Login", body)
        headers: List[Tuple[bytes, bytes]] = [
            (b"content-type", b"text/html; charset=utf-8"),
            (b"content-length", str(len(payload)).encode("ascii")),
            (b"cache-control", b"no-store"),
        ]
        await send({"type": "http.response.start", "status": status, "headers": headers})
        await send({"type": "http.response.body", "body": payload})


def login_enabled() -> bool:
    return env_truthy("MAPKINASE_ENABLE_LOGIN", False)


def maybe_wrap_with_login(inner_app):
    if not login_enabled():
        return inner_app
    return LoginGateASGI(inner_app)
