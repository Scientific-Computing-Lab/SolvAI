#!/usr/bin/env python3
"""Keep the first-party Penpot MCP plugin attached to the SolvAI design file.

The script reads the ignored local credentials file, signs into the loopback-only
Penpot deployment, opens the selected Figure 1 page, and activates Penpot's MCP
plugin. It intentionally prints no credentials or MCP URLs.
"""

from __future__ import annotations

import json
import signal
import subprocess
import time
from hashlib import sha256
from pathlib import Path
from urllib.parse import parse_qsl, urlencode, urlsplit

from selenium import webdriver
from selenium.webdriver.chrome.options import Options

CREDENTIALS = Path(__file__).with_name("credentials.env")
MCP_ENV = Path(__file__).with_name("mcp.env")
DESIGN_URL = (
    "http://localhost:9001/#/workspace?"
    "team-id=ecbaaf6b-3d61-80cb-8008-8f5a0136d0b6&"
    "file-id=ecbaaf6b-3d61-80cb-8008-8f5daf5a7d27&"
    "page-id=cd43df32-78a6-8014-8008-8f618dfd272e"
)


def read_env(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        if not raw or raw.lstrip().startswith("#") or "=" not in raw:
            continue
        key, value = raw.split("=", 1)
        values[key] = value.strip().strip('"').strip("'")
    return values


def evaluate(driver: webdriver.Chrome, expression: str):
    result = driver.execute_cdp_cmd(
        "Runtime.evaluate", {"expression": expression, "returnByValue": True}
    )
    return result["result"].get("value")


def set_input(driver: webdriver.Chrome, selector: str, value: str) -> None:
    expression = (
        "(() => {\n"
        f"  const input = document.querySelector({json.dumps(selector)});\n"
        "  if (!input) return false;\n"
        "  const setter = Object.getOwnPropertyDescriptor(\n"
        "    HTMLInputElement.prototype, 'value'\n"
        "  ).set;\n"
        f"  setter.call(input, {json.dumps(value)});\n"
        "  input.dispatchEvent(new Event('input', {bubbles: true}));\n"
        "  input.dispatchEvent(new Event('change', {bubbles: true}));\n"
        "  return true;\n"
        "})()"
    )
    if not evaluate(driver, expression):
        raise RuntimeError(f"Penpot login input not found: {selector}")


def websocket_token(driver: webdriver.Chrome) -> str:
    """Return the active plugin token without logging or displaying it."""
    for row in driver.get_log("performance"):
        try:
            message = json.loads(row["message"])["message"]
        except (KeyError, TypeError, json.JSONDecodeError):
            continue
        if message.get("method") != "Network.webSocketCreated":
            continue
        url = message.get("params", {}).get("url", "")
        if urlsplit(url).path != "/mcp/ws":
            continue
        query = dict(parse_qsl(urlsplit(url).query))
        access_value = query.get("userToken") or query.get("token")
        if access_value:
            return access_value
    raise RuntimeError("Could not identify the active Penpot MCP plugin token")


def synchronize_local_client(access_value: str) -> None:
    """Align the ignored local MCP URL and Codex registration with the browser."""
    mcp_url = "http://localhost:9001/mcp/stream?" + urlencode(
        {"userToken": access_value}
    )
    MCP_ENV.write_text(f"PENPOT_MCP_URL={mcp_url}\n", encoding="utf-8")
    MCP_ENV.chmod(0o600)
    subprocess.run(
        ["codex", "mcp", "remove", "penpot"],
        check=False,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    subprocess.run(
        ["codex", "mcp", "add", "penpot", "--url", mcp_url],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    # Emit only a non-secret fingerprint so operators can diagnose stale sessions.
    print(
        "Local MCP client synchronized (token fingerprint "
        f"{sha256(access_value.encode()).hexdigest()[:12]}).",
        flush=True,
    )


def main() -> None:
    credentials = read_env(CREDENTIALS)
    options = Options()
    options.add_argument("--headless=new")
    options.add_argument("--no-sandbox")
    options.add_argument("--disable-dev-shm-usage")
    options.add_argument("--enable-unsafe-swiftshader")
    options.add_argument("--window-size=1800,1100")
    options.set_capability("goog:loggingPrefs", {"performance": "ALL"})

    driver = webdriver.Chrome(options=options)
    stopping = False

    def stop(_signum, _frame) -> None:
        nonlocal stopping
        stopping = True

    signal.signal(signal.SIGINT, stop)
    signal.signal(signal.SIGTERM, stop)

    try:
        driver.get("http://localhost:9001/#/auth/login")
        time.sleep(3)
        set_input(driver, "input[name=email]", credentials["PENPOT_PROFILE_EMAIL"])
        set_input(
            driver, "input[name=password]", credentials["PENPOT_PROFILE_PASSWORD"]
        )
        evaluate(driver, "document.querySelector('button[type=submit]').click()")
        time.sleep(5)
        driver.get(DESIGN_URL)
        time.sleep(8)
        connected = evaluate(
            driver,
            """(() => {
              const button = document.querySelector('button[data-testid=mcp-btn]');
              if (!button) return false;
              button.click();
              return true;
            })()""",
        )
        if not connected:
            raise RuntimeError("Penpot MCP button was not available")
        time.sleep(3)
        body_text = evaluate(driver, "document.body.innerText || ''") or ""
        if "MCP connected" not in body_text:
            raise RuntimeError("Penpot MCP plugin did not report a connected state")
        synchronize_local_client(websocket_token(driver))
        print("Penpot MCP browser bridge connected to the SolvAI Figure 1 file.", flush=True)
        while not stopping:
            time.sleep(5)
    finally:
        driver.quit()


if __name__ == "__main__":
    main()
