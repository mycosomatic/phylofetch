"""Persistent project config stored at ~/.phylofetch/config.json."""

import json
from pathlib import Path

CONFIG_PATH = Path.home() / ".phylofetch" / "config.json"


def load_config() -> dict:
    if CONFIG_PATH.exists():
        try:
            return json.loads(CONFIG_PATH.read_text())
        except (json.JSONDecodeError, OSError):
            return {}
    return {}


def save_config(updates: dict) -> dict:
    CONFIG_PATH.parent.mkdir(parents=True, exist_ok=True)
    cfg = load_config()
    cfg.update(updates)
    CONFIG_PATH.write_text(json.dumps(cfg, indent=2, sort_keys=True))
    return cfg
