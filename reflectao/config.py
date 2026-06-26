"""
Loads ReflectAO configuration from a .env file at the repository root.

To configure, copy .env.example to .env and fill in the values for your
machine. System environment variables always take precedence over .env values.
"""

import os
from pathlib import Path


def load_env_file():
    repo_root = Path(__file__).parent.parent
    env_path = repo_root / ".env"
    if not env_path.exists():
        env_path = repo_root / ".env.example"
    if not env_path.exists():
        return
    with open(env_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, _, value = line.partition("=")
            key = key.strip()
            value = value.strip().strip('"').strip("'")
            if key and key not in os.environ:
                os.environ[key] = value


load_env_file()
