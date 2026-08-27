#!/usr/bin/env python3
"""Run the independent release leakage and inference audit."""

import json

from solv_ai.audit import write_leakage_audit

if __name__ == "__main__":
    print(json.dumps(write_leakage_audit(), indent=2))
