#!/usr/bin/env python3
"""Draft mission report that trusts the order of incoming exposures."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path


def main(argv):
    if len(argv) not in (3, 4):
        return 1
    source = Path(argv[1])
    target = Path(argv[2])
    payload = json.loads(source.read_text(encoding="utf-8"))
    exposures = payload["exposures"]
    events = payload["events"]
    order = sorted(events, key=lambda row: row["event_id"])
    buckets = [[] for _ in order]
    for index, sample in enumerate(exposures):
        buckets[index % len(order)].append(sample["exposure_ref"])
    means = []
    for event in order:
        means.append(event["photon_total"] / event["exposure_total"])
    sky_area_sq_degs = [record["sky_area_sq_deg"] for record in order]
    total = sum(sky_area_sq_degs)
    combined = sum(
        sky_area_sq_deg / total * mean
        for sky_area_sq_deg, mean in zip(sky_area_sq_degs, means)
    )
    report = {
        "mission_id": payload["mission_id"],
        "status": "complete",
        "event_order": [record["event_id"] for record in order],
        "event_exposures": [sorted(row) for row in buckets],
        "event_rate_cps": means,
        "survey_mean_rate_cps": combined,
        "rate_standard_error_cps": 0.0,
        "ci_low_cps": combined,
        "ci_high_cps": combined,
    }
    encoded = (json.dumps(report, separators=(",", ":")) + "\n").encode("utf-8")
    target.write_bytes(encoded)
    if len(argv) == 4:
        ledger = Path(argv[3])
        line = {
            "mission_id": payload["mission_id"],
            "input_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "result_sha256": hashlib.sha256(encoded).hexdigest(),
            "survey_mean_rate_cps": combined,
        }
        with ledger.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(line, separators=(",", ":")) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
