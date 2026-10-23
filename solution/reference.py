#!/usr/bin/env python3
"""Recover missing mission-sample events and publish a design-based report."""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import os
import re
import sys
import tempfile
from pathlib import Path


TOP_KEYS = {
    "mission_id",
    "night_id",
    "exposures",
    "events",
    "calibration",
    "analysis_confidence",
}
PANEL_KEYS = {"exposure_ref", "frame_count", "exposure_s", "photon_counts"}
STRATUM_KEYS = {
    "event_id",
    "sky_area_sq_deg",
    "exposure_count",
    "exposure_total",
    "photon_total",
    "photon_sq_total",
}
CALIBRATION_KEYS = {"z_score", "tie_break"}
LEDGER_KEYS = {"mission_id", "input_sha256", "result_sha256", "survey_mean_rate_cps"}
PANEL_PATTERN = re.compile(r"^[A-Z][A-Z0-9_-]{2,15}$")
STRATUM_PATTERN = re.compile(r"^[A-Z][A-Z0-9_-]{0,23}$")
HEX_DIGITS = set("0123456789abcdef")


class Refusal(Exception):
    """A user-facing refusal with a stable business code."""

    def __init__(self, kind: str) -> None:
        super().__init__(kind)
        self.kind = kind


def fail(kind: str) -> None:
    sys.stderr.write(f"运行拒收: {kind}\n")
    raise SystemExit(1)


def finite(value: object) -> bool:
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(float(value))
    )


def integer(value: object, low: int, high: int) -> bool:
    return (
        isinstance(value, int)
        and not isinstance(value, bool)
        and low <= value <= high
    )


def strict_pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise Refusal("BAD_INPUT")
        result[key] = value
    return result


def ledger_pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise Refusal("BAD_LEDGER")
        result[key] = value
    return result


def visible_label(value: object, low: int, high: int) -> bool:
    return (
        isinstance(value, str)
        and low <= len(value) <= high
        and value == value.strip()
        and all(ord(character) >= 32 and ord(character) != 127 for character in value)
    )


def validate(payload):
    if not isinstance(payload, dict) or set(payload) != TOP_KEYS:
        raise Refusal("BAD_INPUT")
    mission_id = payload["mission_id"]
    night_id = payload["night_id"]
    if not visible_label(mission_id, 6, 48) or not mission_id.startswith("OCC-"):
        raise Refusal("BAD_INPUT")
    if not visible_label(night_id, 2, 48):
        raise Refusal("BAD_INPUT")
    exposures = payload["exposures"]
    if not isinstance(exposures, list) or not 12 <= len(exposures) <= 24:
        raise Refusal("BAD_INPUT")
    seen_exposures = set()
    clean_exposures = []
    for row in exposures:
        if not isinstance(row, dict) or set(row) != PANEL_KEYS:
            raise Refusal("BAD_INPUT")
        reference = row["exposure_ref"]
        if (
            not isinstance(reference, str)
            or PANEL_PATTERN.fullmatch(reference) is None
            or reference in seen_exposures
        ):
            raise Refusal("BAD_INPUT")
        seen_exposures.add(reference)
        count = row["frame_count"]
        if not integer(count, 10, 5000):
            raise Refusal("BAD_INPUT")
        weight = row["exposure_s"]
        if not finite(weight) or not 0 < float(weight) <= 100:
            raise Refusal("BAD_INPUT")
        photon_counts = row["photon_counts"]
        if not isinstance(photon_counts, list) or len(photon_counts) != count:
            raise Refusal("BAD_INPUT")
        if any(not integer(photon, 0, 1000) for photon in photon_counts):
            raise Refusal("BAD_INPUT")
        numeric_photon_counts = [int(photon) for photon in photon_counts]
        clean_exposures.append(
            {
                "exposure_ref": reference,
                "frame_count": count,
                "exposure_s": float(weight),
                "photon_counts": numeric_photon_counts,
            }
        )
    events = payload["events"]
    if not isinstance(events, list) or not 3 <= len(events) <= 5:
        raise Refusal("BAD_INPUT")
    seen_events = set()
    clean_events = []
    for row in events:
        if not isinstance(row, dict) or set(row) != STRATUM_KEYS:
            raise Refusal("BAD_INPUT")
        identifier = row["event_id"]
        if (
            not isinstance(identifier, str)
            or STRATUM_PATTERN.fullmatch(identifier) is None
            or identifier in seen_events
        ):
            raise Refusal("BAD_INPUT")
        seen_events.add(identifier)
        sky_area_sq_deg = row["sky_area_sq_deg"]
        if not integer(sky_area_sq_deg, 2, 1_000_000_000):
            raise Refusal("BAD_INPUT")
        exposure_count = row["exposure_count"]
        if not integer(exposure_count, 1, len(clean_exposures)):
            raise Refusal("BAD_INPUT")
        numeric = {}
        for key in (
            "exposure_total",
            "photon_total",
            "photon_sq_total",
        ):
            value = row[key]
            if not finite(value) or float(value) < 0:
                raise Refusal("BAD_INPUT")
            numeric[key] = float(value)
        if numeric["exposure_total"] <= 0 or numeric["photon_total"] <= 0:
            raise Refusal("BAD_INPUT")
        clean_events.append(
            {
                "event_id": identifier,
                "sky_area_sq_deg": sky_area_sq_deg,
                "exposure_count": exposure_count,
                **numeric,
            }
        )
    if sum(row["exposure_count"] for row in clean_events) != len(clean_exposures):
        raise Refusal("BAD_INPUT")
    calibration = payload["calibration"]
    if not isinstance(calibration, dict) or set(calibration) != CALIBRATION_KEYS:
        raise Refusal("BAD_INPUT")
    if not finite(calibration["z_score"]) or not 0 < float(calibration["z_score"]) <= 10:
        raise Refusal("BAD_INPUT")
    if calibration["tie_break"] != "lexicographic":
        raise Refusal("BAD_INPUT")
    confidence = payload["analysis_confidence"]
    if not finite(confidence) or not 0.8 <= float(confidence) <= 0.99:
        raise Refusal("BAD_INPUT")
    return (
        mission_id,
        night_id,
        sorted(clean_exposures, key=lambda item: item["exposure_ref"]),
        sorted(clean_events, key=lambda item: item["event_id"]),
        float(calibration["z_score"]),
        float(confidence),
    )


def sample_signature(sample):
    count = sample["frame_count"]
    weight = sample["exposure_s"]
    photon_total = sum(sample["photon_counts"])
    square_total = sum(value * value for value in sample["photon_counts"])
    return (
        1.0,
        weight * count,
        weight * photon_total,
        weight * square_total,
    )


def event_signature(event):
    return (
        float(event["exposure_count"]),
        event["exposure_total"],
        event["photon_total"],
        event["photon_sq_total"],
    )


def same_signature(left, right, tolerance=1e-7):
    return all(
        abs(left[index] - right[index]) <= tolerance
        for index in range(len(left))
    )


def enumerate_candidates(exposures, event):
    target = event_signature(event)
    count = event["exposure_count"]
    candidates = []
    for indexes in itertools.combinations(range(len(exposures)), count):
        total = [0.0, 0.0, 0.0, 0.0]
        for index in indexes:
            signature = sample_signature(exposures[index])
            for offset in range(4):
                total[offset] += signature[offset]
        if same_signature(total, target):
            candidates.append(indexes)
    candidates.sort(
        key=lambda row: tuple(exposures[index]["exposure_ref"] for index in row)
    )
    return candidates


def recover_partition(exposures, events):
    by_event = {
        row["event_id"]: enumerate_candidates(exposures, row)
        for row in events
    }
    if any(not choices for choices in by_event.values()):
        raise Refusal("BAD_MODEL")

    def search(position, used, assignment):
        if position == len(events):
            return dict(assignment)
        event = events[position]
        for candidate in by_event[event["event_id"]]:
            if any(index in used for index in candidate):
                continue
            updated = set(used)
            updated.update(candidate)
            next_assignment = dict(assignment)
            next_assignment[event["event_id"]] = candidate
            found = search(position + 1, updated, next_assignment)
            if found is not None:
                return found
        return None

    recovered = search(0, set(), {})
    if recovered is None:
        raise Refusal("BAD_MODEL")
    return recovered


def calculate_report(mission_id, exposures, events, z_score):
    assignment = recover_partition(exposures, events)
    total_sky_area_sq_deg = sum(row["sky_area_sq_deg"] for row in events)
    means = []
    variances = []
    references = []
    for event in events:
        chosen = [exposures[index] for index in assignment[event["event_id"]]]
        chosen.sort(key=lambda item: item["exposure_ref"])
        mean = event["photon_total"] / event["exposure_total"]
        exposure_total = event["exposure_total"]
        weight_square_sum = sum(
            sample["frame_count"] * sample["exposure_s"] ** 2
            for sample in chosen
        )
        effective_n = exposure_total * exposure_total / weight_square_sum
        centered = sum(
            sample["exposure_s"] * (photon - mean) ** 2
            for sample in chosen
            for photon in sample["photon_counts"]
        )
        denominator = exposure_total - weight_square_sum / exposure_total
        sample_variance = centered / denominator
        share = event["sky_area_sq_deg"] / total_sky_area_sq_deg
        fpc = max(0.0, 1.0 - effective_n / event["sky_area_sq_deg"])
        variance = share * share * fpc * sample_variance / effective_n
        means.append(mean)
        variances.append(variance)
        references.append([sample["exposure_ref"] for sample in chosen])
    combined = sum(
        event["sky_area_sq_deg"] / total_sky_area_sq_deg * mean
        for event, mean in zip(events, means)
    )
    rate_standard_error_cps = math.sqrt(sum(variances))
    return {
        "mission_id": mission_id,
        "status": "complete",
        "event_order": [row["event_id"] for row in events],
        "event_exposures": references,
        "event_rate_cps": means,
        "survey_mean_rate_cps": combined,
        "rate_standard_error_cps": rate_standard_error_cps,
        "ci_low_cps": combined - z_score * rate_standard_error_cps,
        "ci_high_cps": combined + z_score * rate_standard_error_cps,
    }


def valid_sha(value):
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in HEX_DIGITS for character in value)
    )


def read_ledger(path):
    if not path.exists():
        return []
    try:
        raw = path.read_bytes()
        text = raw.decode("utf-8")
    except (OSError, UnicodeError):
        raise Refusal("BAD_LEDGER")
    if not text:
        return []
    if not text.endswith("\n"):
        raise Refusal("BAD_LEDGER")
    records = []
    for raw_line in text.splitlines():
        try:
            record = json.loads(raw_line, object_pairs_hook=ledger_pairs)
        except (json.JSONDecodeError, UnicodeError, ValueError):
            raise Refusal("BAD_LEDGER")
        if not isinstance(record, dict) or list(record) != [
            "mission_id",
            "input_sha256",
            "result_sha256",
            "survey_mean_rate_cps",
        ]:
            raise Refusal("BAD_LEDGER")
        canonical = json.dumps(record, ensure_ascii=False, separators=(",", ":"))
        if raw_line != canonical:
            raise Refusal("BAD_LEDGER")
        if (
            not visible_label(record["mission_id"], 6, 48)
            or not valid_sha(record["input_sha256"])
            or not valid_sha(record["result_sha256"])
            or not finite(record["survey_mean_rate_cps"])
        ):
            raise Refusal("BAD_LEDGER")
        records.append(record)
    return records


def atomic_write(path, payload):
    descriptor = None
    stage = None
    try:
        descriptor, name = tempfile.mkstemp(prefix=".mission-", dir=path.parent)
        stage = Path(name)
        with os.fdopen(descriptor, "wb") as handle:
            descriptor = None
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(stage, path)
        stage = None
    finally:
        if descriptor is not None:
            os.close(descriptor)
        if stage is not None:
            stage.unlink(missing_ok=True)


def publish_result(target, payload):
    if (
        not target.parent.is_dir()
        or target.is_dir()
        or not os.access(target.parent, os.W_OK)
    ):
        raise Refusal("IO_FAILURE")
    try:
        atomic_write(target, payload)
    except OSError:
        raise Refusal("IO_FAILURE")


def publish_pair(source, target, ledger, input_bytes, output_bytes, combined):
    if (
        not target.parent.is_dir()
        or not ledger.parent.is_dir()
        or target.is_dir()
        or ledger.is_dir()
        or not os.access(target.parent, os.W_OK)
        or not os.access(ledger.parent, os.W_OK)
    ):
        raise Refusal("IO_FAILURE")
    if ledger in (source, target):
        raise Refusal("BAD_LEDGER")
    records = read_ledger(ledger)
    input_sha = hashlib.sha256(input_bytes).hexdigest()
    result_sha = hashlib.sha256(output_bytes).hexdigest()
    for record in records:
        if record["input_sha256"] == input_sha:
            if (
                record["result_sha256"] == result_sha
                and record["mission_id"] == json.loads(output_bytes)["mission_id"]
                and record["survey_mean_rate_cps"] == combined
            ):
                publish_result(target, output_bytes)
                return
            raise Refusal("BAD_LEDGER")
        if (
            record["result_sha256"] == result_sha
            and record["input_sha256"] != input_sha
        ):
            raise Refusal("BAD_LEDGER")
    new_record = {
        "mission_id": json.loads(output_bytes)["mission_id"],
        "input_sha256": input_sha,
        "result_sha256": result_sha,
        "survey_mean_rate_cps": combined,
    }
    new_line = (
        json.dumps(new_record, ensure_ascii=False, separators=(",", ":")) + "\n"
    ).encode("utf-8")
    old_target = target.read_bytes() if target.exists() else None
    old_ledger = ledger.read_bytes() if ledger.exists() else b""
    ledger_bytes = old_ledger + new_line
    try:
        atomic_write(ledger, ledger_bytes)
        try:
            atomic_write(target, output_bytes)
        except OSError:
            if ledger.exists():
                if old_ledger:
                    atomic_write(ledger, old_ledger)
                else:
                    ledger.unlink(missing_ok=True)
            raise
    except OSError:
        if old_target is None:
            target.unlink(missing_ok=True)
        else:
            try:
                atomic_write(target, old_target)
            except OSError:
                pass
        raise Refusal("IO_FAILURE")


def main(argv):
    if len(argv) not in (3, 4):
        fail("BAD_INPUT")
    source = Path(argv[1])
    target = Path(argv[2])
    try:
        input_bytes = source.read_bytes()
        payload = json.loads(
            input_bytes.decode("utf-8"),
            object_pairs_hook=strict_pairs,
        )
        mission_id, _, exposures, events, z_score, _ = validate(payload)
        report = calculate_report(mission_id, exposures, events, z_score)
        encoded = (
            json.dumps(report, ensure_ascii=False, separators=(",", ":")) + "\n"
        ).encode("utf-8")
        if len(argv) == 4:
            publish_pair(
                source,
                target,
                Path(argv[3]),
                input_bytes,
                encoded,
                report["survey_mean_rate_cps"],
            )
        else:
            publish_result(target, encoded)
        return 0
    except Refusal as exc:
        fail(exc.kind)
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError):
        fail("BAD_INPUT")


if __name__ == "__main__":
    sys.exit(main(sys.argv))
