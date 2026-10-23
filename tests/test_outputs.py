"""Acceptance tests for missing-event mission reconstruction."""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys
from pathlib import Path

import pytest


ENTRY = os.environ.get("OCCULTATION_RATE_ENTRY", "/app/occultation_event_rates.py")


def sample(reference, count, weight, start):
    photon_counts = [((start + index * 3) % 91) + 5 for index in range(count)]
    return {
        "exposure_ref": reference,
        "frame_count": count,
        "exposure_s": weight,
        "photon_counts": photon_counts,
    }


def make_document(sample_total=12, event_total=3, shuffle=False):
    if (sample_total, event_total) == (12, 3):
        groups = [
            [
                sample("SED-A01", 10, 1.25, 2),
                sample("SED-A02", 14, 2.5, 13),
                sample("SED-A03", 11, 3.25, 31),
                sample("SED-A04", 17, 1.75, 47),
            ],
            [
                sample("SED-B01", 13, 1.5, 7),
                sample("SED-B02", 16, 2.25, 23),
                sample("SED-B03", 12, 4.0, 39),
                sample("SED-B04", 19, 0.75, 55),
            ],
            [
                sample("SED-C01", 15, 1.1, 11),
                sample("SED-C02", 18, 2.1, 29),
                sample("SED-C03", 10, 5.0, 43),
                sample("SED-C04", 20, 1.4, 61),
            ],
        ]
        sky_area_sq_degs = [1200, 2600, 4100]
        identifiers = ["A", "B", "C"]
    else:
        groups = []
        per_group = sample_total // event_total
        sky_area_sq_degs = []
        identifiers = []
        ordinal = 0
        for group_index in range(event_total):
            rows = []
            for member_index in range(per_group):
                ordinal += 1
                rows.append(
                    sample(
                        f"P-{chr(65 + group_index)}{member_index:02d}",
                        10 + (ordinal % 9),
                        0.5 + ((ordinal * 7) % 19) / 4,
                        3 + ordinal * 11 + group_index * 5,
                    )
                )
            groups.append(rows)
            sky_area_sq_degs.append(900 + group_index * 1100)
            identifiers.append(chr(65 + group_index))
    truth = {
        identifier: [row["exposure_ref"] for row in rows]
        for identifier, rows in zip(identifiers, groups)
    }
    exposures = [row for group in groups for row in group]
    if shuffle:
        exposures = exposures[6:9] + exposures[:6] + exposures[9:]
    events = []
    for identifier, rows, sky_area_sq_deg in zip(identifiers, groups, sky_area_sq_degs):
        exposure_total = sum(row["frame_count"] * row["exposure_s"] for row in rows)
        photon_total = sum(
            row["exposure_s"] * sum(row["photon_counts"]) for row in rows
        )
        square_total = sum(
            row["exposure_s"] * sum(value * value for value in row["photon_counts"])
            for row in rows
        )
        events.append(
            {
                "event_id": identifier,
                "sky_area_sq_deg": sky_area_sq_deg,
                "exposure_count": len(rows),
                "exposure_total": exposure_total,
                "photon_total": photon_total,
                "photon_sq_total": square_total,
            }
        )
    document = {
        "mission_id": "OCC-2026-04",
        "night_id": "NIGHT-AUTUMN",
        "exposures": exposures,
        "events": events,
        "calibration": {"z_score": 1.959963984540054, "tie_break": "lexicographic"},
        "analysis_confidence": 0.95,
    }
    return document, truth


def recalibrate(document):
    groups = {}
    for sample_row in document["exposures"]:
        prefix = sample_row["exposure_ref"].split("-")[1][0]
        groups.setdefault(prefix, []).append(sample_row)
    for event in document["events"]:
        rows = groups[event["event_id"]]
        event["exposure_count"] = len(rows)
        event["exposure_total"] = sum(
            row["frame_count"] * row["exposure_s"] for row in rows
        )
        event["photon_total"] = sum(
            row["exposure_s"] * sum(row["photon_counts"]) for row in rows
        )
        event["photon_sq_total"] = sum(
            row["exposure_s"] * sum(value * value for value in row["photon_counts"])
            for row in rows
        )
    return document


def write_document(path, document, raw=None):
    if raw is None:
        raw = json.dumps(document, ensure_ascii=False, separators=(",", ":")) + "\n"
    path.write_bytes(raw.encode("utf-8"))
    return path


def invoke(tmp_path, document, raw=None, with_ledger=False):
    source = write_document(tmp_path / "mission.json", document, raw)
    target = tmp_path / "report.json"
    ledger = tmp_path / "delivery.jsonl"
    argv = [sys.executable, ENTRY, str(source), str(target)]
    if with_ledger:
        argv.append(str(ledger))
    result = subprocess.run(argv, capture_output=True, text=True, timeout=30)
    return result, source, target, ledger


def expected_values(document, truth):
    by_reference = {row["exposure_ref"]: row for row in document["exposures"]}
    events = sorted(document["events"], key=lambda row: row["event_id"])
    total_sky_area_sq_deg = sum(row["sky_area_sq_deg"] for row in events)
    means = []
    variances = []
    for event in events:
        rows = [by_reference[value] for value in truth[event["event_id"]]]
        mean = event["photon_total"] / event["exposure_total"]
        weight_square_sum = sum(
            row["frame_count"] * row["exposure_s"] ** 2 for row in rows
        )
        effective_n = event["exposure_total"] ** 2 / weight_square_sum
        centered = sum(
            row["exposure_s"] * (photon - mean) ** 2
            for row in rows
            for photon in row["photon_counts"]
        )
        sample_variance = centered / (
            event["exposure_total"] - weight_square_sum / event["exposure_total"]
        )
        share = event["sky_area_sq_deg"] / total_sky_area_sq_deg
        fpc = max(0.0, 1.0 - effective_n / event["sky_area_sq_deg"])
        variances.append(share * share * fpc * sample_variance / effective_n)
        means.append(mean)
    combined = sum(
        row["sky_area_sq_deg"] / total_sky_area_sq_deg * mean
        for row, mean in zip(events, means)
    )
    rate_standard_error_cps = math.sqrt(sum(variances))
    z_score = document["calibration"]["z_score"]
    return means, combined, rate_standard_error_cps, combined - z_score * rate_standard_error_cps, combined + z_score * rate_standard_error_cps


def read_report(path):
    return json.loads(path.read_text(encoding="utf-8"))


def test_success_is_quiet_and_creates_report(tmp_path):
    """题面：「返回码 0 的路径保持 stdout 和 stderr 空白；」；"""
    document, _ = make_document()
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 0 and result.stdout == "" and result.stderr == ""
    assert target.is_file()


def test_existing_compact_json_suffix_is_preserved(tmp_path):
    """题面：「计数率单不是 JSONL：整个文件只有一个 JSON 对象和一个物理行，该行承载上面九个顶层字段」；"""
    document, _ = make_document()
    result, _, target, _ = invoke(tmp_path, document)
    raw = target.read_bytes()
    parsed = json.loads(raw)
    assert result.returncode == 0 and raw.endswith(b"\n") and raw.count(b"\n") == 1
    assert isinstance(parsed, dict) and len(parsed) == 9


def test_report_has_exact_ordered_keys(tmp_path):
    """题面：「计数率单必须是一行紧凑 JSON，末尾单个 LF，成员和顺序固定：」；"""
    document, _ = make_document()
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 0
    assert list(read_report(target)) == [
        "mission_id",
        "status",
        "event_order",
        "event_exposures",
        "event_rate_cps",
        "survey_mean_rate_cps",
        "rate_standard_error_cps",
        "ci_low_cps",
        "ci_high_cps",
    ]


def test_hidden_partition_is_recovered(tmp_path):
    """题面：「每段曝光只能进入一个事件，不能漏掉任何一段；」；"""
    document, truth = make_document()
    result, _, target, _ = invoke(tmp_path, document)
    references = read_report(target)["event_exposures"]
    assert result.returncode == 0
    assert references == [sorted(truth[key]) for key in ("A", "B", "C")]


def test_input_order_does_not_change_business_bytes(tmp_path):
    """题面：「若事件归属存在并列，采用事件编号与事件内曝光编号串联的字典序最小方案；」"""
    first, _ = make_document()
    shuffled, _ = make_document(shuffle=True)
    result_one, _, target_one, _ = invoke(tmp_path, first)
    write_document(tmp_path / "mission.json", shuffled)
    result_two = subprocess.run(
        [sys.executable, ENTRY, str(tmp_path / "mission.json"), str(tmp_path / "second.json")],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result_one.returncode == 0 and result_two.returncode == 0
    assert target_one.read_bytes() == (tmp_path / "second.json").read_bytes()


def test_event_rate_cps_use_published_exposure_totals(tmp_path):
    """题面：「event_mean_h = photon_total_h / exposure_total_h」；"""
    document, truth = make_document()
    result, _, target, _ = invoke(tmp_path, document)
    expected = expected_values(document, truth)[0]
    assert result.returncode == 0
    assert read_report(target)["event_rate_cps"] == pytest.approx(expected)


def test_survey_mean_rate_cps_uses_sky_area_sq_deg_shares(tmp_path):
    """题面：「survey_mean_rate_cps = sum(share_h * event_mean_h)」；"""
    document, truth = make_document()
    result, _, target, _ = invoke(tmp_path, document)
    expected = expected_values(document, truth)[1]
    assert result.returncode == 0 and read_report(target)["survey_mean_rate_cps"] == pytest.approx(expected)


def test_rate_standard_error_cps_uses_effective_sample_size(tmp_path):
    """题面：「rate_standard_error_cps = sqrt(sum(variance_term_h))」；"""
    document, truth = make_document()
    result, _, target, _ = invoke(tmp_path, document)
    expected = expected_values(document, truth)[2]
    assert result.returncode == 0 and read_report(target)["rate_standard_error_cps"] == pytest.approx(expected)


def test_confidence_interval_uses_declared_z_score(tmp_path):
    """题面：「置信区间端点分别是 `survey_mean_rate_cps - z_score*rate_standard_error_cps` 和」；"""
    document, truth = make_document()
    result, _, target, _ = invoke(tmp_path, document)
    expected = expected_values(document, truth)
    report = read_report(target)
    assert result.returncode == 0 and report["ci_low_cps"] == pytest.approx(expected[3])
    assert report["ci_high_cps"] == pytest.approx(expected[4])


def test_analysis_confidence_is_validated_but_not_used_as_z(tmp_path):
    """题面：「业务置信字段不提供计算倍数」。"""
    document, _ = make_document()
    changed = json.loads(json.dumps(document))
    changed["analysis_confidence"] = 0.9
    original, _, original_target, _ = invoke(tmp_path, document)
    altered, _, altered_target, _ = invoke(tmp_path, changed)
    assert original.returncode == 0 and altered.returncode == 0
    assert original_target.read_bytes() == altered_target.read_bytes()


def test_minimum_exposure_boundary_is_accepted(tmp_path):
    """题面：「`exposures` 有 12 到 24 段；」；"""
    document, _ = make_document(12, 3)
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 0 and target.is_file()


def test_maximum_exposure_boundary_is_accepted(tmp_path):
    """题面：「`exposures` 有 12 到 24 段；」；"""
    document, _ = make_document(24, 4)
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 0 and len(read_report(target)["event_order"]) == 4


def test_eleven_exposures_are_refused(tmp_path):
    """题面：「`exposures` 有 12 到 24 段；」；"""
    document, _ = make_document()
    document["exposures"] = document["exposures"][:-1]
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()


def test_extra_top_level_key_is_refused(tmp_path):
    """题面：「只有六个顶层名称可以出现在夜班文档中」。"""
    document, _ = make_document()
    document["operator"] = "night"
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()


def test_mission_identifier_prefix_is_refused(tmp_path):
    """题面：「mission 前缀固定 `OCC-`」"""
    document, _ = make_document()
    document["mission_id"] = "HHS-2026-04"
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()


def test_invalid_event_identifier_is_refused(tmp_path):
    """题面：「`event_id` 匹配 `[A-Z][A-Z0-9_-]{0,23}` 且不可重名；」；"""
    document, _ = make_document()
    document["events"][0]["event_id"] = "a"
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()


def test_night_identifier_with_surrounding_space_is_refused(tmp_path):
    """题面：「`mission_id` 与 `night_id` 均需 2–48 个可打印字符」"""
    document, _ = make_document()
    document["night_id"] = " campaign "
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()
    minimum, _ = make_document()
    minimum["night_id"] = "AB"
    minimum_result, _, minimum_target, _ = invoke(tmp_path, minimum)
    assert minimum_result.returncode == 0 and json.loads(minimum_target.read_text(encoding="utf-8"))["mission_id"] == minimum["mission_id"]
    spaced, _ = make_document()
    spaced["night_id"] = "A B"
    spaced_result, _, spaced_target, _ = invoke(tmp_path, spaced)
    assert spaced_result.returncode == 0 and json.loads(spaced_target.read_text(encoding="utf-8"))["mission_id"] == spaced["mission_id"]
    spaced_bytes = spaced_target.read_bytes()
    control, _ = make_document()
    control["night_id"] = "A\tB"
    control_result, _, control_target, _ = invoke(tmp_path, control)
    assert control_result.returncode == 1 and control_result.stderr == "运行拒收: BAD_INPUT\n" and control_target.read_bytes() == spaced_bytes


def test_duplicate_exposure_reference_is_refused(tmp_path):
    """题面：「`exposure_ref` 匹配 `[A-Z][A-Z0-9_-]{2,15}`，同夜不可重号；」；"""
    document, _ = make_document()
    document["exposures"][1]["exposure_ref"] = document["exposures"][0]["exposure_ref"]
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()


def test_calibration_extra_key_is_refused(tmp_path):
    """题面：「`calibration` 只有以下两个键：」；"""
    document, _ = make_document()
    document["calibration"]["method"] = "manual"
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()


def test_missing_exposure_key_is_refused(tmp_path):
    """题面：「每段只允许四个成员：」；"""
    document, _ = make_document()
    document["exposures"][0].pop("exposure_s")
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()


def test_duplicate_json_key_is_refused(tmp_path):
    """题面：「UTF-8 解码遇到同名键立即退回」。"""
    document, _ = make_document()
    raw = json.dumps(document, ensure_ascii=False, separators=(",", ":"))[:-1] + ',"mission_id":"OCC-X"}\n'
    result, _, target, _ = invoke(tmp_path, document, raw=raw)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()


def test_boolean_frame_count_is_refused(tmp_path):
    """题面：「真/假码不构成这里允许的整数或小数；」；"""
    document, _ = make_document()
    document["exposures"][0]["frame_count"] = True
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()


def test_photon_length_mismatch_is_refused(tmp_path):
    """题面：「每段 `photon_counts` 必须与 `frame_count` 等长」"""
    document, _ = make_document()
    document["exposures"][0]["photon_counts"].pop()
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()


def test_boolean_photon_is_refused(tmp_path):
    """题面：「真/假码不构成这里允许的整数或小数」。"""
    document, _ = make_document()
    document["exposures"][0]["photon_counts"][0] = False
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()


def test_event_exposure_count_must_sum_to_exposure_total(tmp_path):
    """题面：「每个事件的曝光段计数是恢复归属时的硬约束，不允许与实际成员数偏离；」；"""
    document, _ = make_document()
    document["events"][0]["exposure_count"] += 1
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()


def test_inconsistent_aggregate_reports_bad_model(tmp_path):
    """题面：「如果没有任何事件分配能同时满足四项汇总，退出码 1 并只输出 `运行拒收: BAD_MODEL`；」；"""
    document, _ = make_document()
    document["events"][0]["photon_total"] += 9.0
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_MODEL\n" and not target.exists()


def test_wrong_tie_break_is_refused(tmp_path):
    """题面：「`tie_break` 必须是 `lexicographic`；」；"""
    document, _ = make_document()
    document["calibration"]["tie_break"] = "arrival"
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()


def test_ledger_has_exact_four_key_line(tmp_path):
    """题面：「追加流水是一份 UTF-8 JSONL」。"""
    document, _ = make_document()
    result, _, target, ledger = invoke(tmp_path, document, with_ledger=True)
    line = ledger.read_text(encoding="utf-8")
    record = json.loads(line)
    assert result.returncode == 0 and line.endswith("\n") and list(record) == [
        "mission_id",
        "input_sha256",
        "result_sha256",
        "survey_mean_rate_cps",
    ]
    assert record["result_sha256"] == __import__("hashlib").sha256(target.read_bytes()).hexdigest()


def test_repeated_success_keeps_ledger_bytes(tmp_path):
    """题面：「重放同一文件时，追加流水保持在原长度；」；"""
    document, _ = make_document()
    first, _, _, ledger = invoke(tmp_path, document, with_ledger=True)
    before = ledger.read_bytes()
    second = subprocess.run(
        [sys.executable, ENTRY, str(tmp_path / "mission.json"), str(tmp_path / "report.json"), str(ledger)],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert first.returncode == 0 and second.returncode == 0 and ledger.read_bytes() == before


def test_distinct_input_appends_in_success_order(tmp_path):
    """题面：「只要结果字节发生变化，就向流水尾部添加一条新记录；」；"""
    document, _ = make_document()
    first, _, _, ledger = invoke(tmp_path, document, with_ledger=True)
    changed = json.loads(json.dumps(document))
    changed["mission_id"] = "OCC-2026-05"
    changed["night_id"] = "NIGHT-FOLLOWUP"
    changed["exposures"][0]["photon_counts"][0] += 1
    recalibrate(changed)
    second, _, _, _ = invoke(tmp_path, changed, with_ledger=True)
    assert first.returncode == 0 and second.returncode == 0
    assert [json.loads(line)["mission_id"] for line in ledger.read_text(encoding="utf-8").splitlines()] == [
        document["mission_id"],
        changed["mission_id"],
    ]


def test_bad_ledger_is_refused_without_replacing_report(tmp_path):
    """题面：「流水解析或绑定关系不合法时，拒收码为 `BAD_LEDGER`；」；"""
    document, _ = make_document()
    result, source, target, ledger = invoke(tmp_path, document)
    target.write_bytes(b"old bytes\n")
    ledger.write_bytes(b"{}\n")
    attempted = subprocess.run(
        [sys.executable, ENTRY, str(source), str(target), str(ledger)],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0 and attempted.returncode == 1
    assert attempted.stderr == "运行拒收: BAD_LEDGER\n" and target.read_bytes() == b"old bytes\n"


def test_conflicting_ledger_digest_is_refused(tmp_path):
    """题面：「输入摘要重放或结果摘要归属发生冲突时，终止运行并只报告 `运行拒收: BAD_LEDGER`；」"""
    document, _ = make_document()
    first, _, target, ledger = invoke(tmp_path, document, with_ledger=True)
    old_result = target.read_bytes()
    old_ledger = ledger.read_bytes()
    changed = json.loads(json.dumps(document))
    changed["night_id"] = "NIGHT-RETRY"
    second, _, _, _ = invoke(tmp_path, changed, with_ledger=True)
    assert first.returncode == 0 and second.returncode == 1
    assert second.stderr == "运行拒收: BAD_LEDGER\n" and target.read_bytes() == old_result
    assert ledger.read_bytes() == old_ledger


def test_ledger_parent_failure_preserves_old_result(tmp_path):
    """题面：「计数率单和台账应一起成功或都保持旧字节；」；"""
    document, _ = make_document()
    result, source, target, _ = invoke(tmp_path, document)
    target.write_bytes(b"previous report\n")
    absent_ledger = tmp_path / "nested" / "delivery.jsonl"
    attempted = subprocess.run(
        [sys.executable, ENTRY, str(source), str(target), str(absent_ledger)],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0 and attempted.returncode == 1
    assert attempted.stderr == "运行拒收: IO_FAILURE\n" and target.read_bytes() == b"previous report\n"


def test_input_file_is_not_modified(tmp_path):
    """题面：「源夜班数据不得被改写；」；"""
    document, _ = make_document()
    result, source, _, _ = invoke(tmp_path, document, with_ledger=True)
    assert result.returncode == 0 and source.read_bytes() == (
        json.dumps(document, ensure_ascii=False, separators=(",", ":")) + "\n"
    ).encode("utf-8")


def test_success_leaves_no_staging_files(tmp_path):
    """题面：「发布动作结束后目录中不得残留临时产物」；"""
    document, _ = make_document()
    result, _, _, _ = invoke(tmp_path, document, with_ledger=True)
    names = sorted(path.name for path in tmp_path.iterdir())
    assert result.returncode == 0 and names == ["delivery.jsonl", "mission.json", "report.json"]


def test_repeat_without_ledger_is_byte_stable(tmp_path):
    """题面：「计数率单必须是一行紧凑 JSON，末尾单个 LF，成员和顺序固定：」；"""
    document, _ = make_document()
    first, _, target, _ = invoke(tmp_path, document)
    before = target.read_bytes()
    second = subprocess.run(
        [sys.executable, ENTRY, str(tmp_path / "mission.json"), str(target)],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert first.returncode == 0 and second.returncode == 0 and target.read_bytes() == before


def test_missing_input_is_bad_input(tmp_path):
    """题面：「若夜班材料读取或校验失败，退出码为 1，stderr 仅写 `运行拒收: BAD_INPUT`；」；"""
    target = tmp_path / "report.json"
    result = subprocess.run(
        [sys.executable, ENTRY, str(tmp_path / "absent.json"), str(target)],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()


def test_missing_result_parent_reports_io_failure(tmp_path):
    """题面：「路径读写失败时 stderr 中只保留 `运行拒收: IO_FAILURE`；」；"""
    document, _ = make_document()
    source = write_document(tmp_path / "batch.json", document)
    target = tmp_path / "absent" / "report.json"
    result = subprocess.run(
        [sys.executable, ENTRY, str(source), str(target)],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 1 and result.stderr == "运行拒收: IO_FAILURE\n" and not target.exists()


def test_failure_does_not_emit_traceback(tmp_path):
    """题面：「任何失败输出中均不得出现 Python `Traceback`」；"""
    document, _ = make_document()
    document["calibration"]["z_score"] = -1
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and "Traceback" not in result.stderr and not target.exists()


def test_maximum_group_shape_finishes_within_limit(tmp_path):
    """题面：「每段曝光只能进入一个事件」；"""
    document, _ = make_document(24, 4)
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 0 and target.is_file()


def test_three_night_inputs_produce_report(tmp_path):
    """题面：「值班组需要在望远镜交接前恢复事件归属并给出计数率；」"""
    document, _ = make_document()
    result, _, target, _ = invoke(tmp_path, document)
    report = read_report(target)
    assert result.returncode == 0 and report["mission_id"] == document["mission_id"]
    assert report["status"] == "complete" and report["event_order"] == ["A", "B", "C"]


def test_five_events_are_accepted(tmp_path):
    """题面：「事件列表长度可为 3、4 或 5，单条事件的结构键如下：」；"""
    document, _ = make_document(20, 5)
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 0 and read_report(target)["event_order"] == ["A", "B", "C", "D", "E"]


def test_event_order_follows_identifier_order(tmp_path):
    """题面：「`event_order` 按 `event_id` 升序；」；"""
    document, _ = make_document()
    result, _, target, _ = invoke(tmp_path, document)
    report = read_report(target)
    assert result.returncode == 0 and report["event_order"] == sorted(report["event_order"])


def test_lexicographic_normalization_is_stable(tmp_path):
    """题面：「字典序最小方案」。"""
    document, _ = make_document()
    shuffled, _ = make_document(shuffle=True)
    first, _, first_target, _ = invoke(tmp_path, document)
    write_document(tmp_path / "mission.json", shuffled)
    second = subprocess.run(
        [sys.executable, ENTRY, str(tmp_path / "mission.json"), str(tmp_path / "shuffled.json")],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert first.returncode == 0 and second.returncode == 0
    assert first_target.read_bytes() == (tmp_path / "shuffled.json").read_bytes()


def test_event_arrays_follow_identifier_order(tmp_path):
    """题面：「两个事件数组按 `event_order` 对齐」"""
    document, _ = make_document()
    result, _, target, _ = invoke(tmp_path, document)
    report = read_report(target)
    count = len(report["event_order"])
    assert result.returncode == 0 and len(report["event_exposures"]) == count
    assert len(report["event_rate_cps"]) == count


def test_member_references_are_sorted_within_event(tmp_path):
    """题面：「事件内部曝光引用保持升序」。"""
    document, _ = make_document()
    result, _, target, _ = invoke(tmp_path, document)
    report = read_report(target)
    assert result.returncode == 0 and all(row == sorted(row) for row in report["event_exposures"])


def test_ledger_is_append_only_jsonl(tmp_path):
    """题面：「追加流水是一份 UTF-8 JSONL」。"""
    document, _ = make_document()
    first, _, _, ledger = invoke(tmp_path, document, with_ledger=True)
    before = ledger.read_bytes()
    second, _, _, _ = invoke(tmp_path, document, with_ledger=True)
    assert first.returncode == 0 and second.returncode == 0 and ledger.read_bytes() == before
    assert before.decode("utf-8").endswith("\n")
