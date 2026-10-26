#!/usr/bin/env python3
"""Bind supplied substantive review notes to the current advisor source inventory.

Usage: python3 tools/product/review-advisor-semantics.py NOTES.json NEW_OUTPUT_DIR
NOTES maps every variant id to [internal consistency review, unresolved judgments].
This assembles evidence; it does not perform or certify a scientific review.
"""
import csv
import hashlib
import json
from pathlib import Path
import sys

root = Path(__file__).resolve().parents[2]
if len(sys.argv) != 3:
    raise SystemExit(__doc__)
notes_path, out = Path(sys.argv[1]), Path(sys.argv[2])
if out.exists():
    raise SystemExit("Use a new output directory; existing evidence is never overwritten.")
out.mkdir(parents=True)
summary = {"matrix_complete": False, "scientific_acceptance": "not_established",
           "scope": "Supplied substantive notes joined to current source identities; not automatic semantic certification."}

def load(path):
    return json.loads(path.read_text(encoding="utf-8"))

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

try:
    notes = load(notes_path)
    variant_path = root / "inst/advisor/analysis-variants.json"
    rules_path = root / "inst/advisor/rules.json"
    resources_path = root / "inst/metadata/learning-resources.json"
    crosswalk_path = root / "inst/advisor/source-crosswalk.csv"
    variants = load(variant_path)
    rules = load(rules_path)
    resources = load(resources_path)
    with crosswalk_path.open(newline="", encoding="utf-8") as f:
        crosswalk = list(csv.DictReader(f))
    ids = [v["analysis_variant_id"] for v in variants["variants"]]
    if len(ids) != len(set(ids)) or set(notes) != set(ids):
        raise ValueError("Review notes must cover every unique registered variant exactly once.")
    if any(not isinstance(notes[v], list) or len(notes[v]) != 2 or
           any(not isinstance(s, str) or not s.strip() for s in notes[v]) for v in ids):
        raise ValueError("Each variant needs substantive consistency and unresolved-judgment notes.")
    references = {r["ref_id"]: r for r in rules["references"]}
    uzh = load(root / "inst/advisor/uzh-tree.json")
    dt2 = load(root / "inst/advisor/descrtab2-tree.json")
    uzh_nodes = {n["node_id"]: n for n in uzh["nodes"]}
    bindings = {v: [] for v in ids}
    for outcome in rules["outcomes"]:
        associations = [(outcome.get("analysis_variant_id"), "primary")]
        associations += [(a.get("analysis_variant_id"), "override") for a in outcome.get("alternatives", [])]
        associations += [(a.get("analysis_variant_id"), "rejected") for a in outcome.get("rejected_routes", [])]
        for variant, role in associations:
            if variant is None:
                continue
            if variant not in bindings:
                raise ValueError("Unknown rule variant: " + variant)
            bindings[variant].append({"outcome_id": outcome["outcome_id"], "role": role,
                "estimand_en": outcome.get("estimand_en"), "assumptions": outcome.get("assumptions", []),
                "references": [references[r] for r in outcome.get("references", [])],
                "r_template_id": outcome.get("r_template_id"),
                "interpretation": "Assumptions of an original primary route are context, not proof of override assumptions."})
    rows, detail = [], []
    for v in variants["variants"]:
        vid = v["analysis_variant_id"]
        template = root / "inst/advisor/templates" / ("tpl." + vid + ".v1.R")
        routes = bindings[vid]
        permitted = any(b["role"] != "rejected" for b in routes)
        if permitted != template.is_file():
            raise ValueError("Executable route/template mismatch for " + vid)
        learning_ids = v.get("learning_method_ids", [])
        links = []
        resolved = set()
        for r in resources:
            hits = [c for c in r.get("crosswalk", []) if c["method_id"] in learning_ids]
            if hits:
                resolved.update(c["method_id"] for c in hits)
                links.append({k: r.get(k) for k in ("resource_id", "title", "resource_type", "published_url",
                    "repo_url", "commit", "last_checked_utc", "known_issues", "variant_notes", "review_depth")})
                links[-1]["method_matches"] = hits
        if set(learning_ids) != resolved:
            raise ValueError("Unresolved learning method ids for " + vid)
        cw = [c for c in crosswalk if c["analysis_variant_id"] == vid]
        source_links = []
        for c in cw:
            node = uzh_nodes.get(c["uzh_node_id"])
            if node:
                source_links.append({"node_id": c["uzh_node_id"], "url": node.get("source_url"),
                    "snapshot_sha256": node.get("page_sha256"), "historical_soundness": node.get("condition_soundness")})
            if c["descrtab2_node_id"]:
                source_links.append({"node_id": c["descrtab2_node_id"], "url": dt2["source"]["url"],
                    "snapshot_sha256": dt2["source"]["sha256"]})
        route = "primary" if any(b["role"] == "primary" for b in routes) else "override_only" if permitted else "rejected_only"
        row = {"variant": vid, "route": route, "template_sha256": sha(template) if permitted else "",
               "estimand": v["estimand_en"], "substantive_internal_review": notes[vid][0],
               "unresolved_judgments_or_discrepancies": notes[vid][1],
               "bound_outcomes": "; ".join(b["role"] + ":" + b["outcome_id"] for b in routes),
               "source_crosswalk_ids": "; ".join(c["crosswalk_id"] for c in cw),
               "source_relationships": "; ".join(sorted({c["relationship"] for c in cw})),
               "learning_method_ids": "; ".join(learning_ids),
               "learning_resources": "; ".join(r["resource_id"] for r in links),
               "learning_known_issues": "; ".join(sorted({i for r in links for i in (r["known_issues"] or [])})),
               "source_acceptance": "Historical snapshots/crosswalk only; source fidelity and scientific acceptance not certified",
               "learning_acceptance": "Resolved links and retained caveats; current content not comprehensively re-reviewed"}
        rows.append(row)
        detail.append(dict(row, template_path=str(template.relative_to(root)) if permitted else None,
            metadata_notes=v["notes"], bindings=routes, source_crosswalk=cw, source_links=source_links,
            learning_links=links))
    with (out / "matrix.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (out / "matrix.json").write_text(json.dumps(detail, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (out / "review-notes.json").write_text(json.dumps(notes, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    summary.update(matrix_complete=True, variants=len(rows),
        primary_variants=sum(r["route"] == "primary" for r in rows),
        override_only_variants=sum(r["route"] == "override_only" for r in rows),
        rejected_only_variants=sum(r["route"] == "rejected_only" for r in rows),
        executable_templates=sum(bool(r["template_sha256"]) for r in rows),
        variants_without_learning_links=[r["variant"] for r in rows if not r["learning_resources"]],
        variants_without_source_crosswalk=[r["variant"] for r in rows if not r["source_crosswalk_ids"]],
        variants_with_learning_cautions=sum(bool(r["learning_known_issues"]) for r in rows),
        historical_references_without_online_reverification=sum("not re-verified online" in r.get("verification", "") for r in references.values()),
        source_versions={"rules": rules["rules_version"], "variants": variants["variants_version"]},
        source_sha256={str(p.relative_to(root)): sha(p) for p in
            (variant_path, rules_path, resources_path, crosswalk_path,
             root / "inst/advisor/uzh-tree.json", root / "inst/advisor/descrtab2-tree.json",
             Path(__file__).resolve())},
        notes_sha256=sha(notes_path))
except Exception as exc:
    summary["error"] = str(exc)
    raise
finally:
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
print(json.dumps({k: summary[k] for k in ("matrix_complete", "variants", "executable_templates", "scientific_acceptance")}))
