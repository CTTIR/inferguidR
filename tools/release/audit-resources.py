#!/usr/bin/env python3
"""Explicit online availability audit of the learning-resource catalog."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import datetime
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urldefrag, unquote, urlsplit
from urllib.request import Request, urlopen


class Page(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids, self.title, self.in_title = set(), [], False

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if attrs.get("id"):
            self.ids.add(attrs["id"])
        if tag == "a" and attrs.get("name"):
            self.ids.add(attrs["name"])
        if tag == "title":
            self.in_title = True

    def handle_endtag(self, tag):
        if tag == "title":
            self.in_title = False

    def handle_data(self, data):
        if self.in_title:
            self.title.append(data)


def inspect(resource):
    url, anchor = urldefrag(resource["published_url"])
    row = {"resource_id": resource["resource_id"], "url": resource["published_url"],
           "checked_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
           "pinned_commit": resource.get("commit"), "source_url": resource.get("repo_url"),
           "http_status": None, "anchor_present": None, "title": None,
           "previous_title": resource.get("published_page_title"), "error": None}
    try:
        request = Request(url, headers={"User-Agent": "inferguidR-resource-audit/1.0"})
        with urlopen(request, timeout=20) as response:
            payload = response.read(8 * 1024 * 1024 + 1)
            row.update(http_status=response.status, final_url=response.url)
        if len(payload) > 8 * 1024 * 1024:
            raise ValueError("Page exceeds the audit size limit")
        page = Page()
        page.feed(payload.decode("utf-8", errors="replace"))
        row.update(title=" ".join("".join(page.title).split()),
                   content_sha256=hashlib.sha256(payload).hexdigest(),
                   anchor_present=(unquote(anchor) in page.ids) if anchor else None)
    except HTTPError as error:
        row.update(http_status=error.code, error=str(error))
    except (OSError, URLError, ValueError) as error:
        row["error"] = str(error)
    row["available"] = row["http_status"] == 200 and row["error"] is None and row["anchor_present"] is not False
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("catalog", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists() and (not args.output.is_dir() or any(args.output.iterdir())):
        parser.error("output must be a new or empty directory")
    args.output.mkdir(parents=True, exist_ok=True)
    payload = None
    try:
        payload = args.catalog.read_bytes()
        resources = json.loads(payload)
        if not isinstance(resources, list) or not resources:
            raise ValueError("catalog must contain at least one resource")
        for resource in resources:
            if not isinstance(resource, dict) or not all(
                isinstance(resource.get(key), str) and resource[key].strip()
                for key in ("resource_id", "published_url")
            ):
                raise ValueError("each resource needs nonempty resource_id and published_url strings")
            url = urlsplit(resource["published_url"])
            if url.scheme not in {"http", "https"} or not url.netloc:
                raise ValueError("published_url must be an HTTP(S) URL")
        if len({r["resource_id"] for r in resources}) != len(resources):
            raise ValueError("catalog has duplicate resource identifiers")
    except (OSError, ValueError, UnicodeError) as error:
        receipt = {"status": "FAIL", "resources": 0, "error": str(error),
                   "catalog_sha256": hashlib.sha256(payload).hexdigest() if payload is not None else None}
        (args.output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
        print(json.dumps(receipt))
        return 1
    rows = []
    with (args.output / "results.jsonl").open("w") as stream:
        with ThreadPoolExecutor(max_workers=2) as pool:
            for row in pool.map(inspect, resources):
                rows.append(row)
                stream.write(json.dumps(row) + "\n")
                stream.flush()
    receipt = {
        "scope": "current HTTP availability and exact fragment presence; not semantic or source-content certification",
        "catalog_sha256": hashlib.sha256(payload).hexdigest(),
        "resources": len(rows), "available": sum(r["available"] for r in rows),
        "unavailable": sum(not r["available"] for r in rows),
        "title_changes": sum(r["title"] is not None and r["title"] != r["previous_title"] for r in rows),
        "status": "PASS" if all(r["available"] for r in rows) else "FAIL",
    }
    (args.output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt))
    return receipt["unavailable"] > 0


if __name__ == "__main__":
    raise SystemExit(main())
