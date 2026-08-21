#!/usr/bin/env python3
"""Check the generated static site before deployment."""
from __future__ import annotations

import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
from urllib.parse import urlparse

SITE_ROOT = Path(__file__).resolve().parent
SITE = SITE_ROOT / "dist"
CONTENT = SITE_ROOT / "content"
GENERATED_PAGES = ["index.html", "methodology.html", "data.html", "cofolding.html"]
PAGES = GENERATED_PAGES + [
    "classification-submissions.html",
    "regression-submissions.html",
]
CHALLENGE_URL = "https://huggingface.co/spaces/openadmet/cyp-challenge"
SUMMARY_RESOURCES = {
    "regression-submissions.html",
    "classification-submissions.html",
    "evaluation-dashboard/",
}
GENERATED_SOURCES = {
    "index.html": CONTENT / "summary.md",
    "methodology.html": CONTENT / "methodology.md",
    "data.html": CONTENT / "data.md",
    "cofolding.html": CONTENT / "cofolding.md",
}


class Document(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []
        self.stylesheets: list[str] = []
        self.scripts: list[str] = []
        self.active_nav = 0
        self.doi_links = 0
        self.mermaid_blocks = 0
        self.source_items = 0
        self._in_sources = False

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        data = dict(attrs)
        classes = set((data.get("class") or "").split())
        if tag == "a":
            href = data.get("href")
            if not href:
                return
            self.links.append(href)
            if "active" in classes:
                self.active_nav += 1
            if "doi.org" in href:
                self.doi_links += 1
        elif tag == "link" and data.get("rel") == "stylesheet":
            href = data.get("href")
            if href:
                self.stylesheets.append(href)
        elif tag == "script":
            src = data.get("src")
            if src:
                self.scripts.append(src)
        elif tag == "pre" and "mermaid" in classes:
            self.mermaid_blocks += 1
        elif tag == "section" and data.get("id") == "sources":
            self._in_sources = True
        elif tag == "li" and self._in_sources:
            self.source_items += 1

    def handle_endtag(self, tag: str) -> None:
        if tag == "section" and self._in_sources:
            self._in_sources = False


def local_target(page: Path, href: str) -> Path | None:
    parsed = urlparse(href)
    if parsed.scheme or parsed.netloc or href.startswith(("mailto:", "#")):
        return None
    target = (page.parent / parsed.path).resolve()
    if parsed.path.endswith("/"):
        target /= "index.html"
    return target


def expected_source_digest(source: Path) -> str:
    digest = hashlib.sha256()
    digest.update(source.read_bytes())
    digest.update(b"\0")
    digest.update((CONTENT / "site.toml").read_bytes())
    return digest.hexdigest()


def main() -> int:
    failures: list[str] = []
    documents: dict[str, Document] = {}
    for name in PAGES:
        path = SITE / name
        if not path.is_file():
            failures.append(f"missing page: {name}")
            continue
        text = path.read_text()
        doc = Document()
        doc.feed(text)
        documents[name] = doc
        if name in GENERATED_PAGES:
            source = GENERATED_SOURCES[name]
            provenance = (
                f"generated-from: {source.relative_to(SITE_ROOT.parent).as_posix()} "
                f"sha256={expected_source_digest(source)}"
            )
            if provenance not in text:
                failures.append(f"{name}: stale relative to canonical Markdown")
            if '<meta name="robots" content="noindex,nofollow">' not in text:
                failures.append(f"{name}: missing robots exclusion")
            if doc.active_nav != 1:
                failures.append(f"{name}: expected one active nav link, found {doc.active_nav}")
        for href in doc.links + doc.stylesheets + doc.scripts:
            target = local_target(path, href)
            if target is not None and not target.is_file():
                failures.append(f"{name}: broken local link {href}")

    method_text = (SITE / "methodology.html").read_text()
    method = documents.get("methodology.html")
    if method:
        if method.doi_links != 16:
            failures.append(f"methodology.html: expected 16 DOI links, found {method.doi_links}")
        if method.source_items != 16:
            failures.append(f"methodology.html: expected 16 source rows, found {method.source_items}")
        if method.mermaid_blocks != 1:
            failures.append(f"methodology.html: expected one Mermaid block, found {method.mermaid_blocks}")
    if "Syntax error in text" in method_text or "binding affinities. T<" in method_text:
        failures.append("methodology.html: stale rendering or dangling source typo")

    if "cofolding.html" not in documents.get("data.html", Document()).links:
        failures.append("data.html: missing cofolding subpage link")
    if CHALLENGE_URL not in documents.get("index.html", Document()).links:
        failures.append("index.html: missing official OpenADMET challenge link")
    if CHALLENGE_URL not in documents.get("methodology.html", Document()).links:
        failures.append("methodology.html: missing official OpenADMET challenge link")
    summary_links = set(documents.get("index.html", Document()).links)
    for resource in sorted(SUMMARY_RESOURCES - summary_links):
        failures.append(f"index.html: missing challenge resource link {resource}")
    summary_text = (SITE / "index.html").read_text()
    aside_start = summary_text.find('<aside class="submission-box">')
    aside_end = summary_text.find("</aside>", aside_start)
    current_focus = summary_text.find('id="current-focus"')
    if aside_start < 0 or aside_end < 0:
        failures.append("index.html: submission logs are not in the summary aside")
    else:
        if "<h2" in summary_text[aside_start:aside_end]:
            failures.append("index.html: submission box must not be a text-body subheading")
        if current_focus < 0 or aside_start > current_focus:
            failures.append("index.html: submission box is not above the text body")
    hf = "https://huggingface.co/datasets/dargason/ADMET-CYP-cofolding"
    if hf not in documents.get("cofolding.html", Document()).links:
        failures.append("cofolding.html: missing Hugging Face dataset link")

    dashboard = SITE / "evaluation-dashboard"
    for relative in ("index.html", "data.js", "plotly.min.js", "manifest.json"):
        if not (dashboard / relative).is_file():
            failures.append(f"evaluation dashboard: missing {relative}")
    dashboard_manifest = dashboard / "manifest.json"
    if dashboard_manifest.is_file():
        dashboard_data = json.loads(dashboard_manifest.read_text())
        if len(dashboard_data.get("endpoints", {})) != 8:
            failures.append("evaluation dashboard: expected eight endpoints")
        if any(str(value).startswith("/") for value in (
            dashboard_data.get("source", ""),
            dashboard_data.get("challenge_source", ""),
        )):
            failures.append("evaluation dashboard: contains machine-local source path")

    for name in PAGES:
        text = (SITE / name).read_text().lower()
        if "stub" in text or "under development" in text:
            failures.append(f"{name}: contains placeholder wording")

    if failures:
        print("Site check: FAIL")
        for failure in failures:
            print(f"- {failure}")
        return 1
    print("Site check: PASS (6 pages, 16 sources, local links resolved)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
