#!/usr/bin/env python3
"""Build the static OpenADMET CYP project pages."""
from __future__ import annotations

import html
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

SITE = Path(__file__).resolve().parent
PROJECT = SITE.parent
METHODOLOGY_SOURCE = PROJECT / "intro" / "BAYESIAN_OUTCOME_PREREGISTRATION.md"


def updated() -> tuple[str, str]:
    now = datetime.now(timezone.utc).replace(second=0, microsecond=0)
    return now.strftime("%Y-%m-%dT%H:%MZ"), now.strftime("%Y-%m-%d %H:%M UTC")


def nav(active: str) -> str:
    links = [
        ("summary", "index.html", "Summary"),
        ("methodology", "methodology.html", "Methodology"),
        ("data", "data.html", "Data"),
    ]
    rendered = []
    for key, href, label in links:
        active_attr = ' class="active"' if key == active else ""
        rendered.append(f'<a href="{href}"{active_attr}>{label}</a>')
    return "".join(rendered)


def page(*, title: str, description: str, active: str, body: str, footer: str, timestamp: tuple[str, str], page_class: str = "") -> str:
    machine_time, display_time = timestamp
    class_attr = f' class="{page_class}"' if page_class else ""
    return f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta name="robots" content="noindex,nofollow"><meta name="description" content="{html.escape(description, quote=True)}">
<title>{html.escape(title)} · OpenADMET CYP</title><link rel="stylesheet" href="assets/style.css"></head>
<body><header><div class="masthead"><div class="brand">OpenADMET · CYP Challenge</div><p class="updated-at">Updated <time datetime="{machine_time}">{display_time}</time></p></div><nav>{nav(active)}</nav></header>
<main{class_attr}>
{body}
</main><footer>{footer}</footer></body></html>
'''


def methodology_html() -> str:
    source = METHODOLOGY_SOURCE.read_text()
    source = source.replace(
        "Daniel Saltzberg - @dargason\nVersion 0.1, 2026-08-10",
        '<p class="byline">Daniel Saltzberg · @dargason<br>Version 0.1, 2026-08-10</p>',
        1,
    )
    source_head, source_sep, source_tail = source.partition("## Sources\n")
    if not source_sep:
        raise RuntimeError("Methodology source has no Sources section")
    source_lines = [f"- {line}" if re.match(r"^\[\d+\] https?://", line) else line for line in source_tail.splitlines()]
    source = source_head + source_sep + "\n".join(source_lines)

    diagram_match = re.search(r"```mermaid\n(.*?)\n```", source, flags=re.DOTALL)
    if diagram_match is None:
        raise RuntimeError("Methodology source has no Mermaid diagram")
    diagram = html.escape(diagram_match.group(1).strip())
    source = source[: diagram_match.start()] + "MERMAID_DIAGRAM_PLACEHOLDER" + source[diagram_match.end() :]
    result = subprocess.run(
        ["pandoc", "--from=gfm+raw_html", "--to=html5", "--section-divs"],
        input=source,
        text=True,
        capture_output=True,
        check=True,
    )
    placeholder = "<p>MERMAID_DIAGRAM_PLACEHOLDER</p>"
    rendered_diagram = f'<div class="diagram-shell"><pre class="mermaid" aria-label="Bayesian claim-evidence framework flowchart">{diagram}</pre></div>'
    if placeholder not in result.stdout:
        raise RuntimeError("Pandoc did not preserve the Mermaid placeholder")
    return result.stdout.replace(placeholder, rendered_diagram, 1)


def build() -> None:
    timestamp = updated()
    summary_body = '''<div class="eyebrow">Project summary · stub</div>
<h1>OpenADMET CYP challenge modeling</h1>
<p class="lede">Work on the OpenADMET CYP challenge is in progress. The methodology and first data release are linked below.</p>
<h2>Project pages</h2>
<div class="page-grid">
  <a class="page-link" href="methodology.html"><strong>Methodology</strong><span>Read the Bayesian claim-evidence plan.</span></a>
  <a class="page-link" href="data.html"><strong>Data</strong><span>Browse the CYP structure release.</span></a>
</div>'''
    (SITE / "index.html").write_text(page(
        title="Summary",
        description="Summary of the OpenADMET CYP challenge modeling project.",
        active="summary",
        body=summary_body,
        footer='Project summary under development. <a href="methodology.html">Read the methodology</a>.',
        timestamp=timestamp,
    ))

    method_body = '<div class="eyebrow">Methodology · preregistration</div>\n<article class="document">\n' + methodology_html() + '''\n</article>
<script src="assets/mermaid.min.js" defer></script>
<script>window.addEventListener('DOMContentLoaded', function () { mermaid.initialize({startOnLoad:false, theme:'base', securityLevel:'strict', themeVariables:{fontFamily:'system-ui, sans-serif', primaryColor:'#f2f0ea', primaryBorderColor:'#8e8a82', lineColor:'#66635f', clusterBkg:'#fffffb', clusterBorder:'#d8d4cc'}}); mermaid.run({querySelector:'.mermaid'}); });</script>'''
    (SITE / "methodology.html").write_text(page(
        title="Methodology",
        description="Bayesian claim-evidence framework preregistration for the OpenADMET CYP challenge.",
        active="methodology",
        body=method_body,
        footer='<a href="index.html">Back to summary</a> · Source document: <a href="https://github.com/saltzberg/OpenADMET-CYP/blob/main/intro/BAYESIAN_OUTCOME_PREREGISTRATION.md">GitHub</a>.',
        timestamp=timestamp,
        page_class="methodology-page",
    ))

    data_body = '''<div class="eyebrow">Project data</div>
<h1>Data releases</h1>
<p class="lede">This page will keep track of project data. There is one release so far.</p>
<h2>Available now</h2>
<div class="page-grid">
  <a class="page-link" href="cofolding.html"><strong>CYP cofolding structures</strong><span>149 experimental structures and 16,590 predictions, with QC results, checksums, and provenance.</span></a>
</div>'''
    (SITE / "data.html").write_text(page(
        title="Data",
        description="Data releases for the OpenADMET CYP challenge modeling project.",
        active="data",
        body=data_body,
        footer='<a href="index.html">Back to summary</a>.',
        timestamp=timestamp,
    ))

    cofold_body = '''<div class="eyebrow">Data · cofolding structures</div>
<h1>CYP cofolding structures</h1>
<p class="lede">We generated an aligned structure set to support structural feature work.</p>
<div class="metrics">
  <div class="metric"><strong>16,739</strong><span>structure rows</span></div>
  <div class="metric"><strong>149</strong><span>experimental structures</span></div>
  <div class="metric"><strong>16,590</strong><span>predicted structures</span></div>
  <div class="metric"><strong>5</strong><span>prediction methods</span></div>
</div>
<h2>Why we made it</h2>
<p>We want to test whether structural features add useful information to the CYP models. The files are aligned to make those comparisons easier.</p>
<div class="callout pass"><div class="label">Checks</div><p>Every structure passed the release parser. The release table, paths, and hashes also passed their checks.</p></div>
<h2>Get the data</h2>
<div class="page-grid single-card">
  <a class="page-link external-link" href="https://huggingface.co/datasets/dargason/ADMET-CYP-cofolding"><strong>Open the dataset on Hugging Face</strong><span>Dataset card, Parquet table, aligned coordinates, QC fields, checksums, and provenance.</span></a>
</div>
<h2>Interpretation</h2>
<p>These structures are hypotheses, not inhibition measurements. They also do not describe the full mechanism of time-dependent inhibition.</p>'''
    (SITE / "cofolding.html").write_text(page(
        title="Cofolding data",
        description="CYP cofolding structures and provenance for the OpenADMET CYP challenge project.",
        active="data",
        body=cofold_body,
        footer='<a href="data.html">Back to data</a> · <a href="index.html">Project summary</a>.',
        timestamp=timestamp,
    ))


if __name__ == "__main__":
    build()
