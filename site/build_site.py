#!/usr/bin/env python3
"""Build the static ADMET-CYP site from canonical Markdown copy."""
from __future__ import annotations

import hashlib
import html
import re
import shutil
import subprocess
import tomllib
from datetime import datetime, timezone
from pathlib import Path

SITE_ROOT = Path(__file__).resolve().parent
PROJECT = SITE_ROOT.parent
CONTENT = SITE_ROOT / "content"
ASSETS = SITE_ROOT / "assets"
STATIC = SITE_ROOT / "static"
DASHBOARD = PROJECT / "cyp_sar_dashboard" / "dist"
OUTPUT = SITE_ROOT / "dist"
STAGING = SITE_ROOT / ".dist.building"
PREVIOUS = SITE_ROOT / ".dist.previous"
SITE_CONFIG = CONTENT / "site.toml"


class BuildError(RuntimeError):
    pass


def updated() -> tuple[str, str]:
    now = datetime.now(timezone.utc).replace(second=0, microsecond=0)
    return now.strftime("%Y-%m-%dT%H:%MZ"), now.strftime("%Y-%m-%d %H:%M UTC")


def site_config() -> dict:
    return tomllib.loads(SITE_CONFIG.read_text())


def content_page(name: str) -> tuple[dict, str, Path]:
    path = CONTENT / f"{name}.md"
    text = path.read_text()
    if not text.startswith("+++\n"):
        raise BuildError(f"{path}: missing TOML front matter")
    marker = "\n+++\n"
    end = text.find(marker, 4)
    if end < 0:
        raise BuildError(f"{path}: unterminated TOML front matter")
    metadata = tomllib.loads(text[4:end])
    required = {"title", "description", "active"}
    missing = required - metadata.keys()
    if missing:
        raise BuildError(f"{path}: missing front matter: {', '.join(sorted(missing))}")
    return metadata, text[end + len(marker) :], path


def extract_block(source: str, name: str, *, required: bool = True) -> tuple[str, str]:
    pattern = re.compile(
        rf"\n?<!-- {re.escape(name)}-start -->\n(.*?)\n<!-- {re.escape(name)}-end -->\n?",
        flags=re.DOTALL,
    )
    match = pattern.search(source)
    if match is None:
        if required:
            raise BuildError(f"missing Markdown block: {name}")
        return "", source
    return match.group(1).strip() + "\n", source[: match.start()] + source[match.end() :]


def methodology_source(source: str) -> str:
    source = re.sub(
        r"^Daniel Saltzberg - @dargason\n(Published: [^\n]+)\n(Updated: [^\n]+)$",
        r'<p class="byline">Daniel Saltzberg - @dargason<br>\1<br>\2</p>',
        source,
        count=1,
        flags=re.MULTILINE,
    )
    source_head, separator, source_tail = source.partition("## Sources\n")
    if not separator:
        raise BuildError("methodology Markdown has no Sources section")
    source_lines = [
        f"- {line}" if re.match(r"^\[\d+\] https?://", line) else line
        for line in source_tail.splitlines()
    ]
    return source_head + separator + "\n".join(source_lines) + "\n"


def render_markdown(source: str, *, mermaid: bool = False) -> str:
    placeholder = "MERMAID_DIAGRAM_PLACEHOLDER"
    rendered_diagram = ""
    if mermaid:
        diagram_match = re.search(r"```mermaid\n(.*?)\n```", source, flags=re.DOTALL)
        if diagram_match is None:
            raise BuildError("methodology Markdown has no Mermaid diagram")
        diagram = html.escape(diagram_match.group(1).strip())
        source = source[: diagram_match.start()] + placeholder + source[diagram_match.end() :]
        rendered_diagram = (
            '<div class="diagram-shell"><pre class="mermaid" '
            f'aria-label="Bayesian claim-evidence framework flowchart">{diagram}</pre></div>'
        )

    result = subprocess.run(
        ["pandoc", "--from=gfm+raw_html", "--to=html5", "--section-divs"],
        input=source,
        text=True,
        capture_output=True,
        check=True,
    )
    rendered = result.stdout
    if mermaid:
        rendered_placeholder = f"<p>{placeholder}</p>"
        if rendered_placeholder not in rendered:
            raise BuildError("Pandoc did not preserve the Mermaid placeholder")
        rendered = rendered.replace(rendered_placeholder, rendered_diagram, 1)
    return rendered.strip()


def render_inline(source: str) -> str:
    rendered = render_markdown(source).strip()
    match = re.fullmatch(r"<p>(.*)</p>", rendered, flags=re.DOTALL)
    return match.group(1) if match else rendered


def add_lede(rendered: str) -> str:
    return re.sub(r"(</h1>\s*)<p>", r'\1<p class="lede">', rendered, count=1)


def source_digest(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    digest.update(b"\0")
    digest.update(SITE_CONFIG.read_bytes())
    return digest.hexdigest()


def asset_version(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()[:12]


def nav(active: str, config: dict) -> str:
    rendered = []
    for link in config["nav"]:
        active_attr = ' class="active"' if link["key"] == active else ""
        rendered.append(
            f'<a href="{html.escape(link["href"], quote=True)}"{active_attr}>'
            f'{html.escape(link["label"])}</a>'
        )
    return "".join(rendered)


def page(
    *,
    metadata: dict,
    body: str,
    footer: str,
    timestamp: tuple[str, str],
    config: dict,
    source_path: Path,
) -> str:
    machine_time, display_time = timestamp
    page_class = metadata.get("page_class", "")
    class_attr = f' class="{html.escape(page_class, quote=True)}"' if page_class else ""
    source_relative = source_path.relative_to(PROJECT).as_posix()
    provenance = f"generated-from: {source_relative} sha256={source_digest(source_path)}"
    return f'''<!doctype html>
<!-- {provenance} -->
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta name="robots" content="noindex,nofollow"><meta name="description" content="{html.escape(metadata["description"], quote=True)}">
<title>{html.escape(metadata["title"])} · CYP claim-evidence framework</title><link rel="stylesheet" href="assets/style.css?v={asset_version(ASSETS / "style.css")}"></head>
<body><header><div class="masthead"><div class="brand">{html.escape(config["brand"])}</div><p class="updated-at">Updated <time datetime="{machine_time}">{display_time}</time></p></div><nav>{nav(metadata["active"], config)}</nav></header>
<main{class_attr}>
{body}
</main><footer>{footer}</footer></body></html>
'''


def copy_public_support() -> None:
    required_dashboard = ("index.html", "data.js", "plotly.min.js", "manifest.json")
    missing = [name for name in required_dashboard if not (DASHBOARD / name).is_file()]
    if missing:
        raise BuildError("dashboard build is incomplete: " + ", ".join(missing))
    if not ASSETS.is_dir() or not STATIC.is_dir():
        raise BuildError("site assets or static source directory is missing")

    shutil.rmtree(STAGING, ignore_errors=True)
    STAGING.mkdir(parents=True)
    shutil.copytree(ASSETS, STAGING / "assets")
    shutil.copytree(DASHBOARD, STAGING / "evaluation-dashboard")
    for source in STATIC.iterdir():
        destination = STAGING / source.name
        if source.is_dir():
            shutil.copytree(source, destination)
        else:
            shutil.copy2(source, destination)


def publish_staging() -> None:
    shutil.rmtree(PREVIOUS, ignore_errors=True)
    if OUTPUT.exists():
        OUTPUT.rename(PREVIOUS)
    try:
        STAGING.rename(OUTPUT)
    except Exception:
        if PREVIOUS.exists() and not OUTPUT.exists():
            PREVIOUS.rename(OUTPUT)
        raise
    shutil.rmtree(PREVIOUS, ignore_errors=True)


def build() -> None:
    config = site_config()
    timestamp = updated()
    copy_public_support()
    try:
        summary_meta, summary_source, summary_path = content_page("summary")
        hero, summary_source = extract_block(summary_source, "hero")
        submissions, summary_source = extract_block(summary_source, "submission")
        summary_footer, summary_source = extract_block(summary_source, "footer")
        hero_html = add_lede(render_markdown(hero))
        submission_html = render_markdown(submissions)
        summary_body = (
            '<div class="summary-hero">\n'
            f'<div class="summary-copy">{hero_html}</div>\n'
            f'<aside class="submission-box">{submission_html}</aside>\n'
            '</div>\n'
            + render_markdown(summary_source)
        )
        (STAGING / "index.html").write_text(
            page(
                metadata=summary_meta,
                body=summary_body,
                footer=render_inline(summary_footer),
                timestamp=timestamp,
                config=config,
                source_path=summary_path,
            )
        )

        for name, output_name in (("data", "data.html"), ("cofolding", "cofolding.html")):
            metadata, source, source_path = content_page(name)
            footer, source = extract_block(source, "footer")
            body = add_lede(render_markdown(source))
            (STAGING / output_name).write_text(
                page(
                    metadata=metadata,
                    body=body,
                    footer=render_inline(footer),
                    timestamp=timestamp,
                    config=config,
                    source_path=source_path,
                )
            )

        method_meta, method_source, method_path = content_page("methodology")
        method_footer, method_source = extract_block(method_source, "footer")
        method_body = (
            '<div class="eyebrow">In progress · framework proposal</div>\n'
            '<article class="document">\n'
            + render_markdown(methodology_source(method_source), mermaid=True)
            + '\n</article>\n'
            '<script src="assets/mermaid.min.js" defer></script>\n'
            "<script>window.addEventListener('DOMContentLoaded', function () { mermaid.initialize({startOnLoad:false, theme:'base', securityLevel:'strict', themeVariables:{fontFamily:'system-ui, sans-serif', primaryColor:'#f2f0ea', primaryBorderColor:'#8e8a82', lineColor:'#66635f', clusterBkg:'#fffffb', clusterBorder:'#d8d4cc'}}); mermaid.run({querySelector:'.mermaid'}); });</script>"
        )
        (STAGING / "methodology.html").write_text(
            page(
                metadata=method_meta,
                body=method_body,
                footer=render_inline(method_footer),
                timestamp=timestamp,
                config=config,
                source_path=method_path,
            )
        )
        publish_staging()
    except Exception:
        shutil.rmtree(STAGING, ignore_errors=True)
        raise


if __name__ == "__main__":
    build()
