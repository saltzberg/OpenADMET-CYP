#!/usr/bin/env python3
"""Real-browser integration check. Requires a built release and Playwright Chromium."""
import argparse
import asyncio
import hashlib
import json
from pathlib import Path
from playwright.async_api import async_playwright

PREFIX = '/articles/admet-cyp-structure/cofolding/'


async def check(args):
    root = args.runtime_root.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    url = args.origin.rstrip('/') + PREFIX
    models = json.loads((root / 'manifest.json').read_text())['final20']
    if not args.all_models:
        models = models[:1]
    errors, http_errors, records = [], [], []
    async with async_playwright() as pw:
        browser = await pw.chromium.launch(headless=True, args=['--no-sandbox', '--use-gl=angle', '--use-angle=swiftshader', '--enable-unsafe-swiftshader'])
        page = await browser.new_page(viewport={'width': 1500, 'height': 1000})
        page.on('pageerror', lambda error: errors.append(str(error)))
        page.on('response', lambda response: http_errors.append([response.status, response.url]) if response.status >= 400 else None)
        try:
            await page.goto(url)
            await page.wait_for_function('poseExplorer?.viewerDiagnostics?.readiness === "ready" && poseExplorer.viewerDiagnostics.roots.length > 0 && !poseExplorer.viewerDiagnostics.activeCount && !poseExplorer.viewerDiagnostics.queuedCount', timeout=120000)
            assert await page.locator('#showRevealed').evaluate('e => e.parentElement.textContent.trim()') == 'Experimental'
            assert 'Download submission' not in await page.locator('header').inner_text()
            for pid, state, expected in models:
                lid = pid.split('/')[1]
                await page.evaluate('async lid => { if (poseExplorer.state.ligand !== lid) await poseExplorer.setLigand(lid); }', lid)
                await page.evaluate('pid => poseExplorer.selectPose(pid)', pid)
                await page.wait_for_function('sha => poseExplorer.viewerDiagnostics?.roots.some(r => r.sha256 === sha) && poseExplorer.viewerDiagnostics.readiness === "ready" && !poseExplorer.viewerDiagnostics.activeCount && !poseExplorer.viewerDiagnostics.queuedCount', arg=expected, timeout=120000)
                strings = await page.evaluate('() => [...poseExplorer.viewer.plugin.state.data.cells.values()].filter(c => c.obj?.type?.name === "String Data").map(c => c.obj.data)')
                hashes = [hashlib.sha256(text.encode()).hexdigest() for text in strings]
                assert expected in hashes, (pid, hashes)
                assert await page.locator('#molecule svg').count() > 0
                records.append({'id': pid, 'state': state, 'sha256': expected, 'actual_molstar_hashes': hashes})
                print(pid, flush=True)
            await page.screenshot(path=str(args.output / 'dashboard.png'))
            assert not errors and not http_errors, (errors, http_errors)
            result = {'status': 'PASS', 'url': url, 'models': records, 'page_errors': errors, 'http_errors': http_errors}
            (args.output / 'browser-results.json').write_text(json.dumps(result, indent=2) + '\n')
        finally:
            await browser.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtime-root', required=True, type=Path)
    parser.add_argument('--origin', required=True)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--all-models', action='store_true')
    asyncio.run(check(parser.parse_args()))


if __name__ == '__main__':
    main()
