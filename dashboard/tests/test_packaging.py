"""Synthetic packaging/ASGI contracts; these fixtures are not scientific results."""
import asyncio
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import build_runtime as build


class BuildContracts(unittest.TestCase):
    def test_relative_paths_and_hashes_fail_closed(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            for relative in ('../escape', '/absolute'):
                with self.assertRaises(ValueError):
                    build.safe_path(root, relative)
            (root / 'escape').symlink_to(root.parent)
            with self.assertRaises(ValueError):
                build.safe_path(root, 'escape/other')
            target = root / 'input.json'
            target.write_bytes(b'contract fixture')
            row = {'bytes': target.stat().st_size, 'sha256': build.digest(target)}
            build.verify_file(target, row)
            target.write_bytes(b'changed fixture')
            with self.assertRaises(ValueError):
                build.verify_file(target, row)

    def test_coordinate_manifest_hash_count_and_gzip_bytes(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / 'fixture.gz').write_bytes(b'opaque synthetic gzip fixture')
            row = {'path': 'fixture.gz', 'gzip_bytes': (root / 'fixture.gz').stat().st_size,
                   'gzip_sha256': build.digest(root / 'fixture.gz')}
            manifest = root / 'manifest.json'
            manifest.write_text(json.dumps({'files': [row]}))
            expected = {'coordinate_manifest_sha256': build.digest(manifest), 'coordinate_count': 1}
            self.assertEqual(build.verify_coordinates(root, manifest, expected), [row])
            with self.assertRaises(ValueError):
                build.verify_coordinates(root, manifest, dict(expected, coordinate_count=2))
            (root / 'fixture.gz').write_bytes(b'corrupt')
            with self.assertRaises(ValueError):
                build.verify_coordinates(root, manifest, expected)


class PublicationContracts(unittest.IsolatedAsyncioTestCase):
    async def test_versioning_and_source_commit_header(self):
        prefix = '/articles/admet-cyp-structure/cofolding/'
        async def source(scope, receive, send):
            await send({'type': 'http.response.start', 'status': 200,
                        'headers': [(b'etag', b'old'), (b'cache-control', b'old')]})
            await send({'type': 'http.response.body', 'body': ('<script src="' + prefix + 'static/app.js"></script>').encode()})
        fake = types.ModuleType('backend.runtime')
        fake.__dict__.update(app=source, SID='synthetic-sid', PREFIX=prefix)
        spec = importlib.util.spec_from_file_location('publication_contract', ROOT / 'publication_entry.py')
        assert spec is not None and spec.loader is not None
        module = importlib.util.module_from_spec(spec)
        with patch.dict(sys.modules, {'backend.runtime': fake}), patch.dict(os.environ, {'CYP_SOURCE_COMMIT': 'a' * 40}):
            spec.loader.exec_module(module)
        for path in (prefix, prefix + 'api/bootstrap'):
            messages = []
            async def send(message):
                messages.append(message)
            await module.app({'type': 'http', 'path': path}, None, send)
            headers = dict(messages[0]['headers'])
            self.assertEqual(headers[b'x-source-commit'], b'a' * 40)
            if path == prefix:
                self.assertIn(b'?v=synthetic-sid', messages[1]['body'])
                self.assertEqual(headers[b'cache-control'], b'no-store')
                self.assertNotIn(b'etag', headers)
            else:
                self.assertNotIn(b'?v=', messages[1]['body'])
        with patch.dict(sys.modules, {'backend.runtime': fake}), patch.dict(os.environ, {'CYP_SOURCE_COMMIT': 'not-a-commit'}):
            with self.assertRaises(ValueError):
                spec.loader.exec_module(module)


if __name__ == '__main__':
    unittest.main(verbosity=2)
