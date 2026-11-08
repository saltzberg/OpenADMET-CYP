# Cofolding dashboard

[Public viewer](https://aetherark.com/articles/admet-cyp-structure/cofolding/).
This directory contains its Python API, browser code, vendor assets, tests and release tools.
It does not contain the scientific dataset.

## Code and data

`runtime-data.manifest.json` specifies two input groups:

- `code_files`: exact runtime code and assets in this checkout.
- `files`: external scientific inputs, with relative paths, sizes and SHA-256 hashes.

An operator must provide an existing runtime-data directory and the original
coordinate store with its coordinate manifest. No public download is supplied.
The manifest does not grant permission to release those inputs. The builder does
not download, compute, register or change scientific data. It copies verified
inputs into a new release and uses the coordinate store read-only.

The frozen source manifest, score packets, reference packets, 31-residue receptor
map, H-bond shards, depictions and 20 experimental PDBs stay external. Logs,
caches, private source paths and historical receipts are not Git inputs.
External packet bytes are retained, including their historical provenance;
the API removes private paths from responses and serves only allowed routes.
Do not expose the release directory with a general file server.

## Build and start

Use Python 3.12 with `requirements.txt`. Use an existing compatible environment,
or install these requirements in a separate virtual environment. Node.js is
needed for JavaScript tests; Playwright with Chromium is needed for browser tests.

Run from `dashboard/`. Set these paths to operator-provided directories:

```bash
export PYTHON=/path/to/venv/bin/python
export DATA_ROOT=/path/to/existing-runtime-data
export COORDINATE_ROOT=/path/to/existing-coordinate-store
export COORDINATE_MANIFEST=/path/to/coordinate-manifest.json
export RELEASE=/path/to/new-release
export TMPDIR=/path/to/writable-scratch

"$PYTHON" build_runtime.py --data-root "$DATA_ROOT" \
  --coordinate-root "$COORDINATE_ROOT" \
  --coordinate-manifest "$COORDINATE_MANIFEST" --output "$RELEASE"
"$PYTHON" "$RELEASE/run_runtime.py" \
  --coordinate-root "$COORDINATE_ROOT" --port 8896
```

Open `http://127.0.0.1:8896/articles/admet-cyp-structure/cofolding/`.
The launcher binds to loopback and starts one CPU-only worker. Use `--host 0.0.0.0`
only when LAN access is intended. Keep numerical libraries single-threaded.

The output must not exist and must be outside this source directory. The builder
checks all 331 reviewed runtime files and all 42,213 compressed coordinates,
checks coordinate membership against the structure allowlist, then imports the
real API and checks its runtime identity. It publishes the completed directory
only after these checks pass. No response archive is used as an API substitute.
The 31-residue map path is resolved inside the release, not through an app alias.

`build-receipt.json` records file hashes, Git HEAD and dirty status. Only a clean
checkout produces an `X-Source-Commit` response header. The launcher checks the
receipt's file hashes before startup. The receipt is an operator record, not a
cryptographic signature. Runtime SID and source Git commit are distinct:
changing documentation does not change scientific or runtime identity.

Current reviewed SID:
`cc2e161dc06cc5e4d3de1dc2fdce7cc921421e0b01db640d6576de348093574a`.

## Checks

Run HTTP and browser checks one after the other on a small host. Full-population
requests run through a bounded worker queue; concurrent verification can exceed
its deadline. Write test output outside the source and served directories.

```bash
"$PYTHON" tests/test_packaging.py
"$PYTHON" build_assets.py --check
export CYP_TEST_APP="$RELEASE"
export CYP_COORDINATE_ROOT="$COORDINATE_ROOT"
export PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
"$PYTHON" tests/test_runtime_lane.py
node --test tests/test_*.js tests/test_*.cjs
"$PYTHON" verify_runtime.py --runtime-root "$RELEASE" \
  --origin http://127.0.0.1:8896 --all-ligands
/path/to/playwright-python tests/browser_integration.py \
  --runtime-root "$RELEASE" --origin http://127.0.0.1:8896 \
  --output /path/to/private-browser-check --all-models
```

Packaging tests use marked synthetic contracts. Runtime and browser integration
tests require the real external dataset. Worker tests copy runtime files and use
separate caches. The transport fixture also uses an isolated analysis cache.
The HTTP check verifies every served metric and receptor-map coordinate for
20 ligands in raw, relaxed and paired states, plus coordinate response hashes.
The browser check verifies the actual Mol* source hash and a 2D depiction for
all 20 selected originals. Inspect its screenshot to check visible rendering.

## Deployment

A Git push does not deploy this service. Build from the reviewed **clean commit**
on the target host, with its existing data store and Python environment. Do not
copy workstation caches or synchronize the whole repository into the web root.

1. Build a new immutable release with the command above.
2. Start it on an unused loopback port and run the HTTP and browser checks.
3. Stop only the existing dashboard service. Retain its previous release and
   service configuration for rollback.
4. Point its application alias at the new release. Set its working directory
   and executable to the resolved release. The service command is
   `PYTHON RELEASE/run_runtime.py --coordinate-root COORDINATE_ROOT --port 8896`,
   with the real absolute paths substituted. Start the existing service.
5. Verify the unchanged public URL with `verify_runtime.py --all-ligands` and
   `--source-commit FULL_COMMIT_ID`, then run the browser check at that origin.
   If verification fails, restore the previous alias and service configuration.

Retain the existing proxy route, sibling site files, coordinate store, memory
cap and service name. Do not replace global proxy configuration. On the current
small public host, retain one worker and the 750 MiB cap; this is not a claim of
unlimited concurrent capacity. Recheck service restarts, OOM events and memory
peak after the serialized checks.

`publication_entry.py` retains the public ASGI adapter. It versions HTML asset
URLs with the runtime SID and marks entry HTML `no-store`, to avoid stale CDN
assets. Its only added behavior is the optional source-commit header. It does
not enter the scientific SID or change API values.

## Scientific meaning

Each model row retains its model ID and coordinate state. Raw and relaxed
coordinates are distinct. Paired view contains both states only for models
with a relaxed counterpart: 40,426 rows, not all 42,213 state-specific inputs.
There are 22,000 raw rows and 20,213 original relaxed rows.

LDDT-PLI and BiSymRMSD compare the exact displayed coordinates with the
experimental reference. The UI name BiSymRMSD maps to source key `BiSyRMSD`.
All 20,213 original relaxed poses and 21,934 raw poses have score pairs; 66 raw
poses retain protocol no-assignment nulls. Missing values are not zero.
These geometry metrics, native confidence, energy and H-bond evidence do not
measure binding affinity or inhibition. The 31-residue receptor map is a saved
shared projection; viewing or filtering does not fit a new basis.

Submitted markers refer to exact original selected structures. Experimental
markers are separate observed references, not generated model rows. UI labels,
reference identities and coordinate hashes are retained without substitution.

## Vendor licenses

Mol* 4.17.0 and Plotly.js 3.5.0 are vendored byte-for-byte. Their MIT texts are in
`licenses/`; bundled dependency notices remain in
`static/vendor/molstar.js.LICENSE.txt`. `static/vendor/manifest.json` records the
reviewed asset hashes. See [third-party notices](THIRD_PARTY_NOTICES.md).
