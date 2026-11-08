"""Publication-only ASGI adapter: version HTML asset URLs past stale CDN entries.
Reviewed application files and their runtime/scientific identity remain unchanged.
"""
import os
import re
from backend.runtime import app as reviewed_app, SID, PREFIX

SOURCE_COMMIT = os.environ.get('CYP_SOURCE_COMMIT', '')
if SOURCE_COMMIT and not re.fullmatch(r'[0-9a-f]{40}', SOURCE_COMMIT):
    raise ValueError('CYP_SOURCE_COMMIT must be a full Git commit ID')

async def app(scope, receive, upstream_send):
    original_send = upstream_send
    async def send(message):
        if message['type'] == 'http.response.start' and SOURCE_COMMIT:
            message = dict(message, headers=[*message['headers'],
                           (b'x-source-commit', SOURCE_COMMIT.encode('ascii'))])
        await original_send(message)
    if scope['type'] != 'http' or scope.get('path') not in {PREFIX, PREFIX+'index.html'}:
        return await reviewed_app(scope, receive, send)
    start = None
    chunks = []
    async def capture(message):
        nonlocal start
        if message['type'] == 'http.response.start':
            start = message
        elif message['type'] == 'http.response.body':
            chunks.append(message.get('body', b''))
            if not message.get('more_body', False):
                body = b''.join(chunks)
                if start['status'] == 200:
                    text = body.decode('utf-8')
                    text = re.sub(r'((?:src|href)=\")('+re.escape(PREFIX)+r'static/[^\"?]+)(\")', lambda m: m[1]+m[2]+'?v='+SID+m[3], text)
                    body = text.encode('utf-8')
                    headers = [(k,v) for k,v in start['headers'] if k.lower() not in {b'content-length', b'etag', b'cache-control'}]
                    headers += [(b'content-length', str(len(body)).encode()),(b'cache-control',b'no-store')]
                    start = dict(start, headers=headers)
                await send(start)
                await send(dict(type='http.response.body',body=body))
    await reviewed_app(scope, receive, capture)
