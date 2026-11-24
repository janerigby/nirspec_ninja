"""Shared archive operations; entry points retain their scientific policies."""
import json
import os

def table_chunks(table, size):
    for start in range(0, len(table), size):
        yield table[start:start + size]

def configure_cloud_dataset(observations, source, profile=None):
    if source == 'MAST':
        observations.disable_cloud_dataset()
    elif source == 'AWS':
        observations.enable_cloud_dataset(provider='AWS', profile=profile)

def atomic_json(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    os.replace(temporary, path)
