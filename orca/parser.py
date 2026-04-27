# Built-in
import inspect

from . import utils
from . import blocks

class OrcaParser:
    def __init__(self, file_path):
        self.file_path = utils._parse_path(file_path)
        self._lines = utils._read_lines(self.file_path)
    
        self._cache = set()
        self._behavior_overrides = {}
        self._ignore_warnings = False
    
        self._build_maps()

    def parse(self):
        for i, line in enumerate(self._lines):
            for tag, funcs in self._tag_map.items():
                if tag in line:
                    for func in funcs:
                        result = func(self._lines, i=i)
                        if result is not None:
                            self._store_result(func, result)

    def get_property(self, key):
        assert key in self._block_map, 'Invalid property key'
        if key not in self._cache:
            func = self._block_map[key]
            mode = self._behavior_overrides.get(key, func.mode)
            i = 0
            while True:
                i, _ = utils._seek_tag(self._lines, func.tags, start=i)
                if i is None:
                    break
                result = func(self._lines, i=i)
                if result is not None:
                    self._store_result(func, result)
                if not func.multiple or mode == 'first':
                    break
                i += 1  # continue scanning
        return getattr(self, key, None)

    def set_behavior(self, key, mode):
        assert key in self._block_map, f'Unknown key: {key}'
        assert mode in {'first', 'last', 'all'}, f'Invalid mode: {mode}'
        self._behavior_overrides[key] = mode
    
    def _build_maps(self):
        funcs = [f for _, f in inspect.getmembers(blocks, inspect.isfunction)
                 if hasattr(f, 'key')]
        self._block_map = {f.key: f for f in funcs}
        self._tag_map = {}
        for f in funcs:
            for tag in f.tags:
                self._tag_map.setdefault(tag, []).append(f)

    def _store_result(self, func, result):
        key = func.key
        mode = self._behavior_overrides.get(key, func.mode)
        if mode == 'all':
            if not func.multiple:
                raise ValueError(f"Block '{key}' does not support mode='all'")
            if not hasattr(self, key):
                setattr(self, key, [])
            elif not isinstance(getattr(self, key), list):
                setattr(self, key, [getattr(self, key)])
            getattr(self, key).append(result)
        elif mode == 'first':
            if not hasattr(self, key):
                setattr(self, key, result)
        elif mode == 'last':
            setattr(self, key, result)
        else:
            raise ValueError(f"Unknown mode '{mode}'")
        self._cache.add(key)

