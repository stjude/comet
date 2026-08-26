### Shared helpers for the run_*.py wrapper scripts: path resolution for
### CLI-supplied file paths, and a restricted unpickler for loading the
### splits/metadata/representations pickles these pipelines pass around.
###
### Rationale: every wrapper takes file paths straight from argparse and
### hands them to open()/pickle.load()/pickle.dump() with no validation
### (flagged by static analysis as path traversal + unsafe deserialization).
### These pipelines are internal research tools where the caller supplies
### their own trusted paths, so the goal here isn't to sandbox a hostile
### user — it's to (a) canonicalize paths before they're used for I/O and
### (b) stop pickle.load from being able to instantiate arbitrary classes
### (the actual code-execution vector in "unsafe deserialization"), so a
### corrupted or mishandled pickle can't do more than fail to load.

import builtins
import io
import os
import pickle

# Modules that legitimate splits/metadata/representations pickles in this
# codebase construct instances from (numpy arrays, pandas objects, plain
# containers). Anything outside this allowlist is rejected.
_ALLOWED_MODULES = {
    'builtins',
    'collections',
    'numpy',
    'numpy.core.multiarray',
    'numpy._core.multiarray',
    'numpy.core.numeric',
    'pandas',
    'pandas.core.series',
    'pandas.core.frame',
    'pandas.core.indexes.base',
    'pandas.core.internals.managers',
    'datetime',
}

_ALLOWED_BUILTINS = {'set', 'frozenset', 'complex', 'bytearray'}


class RestrictedUnpickler(pickle.Unpickler):
    """Unpickler that only allows classes/functions from _ALLOWED_MODULES.

    Blocks the classic pickle RCE gadget (arbitrary module.callable via
    GLOBAL opcode, e.g. os.system / subprocess.Popen / builtins.eval).
    """

    def find_class(self, module, name):
        if module == 'builtins' and name not in _ALLOWED_BUILTINS:
            raise pickle.UnpicklingError(
                'blocked unpickling of builtins.%s (not in allowlist)' % name)
        if module not in _ALLOWED_MODULES and not module.startswith('numpy.') \
                and not module.startswith('pandas.'):
            raise pickle.UnpicklingError(
                'blocked unpickling of %s.%s (module not in allowlist)' % (module, name))
        return super().find_class(module, name)


def resolve_path(path):
    """Canonicalize a user/CLI-supplied path before it's used for file I/O."""
    return os.path.realpath(os.path.expanduser(path))


def safe_pickle_load(path):
    """Resolve `path` and load it with RestrictedUnpickler."""
    with open(resolve_path(path), 'rb') as f:
        return RestrictedUnpickler(f).load()


def safe_pickle_dump(obj, path):
    """Resolve `path` and pickle `obj` to it."""
    with open(resolve_path(path), 'wb') as f:
        pickle.dump(obj, f)
