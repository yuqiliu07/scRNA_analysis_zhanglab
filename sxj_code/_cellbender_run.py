"""Run `cellbender` with a checkpoint fix; called by cellbender_gpu.py with the usual CLI arguments.

pyro attaches a weak reference (`param.unconstrained = weakref.ref(...)`) to every constrained parameter, and
torch 2.x pickles a Parameter's __dict__, so CellBender 0.3.2's torch.save(model) fails with "cannot pickle
'weakref.ReferenceType' object". Training continues, but the final posterior step reloads the checkpoint and
then aborts. Here those attributes are removed only while the checkpoint is written and restored right after;
nothing else in CellBender changes.

Also puts this interpreter's Scripts/bin directory on PATH: CellBender builds its HTML report by shelling out to
`jupyter nbconvert`, which is not found when the venv is not activated (the report then fails, counts unaffected).
"""
import os
import sys
import weakref

import torch
from cellbender.remove_background import checkpoint as _ck

_orig_save = torch.save


def _strip_weakrefs(obj):
    """Remove weakref attributes from every tensor reachable from obj (the CellBender model is a plain object
    holding nn.Modules as attributes, so walk the object graph rather than obj.parameters())."""
    removed, seen, stack = [], set(), [obj]
    while stack:
        o = stack.pop()
        if id(o) in seen:
            continue
        seen.add(id(o))
        if isinstance(o, torch.Tensor):
            for k, v in list(getattr(o, "__dict__", {}).items()):
                if isinstance(v, weakref.ReferenceType):
                    removed.append((o, k, v))
                    delattr(o, k)
            continue
        if isinstance(o, dict):
            stack.extend(o.values())
        elif isinstance(o, (list, tuple, set)):
            stack.extend(o)
        elif hasattr(o, "__dict__") and not isinstance(o, type):
            stack.extend(vars(o).values())
    return removed


def _save(obj, f, *args, **kwargs):
    removed = _strip_weakrefs(obj)
    try:
        return _orig_save(obj, f, *args, **kwargs)
    finally:
        for p, k, v in removed:
            setattr(p, k, v)


_ck.torch.save = _save
os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ.get("PATH", "")

if __name__ == "__main__":
    from cellbender.base_cli import main
    sys.argv = ["cellbender"] + sys.argv[1:]
    main()
