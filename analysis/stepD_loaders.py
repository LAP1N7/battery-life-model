# -*- coding: utf-8 -*-
"""
Loader registry + validation of the standard long table.

Standard dataset object returned by every loader:
    dict(curves=DataFrame, rest=DataFrame (optional, may be empty), meta=DataFrame, info=dict)

curves (long): cell, checkup_idx, ah_throughput, direction ('charge'|'discharge'), Q_Ah, V, T   [+ thickness_um optional]
    Q_Ah  = Ah since the start of that branch, increasing (discharge counted positive). Include the CV tail in the charge branch
            if the protocol has one (the discharge starts from the post-CV state).
    one charge and one discharge branch per check-up, quasi-OCV rate, ordered by checkup_idx.
rest (long, optional, for D2): cell, checkup_idx, ah_throughput, phase ('top_rest'), t_rest_s (0 = start of rest), V, thickness_um, T
meta: cell + free columns; recommended: T_group (R/C/H or any label), room_mate (cell id of a same-batch room-temperature reference)
info: dict(dataset, vmin, vmax, c_nom_ah, ocp (registry name), ...)

External loaders: --loader path/to/my_loader.py:load_function. The file is executed (it is YOUR code) -> never point it at a
downloaded data folder; only at a file you wrote or reviewed.
"""
import importlib.util
import os

REQUIRED = ["cell", "checkup_idx", "ah_throughput", "direction", "Q_Ah", "V", "T"]


def validate(ds):
    cu = ds["curves"]
    miss = [c for c in REQUIRED if c not in cu.columns]
    if miss:
        raise ValueError(f"curves table lacks columns {miss}")
    bad = set(cu["direction"].unique()) - {"charge", "discharge"}
    if bad:
        raise ValueError(f"direction values must be charge|discharge, got {bad}")
    for need in ("vmin", "vmax", "ocp"):
        if need not in ds["info"]:
            raise ValueError(f"info lacks '{need}'")
    return True


def _mohtat(data_dir, cells=None):
    import stepD_loader_mohtat
    import stepD_ocp_mohtat   # noqa: F401  (registers the OCP pair)
    return stepD_loader_mohtat.load(data_dir, cells)


HERE = os.path.dirname(os.path.abspath(__file__))
BUILTIN = {"mohtat": (_mohtat, os.path.join(os.path.dirname(HERE), "data"))}   # same convention as common.py: analysis/../data


def get_loader(name):
    if ":" in name and name.lower().rsplit(":", 1)[0].endswith(".py"):
        path, func = name.rsplit(":", 1)
        spec = importlib.util.spec_from_file_location("stepD_user_loader", os.path.abspath(path))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return getattr(mod, func), None
    if name not in BUILTIN:
        raise KeyError(f"unknown dataset '{name}'. Built-in: {sorted(BUILTIN)}; or use path.py:function")
    return BUILTIN[name]


def load_dataset(name, data_dir=None, cells=None):
    fn, default_dir = get_loader(name)
    ds = fn(data_dir or default_dir, cells)
    validate(ds)
    return ds
