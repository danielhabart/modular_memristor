from __future__ import annotations
import numpy as np
from .model import MemristorModel, ModelState


def simulate(
    model:MemristorModel,
    t:np.ndarray,
    v:np.ndarray,
    initState: ModelState
    ):
    return model.evaluate(t,v)