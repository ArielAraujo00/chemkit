from dataclasses import replace, fields

from .mechanism import Step
from .styles import StepStyle


_TS_STYLE = StepStyle(label_dy=0.4, energy_dy=-0.8, linestyle='-', energy_fontsize=10)
_IM_STYLE = StepStyle(label_dy=-0.8, energy_dy=0.4, linestyle='-', energy_fontsize=10)


def _make_step(style, **kwargs):
    style_fields = {f.name for f in fields(StepStyle)}
    style_kwargs = {k: v for k, v in kwargs.items()
                    if k in style_fields}
    step_kwargs = {k: v for k, v in kwargs.items()
                   if k not in style_fields}
    style = replace(style, **style_kwargs)
    return style, step_kwargs

def TS(**kwargs):
    style, step_kwargs = _make_step(_TS_STYLE, **kwargs)
    return Step(style=style, label=step_kwargs.pop('label', step_kwargs['id']),
                codes=step_kwargs.pop('codes', [step_kwargs['id']]), **step_kwargs)

def IM(**kwargs):
    style, step_kwargs = _make_step(_IM_STYLE, **kwargs)
    return Step(style=style, label=step_kwargs.pop('label', step_kwargs['id']),
        codes=step_kwargs.pop('codes', [step_kwargs['id']]), **step_kwargs)