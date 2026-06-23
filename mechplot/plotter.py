from dataclasses import fields
from collections.abc import Iterable
#from typing import Literal

import matplotlib.pyplot as plt
import pandas as pd

from chemkit.mechplot.mechanism import *
from chemkit.mechplot.styles import *

class MechPlotter:
    def __init__(self, config: PlotConfig | None = None):
        self.config = config or PlotConfig()
        self._global_ylim = None
    
    def plot(self, mechanism: Mechanism, df: pd.DataFrame,
             energy: str | Iterable[str] = 'gibbs',
             profile_styles: dict[str, ProfileStyle] | None = None,
             reference: str | None = None):
        
        energy_cols = [energy] if isinstance(energy, str) else list(energy)
        fig, ax = plt.subplots(figsize=self.config.figsize, dpi=self.config.dpi)
        self._global_ylim = None
        for energy_col in energy_cols:
            profile_style = (profile_styles.get(energy_col, ProfileStyle())
                             if profile_styles is not None else ProfileStyle())
            
            self._plot_single(ax=ax, mechanism=mechanism, df=df, energy_col=energy_col,
                              reference=reference, profile_style=profile_style)

        self._apply_plot_config(ax)
        plt.tight_layout()
        return fig, ax

    def _plot_single(self, ax, mechanism: Mechanism, df: pd.DataFrame,
                     energy_col: str, reference: str | None,
                     profile_style: ProfileStyle):

        config = self.config
        df_plot = mechanism.compute_energies(df=df, energy_col=energy_col,
                                             reference=reference)

        self._draw_steps(ax=ax, mechanism=mechanism, df_plot=df_plot,
                         profile_style=profile_style)

        self._draw_connections(ax=ax, mechanism=mechanism, df_plot=df_plot,
                               profile_style=profile_style)

        style = profile_style.step_style()
        ax.plot([], [], label=profile_style.label,
                color=style.color,
                linewidth=style.linewidth,
                linestyle=style.linestyle)
        
        ymin = df_plot['relative_energy'].min() - config.energy_padding
        ymax = df_plot['relative_energy'].max() + config.energy_padding

        if self._global_ylim is None:
            self._global_ylim = [ymin, ymax]
        else:
            self._global_ylim[0] = min(self._global_ylim[0], ymin)
            self._global_ylim[1] = max(self._global_ylim[1], ymax)

        if config.auto_ylim:
            ax.set_ylim(*self._global_ylim)
        else:
            ax.set_ylim(config.ymin, config.ymax)

    def _draw_steps(self, ax, mechanism: Mechanism,
                    df_plot: pd.DataFrame, profile_style: ProfileStyle):
        config = self.config
        for row in df_plot.itertuples():
            step = mechanism.get_step(row.step_id)
            if not step.visible:
                continue
            
            block = (mechanism.get_block(step.block) 
                     if step.block is not None else None)
            style = self._resolve_style(StepStyle,
                                        step.style,
                                        block.style.step_style() if block else None,
                                        profile_style.step_style(),
                                        StepStyle())
            x = row.x
            y = row.relative_energy
            ax.hlines(y=y, color=style.color, alpha=style.alpha, zorder=style.zorder,
                      xmin=x-config.plateau_width, xmax=x+config.plateau_width,
                      linewidth=style.linewidth, linestyle=style.linestyle)
            
            if profile_style.show_step_labels and row.label:
                label = (f'{profile_style.step_label_prefix or ""}'
                         f'{row.label}'
                         f'{profile_style.step_label_suffix or ""}')

                ax.text(x + style.label_dx, y + style.label_dy,
                        label, ha='center',
                        fontsize=style.label_fontsize or config.label_fontsize,
                        color=style.label_color or config.fontcolor,
                        alpha=style.alpha, zorder=style.zorder)
                
            if profile_style.show_step_energies if style.show_energy is None else style.show_energy:
                energy = format(row.relative_energy, style.energy_format)
                ax.text(x + style.energy_dx, y + style.energy_dy,
                        energy, ha='center',
                        fontsize=style.energy_fontsize or config.label_fontsize,
                        color=style.energy_color or style.color,
                        alpha=style.alpha, zorder=style.zorder)

    def _draw_connections(self, ax, mechanism: Mechanism, df_plot: pd.DataFrame,
                          profile_style: ProfileStyle):
        config = self.config
        connections = (self._auto_connections(mechanism, df_plot)
                       if mechanism.auto_connect else mechanism.connections)
        energy_map = df_plot.set_index('step_id').to_dict(orient='index')
        for connection in connections:
            source = mechanism.get_step(connection.source)
            target = mechanism.get_step(connection.target)
            step_size = abs(target.rc-source.rc)
            if config.skip_missing_steps and step_size > 1:
                continue
            source_data = energy_map[source.id]
            target_data = energy_map[target.id]
            block = (mechanism.get_block(source.block)
                     if source.block is not None else None)
            style = self._resolve_style(ConnectionStyle,
                                        connection.style,
                                        block.style.connection_style() if block else None,
                                        profile_style.connection_style(),
                                        ConnectionStyle())
            x1 = source_data['x'] + config.plateau_width
            x2 = target_data['x'] - config.plateau_width
            y1 = source_data['relative_energy']
            y2 = target_data['relative_energy']
            ax.plot([x1, x2], [y1, y2], color=style.color,
                    linewidth=style.linewidth, linestyle=style.linestyle,
                    alpha=style.alpha, zorder=style.zorder)

    def _resolve_style(self, style_type, *styles):
        styles = [s for s in styles if s is not None]
        resolved = {}
        for field in fields(style_type):
            name = field.name
            for style in styles:
                # Direct style (StepStyle/ConnectionStyle)
                if isinstance(style, style_type):
                    value = getattr(style, name)
                # Container style (BlockStyle/ProfileStyle)
                elif hasattr(style, f'{style_type.__name__.replace("Style","").lower()}_style'):
                    nested = getattr(style, f'{style_type.__name__.replace("Style","").lower()}_style')()
                    value = getattr(nested, name)
                else:
                    value = None
                if value is not None:
                    resolved[name] = value
                    break
        return style_type(**resolved)

    def _auto_connections(self, mechanism, df_plot) -> list[Connection]:
        idx = df_plot.groupby('rc')['relative_energy'].idxmin()

        rows = df_plot.loc[idx].sort_values('rc').itertuples()
        rows = list(rows)
        return [Connection(source=a.step_id, target=b.step_id)
                for a, b in zip(rows[:-1], rows[1:])]

    def _apply_plot_config(self, ax):
        config = self.config
        if config.show_zero_line:
            ax.axhline(0, color='black', linestyle='--', linewidth=0.8)
        if config.show_grid:
            ax.grid(alpha=0.2)
        if config.rc_ticks is not None:
            ax.set_xticks(config.rc_ticks)
        if ax.get_legend_handles_labels()[1]:
            ax.legend(prop={'family': config.fontname, 'size': config.legend_fontsize,
                            'weight': config.fontweight, 'style': config.fontstyle})
        
        ax.set_xlabel(config.xlabel, fontsize=config.axis_fontsize,
                      fontname=config.fontname, fontstyle=config.fontstyle,
                      fontweight=config.fontweight)
        ax.set_ylabel(config.ylabel, fontsize=config.axis_fontsize,
                      fontname=config.fontname, fontstyle=config.fontstyle,
                      fontweight=config.fontweight)
        ax.set_title(config.title, fontsize=config.title_fontsize, fontname=config.fontname,
                     fontstyle=config.fontstyle, fontweight=config.fontweight)
