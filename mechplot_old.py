from dataclasses import dataclass, field, fields, replace
from collections.abc import Iterable
from typing import Literal

import matplotlib.pyplot as plt
import pandas as pd

# ----- Styles -----
@dataclass(slots=True)
class StepStyle:
    color: str | None = None
    linewidth: float | None = None
    linestyle: str | None = None
    alpha: float | None = None
    zorder: int | None = None

    label_dx: float = 0.0
    label_dy: float = -0.4
    label_fontsize: int | None = None
    label_color: str | None = None


@dataclass(slots=True)
class ConnectionStyle:
    color: str | None = None
    linewidth: float | None = None
    linestyle: str | None = None
    alpha: float | None = None
    zorder: int | None = None

@dataclass(slots=True)
class BlockStyle:
    color: str | None = 'Black'
    linewidth: float | None = None
    linestyle: str | None = None
    alpha: float | None = None
    zorder: int | None = None
    
    step: StepStyle = field(
        default_factory=lambda: StepStyle(
            color='black', linewidth=1.0, linestyle='-', alpha=1.0, zorder=5))

    connection: ConnectionStyle = field(
        default_factory=lambda: ConnectionStyle(
            color='black', linewidth=1.0, linestyle='--', alpha=1.0, zorder=5))
    
    def step_style(self):
        return StepStyle(color=self.color, linewidth=self.linewidth,
                         linestyle=self.linestyle, alpha=self.alpha,
                         zorder=self.zorder)

    def connection_style(self):
        return ConnectionStyle(color=self.color, linewidth=self.linewidth,
                               linestyle=self.linestyle, alpha=self.alpha,
                               zorder=self.zorder)

@dataclass(slots=True)
class ProfileStyle:
    color: str | None = 'Black'
    linewidth: float | None = None
    linestyle: str | None = None
    alpha: float | None = None
    zorder: int | None = None
    
    step: StepStyle = field(
        default_factory=lambda: StepStyle(
            color='black', linewidth=1.0, linestyle='-', alpha=1.0, zorder=5))

    connection: ConnectionStyle = field(
        default_factory=lambda: ConnectionStyle(
            color='black', linewidth=1.0, linestyle='--', alpha=1.0, zorder=5))

    label: str | None = None
    show_step_labels: bool = True
    step_label_prefix: str | None = None
    step_label_suffix: str | None = None

    def step_style(self):
        return StepStyle(color=self.color, linewidth=self.linewidth,
                         linestyle=self.linestyle, alpha=self.alpha,
                         zorder=self.zorder)

    def connection_style(self):
        return ConnectionStyle(color=self.color, linewidth=self.linewidth,
                               linestyle=self.linestyle, alpha=self.alpha,
                               zorder=self.zorder)
    

@dataclass(slots=True)
class PlotConfig:
    figsize: tuple[float, float] | None = (8, 6)
    dpi: int | None = 300
    
    plateau_width: float = 0.32
    energy_padding: float = 2.0
    show_grid: bool = False
    show_zero_line: bool = True
    skip_missing_steps: bool = True
    
    auto_ylim: bool = True
    ymin: float | None = None
    ymax: float | None = None

    fontname: str | None = None
    fontstyle: str | None = None
    fontweight: str | None = None
    fontcolor: str | None = None
    label_fontsize: int = 12
    legend_fontsize: int = 14
    axis_fontsize: int = 14
    title_fontsize: int = 16

    rc_ticks: list[float] | None = None
    xlabel: str | None = 'Reaction Coordinate'
    ylabel: str | None = r'$\Delta\Delta G$ (kcal/mol)'
    title: str | None = 'Mechanism Profile'


# ----- Main Structure -----
@dataclass(slots=True)
class Step:
    id: str
    rc: float
    label: str
    codes: list[str]
    style: StepStyle | None = None
    block: str | None = None
    x: float | None = None
    visible: bool = True
    def get_x(self) -> float:
        return self.rc if self.x is None else self.x

@dataclass(slots=True)
class Block:
    name: str
    step_ids: list[str]
    style: BlockStyle = field(default_factory=BlockStyle)
    @classmethod
    def from_steps(cls, name: str, steps: list[Step], **kwargs):
        return cls(name=name, step_ids=[s.id for s in steps], **kwargs)

    def linear_connections(self, style: ConnectionStyle | None = None) -> list['Connection']:
        return [Connection(source=a, target=b, style=style)
                for a, b in zip(self.step_ids[:-1], self.step_ids[1:])]

@dataclass(slots=True)
class Connection:
    source: str
    target: str
    style: ConnectionStyle | None = None
    curved: bool = False
    
    @classmethod
    def from_blocks(cls, block1: Block, block2: Block, steps: list[Step],
                    style: ConnectionStyle | None = None) -> list['Connection']:
        step_map = {s.id: s for s in steps}
        b1 = [step_map[sid] for sid in block1.step_ids]
        b2 = [step_map[sid] for sid in block2.step_ids]
        b2_rc_map = {s.rc: s for s in b2}
        b1_rcs = {s.rc for s in b1}
        connections = []
        for s1 in b1:
            prev_rc = s1.rc - 1
            if prev_rc in b2_rc_map and prev_rc not in b1_rcs:
                connections.append(cls(source=s1.id,
                                       target=b2_rc_map[prev_rc].id,
                                       style=style))
            next_rc = s1.rc + 1
            if next_rc in b2_rc_map and next_rc not in b1_rcs:
                connections.append(cls(source=s1.id,
                                       target=b2_rc_map[next_rc].id,
                                       style=style))
        return connections

@dataclass
class Mechanism:
    steps: list[Step]
    connections: list[Connection] = field(default_factory=list)
    blocks: list[Block] = field(default_factory=list)
    auto_connect: bool = False
    def __post_init__(self):
        self.step_map = {step.id: step for step in self.steps}
        self.block_map = {block.name: block for block in self.blocks}

    def get_step(self, step_id: str) -> Step:
        return self.step_map[step_id]

    def get_block(self, block_name: str) -> Block:
        return self.block_map[block_name]

    def add_step(self, step: Step | list[Step]):
        if isinstance(step, list):
            self.step_map.update({s.id: s for s in step})
        else:
            self.step_map.update({step.id: step})
    
    def add_block(self, block: Block | list[Block]):
        if isinstance(block, list):
            self.block_map.update({b.name: b for b in block})
        else:
            self.block_map.update({block.name: block})

    def add_connection(self, connection: Connection | list[Connection]):
        if isinstance(connection, list):
            self.connections.extend(connection)
        else:
            self.connections.append(connection)
    
    def compute_energies(self, df: pd.DataFrame,
                         energy_col: str,
                         reference: str | None = None) -> pd.DataFrame:

        # Calculate Absolute
        records = []
        for step in self.steps:
            energy = sum(df.loc[df['code'] == code, energy_col].iloc[0]
                         for code in step.codes)
            records.append({'step_id': step.id,
                            'label': step.label,
                            'rc': step.rc,
                            'x': step.get_x(),
                            'energy': energy,
                            'block': step.block})
        df_plot = pd.DataFrame(records)

        # Calculate Relative
        if reference is None:
            reference_energy = df_plot['energy'].min()
        else:
            reference_energy = df_plot.loc[df_plot['step_id'] == reference, 'energy'].iloc[0]
        df_plot['relative_energy'] = (df_plot['energy'] - reference_energy)

        return df_plot
    
    @classmethod
    def from_blocks(cls, steps: list[Step], blocks: list[Block],
                    connections: list[Connection] | None = None):
        internal_connections = []
        for block in blocks:
            internal_connections.extend(block.linear_connections())
        if connections is not None:
            internal_connections.extend(connections)
        return cls(steps=steps, blocks=blocks,
                   connections=internal_connections)

# ----- PLOTTER -----
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

        ax.plot([], [], color=profile_style.step.color,
                linewidth=profile_style.step.linewidth,
                linestyle=profile_style.step.linestyle,
                label=profile_style.label)
        
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
            style = self._resolve_style(step.style, block.step if block else None,
                                        block.step_style() if block else None,
                                        profile_style.step, profile_style.step_style(),
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
            style = self._resolve_style(connection.style, block.connection if block else None,
                                        block.connection_style() if block else None,
                                        profile_style.connection, profile_style.connection_style(),
                                        ConnectionStyle())

            x1 = source_data['x'] + config.plateau_width
            x2 = target_data['x'] - config.plateau_width
            y1 = source_data['relative_energy']
            y2 = target_data['relative_energy']
            ax.plot([x1, x2], [y1, y2], color=style.color,
                    linewidth=style.linewidth, linestyle=style.linestyle,
                    alpha=style.alpha, zorder=style.zorder)

    def _resolve_style(self, *styles):
        styles = [s for s in styles if s is not None]
        style_type = type(styles[-1])
        resolved = {}
        for field in fields(style_type):
            for style in styles:
                value = getattr(style, field.name)
                if value is not None:
                    resolved[field.name] = value
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
