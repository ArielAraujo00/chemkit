from dataclasses import dataclass, field, replace

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

    show_energy: bool | None = None
    energy_dx: float = 0.0
    energy_dy: float = 0.4
    energy_fontsize: int | None = None
    energy_color: str | None = None
    energy_format: str = '.1f'

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
        return replace(
            self.step, **{k: v for k, v in {'color': self.color,
                                            'linewidth': self.linewidth,
                                            'linestyle': self.linestyle,
                                            'alpha': self.alpha,
                                            'zorder': self.zorder}.items() 
                          if v is not None})
    
    def connection_style(self):
        return replace(
            self.connection, **{k: v for k, v in {'color': self.color,
                                                  'linewidth': self.linewidth,
                                                  'linestyle': self.linestyle,
                                                  'alpha': self.alpha,
                                                  'zorder': self.zorder}.items()
                                if v is not None})

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
    show_step_labels: bool = False
    step_label_prefix: str | None = None
    step_label_suffix: str | None = None
    show_step_energies: bool = False


    def step_style(self):
        return replace(
            self.step, **{k: v for k, v in {'color': self.color,
                                            'linewidth': self.linewidth,
                                            'linestyle': self.linestyle,
                                            'alpha': self.alpha,
                                            'zorder': self.zorder}.items() 
                          if v is not None})
    
    def connection_style(self):
        return replace(
            self.connection, **{k: v for k, v in {'color': self.color,
                                                  'linewidth': self.linewidth,
                                                  'linestyle': self.linestyle,
                                                  'alpha': self.alpha,
                                                  'zorder': self.zorder}.items()
                                if v is not None})
    

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