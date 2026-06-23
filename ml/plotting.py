from dataclasses import dataclass, field
import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text

_metric_map = {'corr': 'r',
                'r2': 'R²',
               'mae': 'MAE',
               'rmse': 'RMSE',
               'rss': 'RSS',
               'tss': 'TSS',
               'roc_auc': 'ROC-AUC',
               'accuracy': 'Accuracy',
               'precision': 'Precision',
               'recall': 'Recall',
               'f1': 'F1',
               'n': 'n',
               'dof': 'DOF'}

@dataclass(slots=True)
class ScatterStyle:
    marker: str = 'o'
    color: str = 'black'
    size: float = 30
    alpha: float = 1.0
    edgecolor: str = 'face'
    linewidth: float = 0

@dataclass(slots=True)
class LineStyle:
    color: str = 'red'
    linewidth: float = 1.2
    linestyle: str = '-'
    alpha: float = 1.0

@dataclass(slots=True)
class LabelStyle:
    fontsize: int = 10
    adjust: bool = True

@dataclass
class PlotConfig:
    figsize: tuple = (10, 8)
    dpi: int = 300
    
    scatter: ScatterStyle = field(default_factory=ScatterStyle)
    line: LineStyle = field(default_factory=LineStyle)
    labels: LabelStyle = field(default_factory=LabelStyle)
    
    title_fontsize: int = 20
    label_fontsize: int = 18
    ticks_fontsize: int = 14
    legend_fontsize: int = 14
    notes_fontsize: int = 10

    confidence_color: str = '#d3d3d3'
    confidence_alpha: float = 0.8

    prediction_color: str = '#5a5a5a'
    prediction_style: str = '--'

    cluster_palette: str = 'tab10'

    show_grid: bool = False
    show_confidence: bool = False

    metrics_x: float = 0.98
    metrics_y: float = 0.98
    metrics_ha: str = 'right'
    metrics_va: str = 'top'
    
    xlabel: str | None = None
    ylabel: str | None = None
    title: str | None = None


class Plotter:
    def __init__(self, config=None):
        self.config = config or PlotConfig()
    def figure(self, ax=None):
        if ax is not None:
            return ax.figure, ax
        return plt.subplots(figsize=self.config.figsize, dpi=self.config.dpi)

    def style(self, ax, title=None, xlabel=None, ylabel=None):
        title = title or self.config.title
        if title:
            ax.set_title(title, fontsize=self.config.title_fontsize)
        
        xlabel = xlabel or self.config.xlabel
        if xlabel:
            ax.set_xlabel(xlabel, fontsize=self.config.label_fontsize)
        
        ylabel = ylabel or self.config.ylabel
        if ylabel:
            ax.set_ylabel(ylabel, fontsize=self.config.label_fontsize)
        
        ax.tick_params(labelsize=self.config.ticks_fontsize)
        if self.config.show_grid:
            ax.grid(True)
        return ax

    def annotate(self, ax, x, y, labels):
        texts = [ax.annotate(label, (xi, yi), fontsize=self.config.labels.fontsize)
                 for xi, yi, label in zip(x, y, labels)]
        if self.config.labels.adjust:
            adjust_text(texts, ax=ax)
        return ax

    def save(self, fig, path, **kwargs):
        fig.savefig(path, bbox_inches='tight',
                    dpi=self.config.dpi, **kwargs)

    def close(self, fig):
        plt.close(fig)
    
    def regression(self, X, y, model, labels=None, metrics=None,
                   title=None, xlabel=None, ylabel=None, ax=None):
        fig, ax = self.figure(ax)
        
        X = np.asarray(X).ravel()
        y = np.asarray(y)
        
        mask = ~(np.isnan(X) | np.isnan(y))
        X = X[mask]
        y = y[mask]
        
        if labels is not None:
            labels = np.asarray(labels)[mask]
        
        X_line = np.linspace(X.min(), X.max(), 100)
        y_line = model.predict(X_line.reshape(-1, 1))
        
        ax.scatter(X, y, marker=self.config.scatter.marker,
                   color=self.config.scatter.color,
                   s=self.config.scatter.size,
                   alpha=self.config.scatter.alpha,
                   edgecolors=self.config.scatter.edgecolor,
                   linewidths=self.config.scatter.linewidth)

        ax.plot(X_line, y_line,
                linestyle=self.config.line.linestyle,
                color=self.config.line.color,
                linewidth=self.config.line.linewidth,
                alpha=self.config.line.alpha)
        if metrics is not None:
            text = '\n'.join(f'{_metric_map.get(k, k)} = {v:.3f}'
                             for k, v in metrics.items()
                             if np.isscalar(v))
        
            ax.text(self.config.metrics_x, self.config.metrics_y, text, transform=ax.transAxes,
                    ha=self.config.metrics_ha, va=self.config.metrics_va,
                    
                    fontsize=self.config.notes_fontsize,
                    bbox={'facecolor': 'white', 'alpha': 0.8,
                          'edgecolor': 'none'})
        
        if labels is not None:
            self.annotate(ax, X, y, labels)
        self.style(ax, title, xlabel, ylabel)
        return fig, ax

    def parity(self, y_true, y_pred, labels=None, title=None,
               xlabel='Experimental', ylabel='Predicted', ax=None):
        fig, ax = self.figure(ax)
        y_true = np.asarray(y_true)
        y_pred = np.asarray(y_pred)
        
        mask = ~(np.isnan(y_true) | np.isnan(y_pred))
        
        y_true = y_true[mask]
        y_pred = y_pred[mask]
        
        if labels is not None:
            labels = np.asarray(labels)[mask]
        
        lim = max(np.max(y_true), np.max(y_pred))
        ax.scatter(y_true, y_pred,
                   marker=self.config.scatter.marker,
                   color=self.config.scatter.color,
                   s=self.config.scatter.size,
                   alpha=self.config.scatter.alpha,
                   edgecolors=self.config.scatter.edgecolor,
                   linewidths=self.config.scatter.linewidth)
        ax.plot([0, lim], [0, lim],
                linestyle=self.config.line.linestyle,
                color=self.config.line.color,
                linewidth=self.config.line.linewidth,
                alpha=self.config.line.alpha)
        
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        if labels is not None:
            self.annotate(ax, y_true, y_pred, labels)
        self.style(ax, title, xlabel, ylabel)
        return fig, ax
    
    def embedding(self, embedding, mask=None, x=0, y=1, title=None, ax=None):
        fig, ax = self.figure(ax)
        embedding = np.asarray(embedding)
        mask_valid = ~np.isnan(embedding[:, [x, y]]).any(axis=1)
        embedding = embedding[mask_valid]
        if mask is not None:
            mask = np.asarray(mask)[mask_valid]
        if mask is None:
            ax.scatter(embedding[:, x], embedding[:, y],
                       marker=self.config.scatter.marker,
                       s=self.config.scatter.size,
                       alpha=self.config.scatter.alpha,
                       edgecolors=self.config.scatter.edgecolor,
                       linewidths=self.config.scatter.linewidth,
                       color=self.config.scatter.color)
        else:
            ax.scatter(embedding[:, x], embedding[:, y], c=mask,
                       cmap=self.config.cluster_palette,
                       s=self.config.scatter.size,
                       alpha=self.config.scatter.alpha,
                       edgecolors=self.config.scatter.edgecolor,
                       linewidths=self.config.scatter.linewidth)
        
        self.style(ax, title, f'Dim{x + 1}', f'Dim{y + 1}')
        return fig, ax

    def overlay(self, ax, embedding, mask, x=0, y=1,
                marker='*', color='black', size=100, label=None):
        embedding = np.asarray(embedding)

        mask_valid = ~np.isnan(embedding[:, [x, y]]).any(axis=1)
        embedding = embedding[mask_valid]
        ax.scatter(embedding[:, x], embedding[:, y],
                   marker=marker, color=color, s=size, label=label)
        if label:
            ax.legend()
        return ax