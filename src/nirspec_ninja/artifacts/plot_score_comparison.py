"""Standalone Matplotlib plot from scores.csv; requires matplotlib and numpy.
Edit the settings below or plot_scores() to customize the figure.
Usage: python plot_score_comparison.py /path/to/scores.csv
Optional: --new-csv another.csv --output comparison.png
Combine IFU and MSA: scores.csv --combine-modes
Compare a run to training: scores.csv --training-csv training/scores.csv
"""
import os,csv,argparse
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).resolve().parent/'mpl_config'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, LogFormatterMathtext, NullFormatter
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Editable appearance settings.
FIGURE_SIZE = (14, 10)
DPI = 160
FONT_SIZE = 14
LEGEND_FONT_SIZE = 12
HISTOGRAM_BINS = 64
HISTOGRAM_PANEL_SIZE = '32%'
# Linear counts below 5, logarithmic above: small populations stay visible.
HISTOGRAM_COUNT_LINTHRESH = 5
HISTOGRAM_HATCHES = {'Dev Bad': '///', 'Dev OK': '\\\\'}
# Scatter areas in points squared. Distinct outline sizes expose coincident groups.
MARKER_SIZES = {'Dev Bad': 100, 'Mystery': 65, 'New': 40, 'Dev OK': 100}
MARKER_LINEWIDTH = 1.4
MARKER_ALPHA = 0.85
ZOOM_X_LIMITS = (0.8, 10)
ZOOM_Y_LIMITS = (-0.12, 1)
REFERENCE_COLUMN_SCORE = 10
REFERENCE_GLOW_SCORE = 10

def plot_scores(baseline, new_csvs=(), output=None, training_csvs=(), show_mystery=True, input_label='Input', show_histograms=False, combine_modes=False):
    """Read CSVs, save the comparison PNG, and return its path."""
    plt.rcParams.update({'font.size': FONT_SIZE,
                         'axes.labelsize': FONT_SIZE + 2,
                         'axes.titlesize': FONT_SIZE + 2,
                         'xtick.labelsize': FONT_SIZE,
                         'ytick.labelsize': FONT_SIZE})
    base=Path(baseline);base=base/'scores.csv' if base.is_dir() else base
    rows=list(csv.DictReader(base.open()))
    # The positional CSV is always Input; labeled comparison samples are opt-in.
    for r in rows:
        r['plot_group'] = 'New'
    if training_csvs:
        training_rows = []
        for path in training_csvs:
            with Path(path).open(newline='') as handle:
                training_rows.extend(csv.DictReader(handle))
    else:
        training_rows = []
    for r in training_rows:
     folder=Path(r['file']).parent.name
     r['plot_group']='Dev OK' if folder in {'IFU_okay','MSA_okay'} else 'Dev Bad' if folder in {'IFU_bad','MSA_bad'} else 'Mystery'
    if training_csvs:
        rows.extend(training_rows)
    for path in new_csvs:
     for r in csv.DictReader(open(path)):
      r['plot_group']='New';rows.append(r)
    hidden_mystery = sum(r['plot_group'] == 'Mystery' for r in rows) if not show_mystery else 0
    included_rows = [r for r in rows if show_mystery or r['plot_group'] != 'Mystery']
    valid=[r for r in included_rows if r.get('status')=='ok' and r.get('glow_status')=='ok' and r.get('glow_score')!='']
    panels = [('IFU + MSA', {'IFU', 'MSA'})] if combine_modes else [('IFU', {'IFU'}), ('MSA', {'MSA'})]
    figure_size = (FIGURE_SIZE[0] / 2, FIGURE_SIZE[1]) if combine_modes else FIGURE_SIZE
    fig,axs=plt.subplots(2,len(panels),figsize=figure_size,squeeze=False,
                         gridspec_kw={'height_ratios':[2,1]})
    styles=[('Dev Bad','tomato','^'),('Mystery','purple','x'),('New','royalblue','s'),('Dev OK','teal','o')]
    sample_labels = {'Dev Bad': 'Visual Bad', 'Mystery': 'Mystery',
                     'New': input_label, 'Dev OK': 'Visual OK'}
    for col,(mode, panel_modes) in enumerate(panels):
     for ax in axs[:,col]:
      for group,color,marker in styles:
       rr=[r for r in valid if r.get('mode') in panel_modes and r['plot_group']==group]
       if not rr:continue
       style=dict(s=MARKER_SIZES[group],alpha=MARKER_ALPHA,marker=marker,
                  linewidths=MARKER_LINEWIDTH,label=f'{sample_labels[group]} ({len(rr)})',
                  zorder=4 if group in {'Dev Bad', 'Dev OK'} else 3)
       # Hollow triangles/squares/circles keep underlying groups visible.
       # An x has no face, so give it a normal color rather than edgecolors.
       if marker == 'x':style.update(color=color)
       else:style.update(facecolors='none',edgecolors=color)
       ax.scatter([float(r['column_score']) for r in rr],[float(r['glow_score']) for r in rr],**style)
      ax.axvline(REFERENCE_COLUMN_SCORE, color='gray', linestyle='--', linewidth=1.2, zorder=1)
      ax.axhline(REFERENCE_GLOW_SCORE, color='gray', linestyle='--', linewidth=1.2, zorder=1)
      ax.set_xscale('log')
      ax.xaxis.set_major_locator(LogLocator(base=100, subs=(1.0,), numticks=100))
      ax.xaxis.set_major_formatter(LogFormatterMathtext(base=10))
      ax.xaxis.set_minor_locator(LogLocator(base=10, subs=(1.0,), numticks=100))
      ax.xaxis.set_minor_formatter(NullFormatter())
      ax.tick_params(axis='x', which='major', length=8, width=1.2)
      ax.tick_params(axis='x', which='minor', length=4, width=1)
      ax.grid(alpha=.2);ax.set_xlabel('Column score')
      ax.set_ylabel('glow score')
     axs[0,col].set_yscale('symlog',linthresh=1,linscale=1)
     axs[0,col].set_title(mode);axs[0,col].legend(fontsize=LEGEND_FONT_SIZE,loc='upper left')
     axs[1,col].set_xlim(*ZOOM_X_LIMITS);axs[1,col].set_ylim(*ZOOM_Y_LIMITS)
     # Label every decade in the zoom panels so both 1 and 10 are marked.
     axs[1,col].xaxis.set_major_locator(LogLocator(base=10, subs=(1.0,), numticks=100))
     axs[1,col].xaxis.set_major_formatter(LogFormatterMathtext(base=10))
     # Also retain two labeled ticks when the editable zoom spans less than a decade.
     zoom_low, zoom_high = ZOOM_X_LIMITS
     decade_ticks = [10.0 ** exponent for exponent in range(
         int(np.ceil(np.log10(zoom_low))), int(np.floor(np.log10(zoom_high))) + 1)]
     if len(decade_ticks) < 2:
      axs[1,col].set_xticks([zoom_low, zoom_high], labels=[f'{zoom_low:g}', f'{zoom_high:g}'])
     axs[1,col].set_title(f'{mode}: enlarged low-score region',fontsize=FONT_SIZE + 1)
    xlimits=(min(a.get_xlim()[0] for a in axs[0]),max(a.get_xlim()[1] for a in axs[0]))
    ymax=max(a.get_ylim()[1] for a in axs[0])
    for ax in axs[0]:ax.set_xlim(xlimits);ax.set_ylim(-.15,ymax)
    if show_histograms:
        # Marginal histograms share score axes and bin edges across modes.
        xbins=np.geomspace(xlimits[0],xlimits[1],HISTOGRAM_BINS + 1)
        transform=axs[0,0].yaxis.get_transform()
        ybins=transform.inverted().transform(np.linspace(transform.transform(np.array([-.15]))[0],transform.transform(np.array([ymax]))[0],HISTOGRAM_BINS + 1))
        marginals=[]
        for col,(mode, panel_modes) in enumerate(panels):
         ax=axs[0,col];divider=make_axes_locatable(ax)
         top=divider.append_axes('top',size=HISTOGRAM_PANEL_SIZE,pad=.12,sharex=ax)
         right=divider.append_axes('right',size=HISTOGRAM_PANEL_SIZE,pad=.12,sharey=ax)
         top.set_title(mode);ax.set_title('')
         max_column_count = max_glow_count = 0
         for group,color,marker in styles:
          rr=[r for r in valid if r.get('mode') in panel_modes and r['plot_group']==group]
          if not rr:continue
          # Training populations use transparent hatched interiors; comparison
          # populations stay above them as plain outlines.
          histogram_style = dict(histtype='step', color=color, lw=1.5,
                                 label=sample_labels[group], zorder=3)
          if group in HISTOGRAM_HATCHES:
           histogram_style.update(histtype='stepfilled', facecolor='none',
                                  edgecolor=color, hatch=HISTOGRAM_HATCHES[group],
                                  zorder=2)
          column_counts, _, _ = top.hist([float(r['column_score']) for r in rr],bins=xbins,**histogram_style)
          glow_counts, _, _ = right.hist([float(r['glow_score']) for r in rr],bins=ybins,orientation='horizontal',**histogram_style)
          max_column_count = max(max_column_count, column_counts.max())
          max_glow_count = max(max_glow_count, glow_counts.max())
         top.set_yscale('symlog', linthresh=HISTOGRAM_COUNT_LINTHRESH)
         right.set_xscale('symlog', linthresh=HISTOGRAM_COUNT_LINTHRESH)
         # Keep at most three count labels readable on the small symlog panels.
         for histogram, count_axis, peak in [(top, 'y', max_column_count), (right, 'x', max_glow_count)]:
          upper = max(2, peak * 1.15)
          positive_ticks = []
          power = 10
          while power <= upper:
           positive_ticks.append(power)
           power *= 10
          if positive_ticks:
           ticks = sorted({0, positive_ticks[0], positive_ticks[-1]})
          else:
           ticks = [0, int(upper)]
          if count_axis == 'y':
           histogram.set_ylim(0, upper)
           histogram.set_yticks(ticks, labels=[str(tick) for tick in ticks])
          else:
           histogram.set_xlim(0, upper)
           histogram.set_xticks(ticks, labels=[str(tick) for tick in ticks])
         top.tick_params(axis='x',labelbottom=False);right.tick_params(axis='y',labelleft=False)
         top.set_ylabel('Count',fontsize=FONT_SIZE);right.set_xlabel('Count',fontsize=FONT_SIZE)
         top.legend(fontsize=LEGEND_FONT_SIZE - 2, loc='upper right', framealpha=.8)
         top.grid(alpha=.15);right.grid(alpha=.15)
         marginals.append((top,right))
    #fig.text(.5,.025,f'Glow scale is linear below 1 and logarithmic above.\n{len(valid)} points shown; {len(included_rows)-len(valid)} rows without valid masked scores omitted; {hidden_mystery} Mystery rows hidden ,ha='center',fontsize=FONT_SIZE - 2)
    fig.tight_layout(rect=(0,.08,1,1))
    out=Path(output) if output else base.parent/'score_comparison.png'
    out.parent.mkdir(parents=True,exist_ok=True);fig.savefig(out,dpi=DPI);plt.close(fig)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline', help='Input results directory or scores.csv')
    parser.add_argument('--new-csv', nargs='+', default=[], help='Additional CSVs shown as Input')
    parser.add_argument('--training-csv', nargs='+', default=[],
                        help='Optional labeled comparison CSVs. None are loaded by default; use scores calculated with the same method as Input')
    parser.add_argument('--input-label', default='Input',
                        help='Name of the input sample in scatter and histogram legends (default: Input)')
    parser.add_argument('--output', help='Output PNG path (default: next to baseline CSV)')
    parser.add_argument('--combine-modes', action='store_true',
                        help='Combine IFU and MSA in one score panel and one low-score zoom (default: separate panels)')
    histograms = parser.add_mutually_exclusive_group()
    histograms.add_argument('--show-histograms', dest='show_histograms', action='store_true',
                            help='Show marginal score histograms')
    histograms.add_argument('--hide-histograms', dest='show_histograms', action='store_false',
                            help='Hide marginal score histograms (default)')
    parser.set_defaults(show_histograms=False)
    mystery = parser.add_mutually_exclusive_group()
    mystery.add_argument('--show-mystery', dest='show_mystery', action='store_true',
                         help='Show Mystery samples (default)')
    mystery.add_argument('--hide-mystery', dest='show_mystery', action='store_false',
                         help='Hide Mystery samples in scatter plots and histograms')
    parser.set_defaults(show_mystery=True)
    args = parser.parse_args()
    print(plot_scores(args.baseline, args.new_csv, args.output, args.training_csv, args.show_mystery, args.input_label, args.show_histograms, args.combine_modes))


if __name__ == '__main__':
    main()
