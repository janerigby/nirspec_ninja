"""Run the reusable gradient diagnostic on the supplied COSMIC.EYE example."""
from pathlib import Path
from nirspec_ninja.diagnostics.fs_gradients import make_fs_gradient_report

if __name__ == '__main__':
    make_fs_gradient_report(Path.home() / 'AI_experiments/Temp/COSMIC.EYE/FS',
                            Path.home() / 'AI_experiments/Temp/FS_gradient_review', 'COSMIC.EYE')
