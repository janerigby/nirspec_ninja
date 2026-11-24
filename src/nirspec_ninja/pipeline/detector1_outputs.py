"""Control optional Detector1 products without changing installed JWST code."""


def configure_rateints_saving(pipeline, save_rateints=False):
    """Skip only the explicit RATEINTS save; retain all fitting and RATE output."""
    if save_rateints:
        return
    save_model = pipeline.save_model

    def save_without_rateints(model, suffix=None, *args, **kwargs):
        if suffix == 'rateints':
            return None
        return save_model(model, suffix, *args, **kwargs)

    pipeline.save_model = save_without_rateints
