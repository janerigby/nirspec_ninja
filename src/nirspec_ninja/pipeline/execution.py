"""Shared pipeline mechanics; caller-supplied parameters remain explicit."""

def apply_pipeline_parameters(pipeline, parameters):
    """Recursively apply a nested parameter dictionary to a pipeline."""
    for name, value in parameters.items():
        if isinstance(value, dict):
            apply_pipeline_parameters(getattr(pipeline, name), value)
        else:
            setattr(pipeline, name, value)



def check_pipeline_step_statuses(result, input_file):
    """Raise if any returned data model records a pipeline step as FAILED."""
    if result is None:
        return
    models = [result] if hasattr(result, 'meta') else result
    failed_steps = set()
    for model in models:
        cal_step = getattr(getattr(model, 'meta', None), 'cal_step', None)
        statuses = getattr(cal_step, '_instance', {})
        failed_steps.update(
            step_name for step_name, status in statuses.items()
            if status == 'FAILED'
        )
    if failed_steps:
        failed_list = ', '.join(sorted(failed_steps))
        raise RuntimeError(
            f'Pipeline step(s) reported FAILED for {input_file}: {failed_list}'
        )



def close_pipeline_result(result):
    """Close a returned data model or model container when possible."""
    if isinstance(result, (list, tuple)):
        for model in result:
            close_pipeline_result(model)
        return
    close = getattr(result, 'close', None)
    if close is not None:
        close()

