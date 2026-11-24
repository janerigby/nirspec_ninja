"""Shared transfer attempts, retry timing, and archive download calls."""

def retry_transfer(operation, *, attempts, exceptions, clock, sleep, delay,
                   success, failure, retry_message, exhausted):
    """Run a transfer policy while callers retain logging and restart decisions."""
    last_error = None
    for attempt in range(1, attempts + 1):
        started = clock()
        try:
            result = operation()
            success(result, started, attempt)
            return result
        except exceptions as error:
            failure(error, started, attempt)
            last_error = error
            if attempt == attempts:
                break
            seconds = delay(attempt)
            retry_message(error, attempt, seconds)
            sleep(seconds)
    assert last_error is not None
    raise exhausted(last_error) from last_error


def download_products(client, product, staging_dir, source):
    """Use the same explicit Astroquery transfer settings in both workflows."""
    return client.download_products(product, curl_flag=False,
                                    cloud_only=(source == 'AWS'),
                                    download_dir=str(staging_dir))
