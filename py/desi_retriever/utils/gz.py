import os
import tempfile
import astropy.io.fits as pyfits
import httpio
import indexed_gzip
import pyarrow.parquet as pq
import fsspec
import aiohttp
from pylru import lrudecorator
from .fetcher import process_spectra_hdus, get_desi_login_password


@lrudecorator(20)
def get_index_for_file(
    filename,
    dataset,
    parquet_url_template="https://data.desi.lbl.gov/desi/users/koposov/spec_gz_db/spec_gz_index_{dataset}.parquet"
):
    """
    Fetch the index bytes for a specific filename from the remote parquet file.
    """
    parquet_url = parquet_url_template.format(dataset=dataset)
    print(f"Looking up index for {filename} in {parquet_url}...")
    user, pwd = get_desi_login_password()
    auth = aiohttp.BasicAuth(user, pwd)
    # fs = fsspec.filesystem("https", client_kwargs={'auth': auth})

    try:
        # Use fsspec open with auth directly
        with fsspec.open(parquet_url, "rb", auth=auth) as f:
            pf = pq.ParquetFile(f)

            # Optimization: Use row group statistics to find the file
            # This avoids reading the entire filename column

            filename_col_idx = -1
            for k in range(pf.metadata.num_columns):
                path = pf.metadata.schema.column(k).path
                if path == 'filename':
                    filename_col_idx = k
                    break

            if filename_col_idx == -1:
                print("Column 'filename' not found in parquet metadata")
                return None

            for i in range(pf.metadata.num_row_groups):
                rg = pf.metadata.row_group(i)
                col_meta = rg.column(filename_col_idx)

                if col_meta.is_stats_set:
                    min_val = col_meta.statistics.min
                    max_val = col_meta.statistics.max

                    # Handle bytes/str mismatch
                    if isinstance(min_val, bytes) and isinstance(
                            filename, str):
                        min_val = min_val.decode('utf-8')
                    if isinstance(max_val, bytes) and isinstance(
                            filename, str):
                        max_val = max_val.decode('utf-8')

                    if min_val <= filename <= max_val:
                        # Found candidate row group
                        t = pf.read_row_group(i, columns=['filename', 'index'])
                        filenames = t['filename'].to_pylist()
                        try:
                            idx = filenames.index(filename)
                            return t['index'][idx].as_py()
                        except ValueError:
                            # Not found in this group
                            pass

    except Exception as e:
        print(f"Error reading parquet: {e}")
        import traceback
        traceback.print_exc()
        return None

    return None


def read_spectra_gz(url,
                    targetid=None,
                    fiber=None,
                    expid=None,
                    mask=False,
                    ivar=False,
                    fibermap=False,
                    user=None,
                    pwd=None,
                    dataset=None):
    """
    Modified read_spectra that handles .gz files using indexed_gzip if index is available.
    """
    filename = url.split('/')[-1]

    # Fetch index specific to this file
    index_bytes = get_index_for_file(filename, dataset=dataset)

    if index_bytes is None:
        print(
            f"No index found for {filename}, falling back to standard read (might be slow or fail for .gz)"
        )
        raise Exception(
            'the gzip file cannot be read remotely without index, and'
            ' the index file was not found')
        pass

    kw = dict(verify=False)
    if user is not None:
        kw['auth'] = (user, pwd)
    kw['block_size'] = 2880 * 10

    # Create a temporary file for the index
    # indexed_gzip needs the index on disk or as a file object.
    # We'll use a named temporary file.

    index_file = None
    if index_bytes:
        index_file = tempfile.NamedTemporaryFile(delete=False)
        index_file.write(index_bytes)
        index_file.close()

    try:
        print(f"Opening {url} with indexed_gzip...")
        with httpio.open(url, **kw) as fp:

            if index_bytes:
                # Wrap with IndexedGzipFile
                # drop_handles=False keeps the underlying httpio stream open
                gz_fp = indexed_gzip.IndexedGzipFile(
                    fileobj=fp, index_file=index_file.name, drop_handles=False)
                fits_fp = gz_fp
            else:
                fits_fp = fp

            try:
                hdus = pyfits.open(fits_fp)
                return process_spectra_hdus(hdus,
                                            targetid=targetid,
                                            fiber=fiber,
                                            expid=expid,
                                            mask=mask,
                                            ivar=ivar,
                                            fibermap=fibermap)

            finally:
                if index_bytes:
                    gz_fp.close()  # Close wrapper
                    pass

    finally:
        if index_file:
            os.unlink(index_file.name)
