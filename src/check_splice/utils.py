import os

import pandas as pd
import pyarrow as pa
from pyarrow import ipc


def jsonl2feather(jsonl_file: os.PathLike, feather_file: os.PathLike):
    writer = None
    # Stream JSON in chunks of 10,000 rows
    for chunk in pd.read_json(jsonl_file, chunksize=100000, lines=True):
        batch = pa.RecordBatch.from_pandas(chunk)
        if writer is None:
            writer = ipc.RecordBatchFileWriter(feather_file, batch.schema)
        writer.write_batch(batch)

    if writer is not None:
        writer.close()

    df = pd.read_feather(feather_file)
    df.to_feather(feather_file)
