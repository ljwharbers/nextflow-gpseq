# Container specifications

The pipeline's custom images were originally hosted on the Sylabs Cloud Library
(`library://ljwharbers/gpseq/...`), which was shut down in October 2026. They were
rebuilt as Seqera community containers (https://seqera.io/containers) from the conda
specifications below, using `format: sif`, `linux/amd64`, `freeze: true`.

| Replaces | Spec | Image |
|---|---|---|
| `fastx-barber:0.0.4`, `gpseq_pyenv:0.0.1` | `python.json` | `oras://community.wave.seqera.io/library/python_pip_biopython_joblib_pruned:3c7e81d08cfa618a` |
| `gpseq_renv:0.0.3`, `gpseq_renv:0.0.4` | `r.json` | `oras://community.wave.seqera.io/library/r-base_r-argparser_r-data.table_r-logging_pruned:de3cd9d98b9a253e` |

`main.nf` references the HTTPS blob URL of each SIF so that any Nextflow version can pull it.

To rebuild (e.g. after adding a package), POST the spec to the Wave API and update `main.nf`:

```bash
curl -s -X POST -H 'Content-Type: application/json' -d @r.json https://wave.seqera.io/v1alpha2/container
```

Note: fastx-barber is only available on PyPI (latest 0.1.5); its pinned dependencies
require Python 3.9.
