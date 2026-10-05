# ToF-SIMS openBIS Parser

Parser for ION-TOF ToF-SIMS SurfaceLab files that creates `Sample` and `SIMS` objects for openBIS using the `bam-masterdata` parser interface.

The parser is designed to be used with the [`openbis-upload-helper`](https://github.com/BAMresearch/openbis-upload-helper) and follows the parser structure of the [`openbis-parser-example`](https://github.com/BAMresearch/openbis-parser-example).

## Overview

The parser reads ION-TOF SurfaceLab container files, extracts sample/business metadata and technical measurement metadata, and creates the corresponding openBIS objects.

For each detected file group, the parser creates:

- one `Sample`
- one `SIMS` experimental step
- datasets for the uploaded source files
- relationships from the SIMS measurement to the sample and the existing ToF-SIMS device

The parser expects ION-TOF SurfaceLab container files.

## Supported files

The following file extensions are considered:

| Extension | Purpose |
|---|---|
| `.itm` | SurfaceLab metadata / instrument file |
| `.itmx` | SurfaceLab metadata / instrument file |
| `.ita` | SurfaceLab metadata / instrument file |
| `.itax` | SurfaceLab metadata / instrument file |
| `.txt` | File belonging to a measurement group |

Metadata files are preferred in this order:

1. `.itm`
2. `.itmx`
3. `.ita`
4. `.itax`

Files with other extensions are ignored.

## File grouping

Files are grouped by their filename stem.

For example:

```text
sample01.itm
sample01.itmx
sample01.txt
sample01.itax
```

are treated as one measurement group:

```text
sample01
```

### Technical metadata

The parser maps the following SurfaceLab properties:

| SurfaceLab property | SIMS property |
|---|---|
| `Analysis.Timestamp` | Date/Time |
| `Instrument.PrimaryGun.Species` | Primary ion |
| `Instrument.SputterGun.Species` | Sputter ion |
| `Instrument.Analyzer.Polarity` | Polarity |
| `Registration.Raster.FieldOfView` | Image size [µm] |
| `Instrument.SputterGun.Energy` | Sputter V [kV] |
| `Analysis.SputterTime` | Sputter time [s] |
| `Profile.CraterSize.X` | Krater size X [µm] |
| `Profile.CraterSize.Y` | Krater size Y [µm] |

`Analysis.Timestamp` is converted from the SurfaceLab timestamp format to:

```text
YYYY-MM-DD HH:MM:SS
```


For each valid measurement group, the parser creates a `Sample` and a `SIMS` experimental step.

### Sample

The sample name is taken from the SurfaceLab metadata and converted into an openBIS code by:

- converting to uppercase
- replacing unsupported characters with `_`
- collapsing repeated underscores
- removing leading/trailing `_.-`
- limiting the code to 200 characters

The sample is added to the collection passed to the parser.

### SIMS experimental step

The SIMS object uses the group name as its name and receives a deterministic code:

```text
<SAMPLE_CODE>_<GROUP_NAME>
```

The following properties are assigned when available:

- operator
- customer
- start date
- primary ion
- sputter ion
- polarity

All files in the group except files ending in `_no-upload` are added as datasets.

### Relationships

The parser creates:

```text
Sample
  └── SIMS measurement
        └── datasets
```

The SIMS measurement is also linked to the pre-existing ToF-SIMS device.

The device is **not created or modified by the parser**. The parser references the existing device by its configured permanent ID.

## Required metadata

A group is skipped if no readable metadata file can be found.

The following business fields are required:

- Sample
- Customer

The `Analysis.Timestamp` is also required because it is used as the SIMS `start_date`.

If one of these required values is missing, the group is skipped and an error is written to the parser log.

## Legacy fallback values

The parser currently contains fallback values for older SurfaceLab files where certain metadata was not recorded:

```python
LEGACY_OPERATOR_FALLBACK = "..."
LEGACY_CUSTOMER_FALLBACK = "..."
LEGACY_SAMPLE_FALLBACK = "..."
```

These values should be reviewed before processing a legacy batch. They are intended for batches where the missing metadata is known to be identical for all files in that batch.

For normal/new data, these fallbacks should be empty.

## Installation

Install the parser together with its dependencies in a Python environment.

The parser requires, among other dependencies:

- `bam-masterdata`
- `pySPM`

The parser implements the `bam-masterdata` `AbstractParser` interface and should be registered through the `bam.parsers` Python entry-point group.

A typical entry-point configuration looks like:

```toml
[project.entry-points."bam.parsers"]
tofsims = "<package-name>:<entry-point-variable>"
```

See the [`openbis-parser-example`](https://github.com/BAMresearch/openbis-parser-example) repository for the expected parser package structure and entry-point configuration.

## Using the parser with openbis-upload-helper

Once installed and registered as a `bam.parsers` entry point, the parser can be discovered by `openbis-upload-helper`.

In the upload helper:

1. Select the target openBIS project.
2. Select the collection into which the measurements should be uploaded.
3. Add the ToF-SIMS files.
4. Assign the ToF-SIMS parser to the files or folder.
5. Start the parsing/upload process.
6. Check the parser log for skipped or failed measurement groups.


## Local testing

The repository provides a local runner/notebook for testing the parser independently of the desktop upload workflow.

A minimal invocation follows the parser interface:

```python
parser = MasterdataParserTofsims()
from bam_masterdata.logger import logger
from bam_masterdata.metadata.entities import CollectionType

collection = CollectionType()
parser.parse(
    files,
    collection,
    logger,
)
```

The local runner needs to provide:

- a list of input file paths
- a collection object
- a logger

For the repository's local parser runner, use the run notebook/script added for this parser. This is useful for testing a batch of SurfaceLab files before using the parser through `openbis-upload-helper`.

