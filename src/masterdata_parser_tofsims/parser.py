from __future__ import annotations

import re
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any

import bam_masterdata.datamodel.vocabulary_types  # noqa: F401 - triggers vocabulary registration
from bam_masterdata.datamodel.mass_spectrometry.object_types import SIMS, MassSpec

# NOTE (2026-08-25): in bam-masterdata 0.11.5, the ToF-SIMS experimental
# step class is called `SIMS` (not `ExperimentalStep_SIMS`) and, together
# with `MassSpec`, lives in the mass_spectrometry submodule - not the
# top-level object_types module. `Sample` is still top-level. Confirmed
# directly against the installed package source.
from bam_masterdata.datamodel.object_types import Sample
from bam_masterdata.parsing import AbstractParser
from pySPM.Block import Block, MissingBlockError

# ------------------------------------------------------------------
# KONFIGURATION
# ------------------------------------------------------------------

GROUP_EXTENSIONS = {".itm", ".ita", ".itax", ".itmx", ".txt"}
METADATA_EXTENSIONS = {".itm", ".ita", ".itax", ".itmx"}
METADATA_PREFERENCE = [".itm", ".itmx", ".ita", ".itax"]
NO_UPLOAD_SUFFIX = "_no-upload"

# All Sample objects go into this collection (inside the project selected in
# the upload helper UI; it is created if it does not exist yet). The
# measurements (SIMS steps) stay in the collection selected in the UI.
SAMPLES_COLLECTION = "SAMPLES"
# True  -> if the runner cannot be hooked, write nothing (so nothing ends up
#          in the wrong collection).
# False -> fall back to the UI-selected collection for samples and only warn.
SAMPLES_COLLECTION_STRICT = True

PRIMARY_ION_CODES = {
    "BI1",
    "BI2",
    "BI3",
    "BI4",
    "BI5",
    "BI6",
    "BI7",
    "BI1PP",
    "BI3PP",
    "BI5PP",
    "BI7PP",
}
SPUTTER_ION_MAP = {"o2": "OXYGEN", "cs": "CESIUM"}
POLARITY_MAP = {
    "positive": "POSITIVE",
    "negative": "NEGATIVE",
    "+": "POSITIVE",
    "-": "NEGATIVE",
}

PROPERTY_PATH_TEMPLATES = [
    "propend/{key}",
    "propstart/{key}",
    "CommonDataObjects/DataViewCollection/*/properties/{key}",
]

# Business-metadata fields we need per group, and the instrument-property
# keys to try for each, in order. All CONFIRMED against a real, filled-in
# .itax file (they live under a "User.*" branch that is entirely absent from
# the property list when never filled in - that's why it didn't show up in
# the first, empty test file). Analysis.Description kept as a fallback for
# "sample" since it's a separate, also-confirmed field with the same value.
#
# CAVEAT (found 2026-08-25): on *derived/processed* exports (e.g. files
# produced by applying a "Stream Filter" or similar operation in SurfaceLab),
# User.SampleName is often not carried over, and Analysis.Description instead
# contains SurfaceLab's auto-generated *operation name* (e.g. "Stream
# Filter") rather than the real sample name. The fallback therefore only
# fires when User.SampleName is not found at all - and its value is checked
# against KNOWN_NON_SAMPLE_DESCRIPTIONS below to catch this case.
BUSINESS_FIELD_KEYS = {
    "sample": ["User.SampleName", "Analysis.Description"],
    "operator": ["User.Operator"],
    "customer": ["User.SampleOrigin"],
    "comment": ["User.Comment"],
}

# Auto-generated SurfaceLab processing/operation names that sometimes end up
# in Analysis.Description on derived exports. If the "sample" fallback
# resolves to one of these (case-insensitive), treat it as NOT a real value -
# i.e. as if no key had matched - rather than silently naming the Sample
# after a processing step. Extend this list as new false positives turn up.
KNOWN_NON_SAMPLE_DESCRIPTIONS = {
    "stream filter",
    "depth profile",
}

# Placeholder strings SurfaceLab writes for "field left empty".
EMPTY_PLACEHOLDERS = {
    "<no sample name>",
    "<no description>",
    "<no comment>",
    "",
    "none",
}

# Technical instrument properties copied onto the ExpSims step.
INSTRUMENT_KEYS = {
    "Analysis.Timestamp": "Date/Time",
    "Instrument.PrimaryGun.Species": "Primary ion",
    "Instrument.SputterGun.Species": "Sputter ion",
    "Instrument.Analyzer.Polarity": "Polarity",
    "Registration.Raster.FieldOfView": "Image size [µm]",
    "Instrument.SputterGun.Energy": "Sputter V [kV]",
    "Analysis.SputterTime": "Sputter time [s]",
    "Profile.CraterSize.X": "Krater size X [µm]",
    "Profile.CraterSize.Y": "Krater size Y [µm]",
    "Storage.File.Name": "Original filename",
    "Storage.File.Path": "Original file path",
}

DEVICE_PERM_ID = "20260624132605462-55463"
BAM_OE = "OE_6.1"
DEVICE_NAME = "ToF SIMS"
DEVICE_MANUFACTURER = "IONTOF"
DEVICE_LOCATION_COMPLETE = "FB/80/0/148"

# For batches of old files where SurfaceLab never recorded these -
# leave both as "" for normal/new batches. Edit before running a batch
# of legacy files, since one upload batch is always the same person/customer.
LEGACY_OPERATOR_FALLBACK = "Elisabeth John"
LEGACY_CUSTOMER_FALLBACK = "Simon Schroeder"
LEGACY_SAMPLE_FALLBACK = "Bitumen+Rejuvenator+Polymer_V8"
# ------------------------------------------------------------------
# MATCH META DATA WITH OPENBIS CODE OF CONTROLLED VOCABS
# ------------------------------------------------------------------


def map_primary_ion(raw: str | None) -> str | None:
    if not raw:
        return None
    m = re.match(r"Bi(\d*)(\+{1,2})$", raw.strip())
    if not m:
        return None
    n = int(m.group(1)) if m.group(1) else 1
    code = f"BI{n}" + ("PP" if len(m.group(2)) == 2 else "")
    return code if code in PRIMARY_ION_CODES else None


# ------------------------------------------------------------------
# LOW-LEVEL: read ION-TOF instrument property files with pySPM
# ------------------------------------------------------------------


def _open_root_block(path: Path):
    f = open(path, "rb")
    magic = f.read(8)
    if magic != b"ITStrF01":
        f.close()
        raise ValueError(f"Keine ION-TOF SurfaceLab-Containerdatei: {path}")
    return f, Block(f)


def _get_property(root: Block, key: str) -> dict | None:
    for template in PROPERTY_PATH_TEMPLATES:
        path = template.format(key=key)
        try:
            return root.goto(path, lazy=True).get_key_value()
        except MissingBlockError:
            continue
    return None


def _first_resolved(
    root: Block, candidate_keys: list[str]
) -> tuple[str | None, str | None]:
    """
    Tries each candidate key in order; returns (value, key_that_matched).
    Business fields are always text, so we use the 'string' entry as-is
    (including when it's an empty string) and never fall back to 'float' -
    doing so previously turned an empty User.Comment ('') into '0.0'.
    """
    for key in candidate_keys:
        kv = _get_property(root, key)
        if kv is not None:
            value = kv.get("string")
            return value, key
    return None, None


def is_empty_value(value: Any) -> bool:
    if value is None:
        return True
    text = str(value).strip().lower()
    return text in EMPTY_PLACEHOLDERS


def is_non_sample_description(value: Any) -> bool:
    """True if `value` looks like a SurfaceLab auto-generated processing/
    operation name (e.g. 'Stream Filter') rather than a real sample name."""
    if value is None:
        return False
    return str(value).strip().lower() in KNOWN_NON_SAMPLE_DESCRIPTIONS


def parse_timestamp(raw: str | None) -> str | None:
    if not raw:
        return None
    try:
        return datetime.strptime(raw, "%a %b %d %H:%M:%S %Y").strftime(
            "%Y-%m-%d %H:%M:%S"
        )
    except ValueError:
        return raw


def extract_metadata(path: Path, logger) -> tuple[dict[str, str], dict[str, Any]]:
    """
    Returns (business_fields, technical_fields) for one instrument file.
    business_fields: {"sample": ..., "operator": ..., "customer": ..., "comment": ...}
      - value is None if no candidate key resolved at all (field truly absent)
      - value is "" if a key resolved but SurfaceLab reports it as empty/placeholder
    technical_fields: flat dict of INSTRUMENT_KEYS values, missing keys omitted.
    """
    f, root = _open_root_block(path)
    try:
        business: dict[str, str | None] = {}
        for field, candidate_keys in BUSINESS_FIELD_KEYS.items():
            value, matched_key = _first_resolved(root, candidate_keys)

            # Guard against Analysis.Description containing a SurfaceLab
            # processing-operation name (e.g. "Stream Filter") instead of a
            # real sample name on derived/processed exports. Treat this the
            # same as "key not found" so it doesn't silently become the
            # Sample's name in openBIS.
            if (
                field == "sample"
                and matched_key == "Analysis.Description"
                and is_non_sample_description(value)
            ):
                logger.error(
                    f"'{path.name}': Feld 'sample' - 'User.SampleName' fehlt, und "
                    f"'Analysis.Description' enthält nur den Namen einer "
                    f"SurfaceLab-Verarbeitungsoperation ('{value}'), keinen echten "
                    "Probennamen. Dies ist vermutlich eine abgeleitete/verarbeitete "
                    "Export-Datei ohne Sample-Metadaten - Feld wird als fehlend behandelt."
                )
                value, matched_key = None, None

            if matched_key is None:
                logger.warning(
                    f"'{path.name}': keiner der bekannten Keys für '{field}' "
                    f"({candidate_keys}) existiert in dieser Datei (oder wurde verworfen)."
                )
                business[field] = None
            elif is_empty_value(value):
                logger.error(
                    f"'{path.name}': Feld '{field}' (Key '{matched_key}') ist im "
                    "Instrument-File leer."
                )
                business[field] = ""
            else:
                business[field] = str(value).strip()

        technical: dict[str, Any] = {}

        for raw_key, mapped_name in INSTRUMENT_KEYS.items():
            kv = _get_property(root, raw_key)
            if kv is None:
                continue
            value = kv.get("string") or kv.get("float")
            technical[mapped_name] = value
        if "Date/Time" in technical:
            technical["Date/Time"] = parse_timestamp(technical["Date/Time"])
        logger.info(
            f"DEBUG technical keys: {list(technical.keys())}"
        )  # ---------------------------------------------------------------------------
        logger.info(f"DEBUG technical values: {technical}")
        return business, technical

    finally:
        f.close()


# ------------------------------------------------------------------
# GROUPING
# ------------------------------------------------------------------


def group_key(path: Path) -> str:
    stem = path.stem
    if stem.endswith(NO_UPLOAD_SUFFIX):
        stem = stem[: -len(NO_UPLOAD_SUFFIX)]
    return stem


def is_no_upload(path: Path) -> bool:
    return path.stem.endswith(NO_UPLOAD_SUFFIX)


def group_files(files: list[str], logger) -> dict[str, list[Path]]:
    groups: dict[str, list[Path]] = defaultdict(list)
    ignored = []
    for raw in files:
        p = Path(raw)
        if p.suffix.lower() not in GROUP_EXTENSIONS:
            ignored.append(p.name)
            continue
        groups[group_key(p)].append(p)
    if ignored:
        logger.info(
            f"{len(ignored)} Datei(en) mit unbekannter Endung ignoriert: {ignored}"
        )
    logger.info(f"{len(groups)} Messgruppe(n) aus {len(files)} Datei(en) gebildet.")
    return groups


# ------------------------------------------------------------------
# CODE SANITIZING
# ------------------------------------------------------------------


def sanitize_code(name: str, max_length: int = 200) -> str:
    code = name.strip().upper()
    code = re.sub(r"[^A-Z0-9_\-\.]", "_", code)
    code = re.sub(r"_+", "_", code)
    code = code.strip("_.-")
    return code[:max_length]


# ------------------------------------------------------------------
# PARSER
# ------------------------------------------------------------------


class MasterdataParserTofsims(AbstractParser):
    """One Sample (in the SAMPLES collection) + one SIMS Experimental Step
    (in the collection selected in the UI, child of the Sample) per file group."""

    def parse(self, files: list[str], collection, logger):
        groups = group_files(files, logger)
        if not groups:
            logger.warning(
                f"Keine Dateien mit bekannter ToF-SIMS-Endung ({sorted(GROUP_EXTENSIONS)}) "
                "in diesem Batch gefunden - nichts zu tun."
            )
            return

        sample_ids: dict[
            str, str
        ] = {}  # sample code -> object id, dedup within this run
        created, skipped, failed = 0, 0, 0

        for key, group in sorted(groups.items()):
            logger.info(f"--- Verarbeite Gruppe '{key}' ({len(group)} Datei(en)) ---")

            business, technical = self._extract_group_metadata(group, logger)
            if not business and not technical:
                logger.error(
                    f"Gruppe '{key}': keine lesbare Metadaten-Datei gefunden - Gruppe wird übersprungen."
                )
                skipped += 1
                continue

            if not business.get("operator") and LEGACY_OPERATOR_FALLBACK:
                business["operator"] = LEGACY_OPERATOR_FALLBACK
                logger.info(
                    f"Gruppe '{key}': Operator aus LEGACY_OPERATOR_FALLBACK übernommen."
                )
            if not business.get("customer") and LEGACY_CUSTOMER_FALLBACK:
                business["customer"] = LEGACY_CUSTOMER_FALLBACK
                logger.info(
                    f"Gruppe '{key}': Customer aus LEGACY_CUSTOMER_FALLBACK übernommen."
                )
            if not business.get("sample") and LEGACY_SAMPLE_FALLBACK:
                business["sample"] = LEGACY_SAMPLE_FALLBACK
                logger.info(
                    f"Gruppe '{key}': Sample aus LEGACY_SAMPLE_FALLBACK übernommen."
                )

            sample_name = business.get("sample")
            customer = business.get("customer")

            missing = [
                label
                for label, value in [("Sample", sample_name), ("Customer", customer)]
                if not value
            ]
            if missing:
                logger.error(
                    f"Gruppe '{key}': folgende Pflichtfelder fehlen oder sind leer: "
                    f"{missing} - Gruppe wird übersprungen."
                )
                skipped += 1
                continue

            try:
                sample_code = sanitize_code(sample_name)
                if sample_code not in sample_ids:
                    sample = Sample(
                        name=sample_name
                    )  # VERIFY: "name" is a real property on Sample
                    sample.code = (
                        sample_code  # deterministic -> re-runs update, not duplicate
                    )
                    sample.bam_oe = BAM_OE
                    sample_id = collection.add(sample)
                    sample_ids[sample_code] = sample_id
                    logger.info(
                        f"Sample '{sample_name}' (Code {sample_code}) angelegt/aktualisiert (Collection 'UI-Auswahl')."
                    )
                sample_id = sample_ids[sample_code]

                roi_label = key  # group base name; adjust if you have a nicer ROI label
                step = SIMS(name=roi_label)
                step.code = sanitize_code(f"{sample_code}_{key}")

                if business.get("operator"):
                    step.operator = business["operator"]
                step.customer = customer

                start_date = technical.get("Date/Time")
                if not start_date:
                    logger.error(
                        f"Gruppe '{key}': Analysis.Timestamp fehlt - 'start_date' ist Pflichtfeld, Gruppe wird übersprungen."
                    )
                    skipped += 1
                    continue
                step.start_date = start_date

                # Confirmed 2026-08-25 directly against the installed
                # bam-masterdata 0.11.5 source: the SIMS class's Python
                # attributes are sims_primary_ion / sims_polarity /
                # sims_sputter_source, matching the openBIS server's codes.
                primary_ion = map_primary_ion(technical.get("Primary ion"))
                if primary_ion:
                    step.sims_primary_ion = primary_ion
                else:
                    logger.warning(
                        f"Gruppe '{key}': Primary ion '{technical.get('Primary ion')}' nicht erkannt - übersprungen."
                    )

                sputter_source = SPUTTER_ION_MAP.get(
                    (technical.get("Sputter ion") or "").strip().lower()
                )
                if sputter_source:
                    step.sims_sputter_source = sputter_source
                else:
                    logger.warning(
                        f"Gruppe '{key}': Sputter ion '{technical.get('Sputter ion')}' nicht erkannt - übersprungen."
                    )

                polarity = POLARITY_MAP.get(
                    (technical.get("Polarity") or "").strip().lower()
                )
                if polarity:
                    step.sims_polarity = polarity
                else:
                    logger.warning(
                        f"Gruppe '{key}': Polarity '{technical.get('Polarity')}' nicht erkannt - übersprungen."
                    )

                upload_files = [str(p) for p in group if not is_no_upload(p)]
                for f in upload_files:
                    step.add_dataset(f)
                if upload_files:
                    logger.info(
                        f"{len(upload_files)} Datei(en) an '{step.code}' angehängt."
                    )
                else:
                    logger.info(
                        f"Gruppe '{key}': alle Dateien '{NO_UPLOAD_SUFFIX}' - keine Datei angehängt."
                    )

                step_id = collection.add(step)
                collection.add_relationship(sample_id, step_id)

                # The device is a fixed, pre-existing instrument in openBIS — it must
                # never be created or updated by this parser. We reference it purely by
                # its known identifier and let openBIS itself reject the relationship
                # if that identifier doesn't actually exist.
                collection.add_relationship({"permId": DEVICE_PERM_ID}, step_id)
                created += 1

            except Exception as exc:
                logger.error(f"Gruppe '{key}': Fehler beim Anlegen der Objekte - {exc}")
                failed += 1

        logger.info(
            f"Parsing abgeschlossen: {created} Messung(en) angelegt, "
            f"{skipped} übersprungen, {failed} fehlgeschlagen."
        )

    def _extract_group_metadata(
        self, group: list[Path], logger
    ) -> tuple[dict[str, str], dict[str, Any]]:
        by_ext = {p.suffix.lower(): p for p in group}
        for ext in METADATA_PREFERENCE:
            path = by_ext.get(ext)
            if path is None:
                continue
            try:
                return extract_metadata(path, logger)
            except Exception as exc:
                logger.warning(f"Konnte Metadaten nicht aus '{path.name}' lesen: {exc}")
        logger.error(
            f"Keine lesbare Metadaten-Datei in Gruppe gefunden: {[p.name for p in group]}"
        )
        return {}, {}
