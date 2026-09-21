"""
test_signature_validation.py -- the 10X signature layout check.

WHY THIS MATTERS. The slice points that pull the barcode and UMI out of a
matched block encode a fixed layout: a 22 nt 5' handle, a 16 nt barcode, a 10 nt
3' handle, and a UMI whose width is DERIVED as len(signature) - 48.

That derivation is a feature -- it is what lets one signature string carry the
whole chemistry, so switching from 10X v2 to v3 is the signature plus a matching
whitelist and nothing else. It is also a load-bearing assumption, and before this
check nothing enforced it. A signature with the wrong handle length was accepted
silently and produced a UMI sliced at the wrong offsets: mostly barcode, still
consistent, so reads still grouped -- just into the wrong groups, collapsing
distinct molecules or splitting one apart. No error, no way to tell downstream.
"""

from __future__ import annotations

import pytest

from anchovy.barcodes import (
    build_barcode_query_blocks,
    extract_umi,
    validate_signature,
)
from anchovy.config import DEFAULT_TENX_SIGNATURE, ExtractConfig
from anchovy.schema import (
    SIGNATURE_BARCODE_LEN,
    SIGNATURE_NON_UMI_LEN,
    SIGNATURE_PREFIX_LEN,
    SIGNATURE_SUFFIX_LEN,
)

PREFIX = "CTACACGACGCTCTTCCGATCT"      # 22 nt
SUFFIX = "TTTCTTATAT"                  # 10 nt


def _signature(umi_len: int, prefix: str = PREFIX, suffix: str = SUFFIX) -> str:
    return prefix + "N" * (SIGNATURE_BARCODE_LEN + umi_len) + suffix


# --------------------------------------------------------------------------- #
# The layout constants describe the real signature
# --------------------------------------------------------------------------- #
def test_constants_describe_the_shipped_signature():
    """The named constants must match the default, or they document a fiction."""
    assert DEFAULT_TENX_SIGNATURE.startswith(PREFIX)
    assert DEFAULT_TENX_SIGNATURE.endswith(SUFFIX)
    assert len(PREFIX) == SIGNATURE_PREFIX_LEN
    assert len(SUFFIX) == SIGNATURE_SUFFIX_LEN
    assert DEFAULT_TENX_SIGNATURE.count("N") == SIGNATURE_BARCODE_LEN + 10   # v2


def test_constants_agree_with_the_extract_config_offsets():
    """schema's layout and config's UMI offsets are two views of one truth."""
    c = ExtractConfig()
    assert c.umi_start_offset == SIGNATURE_PREFIX_LEN + SIGNATURE_BARCODE_LEN
    assert c.umi_end_trim == SIGNATURE_SUFFIX_LEN


# --------------------------------------------------------------------------- #
# Valid signatures pass
# --------------------------------------------------------------------------- #
def test_the_shipped_default_is_valid():
    validate_signature(DEFAULT_TENX_SIGNATURE)          # must not raise


@pytest.mark.parametrize("umi_len,chemistry", [(10, "v2"), (12, "v3")])
def test_both_10x_chemistries_are_valid(umi_len, chemistry):
    validate_signature(_signature(umi_len))


def test_lowercase_is_accepted():
    """run() upper-cases before validating; the check should not depend on that."""
    validate_signature(_signature(10).lower())


@pytest.mark.parametrize("umi_len", [10, 12])
def test_umi_width_follows_the_signature(umi_len):
    """The property the whole design rests on: chemistry rides on the signature.

    Changing the signature must change the UMI width with no other edits, which
    is why no `chemistry: v2|v3` setting is needed.
    """
    signature = _signature(umi_len)
    validate_signature(signature)
    config = ExtractConfig()

    block = build_barcode_query_blocks(signature, ["A" * SIGNATURE_BARCODE_LEN])[0]
    assert block.count("N") == umi_len
    assert len(signature) - SIGNATURE_NON_UMI_LEN == umi_len

    truth = ("ACGT" * 4)[:umi_len]
    matchseq = PREFIX + "A" * SIGNATURE_BARCODE_LEN + truth + SUFFIX
    assert extract_umi(matchseq, config.umi_start_offset, config.umi_end_trim) == truth


# --------------------------------------------------------------------------- #
# Malformed signatures are refused
# --------------------------------------------------------------------------- #
def test_short_prefix_is_refused():
    """THE REGRESSION: silently sliced a mostly-barcode 'UMI' before this check."""
    with pytest.raises(ValueError, match="5' handle is 4 nt"):
        validate_signature("ACGT" + "N" * 26 + SUFFIX)


def test_long_prefix_is_refused():
    with pytest.raises(ValueError, match="5' handle is 30 nt"):
        validate_signature(_signature(10, prefix="A" * 30))


def test_wrong_suffix_length_is_refused():
    with pytest.raises(ValueError, match="3' handle is 4 nt"):
        validate_signature(_signature(10, suffix="ACGT"))


def test_n_run_too_short_for_a_umi_is_refused():
    """16 Ns is barcode and nothing else -- there is no UMI to slice."""
    with pytest.raises(ValueError, match="leaves no UMI"):
        validate_signature(PREFIX + "N" * SIGNATURE_BARCODE_LEN + SUFFIX)


def test_signature_without_any_ns_is_refused():
    with pytest.raises(ValueError, match="no N-run"):
        validate_signature(PREFIX + "ACGTACGT" + SUFFIX)


def test_non_contiguous_ns_are_refused():
    """Ns split across the handles would make the slice points meaningless."""
    with pytest.raises(ValueError, match="not contiguous"):
        validate_signature(PREFIX + "N" * 13 + "ACGT" + "N" * 13 + SUFFIX)


def test_the_error_says_what_was_wrong_and_what_it_would_have_done():
    """A bad signature is a config mistake; the message has to be actionable."""
    with pytest.raises(ValueError) as excinfo:
        validate_signature("ACGT" + "N" * 26 + SUFFIX)
    message = str(excinfo.value)
    assert "5' handle is 4 nt, expected 22" in message
    assert "v2 is 26 Ns" in message and "v3 is 28" in message
    # and it names the consequence, not just the mismatch
    assert "silently mis-group" in message


# --------------------------------------------------------------------------- #
# It actually runs, at the entry point
# --------------------------------------------------------------------------- #
def test_extract_run_rejects_a_bad_signature(data_dir, tmp_path):
    """Fail at the front door, not three stages downstream."""
    from anchovy import extract

    sam = data_dir / "test.sam"
    whitelist = data_dir / "whitelist.txt"
    for p in (sam, whitelist):
        if not p.exists():
            pytest.skip(f"missing {p.name}")

    with pytest.raises(ValueError, match="does not have the layout"):
        extract.run(sam=str(sam), whitelist=str(whitelist),
                    signature="ACGT" + "N" * 26 + SUFFIX)
