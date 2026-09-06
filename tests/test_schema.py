import pandas as pd
import pytest

from contextshift.schema import MEMBERS, Column, SchemaError, TableSchema

TINY = TableSchema(
    name="tiny",
    columns=(Column("id", "string"), Column("v", "float64", nullable=True)),
    key=("id",),
)


def test_missing_required_column():
    with pytest.raises(SchemaError, match="missing required columns"):
        TINY.validate(pd.DataFrame({"v": [1.0]}))


def test_duplicate_key_rejected():
    with pytest.raises(SchemaError, match="duplicate key"):
        TINY.validate(pd.DataFrame({"id": ["a", "a"], "v": [1.0, 2.0]}))


def test_non_nullable_column_rejects_nulls():
    with pytest.raises(SchemaError, match="null values"):
        TINY.validate(pd.DataFrame({"id": ["a", None], "v": [1.0, 2.0]}))


def test_nullable_column_accepts_nulls():
    out = TINY.validate(pd.DataFrame({"id": ["a"], "v": [None]}))
    assert out["v"].isna().all()


def test_optional_columns_may_be_absent():
    df = pd.DataFrame({"member_id": ["a"], "genome_id": ["g"], "family": ["f"]})
    assert len(MEMBERS.validate(df)) == 1


def test_parquet_roundtrip(tmp_path):
    df = pd.DataFrame({"id": ["a", "b"], "v": [1.0, 2.0]})
    path = tmp_path / "t.parquet"
    TINY.write(df, path)
    assert TINY.read(path).equals(TINY.validate(df))
