from reflectao.schema import SCHEMA, new_empty_observation_table, schema_column_names


def test_new_empty_observation_table_has_all_columns():
    tbl = new_empty_observation_table(n_rows=0)
    assert tbl.colnames == schema_column_names(SCHEMA)


def test_new_empty_observation_table_initializes_none_values():
    tbl = new_empty_observation_table(n_rows=2)
    for name in schema_column_names(SCHEMA):
        assert tbl[name].mask[0] == True
        assert tbl[name].mask[1] == True

