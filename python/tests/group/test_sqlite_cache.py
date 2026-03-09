import libcasm.configuration as casmconfig
import libcasm.group as casmgroup
import libcasm.xtal.prims as xtal_prims
from libcasm.group._get_subgroup_orbits import get_all_subgroup_orbits
from libcasm.group.sqlite_cache import LocalCache, UserCache


def _fcc_factor_group():
    prim = casmconfig.Prim(xtal_prim=xtal_prims.FCC(a=4.0, occ_dof=["A", "B"]))
    return prim.factor_group


# --- UserCache tests ---


def test_user_cache_store_and_get():
    cache = UserCache()
    cache.store("_test_cache", "key1", "value1")
    assert cache.get("_test_cache", "key1") == "value1"
    cache.clear_cache("_test_cache")


def test_user_cache_get_missing_key():
    cache = UserCache()
    cache.clear_cache("_test_cache")
    assert cache.get("_test_cache", "nonexistent") is None
    cache.clear_cache("_test_cache")


def test_user_cache_overwrite():
    cache = UserCache()
    cache.store("_test_cache", "key1", "value1")
    cache.store("_test_cache", "key1", "value2")
    assert cache.get("_test_cache", "key1") == "value2"
    cache.clear_cache("_test_cache")


def test_user_cache_delete():
    cache = UserCache()
    cache.store("_test_cache", "key1", "value1")
    cache.delete("_test_cache", "key1")
    assert cache.get("_test_cache", "key1") is None
    cache.clear_cache("_test_cache")


def test_user_cache_clear_cache():
    cache = UserCache()
    cache.store("_test_cache", "key1", "value1")
    cache.store("_test_cache", "key2", "value2")
    cache.clear_cache("_test_cache")
    assert cache.get("_test_cache", "key1") is None
    assert cache.get("_test_cache", "key2") is None


# --- LocalCache tests ---


def test_local_cache_store_and_get(tmp_path):
    cache = LocalCache(tmp_path)
    cache.store("mycache", "key1", "value1")
    assert cache.get("mycache", "key1") == "value1"


def test_local_cache_get_missing_key(tmp_path):
    cache = LocalCache(tmp_path)
    assert cache.get("mycache", "nonexistent") is None


def test_local_cache_overwrite(tmp_path):
    cache = LocalCache(tmp_path)
    cache.store("mycache", "key1", "value1")
    cache.store("mycache", "key1", "value2")
    assert cache.get("mycache", "key1") == "value2"


def test_local_cache_delete(tmp_path):
    cache = LocalCache(tmp_path)
    cache.store("mycache", "key1", "value1")
    cache.delete("mycache", "key1")
    assert cache.get("mycache", "key1") is None


def test_local_cache_clear_cache(tmp_path):
    cache = LocalCache(tmp_path)
    cache.store("mycache", "key1", "value1")
    cache.store("mycache", "key2", "value2")
    cache.clear_cache("mycache")
    assert cache.get("mycache", "key1") is None
    assert cache.get("mycache", "key2") is None


def test_local_cache_db_file_location(tmp_path):
    cache = LocalCache(tmp_path)
    cache.store("mycache", "key1", "value1")
    db_path = tmp_path / "sqlite_cache" / "mycache.db"
    assert db_path.exists()


def test_local_cache_independent_caches(tmp_path):
    cache = LocalCache(tmp_path)
    cache.store("cache_a", "key1", "value_a")
    cache.store("cache_b", "key1", "value_b")
    assert cache.get("cache_a", "key1") == "value_a"
    assert cache.get("cache_b", "key1") == "value_b"
    cache.clear_cache("cache_a")
    assert cache.get("cache_a", "key1") is None
    assert cache.get("cache_b", "key1") == "value_b"


# --- get_all_subgroup_orbits / UserCache integration tests ---


def test_get_all_subgroup_orbits_returns_list():
    fg = _fcc_factor_group()
    orbits = get_all_subgroup_orbits(fg)
    assert isinstance(orbits, list)
    assert len(orbits) > 0
    # Each orbit is a list of subgroups; each subgroup is a list of int indices
    for orbit in orbits:
        assert isinstance(orbit, list)
        for subgroup in orbit:
            assert isinstance(subgroup, list)
            assert all(isinstance(i, int) for i in subgroup)


def test_get_all_subgroup_orbits_cached(monkeypatch):
    """Second call returns cached result without recomputing."""
    fg = _fcc_factor_group()

    # Warm the cache
    first = get_all_subgroup_orbits(fg)

    # Patch all_subgroups to detect if it is called again
    called = []
    original_all_subgroups = casmgroup.Subset.all_subgroups

    def spy_all_subgroups(self):
        called.append(True)
        return original_all_subgroups(self)

    monkeypatch.setattr(casmgroup.Subset, "all_subgroups", spy_all_subgroups)

    second = get_all_subgroup_orbits(fg)

    assert second == first
    assert called == [], "all_subgroups should not be called on a cache hit"


def test_get_all_subgroup_orbits_cache_cleared():
    """After clearing the cache, result is recomputed and equals the original."""
    fg = _fcc_factor_group()

    first = get_all_subgroup_orbits(fg)

    ucache = UserCache()
    ucache.clear_cache("subgroup_orbits")

    second = get_all_subgroup_orbits(fg)
    assert second == first
