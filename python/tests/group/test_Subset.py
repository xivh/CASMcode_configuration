import libcasm.configuration as casmconfig
import libcasm.group as casmgroup
import libcasm.xtal.prims as xtal_prims


def test_subset_construction():
    # Define a simple lattice primitive
    prim = casmconfig.Prim(xtal_prim=xtal_prims.FCC(a=4.0, occ_dof=["A", "B"]))

    assert len(prim.factor_group.elements) == 48

    subset = casmgroup.Subset(group=prim.factor_group, indices=set(range(48)))
    assert isinstance(subset, casmgroup.Subset)

    assert subset.is_group
    assert subset.is_normal

    print("Conjugacy classes:")
    conjugacy_classes = prim.factor_group.conjugacy_classes()
    for i, cc in enumerate(conjugacy_classes):
        print(f"{i}: Indices: {cc}")
    print()

    def find_class(j):
        for i, cc in enumerate(conjugacy_classes):
            if j in cc:
                return i
        return -1

    maximal_cyclic_subgroups = subset.maximal_cyclic_subgroups()
    generators = subset.maximal_cyclic_generators()
    assert len(maximal_cyclic_subgroups) == len(generators)
    for i, s in enumerate(maximal_cyclic_subgroups):
        gen = generators[i]
        classes = [find_class(j) for j in s.indices]

        print(f"{i}: Generator: {gen}, Indices: {s.indices}, Classes: {classes}")
        # for j in s.indices:
        #     op = s.group.elements[j]
        #     info = xtal.SymInfo(op, prim.xtal_prim.lattice())
        #     print(f"- Element {j}: {info.brief_cart()}")
    print()

    generators = subset.minimal_generators()
    print("Minimal generators:", generators)
    print()

    generators = {1, 32}
    subset = casmgroup.Subset.from_generators(
        group=prim.factor_group,
        generators=generators,
    )
    print("Generators:", generators)
    print("Indices:", subset.indices)
    print("Size:", len(subset.indices))
    print()
