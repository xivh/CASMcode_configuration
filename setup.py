from skbuild import setup

setup(
    name="libcasm-configuration",
    version="3.0a1",
    packages=[
        "libcasm",
        "libcasm.clusterography",
        "libcasm.configuration",
        "libcasm.configuration.io",
        "libcasm.enumerate",
        "libcasm.group",
        "libcasm.irreps",
        "libcasm.local_configuration",
        "libcasm.occ_events",
        "libcasm.sym_info",
    ],
    package_dir={"": "python"},
    cmake_install_dir="python/libcasm",
    include_package_data=False,
)
