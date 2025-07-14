from skbuild import setup

setup(
    name="libcasm-configuration",
    version="2.2.0",
    packages=[
        "libcasm",
        "libcasm.clusterography",
        "libcasm.configuration",
        "libcasm.configuration.io",
        "libcasm.enumerate",
        "libcasm.irreps",
        "libcasm.local_configuration",
        "libcasm.occ_events",
        "libcasm.sym_info",
    ],
    package_dir={"": "python"},
    cmake_install_dir="python/libcasm",
    include_package_data=False,
)
