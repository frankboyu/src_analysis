source /group/halld/Software/build_scripts/gluex_env_boot_jlab.csh
gxenv $HALLD_VERSIONS/version.xml

setenv CCDB_USER boyu
setenv JANA_CALIB_URL "sqlite:///ccdb.sqlite"
setenv CCDB_CONNECTION "sqlite:///ccdb.sqlite"
