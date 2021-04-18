CFLAGS    = @CFLAGS@ @PTHREAD_CFLAGS@ @SIMD_CFLAGS@
DEFAULT_INCLUDES=-I@top_srcdir@/@INFERNAL_ESLDIR@ -I@top_srcdir@/@INFERNAL_HMMERDIR@/src -I@top_srcdir@/@INFERNAL_SRCDIR@ -I@top_srcdir@/@VIENNA_DIR@
LDADD=@top_srcdir@/@INFERNAL_SRCDIR@/libinfernal.a @top_srcdir@/@INFERNAL_HMMERDIR@/src/libhmmer.a @top_srcdir@/@INFERNAL_ESLDIR@/libeasel.a @top_srcdir@/@VIENNA_DIR@/libvienna.a @LIBS@ @PTHREAD_LIBS@ -lm
bin_PROGRAMS=cmfinder04
cmfinder04_SOURCES=base_pair_prediction.c cmfinder.c from_cmbuild.c from_cmsearch.c cand.c em.c from_cmcalibrate.c global.c
