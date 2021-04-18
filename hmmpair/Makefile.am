DEFAULT_INCLUDES=-I. -I@top_srcdir@/@INFERNAL07_EASELDIR@ -I@top_srcdir@/@INFERNAL07_SQUIDDIR@ -I@top_srcdir@/@INFERNAL07_SRCDIR@ -I@top_srcdir@/@VIENNA171_DIR@
LDADD=@top_srcdir@/@INFERNAL07_SRCDIR@/libinfernal.a @top_srcdir@/@INFERNAL07_EASELDIR@/libeasel.a @top_srcdir@/@INFERNAL07_SQUIDDIR@/libsquid.a @top_srcdir@/@VIENNA171_DIR@/libRNA.a -lm
bin_PROGRAMS=hmmpair
hmmpair_CPPFLAGS=-DSRE_STRICT_ANSI -DHAS_MALLOC_OR_DIE -DHMMPAIR
hmmpair_SOURCES=AlignmentConsensusAndScoring.cpp Cm2HMM.cpp Cm2HmmOptimize.cpp CommaSepFileReader.cpp CovarianceModel.cpp GSCConsensus.cpp hmmpair.cpp HmmType1.cpp InfernalHmm.cpp MarkovModelStats.cpp MiscExceptions.cpp MLHeuristic.cpp NaryCounter.cpp Optimize.cpp ScanHmm.cpp SymbolicMath.cpp PositionsToIgnore.cpp

